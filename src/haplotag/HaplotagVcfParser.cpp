#include "HaplotagVcfParser.h"

/**
 * @brief Constructor - Initialize VCF parser with sample type
 * @param tagSample Sample type (NORMAL or TUMOR) to process
 */
VcfParser::VcfParser(Genome tagSample): tagSample(tagSample){
    this->tagSample = tagSample;
    this->mode = VCF_PARSER_LOAD_NODE;
    this->resultVcf = nullptr;
    this->commandline = "";
    this->version = "";
    reset();
}

VcfParser::~VcfParser(){

}

/**
 * @brief Reset parser state and clear internal data structures
 * 
 * Called between parsing different file types to ensure clean state
 */
void VcfParser::reset(){
    parseSnpFile = false;
    parseSVFile = false;
    parseMODFile = false;
    integerPS = false;
    psIndex.clear();
    writeCommandline = false;
}

/**
 * @brief Main entry point for parsing VCF files
 * 
 * Automatically detects file format (.vcf or .vcf.gz) and delegates to appropriate parser
 * 
 * @param variantFile Input VCF file path
 * @param Info VCF metadata and sample information
 * @param chrMultiVariants Output container for parsed variants
 */
void VcfParser::parsingVCF(std::string &variantFile, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants){
    mode = VCF_PARSER_LOAD_NODE;
    if( variantFile.find("gz") != std::string::npos ){
        // .vcf.gz 
        compressParser(variantFile, Info, chrMultiVariants);
    }
    else if( variantFile.find("vcf") != std::string::npos ){
        // .vcf
        unCompressParser(variantFile, Info, chrMultiVariants);
    }
    else{
        std::cerr<<"file: "<< variantFile << "\nnot vcf file. please check filename extension\n";
        exit(EXIT_FAILURE);
    }
    return;
}

/**
 * @brief Write processed variant data to output VCF file
 * 
 * Creates output file with suffix "_sc.vcf" and writes processed variants
 * with somatic variant annotations
 * 
 * @param variantFile Input VCF file path (for header information)
 * @param Info VCF metadata and sample information
 * @param chrMultiVariants Input container with processed variants
 * @param outputPrefix Output file prefix
 */
void VcfParser::writingResultVCF(
    std::string &variantFile,
    VCF_Info &Info,
    std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants,
    const std::string &outputPrefix
){
    mode = VCF_PARSER_WRITE_NODE;

    resultVcf = new std::ofstream(outputPrefix + "_sc.vcf");
    if(!resultVcf->is_open()){
        std::cerr<<"Fail to open output file: " << outputPrefix << "_sc.vcf\n";
        exit(EXIT_FAILURE);
    }
    
    if( variantFile.find("gz") != std::string::npos ){
        // .vcf.gz 
        compressParser(variantFile, Info, chrMultiVariants);
    }
    else if( variantFile.find("vcf") != std::string::npos ){
        // .vcf
        unCompressParser(variantFile, Info, chrMultiVariants);
    }
    else{
        std::cerr<<"file: "<< variantFile << "\nnot vcf file. please check filename extension\n";
        exit(EXIT_FAILURE);
    }
    return;

    resultVcf->close();
    delete resultVcf;
    resultVcf = nullptr;
}

/**
 * @brief Parse compressed VCF file (.vcf.gz) using zlib
 * 
 * Uses gzFile interface to read compressed VCF files with dynamic buffer allocation
 * for efficient memory usage with large files
 * 
 * @param variantFile Input compressed VCF file path
 * @param Info VCF metadata and sample information
 * @param chrMultiVariants Output container for parsed variants
 */
void VcfParser::compressParser(std::string &variantFile, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants){
    gzFile file = gzopen(variantFile.c_str(), "rb");
    if(variantFile=="")
        return;
    if(!file){
        std::cerr<< "Fail to open vcf: " << variantFile << "\n";
    }
    else{  
        int buffer_size = 1048576; // 1M
        char* buffer = (char*) malloc(buffer_size);
        if(!buffer){
            std::cerr<<"Failed to allocate buffer\n";
            exit(EXIT_FAILURE);
        }
        char* offset = buffer;
            
        while(true) {
            int len = buffer_size - (offset - buffer);
            if (len == 0){
                buffer_size *= 2; // Double the buffer size
                char* new_buffer = (char*) realloc(buffer, buffer_size);
                if(!new_buffer){
                    std::cerr<<"Failed to allocate buffer\n";
                    free(buffer);
                    exit(EXIT_FAILURE);
                }
                buffer = new_buffer;
                offset = buffer + buffer_size / 2; // Update the offset pointer to the end of the old buffer
                len = buffer_size - (offset - buffer);
            }

            len = gzread(file, offset, len);
            if (len == 0) break;    
            if (len <  0){ 
                int err;
                fprintf (stderr, "Error: %s.\n", gzerror(file, &err));
                exit(EXIT_FAILURE);
            }

            char* cur = buffer;
            char* end = offset+len;
            for (char* eol; (cur<end) && (eol = std::find(cur, end, '\n')) < end; cur = eol + 1)
            {
                std::string input = std::string(cur, eol);
                processLine(input, Info, chrMultiVariants);
            }
            // any trailing data in [eol, end) now is a partial line
            offset = std::copy(cur, end, buffer);
        }
        gzclose (file);
        free(buffer);
    }    
}

/**
 * @brief Parse uncompressed VCF file (.vcf)
 * 
 * Reads uncompressed VCF files line by line using standard file I/O
 * 
 * @param variantFile Input uncompressed VCF file path
 * @param Info VCF metadata and sample information
 * @param chrMultiVariants Output container for parsed variants
 */
void VcfParser::unCompressParser(std::string &variantFile, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants){
    std::ifstream originVcf(variantFile);
    if(variantFile=="")
        return;
    if(!originVcf.is_open()){
        std::cerr<< "Fail to open vcf: " << variantFile << "\n";
        exit(1);
    }
    else{
        std::string input;
        while(! originVcf.eof() ){
            std::getline(originVcf, input);
            processLine(input, Info, chrMultiVariants);
        }
    }
}

/**
 * @brief Route VCF line processing based on current operation mode
 * 
 * Delegates line processing to either parserProcess (for loading) or writeProcess (for writing)
 * 
 * @param input VCF line content
 * @param Info VCF metadata and sample information
 * @param chrMultiVariants Variant data container
 */
void VcfParser::processLine(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants){
    switch(mode){
        case VCF_PARSER_LOAD_NODE:
            parserProcess(input, Info, chrMultiVariants);
            break;
        case VCF_PARSER_WRITE_NODE:
            writeProcess(input, Info, chrMultiVariants);
            break;
        default:
            std::cerr<< "[ERROR](VcfParser::processLine): invalid mode\n";
            exit(EXIT_FAILURE);
    }
}

/**
 * @brief Parse and load variant data from VCF line
 * 
 * Handles different types of VCF lines:
 * - Header lines (##): Extract contig information and PS field type
 * - Comment lines (#): Skip
 * - Variant lines: Parse phased SNPs, SVs, and MODs based on enabled flags
 * 
 * For phased SNPs, extracts:
 * - Haplotype information (HP1, HP2)
 * - Phase set (PS) information
 * - Genotype and variant type
 * 
 * For SVs and MODs, tracks read-haplotype associations
 * 
 * @param input VCF line content
 * @param Info VCF metadata and sample information
 * @param chrMultiVariants Output container for parsed variants
 */
void VcfParser::parserProcess(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants){
    if( input.substr(0, 2) == "##" && parseSnpFile){
        if( input.find("contig=")!= std::string::npos ){
            // Extract chromosome information from contig header
            int id_start  = input.find("ID=")+3;
            int id_end    = input.find(",length=");
            int len_start = id_end+8;
            int len_end   = input.find(">");
            
            std::string chr = input.substr(id_start,id_end-id_start);
            int chrLen = std::stoi( input.substr(len_start,len_end-len_start) );

            Info.chrVec.push_back(chr);
            Info.chrLength[chr]=chrLen;                
        }
        if( input.substr(0, 16) == "##FORMAT=<ID=PS," ){
            // Determine PS field type (Integer or String)
            if( input.find("Type=Integer")!= std::string::npos ){
                integerPS = true;
            }
            else if( input.find("Type=String")!= std::string::npos ){
                integerPS = false;
                std::cerr<< "PS type is String. Auto index to integer ... ";
            }
            else{
                std::cerr<< "[ERROR](VcfParser::processLine) => not found PS type (Type=Integer or Type=String).\n"; 
                exit(EXIT_SUCCESS);
            }
        }
    }
    else if ( input.substr(0, 1) == "#" ){
        // Skip comment lines
    }
    else{
        // Parse variant data line
        std::istringstream iss(input);
        std::vector<std::string> fields((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());

        if( fields.size() == 0 )
            return;
        else if(fields.size() < 10){
            std::cerr << "[ERROR](VcfParser::parserProcess) => VCF file format not supported: " << input << std::endl;
            exit(EXIT_FAILURE);
        }

        // Convert to 0-based position
        int pos = std::stoi( fields[1] ) - 1;
        std::string chr = fields[0];
                
        // find GT flag colon
        int colon_pos = 0;
        int gt_pos = fields[8].find("GT");
        for(int i =0 ; i< gt_pos ; i++){
            if(fields[8][i]==':')
                colon_pos++;
        }
        // find GT value start
        int current_colon = 0;
        int modifu_start = 0;
        for(unsigned int i =0; i < fields[9].length() ; i++){
            if( current_colon >= colon_pos )
                break;
            if(fields[9][i]==':')
                current_colon++;  
            modifu_start++;
        }

        // phased hetero GT
        if( (fields[9][modifu_start] != fields[9][modifu_start+2]) && fields[9][modifu_start+1] == '|' ){
            // find PS flag
            colon_pos = 0;
            int ps_pos = fields[8].find("PS");
            for(int i =0 ; i< ps_pos ; i++){
                if(fields[8][i]==':')
                    colon_pos++;
            }

            // find PS value start
            current_colon = 0;
            int ps_start = 0;
            for(unsigned int i =0; i < fields[9].length() ; i++){
                if( current_colon >= colon_pos )
                    break;
                if(fields[9][i]==':')
                    current_colon++;  
                ps_start++;
            }
            
            std::string psValue;
            // get PS value
            if( fields[9].find(":",ps_start+1) != std::string::npos ){
                int ps_end_pos = fields[9].find(":",ps_start+1);
                psValue = fields[9].substr(ps_start, ps_end_pos - ps_start);
            }
            else{
                psValue = fields[9].substr(ps_start, fields[9].length() - ps_start);
            }
            
            // snp file
            if( parseSnpFile ){

                VarData varData;
                varData.allele.Ref = fields[3];
                varData.allele.Alt = fields[4];
                varData.GT = GenomeType::PHASED_HETERO;
                varData.setVariantType();
                
                // Handle PS value (integer or string)
                if(integerPS){
                    varData.PhasedSet = std::stoi(psValue);
                }
                else{
                    // Map string PS values to integer indices
                    std::map<std::string, int>::iterator psIter = psIndex.find(psValue);
                    
                    if( psIter == psIndex.end() ){
                        psIndex[psValue] = psIndex.size();
                    }
                    varData.PhasedSet = psIndex[psValue];
                }
                
                // record haplotype allele
                if( fields[9][modifu_start] == '0' && fields[9][modifu_start+2] == '1' ){
                    varData.HP1 = fields[3];
                    varData.HP2 = fields[4];
                }
                else if( fields[9][modifu_start] == '1' && fields[9][modifu_start+2] == '0' ){
                    varData.HP1 = fields[4];
                    varData.HP2 = fields[3];
                }

                // Store variant data for appropriate sample
                if(Info.sample == NORMAL){
                    chrMultiVariants[chr][pos].Variant[NORMAL] = varData;
                }else if(Info.sample == TUMOR){
                    chrMultiVariants[chr][pos].Variant[TUMOR] = varData;
                }
            }
            // sv file
            if( parseSVFile ){
                // get read INFO
                int read_pos = fields[7].find("RNAMES=");
                read_pos = fields[7].find("=",read_pos);
                read_pos++;
                        
                int next_field = fields[7].find(";",read_pos);
                std::string totalRead = fields[7].substr(read_pos,next_field-read_pos);
                std::stringstream totalReadStream(totalRead);
                
                int svHaplotype;
                // In which haplotype does SV occur
                if( fields[9][modifu_start] == '0' && fields[9][modifu_start+2] == '1' ){
                    svHaplotype = 1;
                }
                else if( fields[9][modifu_start] == '1' && fields[9][modifu_start+2] == '0' ){
                    svHaplotype = 0;
                }
                
                // Track read-haplotype associations for SVs
                std::string read;
                while(std::getline(totalReadStream, read, ','))
                {
                   auto readIter = Info.readSVHapCount.find(read);
                   if(readIter==Info.readSVHapCount.end()){
                       Info.readSVHapCount[read][0]=0;
                       Info.readSVHapCount[read][1]=0;
                   }
                   Info.readSVHapCount[read][svHaplotype]++;
                }
                
            }
            // mod file
            if( parseMODFile ){
                // get read INFO
                int read_pos = fields[7].find("MR=");
                read_pos = fields[7].find("=",read_pos);
                read_pos++;
                        
                int next_field = fields[7].find(";",read_pos);
                std::string totalRead = fields[7].substr(read_pos,next_field-read_pos);
                std::stringstream totalReadStream(totalRead);
                
                int modHaplotype;
                // Determine which haplotype the modification occurs in
                if( fields[9][modifu_start] == '0' && fields[9][modifu_start+2] == '1' ){
                    modHaplotype = 1;
                }
                else if( fields[9][modifu_start] == '1' && fields[9][modifu_start+2] == '0' ){
                    modHaplotype = 0;
                }
                
                // Track read-haplotype associations for modifications
                std::string read;
                while(std::getline(totalReadStream, read, ','))
                {
                   auto readIter = Info.readSVHapCount.find(read);
                   if(readIter==Info.readSVHapCount.end()){
                       Info.readSVHapCount[read][0]=0;
                       Info.readSVHapCount[read][1]=0;
                   }
                   Info.readSVHapCount[read][modHaplotype]++;
                }
            }
        }
        // record unphased tumor SNPs
        else if((tagSample == Genome::TUMOR)){
            //homozygous SNPs
            if( fields[9][modifu_start] == '1' && fields[9][modifu_start+1] == '/' && fields[9][modifu_start+2] == '1' ){
                if(parseSnpFile){

                    VarData varData;
                    varData.allele.Ref = fields[3];
                    varData.allele.Alt = fields[4];
                    varData.GT = GenomeType::UNPHASED_HOMO;
                    varData.setVariantType();

                    if(Info.sample == NORMAL){
                        chrMultiVariants[chr][pos].Variant[NORMAL] = varData;
                    }else if(Info.sample == TUMOR){
                        chrMultiVariants[chr][pos].Variant[TUMOR] = varData;
                    }
                }
            //unphased heterozygous
            }else if( fields[9][modifu_start] == '0' && fields[9][modifu_start+1] == '/' && fields[9][modifu_start+2] == '1' ){
                if(parseSnpFile){

                    VarData varData;
                    varData.allele.Ref = fields[3];
                    varData.allele.Alt = fields[4];

                    varData.GT = GenomeType::UNPHASED_HETERO;
                    varData.setVariantType();

                    if(Info.sample == NORMAL){
                        chrMultiVariants[chr][pos].Variant[NORMAL] = varData;
                    }else if(Info.sample == TUMOR){
                        chrMultiVariants[chr][pos].Variant[TUMOR] = varData;
                    }
                }
            }
        }
    }
}

/**
 * @brief Write processed variant data to output VCF file
 * 
 * Handles different types of VCF lines during writing:
 * - Header lines (##): Copy as-is
 * - Column header (#CHROM): Add version and command line headers
 * - Variant lines: Apply somatic variant filtering and quality updates
 * 
 * For somatic variants, updates FILTER field based on somatic status
 * 
 * @param input VCF line content
 * @param Info VCF metadata and sample information
 * @param chrMultiVariants Input container with processed variants
 */
void VcfParser::writeProcess(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants){
    if(resultVcf == nullptr){
        std::cerr<< "[ERROR](VcfParser::writeProcess): resultVcf is nullptr\n";
        exit(EXIT_FAILURE);
    }

    if(input.length() >= 2 && input.substr(0, 2) == "##"){
        // Copy header lines as-is
        *resultVcf << input << std::endl;
    }else if (input.length() >= 6 && (input.substr(0, 6) == "#CHROM" || input.substr(0, 6) == "#chrom")){
        // Add version and command line headers before column header
        if( writeCommandline == false ){
            *resultVcf << "##longphase_s_version=" << version << std::endl;
            *resultVcf << "##commandline=" << commandline << std::endl;
            writeCommandline = true;
        }    
        *resultVcf << input << std::endl;
    }else{
        // Process variant data lines
        std::istringstream iss(input);
        std::vector<std::string> fields((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());

        if( fields.size() == 0 )
            return;

        if (fields.size() >= 7){
            // Convert to 0-based position
            int pos = std::stoi( fields[1] ) - 1;
            std::string chr = fields[0];

            // Check if variant exists in processed data
            auto chrIt = chrMultiVariants.find(chr);
            if (chrIt != chrMultiVariants.end()) {

                auto posIt = chrIt->second.find(pos);
                if (posIt != chrIt->second.end()) {
                    
                    auto varIt = posIt->second.Variant.find(TUMOR);
                    if (varIt != posIt->second.Variant.end()) {
                        // Only process SNPs for somatic variant calling
                        if(varIt->second.variantType == HaplotagVariantType::SNP || varIt->second.variantType == HaplotagVariantType::INSERTION || varIt->second.variantType == HaplotagVariantType::DELETION){
                            // Update FILTER field based on somatic status
                            if (posIt->second.isSomaticVariant) {
                                if (fields[6] != "PASS") {
                                    fields[6] = "PASS";
                                }
                            } else {
                                if (fields[6] == "PASS") {
                                    fields[6] = "LowQual";
                                }
                            }
                            // Write updated variant line
                            std::string output = fields[0];
                            for (size_t i = 1; i < fields.size(); ++i) {
                                output += "\t" + fields[i];
                            }
                            *resultVcf << output << std::endl;
                        }
                    }
                }
            }
        }else{
            std::cerr << "[ERROR](VcfParser::writeProcess) => VCF file format error: " << input << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}

/**
 * @brief Set parser operation mode
 * @param mode VCF_PARSER_LOAD_NODE or VCF_PARSER_WRITE_NODE
 */
void VcfParser::setMode(VcfParserMode mode){
    this->mode = mode;
}

/**
 * @brief Enable/disable SNP variant parsing
 * @param parseSnpFile True to enable SNP parsing
 */
void VcfParser::setParseSnpFile(bool parseSnpFile){
    this->parseSnpFile = parseSnpFile;
}

/**
 * @brief Enable/disable structural variant parsing
 * @param parseSVFile True to enable SV parsing
 */
void VcfParser::setParseSVFile(bool parseSVFile){
    this->parseSVFile = parseSVFile;
} 

/**
 * @brief Enable/disable modification variant parsing
 * @param parseMODFile True to enable MOD parsing
 */
void VcfParser::setParseMODFile(bool parseMODFile){
    this->parseMODFile = parseMODFile;
}

/**
 * @brief Get current SNP parsing status
 * @return True if SNP parsing is enabled
 */
bool VcfParser::getParseSnpFile(){
    return this->parseSnpFile;
}

/**
 * @brief Set command line string for VCF header
 * @param commandline Command line string to include in output VCF
 */
void VcfParser::setCommandLine(std::string &commandline){
    this->commandline = commandline;
}

/**
 * @brief Set version string for VCF header
 * @param version Version string to include in output VCF
 */
void VcfParser::setVersion(std::string &version){
    this->version = version;
}