#include "SomaticBenchmark.h"

/**
 * @brief Constructor - Initialize SomaticReadVerifier with testing flag and metrics
 * @param openTestingFunc Flag to enable testing functionality
 * @param metrics Pointer to metrics for data collection
 */
SomaticReadVerifier::SomaticReadVerifier(bool openTestingFunc, SomaticReadMetrics *metrics): 
openTestingFunc(openTestingFunc),
metrics(metrics)
{}

/**
 * @brief Destructor - Clean up resources
 */
SomaticReadVerifier::~SomaticReadVerifier(){

}

/**
 * @brief Record deletion read count for a variant position
 * 
 * Records deletion count for truth somatic variants and adds position
 * to truth somatic position vector for tracking reads crossing truth somatic SNPs
 * 
 * @param chr Chromosome name
 * @param currentVariantIter Iterator to current variant
 */
void SomaticReadVerifier::recordDelReadCount(const std::string &chr, std::map<int, MultiGenomeVar>::iterator &currentVariantIter){
    if(!openTestingFunc) return;
    
    if(currentVariantIter->second.isExists(TRUTH_SOMATIC)){
        int pos = (*currentVariantIter).first;
        metrics->posAltRefDelCount[pos].delCount++;

        //record somatic position for record crossing truth somatic snp read
        metrics->truthSomaticPosVec.push_back(std::make_pair(pos, SnpHP::NONE_SNP));
    }
}

/**
 * @brief Record reference/alternate allele count for a variant position
 * 
 * Records reference and alternate allele counts for truth somatic variants.
 * Determines base haplotype based on whether the base matches reference or alternate allele.
 * 
 * @param chr Chromosome name
 * @param base Base at the position
 * @param currentVariantIter Iterator to current variant
 */
void SomaticReadVerifier::recordRefAltAlleleCount(const std::string &chr, std::string &base, std::map<int, MultiGenomeVar>::iterator &currentVariantIter){
    if(!openTestingFunc) return;

    if(currentVariantIter->second.isExists(TRUTH_SOMATIC)){
        int pos = currentVariantIter->first;
        std::string& refAllele = currentVariantIter->second.Variant[TRUTH_SOMATIC].allele.Ref;
        std::string& altAllele = currentVariantIter->second.Variant[TRUTH_SOMATIC].allele.Alt;

        int baseHP = SnpHP::NONE_SNP;

        if(base == refAllele){
            metrics->posAltRefDelCount[pos].refCount++;
        }else if(base == altAllele){
            metrics->posAltRefDelCount[pos].altCount++;
            baseHP = SnpHP::SOMATIC_H3;
        }
        // record somatic position for tracking reads crossing truth somatic SNPs
        metrics->truthSomaticPosVec.push_back(std::make_pair(pos, baseHP));
    }
}

/**
 * @brief Create basic somatic read log entry
 * 
 * Creates a basic SomaticReadLog object with chromosome, read ID, haplotype result,
 * similarity scores, and SNP counts for germline and tumor variants.
 * 
 * @param chr Chromosome name
 * @param readID Read identifier
 * @param hpResult Haplotype result
 * @param norHPsimilarity Normal haplotype similarity
 * @param deriveByHpSimilarity Derived haplotype similarity
 * @param hpCount Haplotype count map
 * @return SomaticReadLog object
 */
SomaticReadLog SomaticReadVerifier::createBasicSomaticReadLog(const std::string &chr, std::string &readID, int &hpResult, double &norHPsimilarity, float &deriveByHpSimilarity, std::map<int, int> &hpCount){
    SomaticReadLog tmp;
    tmp.chr = chr;
    tmp.readID = readID;
    tmp.hpResult = hpResult;
    tmp.germlineVarSimilarity = norHPsimilarity;
    tmp.deriveByHpSimilarity = deriveByHpSimilarity;
    tmp.germlineSnpCount = hpCount[SnpHP::GERMLINE_H1] + hpCount[SnpHP::GERMLINE_H2];
    tmp.tumorSnpCount = hpCount[SnpHP::SOMATIC_H3];

    return tmp;
}

/**
 * @brief Record reads crossing truth somatic somatic SNPs
 * 
 * Records reads that cross truth somatic somatic variants and corrects
 * haplotype assignments based on the presence of somatic variants.
 * 
 * For reads with somatic H3 variants:
 * - H1 -> H1_1, H2 -> H2_1, unTag -> H3
 * 
 * For reads without somatic H3 variants:
 * - H2_1 -> H2, H1_1 -> H1, H3 -> unTag
 * 
 * @param chr Chromosome name
 * @param readID Read identifier
 * @param hpResult Haplotype result
 * @param variantsHP Map of variant positions to haplotypes
 * @param hpCount Haplotype count map
 * @param norHPsimilarity Normal haplotype similarity
 * @param deriveByHpSimilarity Derived haplotype similarity
 * @param currentChrVariants Current chromosome variants
 */
void SomaticReadVerifier::recordCrossingTruthSomaticSnpRead(const std::string &chr, std::string &readID, int &hpResult, std::map<int, int> &variantsHP, std::map<int, int> &hpCount, double &norHPsimilarity, float &deriveByHpSimilarity, std::map<int, MultiGenomeVar> &currentChrVariants){
    // If testing function is not enabled, return
    if(!openTestingFunc) return;

    SomaticReadLog tmp = createBasicSomaticReadLog(chr, readID, hpResult, norHPsimilarity, deriveByHpSimilarity, hpCount);

    bool isCrossHighConSomatic = false;
    bool existHighConVaraints = false;

    // Process truth somatic positions
    for(auto varIter : metrics->truthSomaticPosVec){
        int pos = varIter.first;
        int baseHP = varIter.second;

        // check if this is a somatic H3 variant
        if(baseHP == SnpHP::SOMATIC_H3){
            existHighConVaraints = true;
        }

        tmp.somaticSnpHp[pos] = baseHP;
        isCrossHighConSomatic = true;
    }

    if(isCrossHighConSomatic){
        // correct haplotype result based on presence of truth somatic variants
        if(existHighConVaraints){
            // correct haplotype result for reads with truth somatic variants
            if(hpResult == ReadHP::H1){
                tmp.hpResult = ReadHP::H1_1;
            }else if(hpResult == ReadHP::H2){
                tmp.hpResult = ReadHP::H2_1;
            }else if(hpResult == ReadHP::unTag){
                tmp.hpResult = ReadHP::H3;
            }
        }else{
            // correct haplotype result for reads without truth somatic variants
            if(hpResult == ReadHP::H2_1){
                tmp.hpResult = ReadHP::H2;
            }else if(hpResult == ReadHP::H1_1){
                tmp.hpResult = ReadHP::H1;
            }else if(hpResult == ReadHP::H3){
                tmp.hpResult = ReadHP::unTag;
            }
        }
    }

    if(isCrossHighConSomatic){
        metrics->coverTruthSomaticPosReadVec.push_back(tmp);
    }
    
    // clear truth somatic somatic positions for next read
    if(!metrics->truthSomaticPosVec.empty()){
        metrics->truthSomaticPosVec.clear();
    }
}

/**
 * @brief Record tagged reads for performance evaluation
 * 
 * Records tagged reads and checks if they contain somatic variants.
 * Only records reads that are tagged (not unTag) and contain somatic variants.
 * 
 * @param chr Chromosome name
 * @param readID Read identifier
 * @param hpResult Haplotype result
 * @param variantsHP Map of variant positions to haplotypes
 * @param hpCount Haplotype count map
 * @param norHPsimilarity Normal haplotype similarity
 * @param deriveByHpSimilarity Derived haplotype similarity
 * @param currentChrVariants Current chromosome variants
 */
void SomaticReadVerifier::recordTaggedRead(const std::string &chr, std::string &readID, int &hpResult, std::map<int, int> &variantsHP, std::map<int, int> &hpCount, double &norHPsimilarity, float &deriveByHpSimilarity, std::map<int, MultiGenomeVar> &currentChrVariants){
    // If testing function is not enabled or read is untagged, return
    if(!openTestingFunc || hpResult == ReadHP::unTag) return;

    SomaticReadLog tmp = createBasicSomaticReadLog(chr, readID, hpResult, norHPsimilarity, deriveByHpSimilarity, hpCount);

    bool readExistHighConSomatic = false;

    // check if read contains somatic variants
    auto varIter = variantsHP.begin();
    while(varIter != variantsHP.end()){
        int pos = varIter->first;
        int snpHP = varIter->second;
        if(currentChrVariants.find(pos) != currentChrVariants.end()){
            if(currentChrVariants[pos].isExists(TRUTH_SOMATIC) && (snpHP == SnpHP::SOMATIC_H3)){
                tmp.somaticSnpHp[pos] = snpHP;
                readExistHighConSomatic = true;
            }
        }
        varIter++;
    }

    if(readExistHighConSomatic){
        metrics->taggedSomaticReadVec.push_back(tmp);
    }

    metrics->totalReadVec.push_back(tmp);
}

/**
 * @brief Constructor - Initialize SomaticReadBenchmark with benchmark files and threshold
 * @param benchmarkVcf Path to benchmark VCF file
 * @param benchmarkBed Path to benchmark BED file
 * @param mappingQualityThreshold Mapping quality threshold
 */
SomaticReadBenchmark::SomaticReadBenchmark(std::string benchmarkVcf, std::string benchmarkBed, int mappingQualityThreshold){
    setParseSnpFile(true);
    openTestingFunc = false;
    loadedBedFile = false;
    this->benchmarkVcf = benchmarkVcf;
    this->benchmarkBed = benchmarkBed;
    this->mappingQualityThreshold = mappingQualityThreshold;
}

/**
 * @brief Destructor - Clean up resources
 */
SomaticReadBenchmark::~SomaticReadBenchmark(){

}

/**
 * @brief Enable or disable testing functionality
 * @param openTestingFunc True to enable testing
 */
void SomaticReadBenchmark::setEnabled(bool openTestingFunc){
    this->openTestingFunc = openTestingFunc;
}

/**
 * @brief Check if testing functionality is enabled
 * @return True if testing is enabled
 */
bool SomaticReadBenchmark::isEnabled(){
    return openTestingFunc;
}

/**
 * @brief Check if BED file is loaded
 * @return True if BED file is loaded
 */
bool SomaticReadBenchmark::isLoadBedFile(){
    return loadedBedFile;
}

/**
 * @brief Initialize chromosome key for multi-threaded access
 * @param chr Chromosome name
 */
void SomaticReadBenchmark::loadChrKey(const std::string &chr){
    chrMetrics[chr] = SomaticReadMetrics();
}

/**
 * @brief Load truth somatic VCF file
 * @param input Input VCF file path
 * @param Info VCF metadata and sample information
 * @param chrMultiVariants Output container for parsed variants
 */
void SomaticReadBenchmark::loadTruthSomaticVCF(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants){
    if(!openTestingFunc) return;
    parsingVCF(input, Info, chrMultiVariants);
}

/**
 * @brief Override parserProcess to handle truth somatic VCF parsing
 * 
 * Parses truth somatic VCF files and extracts variant information.
 * Only processes SNP variants and stores them as TRUTH_SOMATIC variants.
 * 
 * @param input VCF line content
 * @param Info VCF metadata and sample information
 * @param chrMultiVariants Output container for parsed variants
 */
void SomaticReadBenchmark::parserProcess(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants){
    if(!openTestingFunc) return;
    
    if( input.substr(0, 2) == "##" && getParseSnpFile()){
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
    }
    else if ( input.substr(0, 1) == "#" ){
        // Skip comment lines
    }
    else{
        // Parse variant data line
        std::istringstream iss(input);
        std::vector<std::string> fields((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());
        if( fields.size() == 0 ){
            return;
        }else if(fields.size() < 5){
            std::cerr << "[ERROR] (SomaticReadBenchmark::parserProcess) => VCF file format not supported: " << input << std::endl;
            exit(EXIT_FAILURE);
        }
            
        // Convert to 0-based position
        int pos = std::stoi( fields[1] ) - 1;
        std::string chr = fields[0];

        // Create variant data and store as TRUTH_SOMATIC
        VarData varData;
        varData.allele.Ref = fields[3];
        varData.allele.Alt = fields[4];
        chrMultiVariants[chr][pos].Variant[Genome::TRUTH_SOMATIC] = varData;
    }
}

/**
 * @brief Parse benchmark BED file
 * 
 * Parses BED files containing regions of interest for benchmarking.
 * Currently only supports uncompressed BED files (.bed).
 * 
 * @param bedFile Path to BED file
 */
void SomaticReadBenchmark::parseBedFile(const std::string& bedFile) {
    if(!openTestingFunc) return;

    if( bedFile.find("bed.gz") != std::string::npos){
        std::cerr << "[WARNING] BED .gz files are not supported. Please use an uncompressed .bed file: " << bedFile << std::endl;
        return;
    }else if( bedFile.find("bed") != std::string::npos){
        std::ifstream file(bedFile);
        if (!file.is_open()) {
            std::cerr << "[ERROR] Failed to open BED file: " << bedFile << std::endl;
            return;
        }
        
        std::string line;
        bool valid = true;
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;

            valid = processBedLine(line);
            
            if(!valid) break;
        }

        if(valid) {
            loadedBedFile = true;
        }else {
            std::cerr << "[WARNING] Failed to load BED file due to invalid format" << std::endl;
        }
    }
}

/**
 * @brief Process a single BED file line
 * 
 * Validates and processes a single BED file line, extracting chromosome,
 * start position, and end position. Validates coordinate ranges.
 * 
 * @param line BED file line content
 * @return True if line is valid, false otherwise
 */
bool SomaticReadBenchmark::processBedLine(const std::string& line) {
    std::istringstream iss(line);
    std::string chr;
    int start, end;

    // Check if the line has valid format (chr start end)
    if (!(iss >> chr >> start >> end)) {
        std::cerr << "[WARNING] Invalid BED line format: " << line << std::endl;
        return false;
    }
    
    // Check if the coordinates are valid (start >= 0, end > start)
    if (start < 0 || end <= start) {
        std::cerr << "[WARNING] Invalid BED coordinates (start=" << start << ", end=" << end << "): " << line << std::endl;
        return false;
    }
    
    // Store BED region (convert to 0-based end position)
    bedRegions[chr].push_back(BedRegion{start, end - 1});
    return true;
}

/**
 * @brief Mark variants in BED regions
 * 
 * Marks variants as being inside or outside BED regions for each chromosome.
 * Uses efficient algorithm to traverse bed regions and variants simultaneously.
 * 
 * @param chrVec Vector of chromosome names
 * @param chrMultiVariants Variant data container
 */
void SomaticReadBenchmark::markVariantsInBedRegions(std::vector<std::string> &chrVec, std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants) {
    if(!openTestingFunc || !loadedBedFile) return;
    time_t begin = time(NULL);
    std::cerr<< "[Benchmark] marking variants in bed regions ... ";
    
    for (auto& chr : chrVec) {
        auto& chrPosVariants = chrMultiVariants[chr];
        
        // check if this chromosome has bed regions
        auto bedIt = bedRegions.find(chr);
        if (bedIt == bedRegions.end()) {
            // if this chromosome has no bed regions, mark all points as false
            for (auto& curVar : chrPosVariants) {
                curVar.second.isInBedRegion = false;
                if(curVar.second.isExists(Genome::TUMOR)) {
                    variantOutBedRegionCount[Genome::TUMOR]++;
                }
                if(curVar.second.isExists(Genome::NORMAL)) {
                    variantOutBedRegionCount[Genome::NORMAL]++;
                }
                if(curVar.second.isExists(Genome::TRUTH_SOMATIC)){
                    variantOutBedRegionCount[Genome::TRUTH_SOMATIC]++;
                }
            }
            continue;
        }

        const auto& regions = bedIt->second;
        if (regions.empty()) {
            // If there are no bed regions, mark all variants as outside
            for (auto& curVar : chrPosVariants) {
                curVar.second.isInBedRegion = false;
                if(curVar.second.isExists(Genome::TUMOR)) {
                    variantOutBedRegionCount[Genome::TUMOR]++;
                }
                if(curVar.second.isExists(Genome::NORMAL)) {
                    variantOutBedRegionCount[Genome::NORMAL]++;
                }
                if(curVar.second.isExists(Genome::TRUTH_SOMATIC)){
                    variantOutBedRegionCount[Genome::TRUTH_SOMATIC]++;
                }
            }
            continue;
        }

        auto regionIter = regions.begin();
        auto varPosIter = chrPosVariants.begin();
        
        // traverse bed regions and variants simultaneously for efficiency
        while (varPosIter != chrPosVariants.end()) {
            int variantPos = varPosIter->first;
            
            // skip bed regions that end before current variant position
            while (regionIter != regions.end() && variantPos > regionIter->end) {
                ++regionIter;
            }
            
            // check if variant is in current bed region
            if (regionIter != regions.end() && 
                variantPos >= regionIter->start && 
                variantPos <= regionIter->end) {
                varPosIter->second.isInBedRegion = true;
                if(varPosIter->second.isExists(Genome::TUMOR)) {
                    variantInBedRegionCount[Genome::TUMOR]++;
                }
                if(varPosIter->second.isExists(Genome::NORMAL)) {
                    variantInBedRegionCount[Genome::NORMAL]++;
                }
                if(varPosIter->second.isExists(Genome::TRUTH_SOMATIC)){
                    variantInBedRegionCount[Genome::TRUTH_SOMATIC]++;
                }
            } else {
                varPosIter->second.isInBedRegion = false;
                if(varPosIter->second.isExists(Genome::TUMOR)) {
                    variantOutBedRegionCount[Genome::TUMOR]++;
                }
                if(varPosIter->second.isExists(Genome::NORMAL)) {
                    variantOutBedRegionCount[Genome::NORMAL]++;
                }
                if(varPosIter->second.isExists(Genome::TRUTH_SOMATIC)){
                    variantOutBedRegionCount[Genome::TRUTH_SOMATIC]++;
                }
            }
            
            ++varPosIter;
        }
    }
    std::cerr<< difftime(time(NULL), begin) << "s\n";
}

/**
 * @brief Remove variants outside BED regions
 * 
 * Removes variants that are outside BED regions from the analysis.
 * If a position has no NORMAL data, the entire position is removed.
 * Otherwise, only TUMOR and TRUTH_SOMATIC data are removed.
 * 
 * @param chrMultiVariants Variant data container
 */
void SomaticReadBenchmark::removeVariantsOutBedRegion(std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants){
    if(!openTestingFunc || !loadedBedFile) return;

    // remove variants outside bed regions
    for (auto& chrPair : chrMultiVariants) {
        auto& chrPosVariants = chrPair.second;
        
        auto varPosIter = chrPosVariants.begin();
        while (varPosIter != chrPosVariants.end()) {
            if (!varPosIter->second.isInBedRegion) {
                bool hasTumor = varPosIter->second.isExists(Genome::TUMOR);
                bool hasHighCon = varPosIter->second.isExists(Genome::TRUTH_SOMATIC);
                
                if (hasTumor || hasHighCon) {
                    // If no NORMAL data exists, remove the entire position
                    if (!varPosIter->second.isExists(Genome::NORMAL)) {
                        varPosIter = chrPosVariants.erase(varPosIter);
                        continue;
                    } else {
                        // If NORMAL data exists, only remove TUMOR and TRUTH_SOMATIC data
                        if (hasTumor) {
                            varPosIter->second.Variant.erase(Genome::TUMOR);
                        }
                        if (hasHighCon) {
                            varPosIter->second.Variant.erase(Genome::TRUTH_SOMATIC);
                        }
                        ++varPosIter;
                    }
                } else {
                    ++varPosIter;
                }
            } else {
                ++varPosIter;
            }
        }
    }
}

/**
 * @brief Write BED region log
 * 
 * Generates separate log files for variants inside and outside BED regions.
 * Creates two output files: _var_in_bed.out and _var_out_bed.out
 * 
 * @param chrVec Vector of chromosome names
 * @param chrMultiVariants Variant data container
 * @param outPrefix Output file prefix
 */
void SomaticReadBenchmark::writeBedRegionLog(const std::vector<std::string>& chrVec,
                                           const std::map<std::string, std::map<int, MultiGenomeVar>>& chrMultiVariants,
                                           const std::string& outPrefix) {
    if(!openTestingFunc || !loadedBedFile) return;
    std::ofstream inBedLog(outPrefix + "_var_in_bed.out");
    std::ofstream outBedLog(outPrefix + "_var_out_bed.out");
    
    std::string header = "#Chr\tPosition\tRef\tAlt\tVariant_Type\n";
    inBedLog << header;
    outBedLog << header;

    for (const auto& chr : chrVec) {
        auto chrPosVariants = chrMultiVariants.at(chr);

        auto varPosIter = chrPosVariants.begin();

        while(varPosIter != chrPosVariants.end()){
            int pos = varPosIter->first + 1;  // Convert to 1-based position
            
            std::string baseInfo = chr + "\t" + std::to_string(pos) + "\t";
            
            if(varPosIter->second.isExists(Genome::TUMOR)){
                auto& tumorVar = varPosIter->second.Variant.at(TUMOR);
                std::string varInfo = baseInfo + 
                                    tumorVar.allele.Ref + "\t" + 
                                    tumorVar.allele.Alt + "\t" +
                                    "TUMOR\n";

                if (varPosIter->second.isInBedRegion) {
                    inBedLog << varInfo;
                } else {
                    outBedLog << varInfo;
                }
            }

            varPosIter++;
        }
    }
    
    inBedLog.close();
    outBedLog.close();
}

/**
 * @brief Get metrics pointer for multi-threaded parallel processing
 * @param chr Chromosome name
 * @return Pointer to chromosome metrics
 */
SomaticReadMetrics* SomaticReadBenchmark::getMetricsPtr(const std::string &chr){
    return &(chrMetrics[chr]);
}

/**
 * @brief Write position allele count log
 * 
 * Generates a report showing allele counts (reference, alternate, deletion)
 * for each somatic variant position across all chromosomes.
 * 
 * @param chrVec Vector of chromosome names
 * @param outputFileName Output file name
 * @param chrMultiVariants Variant data container
 */
void SomaticReadBenchmark::writePosAlleleCountLog(
    std::vector<std::string>& chrVec,
    std::string outputFileName,
    std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants
){
    // If testing function is not enabled, return
    if(!openTestingFunc) return;

    std::ofstream *refAltCountLog=NULL;
    refAltCountLog=new std::ofstream(outputFileName);
    int totalVariantCount = 0;

    // Calculate total variant count across all chromosomes
    for(auto chr: chrVec){
        totalVariantCount += chrMetrics[chr].posAltRefDelCount.size();
    }

    if(!refAltCountLog->is_open()){
        std::cerr<< "Fail to open write file: " << outputFileName << "\n";
        exit(1);
    }else{
        (*refAltCountLog) << "#############################\n";
        (*refAltCountLog) << "# Somatic SNP allele count #\n";
        (*refAltCountLog) << "#############################\n";
        (*refAltCountLog) << "##Benchmark VCF:"  << benchmarkVcf << "\n";
        (*refAltCountLog) << "##MappingQualityThreshold:"  << mappingQualityThreshold << "\n";
        (*refAltCountLog) << "##Tatal variants:"  << totalVariantCount << "\n";
        (*refAltCountLog) << "#CHROM\t"
                          << "POS\t"
                          << "REF\t"
                          << "ALT\t"
                          << "REF_COUNT\t"
                          << "ALT_COUNT\t"
                          << "DEL_COUNT\n";
    }

    // Write allele counts for each variant position
    for(auto chr: chrVec){
        for(auto &posIter : chrMetrics[chr].posAltRefDelCount){
            (*refAltCountLog) << chr << "\t"
                              << posIter.first << "\t"
                              << chrMultiVariants[chr][posIter.first].Variant[TRUTH_SOMATIC].allele.Ref << "\t"
                              << chrMultiVariants[chr][posIter.first].Variant[TRUTH_SOMATIC].allele.Alt << "\t"
                              << posIter.second.refCount << "\t"
                              << posIter.second.altCount << "\t"
                              << posIter.second.delCount << "\n";
        }
    }

    refAltCountLog->close();
    delete refAltCountLog;
    refAltCountLog = nullptr;
}

/**
 * @brief Write tagged read report
 * 
 * Generates a report for all tagged reads across all chromosomes.
 * 
 * @param chrVec Vector of chromosome names
 * @param outputFileName Output file name
 */
void SomaticReadBenchmark::writeTaggedReadReport(
    const std::vector<std::string>& chrVec,
    std::string outputFileName
){
    // if testing function is not enabled, return
    if(!openTestingFunc) return;

    std::map<std::string, std::vector<SomaticReadLog>*> somaticReadVecMap;

    // set up chromosome somatic read vector pointers
    auto metricsIter = chrMetrics.begin();
    while(metricsIter != chrMetrics.end()){
        setChrSomaticReadVecPtr(metricsIter->first, somaticReadVecMap, metricsIter->second.totalReadVec);
        metricsIter++;
    }

    writeReadLog(chrVec, outputFileName, somaticReadVecMap);
}

/**
 * @brief Write tagged somatic read report
 * 
 * Generates a comprehensive report for tagged somatic reads including
 * performance metrics (precision, recall, F1-score) for each haplotype.
 * 
 * @param chrVec Vector of chromosome names
 * @param outputFileName Output file name
 */
void SomaticReadBenchmark::writeTaggedSomaticReadReport(
    const std::vector<std::string>& chrVec,
    std::string outputFileName
){
    // If testing function is not enabled, return
    if(!openTestingFunc) return;

    time_t begin = time(NULL);
    std::cerr << "[Benchmark] writing somatic haplotagging metrics report ... ";

    std::map<std::string, std::vector<SomaticReadLog>*> somaticReadVecMap;

    // Set up chromosome somatic read vector pointers
    auto metricsIter = chrMetrics.begin();
    while(metricsIter != chrMetrics.end()){
        setChrSomaticReadVecPtr(metricsIter->first, somaticReadVecMap, metricsIter->second.taggedSomaticReadVec);
        metricsIter++;
    }

    writeReadLog(chrVec, outputFileName, somaticReadVecMap);
    std::cerr<< difftime(time(NULL), begin) << "s\n";
}

/**
 * @brief Write total truth somatic read report
 * 
 * Generates a report for all reads that cover truth somatic positions.
 * 
 * @param chrVec Vector of chromosome names
 * @param outputFileName Output file name
 */
void SomaticReadBenchmark::writeTotalTruthSomaticReadReport(
    const std::vector<std::string>& chrVec,
    std::string outputFileName
){
    // If testing function is not enabled, return
    if(!openTestingFunc) return;

    std::map<std::string, std::vector<SomaticReadLog>*> somaticReadVecMap;

    // Set up chromosome somatic read vector pointers
    auto metricsIter = chrMetrics.begin();
    while(metricsIter != chrMetrics.end()){
        setChrSomaticReadVecPtr(metricsIter->first, somaticReadVecMap, metricsIter->second.coverTruthSomaticPosReadVec);
        metricsIter++;
    }

    writeReadLog(chrVec, outputFileName, somaticReadVecMap);
}

/**
 * @brief Set chromosome somatic read vector pointer for multi-threaded access
 * @param chr Chromosome name
 * @param somaticReadVecMap Map of chromosome names to read vector pointers
 * @param somaticReadVec Somatic read vector
 */
void SomaticReadBenchmark::setChrSomaticReadVecPtr(
    const std::string &chr,
    std::map<std::string, std::vector<SomaticReadLog>*> &somaticReadVecMap,
    std::vector<SomaticReadLog> &somaticReadVec
){
    somaticReadVecMap[chr] = &somaticReadVec;
}

/**
 * @brief Write read log to output file
 * 
 * Generates comprehensive read log with performance metrics including
 * precision, recall, and F1-score for each haplotype type.
 * 
 * @param chrVec Vector of chromosome names
 * @param outputFileName Output file name
 * @param somaticReadVecMap Map of chromosome names to read vector pointers
 */
void SomaticReadBenchmark::writeReadLog(
    const std::vector<std::string>& chrVec,
    std::string outputFileName,
    std::map<std::string, std::vector<SomaticReadLog>*> &somaticReadVecMap
){
    // if testing function is not enabled, return
    if(!openTestingFunc) return;

    std::ofstream *somaticReadLog=NULL;
    somaticReadLog=new std::ofstream(outputFileName);

    int totalReads = 0;

    // Count all truth somatic reads
    int totalTruthSomaticReads = 0;
    std::map<int, int> truthSomaticReadsMap;

    // Count all tagged truth somatic reads (TP)
    int totalTaggedTruthSomaticReads = 0;
    std::map<int, int> taggedTruthSomaticReadsMap;

    // Count all tagged somatic reads (TP+FP)
    int totalTaggedSomaticReads = 0;
    std::map<int, int> toatlTaggedSomaticReadsMap;

    // Calculate statistics across all chromosomes
    for(auto chr: chrVec){
        for(auto readIter :chrMetrics[chr].coverTruthSomaticPosReadVec){
            if(readIter.hpResult == ReadHP::H1_1 || readIter.hpResult == ReadHP::H2_1 || readIter.hpResult == ReadHP::H3){
                truthSomaticReadsMap[readIter.hpResult]++;
                totalTruthSomaticReads++;
            }
        }

        for(auto readIter :chrMetrics[chr].taggedSomaticReadVec){
            if(readIter.hpResult == ReadHP::H1_1 || readIter.hpResult == ReadHP::H2_1 || readIter.hpResult == ReadHP::H3){
                totalTaggedTruthSomaticReads++;
                taggedTruthSomaticReadsMap[readIter.hpResult]++;
            }
        }
        for(auto readIter :chrMetrics[chr].totalReadVec){
            if(readIter.hpResult == ReadHP::H1_1 || readIter.hpResult == ReadHP::H2_1 || readIter.hpResult == ReadHP::H3){
                toatlTaggedSomaticReadsMap[readIter.hpResult]++;
                totalTaggedSomaticReads++;
            }
            totalReads++;
        }
    }

    if(!somaticReadLog->is_open()){
        std::cerr<< "Fail to open write file: " << outputFileName << "\n";
        exit(1);
    }
    
    // Calculate overall performance metrics
    float recall = calculateRecall(totalTaggedTruthSomaticReads, totalTruthSomaticReads);
    float precision = calculatePrecision(totalTaggedTruthSomaticReads, totalTaggedSomaticReads);
    float f1_score = calculateF1Score(recall, precision);


    (*somaticReadLog) << "############################\n";
    (*somaticReadLog) << "# Somatic Haplotag Metrics #\n";
    (*somaticReadLog) << "############################\n";
    (*somaticReadLog) << "##Truth VCF File: "  << benchmarkVcf << "\n";
    (*somaticReadLog) << "##Truth BED File: "  << benchmarkBed << "\n";
    (*somaticReadLog) << "##MappingQualityThreshold: "  << mappingQualityThreshold << "\n";
    (*somaticReadLog) << "##Total reads: "  << totalReads << "\n";
    (*somaticReadLog) << "##Total truth somatic reads: "  << totalTruthSomaticReads << "\n";
    (*somaticReadLog) << "##Total truth HP1-1: "  << truthSomaticReadsMap[ReadHP::H1_1] << "\n";
    (*somaticReadLog) << "##Total truth HP2-1: "  << truthSomaticReadsMap[ReadHP::H2_1] << "\n";
    (*somaticReadLog) << "##Total truth HP3: "  << truthSomaticReadsMap[ReadHP::H3] << "\n";
    
    int separatorLength = 95;
    int columnWidth = 15;

    // Write performance metrics table header
    (*somaticReadLog) << std::left << std::setw(columnWidth) << "## Haplotype" 
                      << std::setw(columnWidth) << "Precision" 
                      << std::setw(columnWidth) << "Recall" 
                      << std::setw(columnWidth) << "F1-Score"
                      << std::setw(columnWidth) << "TP" 
                      << std::setw(columnWidth) << "FP" 
                      << std::setw(columnWidth) << "FN" << "\n";
    (*somaticReadLog) << "##" << std::string(separatorLength, '-') << "\n";
    
    // Calculate and write metrics for each haplotype
    std::vector<ReadHP> haplotypes = {ReadHP::H1_1, ReadHP::H2_1, ReadHP::H3};
    for(const auto& hp : haplotypes) {
        int tp = taggedTruthSomaticReadsMap[hp];
        int fp = toatlTaggedSomaticReadsMap[hp] - taggedTruthSomaticReadsMap[hp];
        int fn = truthSomaticReadsMap[hp] - taggedTruthSomaticReadsMap[hp];
        
        float precision = calculatePrecision(tp, tp + fp);
        float recall = calculateRecall(tp, tp + fn);
        float f1 = calculateF1Score(recall, precision);

        std::string hpStr = ReadHapUtil::readHapIntToString(hp);

        (*somaticReadLog) << std::left << std::setw(columnWidth) << ("## HP" + hpStr)
                        << std::setw(columnWidth) << std::fixed << std::setprecision(4) << precision
                        << std::setw(columnWidth) << std::fixed << std::setprecision(4) << recall
                        << std::setw(columnWidth) << std::fixed << std::setprecision(4) << f1 
                        << std::setw(columnWidth) << tp
                        << std::setw(columnWidth) << fp  
                        << std::setw(columnWidth) << fn << "\n";
    }

    // Write overall statistics
    (*somaticReadLog) << "##" << std::string(separatorLength, '-') << "\n";
    (*somaticReadLog) << std::left << std::setw(columnWidth) << "## Overall"
                    << std::setw(columnWidth) << std::fixed << std::setprecision(4) << precision
                    << std::setw(columnWidth) << std::fixed << std::setprecision(4) << recall
                    << std::setw(columnWidth) << std::fixed << std::setprecision(4) << f1_score 
                    << std::setw(columnWidth) << totalTaggedTruthSomaticReads
                    << std::setw(columnWidth) << totalTaggedSomaticReads - totalTaggedTruthSomaticReads
                    << std::setw(columnWidth) << totalTruthSomaticReads - totalTaggedTruthSomaticReads << "\n";
    
    (*somaticReadLog) << "##\n";

    // Write detailed read information
    (*somaticReadLog) << "#CHROM\t"
                      << "READID\t"
                      << "GERMLINE_VAR_SIMILARITY\t"
                      << "DERIVE_BY_HP_SIMILARITY\t"
                      << "GERMLINE_SNP_COUNT\t"
                      << "TUMOR_SNP_COUNT\t"
                      << "HAPLOTYPE\t"
                      << "TRUTH_VARIANT_POS,HP\n";

    for(auto chr: chrVec){
        for(auto somaticRead: *somaticReadVecMap[chr]){
            (*somaticReadLog) << somaticRead.chr << "\t"
                            << somaticRead.readID << "\t"
                            << somaticRead.germlineVarSimilarity << "\t"
                            << somaticRead.deriveByHpSimilarity << "\t"
                            << somaticRead.germlineSnpCount << "\t"
                            << somaticRead.tumorSnpCount << "\t"
                            << "H" << ReadHapUtil::readHapIntToString(somaticRead.hpResult) << "\t";

            for(auto SnpHp: somaticRead.somaticSnpHp){
                (*somaticReadLog) << SnpHp.first+1 << "," << SnpHp.second << "\t";
            }
                (*somaticReadLog) << "\n";
        }
    }
    (*somaticReadLog).close();
    delete somaticReadLog;
    somaticReadLog = nullptr;
}

/**
 * @brief Display somatic variant count
 * 
 * Prints the total number of somatic variants across all chromosomes.
 * 
 * @param chrVec Vector of chromosome names
 * @param chrMultiVariants Variant data container
 */
void SomaticReadBenchmark::displaySomaticVarCount(std::vector<std::string> &chrVec, std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants){
    if(!openTestingFunc) return;
    int totalVariantCount = 0;
    
    for(auto &chr : chrVec){
        std::map<int, MultiGenomeVar>::iterator chrVariantIter = chrMultiVariants[chr].begin();
        while(chrVariantIter != chrMultiVariants[chr].end()){
            if(chrVariantIter->second.isExists(Genome::TRUTH_SOMATIC)){
                totalVariantCount++;
            }
            chrVariantIter++;
        }
    }

    std::cout << "Total somatic variants: " << totalVariantCount << "\n";
}

/**
 * @brief Display BED region count
 * 
 * Prints statistics about BED regions and variant counts inside/outside
 * BED regions for each genome type (TUMOR, NORMAL, TRUTH_SOMATIC).
 * 
 * @param chrVec Vector of chromosome names
 */
void SomaticReadBenchmark::displayBedRegionCount(std::vector<std::string> &chrVec){
    if(!openTestingFunc || !loadedBedFile) return;
    int totalBedRegionCount = 0;
    for(auto &chr : chrVec){
        totalBedRegionCount += bedRegions[chr].size();
    }
    std::cout << "Total bed regions: " << totalBedRegionCount << "\n";
    std::cout << "--------------------------------" << "\n";
    std::cout << "Variant in bed region count: " << "\n";
    std::cout << "  -Tumor: " << variantInBedRegionCount[Genome::TUMOR] << "\n";
    std::cout << "  -Normal: " << variantInBedRegionCount[Genome::NORMAL] << "\n";
    std::cout << "  -Truth Somatic: " << variantInBedRegionCount[Genome::TRUTH_SOMATIC] << "\n";
    std::cout << "\n";
    std::cout << "Variant out bed region count: " << "\n";
    std::cout << "  -Tumor: " << variantOutBedRegionCount[Genome::TUMOR] << "\n";
    std::cout << "  -Normal: " << variantOutBedRegionCount[Genome::NORMAL] << "\n";
    std::cout << "  -Truth Somatic: " << variantOutBedRegionCount[Genome::TRUTH_SOMATIC] << "\n";
    std::cout << "--------------------------------" << "\n";
}