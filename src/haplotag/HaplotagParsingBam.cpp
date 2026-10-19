#include "HaplotagParsingBam.h"

// bgzf jobs (128 KiB each) queued for the output bam; htslib defaults to 2 x threads
static const int OUTPUT_BGZF_QUEUE_SIZE = 512;

/**
 * @brief Constructor for BamFileRAII class
 * @param BamFile Path to the input BAM file
 * @param fastaFile Path to the reference FASTA file
 * @param threadPool Thread pool for parallel processing
 * @param version Program version string
 * @param command Command line string
 * @param resultPrefix Output file prefix
 * @param outputFormat Output format (bam/cram)
 * @param writeOutputBam Whether to write output BAM file
 * 
 * Initializes BAM file resources and sets up input/output streams.
 * Handles BAM file opening, header reading, index loading, and thread pool setup.
 * Supports both input and output BAM file management with proper error handling.
 */
BamFileRAII::BamFileRAII(
    const std::string& BamFile,
    const std::string& fastaFile,
    htsThreadPool &threadPool,
    const std::string& version,
    const std::string& command,
    const std::string& resultPrefix,
    const std::string& outputFormat,
    const bool writeOutputBam
):
writeOutputBam(writeOutputBam), isReleased(false), in(nullptr), out(nullptr), bamHdr(nullptr), idx(nullptr), aln(nullptr)
{
    try{
        // open bam file
        in = hts_open(BamFile.c_str(), "r");
        checkNullPointer(in, "Cannot open bam file " + BamFile);

        // load reference file
        if (hts_set_fai_filename(in, fastaFile.c_str()) != 0) {
            throw std::runtime_error("Cannot set FASTA index file for " + fastaFile);
        }

        // input reader
        bamHdr = sam_hdr_read(in);
        checkNullPointer(bamHdr, "Cannot read header from bam file " + BamFile);

        // header add pg tag
        sam_hdr_add_pg(bamHdr, "longphase-s", "VN", version.c_str(), "CL", command.c_str(), NULL);

        // check bam file index
        idx = sam_index_load(in, BamFile.c_str());
        checkNullPointer(idx, "Cannot open index for bam file " + BamFile);

        // set thread
        if (hts_set_opt(in, HTS_OPT_THREAD_POOL, &threadPool) != 0) {
            throw std::runtime_error("Cannot set thread pool for input bam file " + BamFile);
        }

        if (writeOutputBam) {
            // output file mangement
            std::string writeBamFile = resultPrefix + "." + outputFormat;
            // Open output BAM file with appropriate mode (bam/cram)
            out = hts_open(writeBamFile.c_str(), (outputFormat == "bam" ? "wb" : "wc" ));
            // load reference file
            hts_set_fai_filename(out, fastaFile.c_str());
            // output writer
            int result = sam_hdr_write(out, bamHdr);
            if (result < 0) {
                throw std::runtime_error("Cannot write header to output bam file " + writeBamFile);
            }
            // set thread pool for output bam file
            // for bam, a deeper output queue lets compressed blocks pile up while the writer thread
            // waits on its periodic fdatasync, instead of blocking the main thread
            htsThreadPool outThreadPool = threadPool;
            if (outputFormat == "bam") {
                outThreadPool.qsize = OUTPUT_BGZF_QUEUE_SIZE;
            }
            if (hts_set_opt(out, HTS_OPT_THREAD_POOL, &outThreadPool) != 0) {
                throw std::runtime_error("Cannot set thread pool for output bam file " + BamFile);
            }
        }

        // initialize an alignment
        aln = bam_init1();
        checkNullPointer(aln, "Cannot initialize alignment for bam file " + BamFile);
    }catch(const std::exception& e){
        std::cerr << "[ERROR](BamFileRAII): " << e.what() << std::endl;
        exit(1);
    }
}

/**
 * @brief Template function to check for null pointers
 * @param ptr Pointer to check
 * @param errorMessage Error message to display if pointer is null
 * 
 * Throws runtime_error if the pointer is null
 */
template<typename T>
void BamFileRAII::checkNullPointer(const T* ptr, const std::string& errorMessage) const {
    if (ptr == nullptr) {
        throw std::runtime_error(errorMessage);
    }
}

/**
 * @brief Validates the state of BAM file resources
 * @return true if all resources are valid, exits program if not
 * 
 * Checks that all BAM file pointers are valid and not null
 */
bool BamFileRAII::validateState(){
    try{
        checkNullPointer(in, "(validateState) in ptr is nullptr");
        checkNullPointer(bamHdr, "(validateState) bamHdr ptr is nullptr");
        checkNullPointer(idx, "(validateState) idx ptr is nullptr");
        checkNullPointer(aln, "(validateState) aln ptr is nullptr");
        if(writeOutputBam){
            checkNullPointer(out, "(validateState) out ptr is nullptr");
        }
        return true;
    }catch(const std::runtime_error& e){
        std::cerr << "[ERROR](BamFileRAII): " << e.what() << std::endl;
        exit(1);
    }
}

/**
 * @brief Writes the current alignment to the output BAM file
 * 
 * Writes the alignment record to the output file and handles errors
 */
void BamFileRAII::samWriteBam(){
    try{
        int result = sam_write1(out, bamHdr, aln);
        if (result < 0) {
            throw std::runtime_error("write output bam file failed");
        }
    }catch(const std::exception& e){
        std::cerr << "[ERROR](BamFileRAII): " << e.what() << std::endl;
        exit(1);
    }
}

/**
 * @brief Destructor for BamFileRAII class
 * 
 * Ensures proper cleanup of BAM file resources.
 * Only destroys resources if they haven't been manually released.
 */
BamFileRAII::~BamFileRAII(){
    if (!isReleased) {
        destroy();
    }
}

/**
 * @brief Destroys all BAM file resources
 * 
 * Frees memory and closes file handles in the correct order.
 * Sets isReleased flag to prevent double cleanup.
 */
void BamFileRAII::destroy(){
    if (aln) bam_destroy1(aln);
    if (idx) hts_idx_destroy(idx);
    if (bamHdr) bam_hdr_destroy(bamHdr);
    if (in) sam_close(in);
    if (out) sam_close(out);
    isReleased = true;
}

/**
 * @brief Constructor for HaplotagBamParser class
 * @param config Configuration parameters for BAM parsing
 * @param control Control parameters for processing behavior
 * 
 * Initializes the BAM parser with configuration and control settings.
 */
HaplotagBamParser::HaplotagBamParser(
    const ParsingBamConfig &config,
    const ParsingBamControl &control
): config(config), control(control)
{

}

/**
 * @brief Destructor for HaplotagBamParser class
 * 
 * Virtual destructor for proper cleanup of derived classes.
 */
HaplotagBamParser::~HaplotagBamParser(){

}

/**
 * @brief Main function for parsing BAM files and performing haplotagging
 * @param ctx Context containing all necessary data structures and files
 * 
 * Orchestrates the entire BAM parsing and haplotagging process.
 * Validates input parameters, sets up thread pool, and delegates processing
 * to appropriate methods based on control mode.
 * 
 * Processing flow:
 * 1. Validate input parameters and context
 * 2. Determine last variant positions for efficient processing
 * 3. Initialize FASTA parser and thread pool
 * 4. Process BAM file based on mode (single/multi-thread)
 * 5. Clean up resources
 */
void HaplotagBamParser::parsingBam(BamParserContext& ctx){
    try{

        // Validate that multi-thread mode is not used with output BAM
        if(control.mode == ParsingBamMode::MULTI_THREAD && control.writeOutputBam) {
            throw std::runtime_error("Cannot set Multi Thread mode when writeOutputBam is true");
        }

        // Validate chromosome vector is not empty
        if(ctx.chrVec.size() == 0){
            throw std::runtime_error("chrVec is empty");
        }

        // Validate chromosome length map is not empty
        if(ctx.chrLength.size() == 0){
            throw std::runtime_error("chrLength is empty");
        }
        
        std::vector<int> last_pos;
        // get the last variant position of the reference
        getLastVarPos(last_pos, ctx.chrVec, ctx.chrMultiVariants, ctx.genomeSample);
        // reference fasta parser
        FastaParser fastaParser(ctx.fastaFile, ctx.chrVec, last_pos, config.numThreads);

        // Initialize thread pool for parallel processing
        htsThreadPool threadPool = {NULL, 0};
        // Create thread pool with specified number of threads
        if (!(threadPool.pool = hts_tpool_init(config.numThreads))) {
            throw std::runtime_error("Error creating thread pool");
        }

        // Process BAM file based on control mode
        switch(control.mode){
            case ParsingBamMode::SINGLE_THREAD:
                // Process the BAM file using a single thread
                // std::cerr << "[single thread]"; // [debug]
                processBamWithOutput(ctx, fastaParser, threadPool);
                break;

            case ParsingBamMode::MULTI_THREAD:
                // Process the BAM file using multiple threads
                // std::cerr << "[multi thread]"; // [debug]
                processBamParallel(ctx, fastaParser, threadPool);
                break;

            default:
                throw std::runtime_error("Unsupported parsing bam mode: " + std::to_string(control.mode));
        }

        // Clean up thread pool
        hts_tpool_destroy(threadPool.pool);
    }catch(const std::exception& e){
        std::cerr << "[ERROR](HaplotagBamParser): " << e.what() << std::endl;
        exit(1);
    }

    return;
}

/**
 * @brief Processes BAM file using parallel processing
 * @param ctx BAM parser context
 * @param fastaParser Reference sequence parser
 * @param threadPool Thread pool for parallel execution
 * 
 * Uses OpenMP to process chromosomes in parallel.
 * Each chromosome is processed independently with its own BAM file resources.
 * This mode does not support output BAM generation due to thread safety concerns.
 */
void HaplotagBamParser::processBamParallel(
    BamParserContext& ctx,
    const FastaParser &fastaParser,
    htsThreadPool &threadPool
){
    // loop all chromosome
    #pragma omp parallel for schedule(dynamic) num_threads(config.numThreads) 
    for(auto chr : ctx.chrVec ){
        // Allocate BAM file resources for each thread
        BamFileRAII bamRAII(ctx.BamFile, ctx.fastaFile, threadPool, config.version, config.command, config.resultPrefix, config.outputFormat, control.writeOutputBam);
        // CRAM only decodes the listed fields; QUAL and aux are left empty (no effect on BAM input)
        if (control.skipQualAndAux &&
            hts_set_opt(bamRAII.in, CRAM_OPT_REQUIRED_FIELDS, SAM_QNAME | SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR | SAM_SEQ) != 0) {
            std::cerr << "[ERROR](HaplotagBamParser): Cannot set required fields for " << ctx.BamFile << std::endl;
            exit(1);
        }

        // Create chromosome processing context
        ChrProcContext commonCtx(chr, ctx.chrLength.at(chr), config, ctx.genomeSample, ctx.vcfSet);
        // Create the chromosome processor using factory method
        auto chrProcessor = createProcessor(chr);
        // Process the chromosome
        chrProcessor->processSingleChrom(commonCtx, bamRAII, fastaParser, ctx.chrMultiVariants);
    }
}

/**
 * @brief Processes BAM file using single-threaded processing with output
 * @param ctx BAM parser context
 * @param fastaParser Reference sequence parser
 * @param threadPool Thread pool for parallel execution
 * 
 * Processes chromosomes sequentially and writes output BAM file.
 * This mode supports output BAM generation and provides progress information.
 * Each chromosome is processed one at a time to maintain output order.
 */
void HaplotagBamParser::processBamWithOutput(
    BamParserContext& ctx,
    const FastaParser &fastaParser,
    htsThreadPool &threadPool
){
    // Allocate BAM file resources for single-threaded processing
    BamFileRAII bamRAII(ctx.BamFile, ctx.fastaFile, threadPool, config.version, config.command, config.resultPrefix, config.outputFormat, control.writeOutputBam);
    // loop all chromosome
    for(auto chr : ctx.chrVec ){
        // Create the chromosome processor
        std::time_t begin = time(NULL);
        std::cerr<<"chr: " << chr << " ... " ;

        // Create chromosome processing context
        ChrProcContext commonCtx(chr, ctx.chrLength.at(chr), config, ctx.genomeSample, ctx.vcfSet);
        // Create the chromosome processor using factory method
        auto chrProcessor = createProcessor(chr);
        // Process the chromosome
        chrProcessor->processSingleChrom(commonCtx, bamRAII, fastaParser, ctx.chrMultiVariants);
        std::cerr<< difftime(time(NULL), begin) << "s\n";
    }
}

/**
 * @brief Gets the last variant position for each chromosome
 * @param last_pos Vector to store last variant positions
 * @param chrVec Vector of chromosome names
 * @param chrMultiVariants Map of variants by chromosome
 * @param genomeSample Genome type (NORMAL/TUMOR)
 * 
 * Determines the last variant position for efficient processing.
 * For NORMAL genome: finds last phased variant
 * For TUMOR genome: finds last variant (tumor or phased normal)
 * If no variants exist, sets position to 0.
 */
void HaplotagBamParser::getLastVarPos(
    std::vector<int>& last_pos, 
    const std::vector<std::string>& chrVec, 
    std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants,
    const Genome& genomeSample
){
    for( auto chr : chrVec ){
        bool existLastPos = false;

        // Search for last variant from end of chromosome
        for (auto lastVariantIter = chrMultiVariants[chr].rbegin(); lastVariantIter != chrMultiVariants[chr].rend(); ++lastVariantIter) {
            if (genomeSample == NORMAL) {
                // For normal genome, find last phased variant
                if ((*lastVariantIter).second.isExists(NORMAL) && (*lastVariantIter).second.Variant[NORMAL].isExistPhasedSet()) {
                    last_pos.push_back((*lastVariantIter).first);
                        existLastPos = true;
                        break;
                    }
            }
            else if (genomeSample == TUMOR) {
                // For tumor genome, find last variant (tumor or phased normal)
                if ((*lastVariantIter).second.isExists(TUMOR) || 
                    ((*lastVariantIter).second.isExists(NORMAL) && (*lastVariantIter).second.Variant[NORMAL].isExistPhasedSet())) {
                    last_pos.push_back((*lastVariantIter).first);
                    existLastPos = true;
                    break;
                }
            }else{
                std::cerr << "[ERROR] (getRefLastVarPos) => unsupported genome sample: " << genomeSample << std::endl;
                exit(EXIT_SUCCESS);
            }
        }
        
        // If no variants found, set position to 0
        if(!existLastPos){
            last_pos.push_back(0);
        }
    }
}

/**
 * @brief Constructor for ChromosomeProcessor class
 * @param writeOutputBam Whether to write output BAM file
 * @param mappingQualityFilter Whether to filter by mapping quality
 * 
 * Initializes chromosome processor with output and filtering settings.
 */
ChromosomeProcessor::ChromosomeProcessor(bool writeOutputBam, bool mappingQualityFilter)
: writeOutputBam(writeOutputBam), mappingQualityFilter(mappingQualityFilter)
{

}

/**
 * @brief Destructor for ChromosomeProcessor class
 * 
 * Virtual destructor for proper cleanup of derived classes.
 */
ChromosomeProcessor::~ChromosomeProcessor(){

}

/**
 * @brief Processes a single chromosome
 * @param ctx Chromosome processing context
 * @param bamRAII BAM file RAII wrapper
 * @param fastaParser Reference sequence parser
 * @param chrMultiVariants Map of variants by chromosome
 * 
 * Main function for processing all reads in a chromosome.
 * 
 * Processing flow:
 * 1. Extract chromosome variants and set up iterators
 * 2. Fetch reference sequence for the chromosome
 * 3. Set up region iterator for chromosome-specific reads
 * 4. Process each read based on its characteristics:
 *    - Low mapping quality reads
 *    - Unmapped reads
 *    - Secondary alignments
 *    - Supplementary alignments
 *    - Reads with no variants
 *    - Normal reads for processing
 * 5. Write reads to output BAM if enabled
 * 6. Perform post-processing operations
 */
void ChromosomeProcessor::processSingleChrom(
    ChrProcContext& ctx,
    BamFileRAII& bamRAII,
    const FastaParser& fastaParser,
    std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants  
){
    const std::string& chr = ctx.chrName;
    const int& chrLength = ctx.chrLength;
    const ParsingBamConfig& params = ctx.params; 

    // Extract all variants within this chromosome for processing
    std::map<int, MultiGenomeVar> currentVariants;

    #pragma omp critical
    {
        currentVariants = chrMultiVariants[chr];
    }
    // since each read is sorted based on the start coordinates, to save time, 
    // firstVariantIter keeps track of the first variant that each read needs to check.
    std::map<int, MultiGenomeVar>::iterator firstVariantIter = currentVariants.begin();
    // get the coordinates of the last variant
    // the tagging process will not be perform if the read's start coordinate are over than last variant.
    std::map<int, MultiGenomeVar>::reverse_iterator last = currentVariants.rbegin();

    // fetch chromosome string
    const std::string& ref_string = fastaParser.chrString.at(chr);

    // Set up region for chromosome-specific read iteration
    std::string region = !params.region.empty() ? params.region : std::string(chr + ":1-" + std::to_string(chrLength));
    hts_itr_t* iter = sam_itr_querys(bamRAII.idx, bamRAII.bamHdr, region.c_str());

    // Process each read in the chromosome
    while (sam_itr_multi_next(bamRAII.in, iter, bamRAII.aln) >= 0) {
        
        int flag = bamRAII.aln->core.flag;

        if ( bamRAII.aln->core.qual < params.qualityThreshold && mappingQualityFilter){
           // Mapping quality is lower than threshold
           processLowMappingQuality();
        }
        else if( (flag & 0x4) != 0 ){
            // Read is unmapped
            processUnmappedRead();
        }
        else if( (flag & 0x100) != 0 ){
            // Secondary alignment
            // A secondary alignment occurs when a given read could align reasonably well to more than one place
            processSecondaryAlignment();
        }
        else if( (flag & 0x800) != 0 && params.tagSupplementary == false ){
            // supplementary alignment
            // A chimeric alignment is represented as a set of linear alignments that do not have large overlaps.
            processSupplementaryAlignment();
        }
        // Check if no variants exist for this chromosome
        else if(last == currentVariants.rend()){ 
            //skip this read
            processEmptyVariants();
        }
        else if(int(bamRAII.aln->core.pos) <= (*last).first){
            //read start coordinate is less than or equal to the last variant coordinate
            processRead(*bamRAII.aln, *bamRAII.bamHdr, ref_string, currentVariants, firstVariantIter, ctx);
        }
        else{
            processOtherCase();
        }

        // Write all reads to output BAM file if enabled
        if(writeOutputBam){
            bamRAII.samWriteBam();
        }
    }

    // Perform post-processing operations for the chromosome
    postProcess(chr, currentVariants);
    
    // Clean up region iterator
    hts_itr_destroy(iter);
}

/**
 * @brief Constructor for CigarParser class
 * @param ctx CIGAR parsing context
 * @param ref_pos Reference position (will be modified)
 * @param query_pos Query position (will be modified)
 * 
 * Initializes CIGAR parser with context and position references.
 */
CigarParser::CigarParser(CigarParserContext& ctx, int& ref_pos, int& query_pos)
: ctx(ctx), ref_pos(ref_pos), query_pos(query_pos){}

/**
 * @brief Destructor for CigarParser class
 * 
 * Virtual destructor for proper cleanup of derived classes.
 */
CigarParser::~CigarParser(){

};

/**
 * @brief Parses CIGAR string and processes alignment operations
 * @param hpCount Map to count haplotype assignments
 * @param variantsHP Map to record variant haplotype assignments
 * @param norCountPS Map to count phase set assignments
 * 
 * Main function for parsing CIGAR operations and determining haplotype assignments.
 * 
 * CIGAR parsing flow:
 * 1. Skip variants that are to the left of the read start position
 * 2. Initialize position trackers (reference and query positions)
 * 3. Iterate through each CIGAR operation:
 *    - Match operations (M, =, X): process variants within match region
 *    - Insertion operations (I): skip reference position, advance query position
 *    - Deletion operations (D): process variants within deletion region
 *    - Skipped operations (N): advance reference position only
 *    - Soft clipping (S): advance query position only
 *    - Hard clipping (H): no position advancement
 *    - Padding (P): no position advancement
 */
void CigarParser::parsingCigar(
    std::map<int, int>& hpCount,
    std::map<int, int>& variantsHP,
    std::map<int, int>& norCountPS
){
    this->hpCount = &hpCount;
    this->variantsHP = &variantsHP;
    this->norCountPS = &norCountPS;

    // Skip variants that are to the left of this read
    while (ctx.firstVariantIter != ctx.currentVariants.end() && (*ctx.firstVariantIter).first < ctx.aln.core.pos) {
        ctx.firstVariantIter++;
    }

    if (ctx.firstVariantIter == ctx.currentVariants.end()) {
        return;
    }

    currentVariantIter = ctx.firstVariantIter;

    // position relative to reference
    ref_pos = ctx.aln.core.pos;
    // position relative to read
    query_pos = 0;

    // Parse CIGAR string to detect variants on this read
    int aln_core_n_cigar = int(ctx.aln.core.n_cigar);
    for (int i = 0; i < aln_core_n_cigar; i++) {
        uint32_t* cigar = bam_get_cigar(&ctx.aln);
        int cigar_op = bam_cigar_op(cigar[i]);
        int length = bam_cigar_oplen(cigar[i]);
        // iterator next variant
        while (currentVariantIter != ctx.currentVariants.end() && (*currentVariantIter).first < ref_pos) {
            currentVariantIter++;
        }

        // CIGAR operators: MIDNSHP=X correspond 012345678
        // 0: alignment match (can be a sequence match or mismatch)
        // 7: sequence match
        // 8: sequence mismatch
        if( cigar_op == 0 || cigar_op == 7 || cigar_op == 8 ){ 
                while( currentVariantIter != ctx.currentVariants.end() && (*currentVariantIter).first < ref_pos + length){
                    int offset = (*currentVariantIter).first - ref_pos;
                    if( offset < 0){
                    }
                    else{
                        uint8_t *q = bam_get_seq(&ctx.aln);
                        char base_chr = seq_nt16_str[bam_seqi(q,query_pos + offset)];
                        std::string base(1, base_chr);
                        HaplotagVariantType::VariantType variantType = HaplotagVariantType::NONE_VAR;
                        if((*currentVariantIter).second.isExists(NORMAL)){
                            variantType = (*currentVariantIter).second.Variant[NORMAL].variantType;
                        }else if((*currentVariantIter).second.isExists(TUMOR)){
                            variantType = (*currentVariantIter).second.Variant[TUMOR].variantType;
                        }
                        bool isAlt = false;

                        // [TODO] Check if the variant is both exist in normal and tumor
                        if((*currentVariantIter).second.isExists(NORMAL)){
                            isAlt = IsAltIndel(ref_pos, length, i, aln_core_n_cigar, cigar, currentVariantIter, base, variantType, NORMAL);
                        }else if((*currentVariantIter).second.isExists(TUMOR)){
                            isAlt = IsAltIndel(ref_pos, length, i, aln_core_n_cigar, cigar, currentVariantIter, base, variantType, TUMOR);
                        }    
                        processMatchOperation(length, cigar, i, aln_core_n_cigar, base, isAlt, offset);                        
                    }
                    currentVariantIter++;
                }
                query_pos += length;
                ref_pos += length;
        }
        // 1: insertion to the reference
        else if( cigar_op == 1 ){            
                processInsertionOperation(length);
                query_pos += length;
        }
        // 2: deletion from the reference
        else if( cigar_op == 2 ){
                // only execute at the first phased normal snp
                bool alreadyJudgeDel = false;
                while( currentVariantIter != ctx.currentVariants.end() && (*currentVariantIter).first < ref_pos + length){
                    // INS variant whose pos is consumed by a D op: record CIGAR flanking
                    // (read carries D where INS was expected; not counted as ALT)
                    if((*currentVariantIter).second.isExists(TUMOR) &&
                       (*currentVariantIter).second.Variant[TUMOR].variantType == HaplotagVariantType::INSERTION){
                        int variantPos = (*currentVariantIter).first;
                        const int FLANK_BP = 10;
                        int insLengthTumor = static_cast<int>((*currentVariantIter).second.Variant[TUMOR].allele.Alt.length()) - 1;
                        int currentInsLength = 0;
                        int lengthDiff = abs(insLengthTumor - currentInsLength);

                        std::string insRegionCigar = buildFlankingCigarString(
                            cigar, aln_core_n_cigar, ctx.aln.core.pos, variantPos, FLANK_BP, FLANK_BP);

                        if(ctx.genomeSample == TUMOR){
                            ctx.tumorInsInfo[variantPos].push_back(
                                TumorInsInfo(bam_get_qname(&ctx.aln), insRegionCigar,
                                             insLengthTumor, currentInsLength, lengthDiff));
                        } else if(ctx.genomeSample == NORMAL){
                            // D doesn't consume read; query_pos is at the read base just before the D
                            int normalInsPos = query_pos;
                            ctx.normalInsInfo[variantPos].push_back(
                                NormalInsInfo(bam_get_qname(&ctx.aln), insRegionCigar,
                                              currentInsLength, normalInsPos));
                        }
                    }
                    // Check deletion length difference when processing tumor BAM
                    if(ctx.genomeSample == TUMOR){
                        // Check if this variant is a DELETION
                        if((*currentVariantIter).second.isExists(TUMOR) && 
                           (*currentVariantIter).second.Variant[TUMOR].variantType == HaplotagVariantType::DELETION){
                            int delLengthTumor = static_cast<int>((*currentVariantIter).second.Variant[TUMOR].allele.Ref.length());
                            int currentDelLength = length; // CIGAR deletion length
                            int lengthDiff = abs(delLengthTumor - currentDelLength);
                            int shiftLength = (*currentVariantIter).first - ref_pos + 1;
                            // Build simplified CIGAR string for the deletion region (start from current D operation)
                            std::string delRegionCigar = buildSimplifiedCigarString(cigar, i, aln_core_n_cigar, delLengthTumor);
                            if(lengthDiff < 100){
                                ctx.tumorDelShiftInfo[(*currentVariantIter).first].push_back(TumorDelShiftInfo(bam_get_qname(&ctx.aln), delRegionCigar,shiftLength, delLengthTumor, currentDelLength, lengthDiff));
                                // std::cerr << "[DELETION_LENGTH_DIFF] Read"  <<  ctx.tumorDelShiftInfo[(*currentVariantIter).first].back().readid 
                                // << " Pos: " << (*currentVariantIter).first 
                                // << " Sample: " << TUMOR
                                // << " Sample_Del_Length: " << delLengthTumor 
                                // << " Shift_Length: " << ctx.tumorDelShiftInfo[(*currentVariantIter).first].back().shiftLength
                                // << " TUMOR_Del_Length: " << ctx.tumorDelShiftInfo[(*currentVariantIter).first].back().delLengthTumor
                                // << " CIGAR_Del_Length: " << ctx.tumorDelShiftInfo[(*currentVariantIter).first].back().CIGAR_Del_Length 
                                // << " CIGAR_string: " << ctx.tumorDelShiftInfo[(*currentVariantIter).first].back().cigarstring
                                // << " Diff: " << ctx.tumorDelShiftInfo[(*currentVariantIter).first].back().lengthDiff << std::endl;
                            }
                            
                        }
                    }
                    processDeletionOperation(length, cigar, i, aln_core_n_cigar, alreadyJudgeDel);
                    currentVariantIter++;
                }
                ref_pos += length;
        }
        // 3: skipped region from the reference
        else if( cigar_op == 3 ){
            processSkippedOperation(length);
            ref_pos += length;
        }
        // 4: soft clipping (clipped sequences present in SEQ)
        else if( cigar_op == 4 ){
            processSoftClippingOperation(length);
            query_pos += length;
        }
        // 5: hard clipping (clipped sequences NOT present in SEQ)
        // 6: padding (silent deletion from padded reference)
        else if( cigar_op == 5 || cigar_op == 6 ){
            // do nothing
        }
        else{
            std::cerr << "Alignment find unsupported CIGAR operation from read: " << bam_get_qname(&ctx.aln) << "\n";
            exit(1);
        }
    }
}

std::string CigarParser::buildSimplifiedCigarString(uint32_t* cigar, int startIdx, int aln_core_n_cigar, int refLengthLimit, int firstOpPartialLen){
    std::string delRegionCigar;
    static const char CIGAR_OPS[] = "MIDNSHP=X";
    int refConsumed = 0;
    char lastOp = '\0';
    int lastOpLen = 0;
    bool stopAfterThisOp = false;  // Flag to stop after current operation (for D exceeding region)
    
    for(int j = startIdx; j < aln_core_n_cigar && refConsumed < refLengthLimit && !stopAfterThisOp; j++){
        int opLen = bam_cigar_oplen(cigar[j]);
        int op = bam_cigar_op(cigar[j]);
        // Convert = (7) and X (8) to M
        char opChar = (op == 7 || op == 8) ? 'M' : CIGAR_OPS[op];
        
        // Check if this operation consumes reference (M, D, N, =, X)
        bool consumesRef = (op == 0 || op == 2 || op == 3 || op == 7 || op == 8);
        // Check if this is a D operation
        bool isDeletion = (op == 2);
        
        // For the first operation, use partial length if specified (but not for D)
        if(j == startIdx && firstOpPartialLen > 0 && consumesRef && !isDeletion){
            opLen = firstOpPartialLen;
        }
        
        // Calculate effective length within the deletion region
        int effectiveLen = opLen;
        int refToConsume = opLen;  // Amount to add to refConsumed
        if(consumesRef){
            int remaining = refLengthLimit - refConsumed;
            if(opLen > remaining && remaining > 0){
                if(isDeletion){
                    // D: keep full length for output, record D and stop
                    effectiveLen = opLen;
                    refToConsume = remaining;
                    stopAfterThisOp = true;  // Stop after recording this D
                } else {
                    // M: truncate both
                    effectiveLen = remaining;
                    refToConsume = remaining;
                }
            } else if(remaining <= 0){
                // No more room in region
                break;
            }
            refConsumed += refToConsume;
        }
        
        // Skip if effective length is 0
        if(effectiveLen <= 0) continue;
        
        // Merge consecutive same operations
        if(opChar == lastOp){
            lastOpLen += effectiveLen;
        } else {
            if(lastOp != '\0'){
                delRegionCigar += std::to_string(lastOpLen) + lastOp;
            }
            lastOp = opChar;
            lastOpLen = effectiveLen;
        }
    }
    // Append the last operation
    if(lastOp != '\0'){
        delRegionCigar += std::to_string(lastOpLen) + lastOp;
    }
    return delRegionCigar;
}

std::string CigarParser::buildFlankingCigarString(uint32_t* cigar, int aln_core_n_cigar, int alnStartRefPos, int variantPos, int leftFlank, int rightFlank){
    static const char CIGAR_OPS[] = "MIDNSHP=X";
    // Inclusive reference window [regionStart, regionEnd] anchored on variantPos
    // leftFlank covers up to and including variantPos; rightFlank covers strictly to the right
    int regionStart = variantPos - leftFlank + 1;
    int regionEnd   = variantPos + rightFlank;

    int currentRefPos = alnStartRefPos;
    std::string result;
    char lastOp = '\0';
    int lastOpLen = 0;

    auto flushLastOp = [&](){
        if(lastOp != '\0'){
            result += std::to_string(lastOpLen) + lastOp;
            lastOp = '\0';
            lastOpLen = 0;
        }
    };
    auto appendOp = [&](char opChar, int len){
        if(len <= 0) return;
        if(opChar == lastOp){
            lastOpLen += len;
        } else {
            flushLastOp();
            lastOp = opChar;
            lastOpLen = len;
        }
    };

    for(int j = 0; j < aln_core_n_cigar; j++){
        int opLen = bam_cigar_oplen(cigar[j]);
        int op    = bam_cigar_op(cigar[j]);
        char opChar = (op == 7 || op == 8) ? 'M' : CIGAR_OPS[op];
        bool consumesRef = (op == 0 || op == 2 || op == 3 || op == 7 || op == 8);

        if(consumesRef){
            int opStart = currentRefPos;             // inclusive
            int opEnd   = currentRefPos + opLen;     // exclusive
            int effStart = std::max(opStart, regionStart);
            int effEnd   = std::min(opEnd, regionEnd + 1);
            if(effStart < effEnd){
                appendOp(opChar, effEnd - effStart);
            }
            currentRefPos += opLen;
            if(currentRefPos > regionEnd) break;
        } else if(op == 1){
            // Insertion: keep it if its anchor sits inside the window.
            // currentRefPos is the next reference base after the I; an I anchored
            // at variantPos sits between variantPos and variantPos+1, i.e.
            // currentRefPos == variantPos + 1, which we want to include.
            if(currentRefPos >= regionStart && currentRefPos <= regionEnd + 1){
                appendOp('I', opLen);
            }
        }
        // Skip S(4)/H(5)/P(6)/N(3 already handled as ref-consuming) for output
    }
    flushLastOp();
    return result;
}

bool CigarParser::IsAltIndel(
    int& ref_pos, 
    int& length, 
    int& i, 
    int& aln_core_n_cigar, 
    uint32_t* cigar,
    std::map<int, MultiGenomeVar>::iterator& currentVariantIter, 
    std::string& base, 
    HaplotagVariantType::VariantType variantType,
    Genome sample
){
    bool isAlt = false;
    if(variantType == HaplotagVariantType::SNP){
        isAlt = base == (*currentVariantIter).second.Variant[sample].allele.Alt;
    }else if(variantType == HaplotagVariantType::INSERTION){
        int variantPos = (*currentVariantIter).first;
        bool tumorExists = (*currentVariantIter).second.isExists(TUMOR);

        // ALT detection: M ends exactly at variant pos and next op is I (no length/seq check).
        // Loose criterion is intentional: long-read aligners often emit slightly shifted/sized
        // I ops in indel-rich regions; counting them as ALT keeps norVAF accurate so that
        // germline INS get filtered out via the somatic norVAF threshold.
        if(i+1 < aln_core_n_cigar){
            isAlt = ref_pos + length - 1 == variantPos && bam_cigar_op(cigar[i+1]) == 1;
        }
        int observedInsLen = isAlt ? bam_cigar_oplen(cigar[i+1]) : 0;

        // Record CIGAR flanking window for every read whose M op covers this INS variant pos,
        // regardless of whether the read carries the expected I (ALT), a D, or just plain M.
        // ALT count is not affected (controlled by isAlt above).
        if(tumorExists){
            const int FLANK_BP = 10;
            int insLengthTumor = static_cast<int>((*currentVariantIter).second.Variant[TUMOR].allele.Alt.length()) - 1;
            int lengthDiff = abs(insLengthTumor - observedInsLen);

            std::string insRegionCigar = buildFlankingCigarString(
                cigar, aln_core_n_cigar, ctx.aln.core.pos, variantPos, FLANK_BP, FLANK_BP);

            if(ctx.genomeSample == TUMOR){
                ctx.tumorInsInfo[variantPos].push_back(
                    TumorInsInfo(bam_get_qname(&ctx.aln), insRegionCigar,
                                 insLengthTumor, observedInsLen, lengthDiff));
            } else if(ctx.genomeSample == NORMAL){
                int normalInsPos = query_pos + (variantPos - ref_pos) + 1;
                ctx.normalInsInfo[variantPos].push_back(
                    NormalInsInfo(bam_get_qname(&ctx.aln), insRegionCigar,
                                  observedInsLen, normalInsPos));
            }
        }
    }else if(variantType == HaplotagVariantType::DELETION && i+1 < aln_core_n_cigar){
        if(ref_pos + length - 1 == (*currentVariantIter).first && bam_cigar_op(cigar[i+1]) == 2){
            // Get deletion length from sample variant
            int delLengthSample = static_cast<int>((*currentVariantIter).second.Variant[sample].allele.Ref.length()) - 1;
            int currentDelLength = bam_cigar_oplen(cigar[i+1]);
            int shiftLength = abs(currentDelLength - delLengthSample);
            // Check if TUMOR variant exists
            bool tumorExists = (*currentVariantIter).second.isExists(TUMOR);
            int delLengthTumor = 0;
            if(tumorExists){
                delLengthTumor = static_cast<int>((*currentVariantIter).second.Variant[TUMOR].allele.Ref.length()) - 1;
            }
            
            int lengthDiff = shiftLength ;
            
            // Only store when processing tumor BAM and tumor variant exists and difference is small
            if(ctx.genomeSample == TUMOR && tumorExists && lengthDiff < 100){
                int firstMInRegion = ref_pos + length - (*currentVariantIter).first;
                std::string delRegionCigar = buildSimplifiedCigarString(cigar, i + 1 , aln_core_n_cigar, delLengthTumor , firstMInRegion);
                ctx.tumorDelShiftInfo[(*currentVariantIter).first].push_back(TumorDelShiftInfo(bam_get_qname(&ctx.aln), delRegionCigar, shiftLength, delLengthTumor, currentDelLength, lengthDiff));
                isAlt = true;
            }else if(ctx.genomeSample == NORMAL && tumorExists && lengthDiff < 100){
                int variantPos = (*currentVariantIter).first;
                int normalDelPos = ref_pos - 1;  // Position where D starts in read
                // Calculate how much of current M is in deletion region
                int firstMInRegion = ref_pos + length - variantPos;
                std::string delRegionCigar = buildSimplifiedCigarString(cigar, i + 1 , aln_core_n_cigar, delLengthTumor, firstMInRegion);
                ctx.normalDelInfo[variantPos].push_back(NormalDelInfo(bam_get_qname(&ctx.aln), delRegionCigar, currentDelLength, normalDelPos));
            }
            // Read has deletion at variant position - mark as ALT if length is similar
            // This ensures correct germline haplotype voting for normal reads
            // Right-shift case handles TUMOR-only restriction separately
            if(lengthDiff < 100){
                isAlt = true;
            }
        }else if(ref_pos + length - 1 > (*currentVariantIter).first && bam_cigar_op(cigar[i+1]) == 2 &&
        ref_pos + length - 1 <= (*currentVariantIter).first + static_cast<int>((*currentVariantIter).second.Variant[sample].allele.Ref.length()) - 1 ){//&& ctx.genomeSample == TUMOR
            // Get deletion length from sample variant
            // int delLengthSample = static_cast<int>((*currentVariantIter).second.Variant[sample].allele.Ref.length()) - 1;
            int currentDelLength = bam_cigar_oplen(cigar[i+1]);
            int shiftLength = ref_pos + length - 1 - (*currentVariantIter).first;
            // Check if TUMOR variant exists
            bool tumorExists = (*currentVariantIter).second.isExists(TUMOR);
            int delLengthTumor = 0;
            if(tumorExists){
                delLengthTumor = static_cast<int>((*currentVariantIter).second.Variant[TUMOR].allele.Ref.length()) - 1;
            }
            
            int lengthDiff = abs(delLengthTumor - currentDelLength);
            
            // Only store when processing tumor BAM and tumor variant exists and difference is small
            if(ctx.genomeSample == TUMOR && tumorExists && lengthDiff < 100){
                int firstMInRegion = ref_pos + length - (*currentVariantIter).first;
                std::string delRegionCigar = buildSimplifiedCigarString(cigar, i , aln_core_n_cigar, delLengthTumor , firstMInRegion);
                ctx.tumorDelShiftInfo[(*currentVariantIter).first].push_back(TumorDelShiftInfo(bam_get_qname(&ctx.aln), delRegionCigar, shiftLength, delLengthTumor, currentDelLength, lengthDiff));
                isAlt = true;
            }else if(ctx.genomeSample == NORMAL && tumorExists && lengthDiff < 100){
                int variantPos = (*currentVariantIter).first;
                int normalDelPos = ref_pos - 1;  // Position where D starts in read
                // Calculate how much of current M is in deletion region
                int firstMInRegion = ref_pos + length - variantPos;
                std::string delRegionCigar = buildSimplifiedCigarString(cigar, i, aln_core_n_cigar, delLengthTumor, firstMInRegion);
                ctx.normalDelInfo[variantPos].push_back(NormalDelInfo(bam_get_qname(&ctx.aln), delRegionCigar, currentDelLength, normalDelPos));
            }
        }
    }
    return isAlt;
}


/**
 * @brief Counts base nucleotides and applies mapping quality filtering
 * @param posBase Position base information structure
 * @param base Base nucleotide
 * @param aln BAM alignment record
 * @param mpqThreshold Mapping quality threshold
 * 
 * Updates base counts and filtered depth based on mapping quality
 */
void CigarParser::countBaseNucleotide(PosBase& posBase, std::string& base, const bam1_t& aln, const float& mpqThreshold, bool isAlt, HaplotagVariantType::VariantType variantType){
    // mapping quality is higher than threshold
    if ( aln.core.qual >= mpqThreshold ){
        if(base == "A"){
            posBase.MPQ_A_count++;
        }else if(base == "C"){
            posBase.MPQ_C_count++;
        }else if(base == "G"){
            posBase.MPQ_G_count++;
        }else if(base == "T"){
            posBase.MPQ_T_count++;
        }else{
            posBase.MPQ_unknow++;
        }
        if(isAlt){
            posBase.MPQ_altCount++;
        }
        posBase.filteredMpqDepth++;
    }
    if(base == "A"){
        posBase.A_count++;
    }else if(base == "C"){
        posBase.C_count++;
    }else if(base == "G"){
        posBase.G_count++;
    }else if(base == "T"){
        posBase.T_count++;
    }else{
        posBase.unknow++;
    }
    if(isAlt){
        if(variantType == HaplotagVariantType::DELETION){
            posBase.delCount++;
        }
        posBase.altCount++;
    }
    posBase.depth++;  
}

/**
 * @brief Counts deletion bases
 * @param posBase Position base information structure
 * 
 * Updates deletion count and total depth
 */
void CigarParser::countDeletionBase(PosBase& posBase){
    posBase.delCount++;
    posBase.depth++;
}