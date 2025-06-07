#include "HaplotagParsingBam.h"

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
            // open output bam file
            out = hts_open(writeBamFile.c_str(), (outputFormat == "bam" ? "wb" : "wc" ));
            // load reference file
            hts_set_fai_filename(out, fastaFile.c_str());
            // output writer
            int result = sam_hdr_write(out, bamHdr);
            if (result < 0) {
                throw std::runtime_error("Cannot write header to output bam file " + writeBamFile);
            }
            // set thread pool for output bam file
            if (hts_set_opt(out, HTS_OPT_THREAD_POOL, &threadPool) != 0) {
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

template<typename T>
void BamFileRAII::checkNullPointer(const T* ptr, const std::string& errorMessage) const {
    if (ptr == nullptr) {
        throw std::runtime_error(errorMessage);
    }
}

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

BamFileRAII::~BamFileRAII(){
    if (!isReleased) {
        destroy();
    }
}

void BamFileRAII::destroy(){
    if (aln) bam_destroy1(aln);
    if (idx) hts_idx_destroy(idx);
    if (bamHdr) bam_hdr_destroy(bamHdr);
    if (in) sam_close(in);
    if (out) sam_close(out);
    isReleased = true;
}


HaplotagBamParser::HaplotagBamParser(
    const ParsingBamConfig &config,
    const ParsingBamControl &control
): config(config), control(control)
{

}

HaplotagBamParser::~HaplotagBamParser(){

}

void HaplotagBamParser::parsingBam(BamParserContext& ctx){
    try{

        if(control.mode == ParsingBamMode::MULTI_THREAD && control.writeOutputBam) {
            throw std::runtime_error("Cannot set Multi Thread mode when writeOutputBam is true");
        }

        if(ctx.chrVec.size() == 0){
            throw std::runtime_error("chrVec is empty");
        }

        if(ctx.chrLength.size() == 0){
            throw std::runtime_error("chrLength is empty");
        }
        
        std::vector<int> last_pos;
        // get the last variant position of the reference
        getLastVarPos(last_pos, ctx.chrVec, ctx.mergedChrVarinat, ctx.genomeSample);
        // reference fasta parser
        FastaParser fastaParser(ctx.fastaFile, ctx.chrVec, last_pos, config.numThreads);

        // init data structure and get core n
        htsThreadPool threadPool = {NULL, 0};
        // creat thread pool
        if (!(threadPool.pool = hts_tpool_init(config.numThreads))) {
            throw std::runtime_error("Error creating thread pool");
        }

        switch(control.mode){
            case ParsingBamMode::SINGLE_THREAD:
                // Process the BAM file using a single thread
                std::cerr<< "[single thread]";
                processBamWithOutput(ctx, fastaParser, threadPool);
                break;

            case ParsingBamMode::MULTI_THREAD:
                // Process the BAM file using multiple threads
                std::cerr<< "[multi thread]";
                processBamParallel(ctx, fastaParser, threadPool);
                break;

            default:
                throw std::runtime_error("Unsupported parsing bam mode: " + std::to_string(control.mode));
        }

        hts_tpool_destroy(threadPool.pool);
    }catch(const std::exception& e){
        std::cerr << "[ERROR](HaplotagBamParser): " << e.what() << std::endl;
        exit(1);
    }

    return;
}

void HaplotagBamParser::processBamParallel(
    BamParserContext& ctx,
    const FastaParser &fastaParser,
    htsThreadPool &threadPool
){
    // loop all chromosome
    #pragma omp parallel for schedule(dynamic) num_threads(config.numThreads) 
    for(auto chr : ctx.chrVec ){
        // bam file resource allocation
        BamFileRAII bamRAII(ctx.BamFile, ctx.fastaFile, threadPool, config.version, config.command, config.resultPrefix, config.outputFormat, control.writeOutputBam);

        ChrProcContext commonCtx(chr, ctx.chrLength.at(chr), config, ctx.genomeSample, ctx.vcfSet);
        //create the chromosome processor
        auto chrProcessor = createProcessor(chr);
        chrProcessor->processSingleChromosome(commonCtx, bamRAII, fastaParser, ctx.mergedChrVarinat);
    }
}

void HaplotagBamParser::processBamWithOutput(
    BamParserContext& ctx,
    const FastaParser &fastaParser,
    htsThreadPool &threadPool
){
    // bam file resource allocation
    BamFileRAII bamRAII(ctx.BamFile, ctx.fastaFile, threadPool, config.version, config.command, config.resultPrefix, config.outputFormat, control.writeOutputBam);
    // loop all chromosome
    for(auto chr : ctx.chrVec ){
        //create the chromosome processor
        std::time_t begin = time(NULL);
        std::cerr<<"chr: " << chr << " ... " ;

        ChrProcContext commonCtx(chr, ctx.chrLength.at(chr), config, ctx.genomeSample, ctx.vcfSet);
        //create the chromosome processor
        auto chrProcessor = createProcessor(chr);
        chrProcessor->processSingleChromosome(commonCtx, bamRAII, fastaParser, ctx.mergedChrVarinat);
        std::cerr<< difftime(time(NULL), begin) << "s\n";
    }
}

void HaplotagBamParser::getLastVarPos(
    std::vector<int>& last_pos, 
    const std::vector<std::string>& chrVec, 
    std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat,
    const Genome& genomeSample
){
    for( auto chr : chrVec ){
        bool existLastPos = false;

        for (auto lastVariantIter = mergedChrVarinat[chr].rbegin(); lastVariantIter != mergedChrVarinat[chr].rend(); ++lastVariantIter) {
            if (genomeSample == NORMAL) {
                if ((*lastVariantIter).second.isExists(NORMAL) && (*lastVariantIter).second.Variant[NORMAL].isExistPhasedSet()) {
                    last_pos.push_back((*lastVariantIter).first);
                        existLastPos = true;
                        break;
                    }
            }
            else if (genomeSample == TUMOR) {
                if ((*lastVariantIter).second.isExists(TUMOR) || 
                    ((*lastVariantIter).second.isExists(NORMAL) && (*lastVariantIter).second.Variant[NORMAL].isExistPhasedSet())) {
                    last_pos.push_back((*lastVariantIter).first);
                    existLastPos = true;
                    break;
                }
            }else{
                std::cerr << "ERROR (getRefLastVarPos) => unsupported genome sample: " << genomeSample << std::endl;
                exit(EXIT_SUCCESS);
            }
        }
        
        if(!existLastPos){
            last_pos.push_back(0);
        }
    }
}

ChromosomeProcessor::ChromosomeProcessor(bool writeOutputBam, bool mappingQualityFilter)
: writeOutputBam(writeOutputBam), mappingQualityFilter(mappingQualityFilter)
{

}

ChromosomeProcessor::~ChromosomeProcessor(){

}

void ChromosomeProcessor::processSingleChromosome(
    ChrProcContext& ctx,
    BamFileRAII& bamRAII,
    const FastaParser& fastaParser,
    std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat  
){
    const std::string& chr = ctx.chrName;
    const int& chrLength = ctx.chrLength;
    const ParsingBamConfig& params = ctx.params; 

    // records all variants within this chromosome.
    std::map<int, MultiGenomeVar> currentVariants;

    #pragma omp critical
    {
        currentVariants = mergedChrVarinat[chr];
    }
    // since each read is sorted based on the start coordinates, to save time, 
    // firstVariantIter keeps track of the first variant that each read needs to check.
    std::map<int, MultiGenomeVar>::iterator firstVariantIter = currentVariants.begin();
    // get the coordinates of the last variant
    // the tagging process will not be perform if the read's start coordinate are over than last variant.
    std::map<int, MultiGenomeVar>::reverse_iterator last = currentVariants.rbegin();

    // fetch chromosome string
    std::string ref_string = fastaParser.chrString.at(chr);

    std::string region = !params.region.empty() ? params.region : std::string(chr + ":1-" + std::to_string(chrLength));
    hts_itr_t* iter = sam_itr_querys(bamRAII.idx, bamRAII.bamHdr, region.c_str());

    while (sam_itr_multi_next(bamRAII.in, iter, bamRAII.aln) >= 0) {
        
        int flag = bamRAII.aln->core.flag;

        if ( bamRAII.aln->core.qual < params.qualityThreshold && mappingQualityFilter){
           // mapping quality is lower than threshold
           processLowMappingQuality();
        }
        else if( (flag & 0x4) != 0 ){
            // read unmapped
            processUnmappedRead();
        }
        else if( (flag & 0x100) != 0 ){
            // secondary alignment. repeat.
            // A secondary alignment occurs when a given read could align reasonably well to more than one place.
            processSecondaryAlignment();
        }
        else if( (flag & 0x800) != 0 && params.tagSupplementary == false ){
            // supplementary alignment
            // A chimeric alignment is represented as a set of linear alignments that do not have large overlaps.
            processSupplementaryAlignment();
        }
        //currentVariants rbegin == rend
        else if(last == currentVariants.rend()){ 
            //skip this read
            processEmptyVariants();
        }
        else if(int(bamRAII.aln->core.pos) <= (*last).first){
            processRead(*bamRAII.aln, *bamRAII.bamHdr, ref_string, currentVariants, firstVariantIter, ctx);
        }
        else{
            processOtherCase();
        }

        //for write all reads in bam file
        if(writeOutputBam){
            bamRAII.samWriteBam();
        }
    }

    postProcess(chr, currentVariants);
    
    hts_itr_destroy(iter);
}

void ChromosomeProcessor::calculateBaseCommonInfo(PosBase& baseInfo, std::string& tumorAltBase){
    int &depth = baseInfo.depth;
    int AltCount = baseInfo.getBaseCount(tumorAltBase);

    int filteredMpqDepth = baseInfo.filteredMpqDepth;
    int filteredMpqAltCount = baseInfo.getMpqBaseCount(tumorAltBase);

    baseInfo.VAF = calculateVAF(AltCount, depth);
    baseInfo.filteredMpqVAF = calculateVAF(filteredMpqAltCount, filteredMpqDepth);
    baseInfo.nonDelVAF = calculateVAF(AltCount, (depth - baseInfo.delCount));

    baseInfo.lowMpqReadRatio = calculateLowMpqReadRatio(depth, filteredMpqDepth);
    baseInfo.delRatio = calculateDelRatio(baseInfo.delCount, depth);

    //read hp count in the normal bam
    int H1readCount = baseInfo.ReadHpCount[ReadHP::H1];
    int H2readCount = baseInfo.ReadHpCount[ReadHP::H2];
    int germlineReadHpCount = H1readCount + H2readCount;

    baseInfo.germlineHaplotypeImbalanceRatio = calculateHaplotypeImbalanceRatio(H1readCount, H2readCount, germlineReadHpCount);

    baseInfo.percentageOfGermlineHp = calculatePercentageOfGermlineHp(germlineReadHpCount, depth);
}

float ChromosomeProcessor::calculateVAF(int altCount, int depth){
    return (depth == 0 || altCount == 0) ? 0.0 : (float)altCount / (float)depth;
}

float ChromosomeProcessor::calculateLowMpqReadRatio(int depth, int filteredMpqDepth){
    return depth == 0 ? 0.0 : (float)(depth - filteredMpqDepth) / (float)depth;
}

float ChromosomeProcessor::calculateDelRatio(int delCount, int depth){
    return (depth == 0 || delCount == 0) ? 0.0 : (float)delCount / (float)depth;
}

double ChromosomeProcessor::calculateHaplotypeImbalanceRatio(int& H1readCount, int& H2readCount, int& totalReadCount) {
    if(H1readCount > 0 && H2readCount > 0) {
        return (H1readCount > H2readCount) ? 
               ((double)H1readCount / (double)totalReadCount) : 
               ((double)H2readCount / (double)totalReadCount);
    } else if(H1readCount == 0 && H2readCount == 0) {
        return 0.0;
    } else {
        return 1.0;
    }
}

double ChromosomeProcessor::calculatePercentageOfGermlineHp(int& totalGermlineReadCount, int& depth) {
    return (depth == 0 || totalGermlineReadCount == 0) ? 0.0 : (double)totalGermlineReadCount / (double)depth;
}

CigarParser::CigarParser(CigarParserContext& ctx, int& ref_pos, int& query_pos)
: ctx(ctx), ref_pos(ref_pos), query_pos(query_pos){}

CigarParser::~CigarParser(){

};


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

    // reading cigar to detect snp on this read
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
                        processMatchOperation(length, cigar, i, aln_core_n_cigar, base);
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

void CigarParser::countBaseNucleotide(PosBase& posBase, std::string& base, const bam1_t& aln, const float& mpqThreshold){
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
    posBase.depth++;  
}

void CigarParser::countDeletionBase(PosBase& posBase){
    posBase.delCount++;
    posBase.depth++;
}



void GermlineHaplotagStrategy::judgeSnpHap(
    const std::string& chrName,
    VarData& norVar,
    const std::string& base,
    int& ref_pos,
    int& length,
    int& i,
    int& aln_core_n_cigar,
    uint32_t* cigar,
    std::map<int, MultiGenomeVar>::iterator &currentVariantIter,
    std::map<int, int>& hpCount,
    std::map<int, int>& variantsHP,
    std::map<int, int>& countPS
){
    int curPos = (*currentVariantIter).first;

    // currentVariant is SNP
    if( norVar.variantType == VariantType::SNP ){
        // Detected that the base of the read is either REF or ALT. 
        if( (base == norVar.allele.Ref) || (base == norVar.allele.Alt) ){


            if(!norVar.isExistPhasedSet()){
                std::cerr << "ERROR (judgeSnpHap) => can't find the position:" 
                          << " chr: " << chrName << "\t"
                          << " pos: " << curPos << "\t"
                          << " ref: " << norVar.allele.Ref << "\t"
                          << " alt: " << norVar.allele.Alt << "\n";
                exit(EXIT_SUCCESS);
            }
            else{
                if( base == norVar.HP1){
                    hpCount[SnpHP::GERMLINE_H1]++;
                    variantsHP[curPos]=0;
                }
                if( base == norVar.HP2){
                    hpCount[SnpHP::GERMLINE_H2]++;
                    variantsHP[curPos]=1;
                }
                countPS[norVar.PhasedSet]++;
            }
            
        }
    }
    // currentVariant is insertion
    else if( norVar.variantType == VariantType::INSERTION && i+1 < aln_core_n_cigar){
        
        int hp1Length = norVar.HP1.length();
        int hp2Length = norVar.HP2.length();
        
        if ( ref_pos + length - 1 == (*currentVariantIter).first && bam_cigar_op(cigar[i+1]) == 1 ) {
            // hp1 occur insertion
            if( hp1Length != 1 && hp2Length == 1 ){
                hpCount[SnpHP::GERMLINE_H1]++;
                variantsHP[curPos]=0;
            }
            // hp2 occur insertion
            else if( hp1Length == 1 && hp2Length != 1 ){
                hpCount[SnpHP::GERMLINE_H2]++;
                variantsHP[curPos]=1;
            }
        }
        else {
            // hp1 occur insertion
            if( hp1Length != 1 && hp2Length == 1 ){
                hpCount[SnpHP::GERMLINE_H2]++;
                variantsHP[curPos]=1;
            }
            // hp2 occur insertion
            else if( hp1Length == 1 && hp2Length != 1 ){
                hpCount[SnpHP::GERMLINE_H1]++;
                variantsHP[curPos]=0;
            }
        }
        countPS[norVar.PhasedSet]++;
    } 
    // currentVariant is deletion
    else if( norVar.variantType == VariantType::DELETION && i+1 < aln_core_n_cigar) {

        int hp1Length = norVar.HP1.length();
        int hp2Length = norVar.HP2.length();
        
        if ( ref_pos + length - 1 == (*currentVariantIter).first && bam_cigar_op(cigar[i+1]) == 2 ) {
            // hp1 occur deletion
            if( hp1Length != 1 && hp2Length == 1 ){
                hpCount[SnpHP::GERMLINE_H1]++;
                variantsHP[curPos]=0;
            }
            // hp2 occur deletion
            else if( hp1Length == 1 && hp2Length != 1 ){
                hpCount[SnpHP::GERMLINE_H2]++;
                variantsHP[curPos]=1;
            }
        }
        else {
            // hp2 occur deletion
            if( hp1Length != 1 && hp2Length == 1 ){
                hpCount[SnpHP::GERMLINE_H2]++;
                variantsHP[curPos]=1;
            }
            // hp1 occur deletion
            else if( hp1Length == 1 && hp2Length != 1 ){
                hpCount[SnpHP::GERMLINE_H1]++;
                variantsHP[curPos]=0;
            }
        }
        countPS[norVar.PhasedSet]++;
    } 
}


void GermlineHaplotagStrategy::judgeDeletionHap(
    const std::string& chrName,
    const std::string& ref_string,
    int& ref_pos,
    int& length,
    int& query_pos,
    std::map<int, MultiGenomeVar>::iterator &currentVariantIter,
    const bam1_t* aln,
    std::map<int, int>& hpCount,
    std::map<int, int>& variantsHP,
    std::map<int, int>& countPS
) {
    if (ref_string != "") {
        int del_len = length;
        if (ref_pos + del_len + 1 == (*currentVariantIter).first) {
            //if (homopolymerLength((*currentVariantIter).first, ref_string) >= 3) {
                // special case
            //}
        } else if ((*currentVariantIter).first >= ref_pos && (*currentVariantIter).first < ref_pos + del_len) {
            // check variant in homopolymer
            if (homopolymerLength((*currentVariantIter).first, ref_string) >= 3) {
                
                int curPos = (*currentVariantIter).first;
                auto norVar = (*currentVariantIter).second.Variant[NORMAL];
                
                // SNP
                if (norVar.variantType == VariantType::SNP) {
                    // get the next match
                    char base_chr = seq_nt16_str[bam_seqi(bam_get_seq(aln), query_pos)];
                    std::string base(1, base_chr);

                    if (base == norVar.HP1) {
                        hpCount[SnpHP::GERMLINE_H1]++;
                        variantsHP[curPos] = 0;
                    }
                    if (base == norVar.HP2) {
                        hpCount[SnpHP::GERMLINE_H2]++;
                        variantsHP[curPos] = 1;
                    }
                    countPS[norVar.PhasedSet]++;
                }
                
                // the read deletion contain VCF's deletion
                else if (norVar.variantType == VariantType::DELETION) {

                    int hp1Length = norVar.HP1.length();
                    int hp2Length = norVar.HP2.length();
                    // hp1 occur deletion
                    if (hp1Length != 1 && hp2Length == 1) {
                        hpCount[SnpHP::GERMLINE_H1]++;
                        variantsHP[curPos] = 0;
                    }
                    // hp2 occur deletion
                    else if (hp1Length == 1 && hp2Length != 1) {
                        hpCount[SnpHP::GERMLINE_H2]++;
                        variantsHP[curPos] = 1;
                    }
                    countPS[norVar.PhasedSet]++;
                }
            }
        }
    }
}

void GermlineHaplotagStrategy::judgeSVHap(const bam1_t &aln, std::map<Genome, VCF_Info> &vcfSet, std::map<int, int>& hpCount, const int& genomeSample){
    auto readIter = vcfSet[Genome::NORMAL].readSVHapCount.find(bam_get_qname(&aln));
    if( readIter != vcfSet[Genome::NORMAL].readSVHapCount.end()){
        hpCount[SnpHP::GERMLINE_H1] += vcfSet[Genome::NORMAL].readSVHapCount[bam_get_qname(&aln)][0];
        hpCount[SnpHP::GERMLINE_H2] += vcfSet[Genome::NORMAL].readSVHapCount[bam_get_qname(&aln)][1];
    }
}

int GermlineHaplotagStrategy::judgeReadHap(
    std::map<int, int>& hpCount, 
    double& min, 
    double& max, 
    double& percentageThreshold, 
    int& pqValue, 
    int& psValue, 
    std::map<int, int>& countPS, 
    int* totalHighSimilarity, 
    int* totalWithOutVaraint
){
    int hpResult = ReadHP::unTag;

    if(hpCount[SnpHP::GERMLINE_H1] > hpCount[SnpHP::GERMLINE_H2]){
        min = hpCount[SnpHP::GERMLINE_H2];
        max = hpCount[SnpHP::GERMLINE_H1];
    }
    else{
        min = hpCount[SnpHP::GERMLINE_H1];
        max = hpCount[SnpHP::GERMLINE_H2];
    }

    if( max/(max+min) < percentageThreshold){
        // no tag
        pqValue = 0;
        if(totalHighSimilarity != nullptr) (*totalHighSimilarity)++;
    }
    else{
        if(hpCount[SnpHP::GERMLINE_H1] > hpCount[SnpHP::GERMLINE_H2]){
            hpResult = ReadHP::H1;
        }
        if(hpCount[SnpHP::GERMLINE_H1] < hpCount[SnpHP::GERMLINE_H2]){
            hpResult = ReadHP::H2;
        }
    }

    if( max == 0 ){
        pqValue=0;
        if(totalWithOutVaraint != nullptr) (*totalWithOutVaraint)++;
    }
    else if( max == ( max + min ) ){
        pqValue=40;
    }
    else{
        pqValue=-10*(std::log10((double)min/double(max+min)));
    }
    
    // cross two block
    if( countPS.size() > 1  ){
        hpResult = ReadHP::unTag;
    }
    //set psValue
    if(hpResult != ReadHP::unTag){
        auto psIter = countPS.begin();
        psValue = (*psIter).first;
    }
    return hpResult;
}

void SomaticJudgeHapStrategy::judgeSomaticSnpHap(std::map<int, MultiGenomeVar>::iterator &currentVariantIter, std::string chrName, std::string base, std::map<int, int> &hpCount, std::map<int, int> &norCountPS, std::map<int, int> &tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *tumorAllelePosVec){
    int curPos = (*currentVariantIter).first;
    auto& curVar = (*currentVariantIter).second;

    // normal & tumor SNP at the current position (base on normal phased SNPs)
    // both normal and tumor samples that do not exist in the high-confidence set
    if(curVar.isExists(NORMAL) && curVar.isExists(TUMOR)){
        // the tumor & normal SNP GT are phased heterozygous 
        if((curVar.Variant[NORMAL].GT == GenomeType::PHASED_HETERO) && 
            (curVar.Variant[TUMOR].GT == GenomeType::PHASED_HETERO)){ 
            judgeNormalSnpHap(chrName, curPos, curVar, base, hpCount, norCountPS, variantsHP);

        }
        //the normal SNP GT is phased heterozgous & the tumor SNP GT is unphased heterozgous 
        else if((curVar.Variant[NORMAL].GT == GenomeType::PHASED_HETERO) && 
                (curVar.Variant[TUMOR].GT == GenomeType::UNPHASED_HETERO)){   
            judgeNormalSnpHap(chrName, curPos, curVar, base, hpCount, norCountPS, variantsHP);

        }
        //the normal SNP GT is phased heterozgous & the tumor SNP GT is homozygous 
        else if((curVar.Variant[NORMAL].GT == GenomeType::PHASED_HETERO) && 
                (curVar.Variant[TUMOR].GT == GenomeType::UNPHASED_HOMO)){ 
            judgeNormalSnpHap(chrName, curPos, curVar, base, hpCount, norCountPS, variantsHP);
        }
    // only normal SNP at the current position
    }else if(curVar.isExists(NORMAL)){
        // the normal SNP GT is phased heterozgous SNP
        if((curVar.Variant[NORMAL].GT == GenomeType::PHASED_HETERO)){
            judgeNormalSnpHap(chrName, curPos, curVar, base, hpCount, norCountPS, variantsHP);
        }
    // only tumor SNP at the current position
    }else if(curVar.isExists(TUMOR)){
        //the tumor SNP GT is phased heterozygous
        if(curVar.Variant[TUMOR].GT == GenomeType::PHASED_HETERO){
            if(curVar.Variant[TUMOR].allele.Ref == base || curVar.Variant[TUMOR].allele.Alt == base){
                if(!curVar.Variant[TUMOR].isExistPhasedSet()){
                    std::cerr<< curPos << "\t"
                             << curVar.Variant[TUMOR].allele.Ref << "\t"
                             << curVar.Variant[TUMOR].allele.Alt << "\n";
                    exit(EXIT_SUCCESS);
                }else{
                    judgeTumorOnlySnpHap(chrName, curPos, curVar, base, hpCount, &tumCountPS, variantsHP, tumorAllelePosVec);
                }
            }
        //the tumor SNP GT is unphased heterozygous
        }else if(curVar.Variant[TUMOR].GT == GenomeType::UNPHASED_HETERO){
            if(curVar.Variant[TUMOR].allele.Ref == base || curVar.Variant[TUMOR].allele.Alt == base){
                judgeTumorOnlySnpHap(chrName, curPos, curVar, base, hpCount, nullptr, variantsHP, tumorAllelePosVec);
            }           
        //the tumor SNP GT is homozygous
        }else if(curVar.Variant[TUMOR].GT == GenomeType::UNPHASED_HOMO){
            if(curVar.Variant[TUMOR].allele.Ref == base || curVar.Variant[TUMOR].allele.Alt == base){
                judgeTumorOnlySnpHap(chrName, curPos, curVar, base, hpCount, nullptr, variantsHP, tumorAllelePosVec);
            }
        }
    }
}

void SomaticJudgeHapStrategy::judgeNormalSnpHap(
    const std::string& chrName, 
    int& curPos,
    MultiGenomeVar& curVar,
    std::string& base,
    std::map<int, int>& hpCount, 
    std::map<int, int>& norCountPS,
    std::map<int, int> *variantsHP
){
    if(curVar.Variant[NORMAL].allele.Ref == base || curVar.Variant[NORMAL].allele.Alt == base){

        if(!curVar.Variant[NORMAL].isExistPhasedSet()){
            std::cerr<< "Unable to locate the phase set of the current normal SNP\n"
                        << curPos << "\t"
                        << curVar.Variant[NORMAL].allele.Ref << "\t"
                        << curVar.Variant[NORMAL].allele.Alt  << "\n";
            exit(EXIT_SUCCESS);
        }

        std::string& norHP1 = curVar.Variant[NORMAL].HP1;
        std::string& norHP2 = curVar.Variant[NORMAL].HP2;

        if( base == norHP1){
            hpCount[1]++;
            if(variantsHP != nullptr) (*variantsHP)[curPos] = SnpHP::GERMLINE_H1;
        }
        if(base == norHP2){
            hpCount[2]++;
            if(variantsHP != nullptr) (*variantsHP)[curPos] = SnpHP::GERMLINE_H2;
        }
        norCountPS[curVar.Variant[NORMAL].PhasedSet]++;
    }
}

int SomaticJudgeHapStrategy::judgeSomaticReadHap(
    std::map<int, int> &hpCount,
    int &pqValue,
    std::map<int, int> &norCountPS,
    double &norHPsimilarity,
    double &tumHPsimilarity,
    double percentageThreshold,
    int *totalHighSimilarity,
    int *totalCrossTwoBlock,
    int *totalWithOutVaraint
){
    double normalMinHPcount = 0;
    double normalMaxHPcount = 0;
    int maxNormalHP = 0;

    double tumorMinHPcount = 0;
    double tumorMaxHPcount = 0;
    int maxTumorHP = 0;

    // determine max and min
    //tumor HP count
    if(hpCount[3] > hpCount[4]){
        tumorMinHPcount = hpCount[4];
        tumorMaxHPcount = hpCount[3];
        maxTumorHP = SnpHP::SOMATIC_H3;
    }
    else{
        tumorMinHPcount = hpCount[3];
        tumorMaxHPcount = hpCount[4];
        maxTumorHP = SnpHP::SOMATIC_H4;
    }

    //normal HP count
    if(hpCount[1] > hpCount[2]){
        normalMinHPcount = hpCount[2];
        normalMaxHPcount = hpCount[1];
        maxNormalHP = SnpHP::GERMLINE_H1;
    }
    else{
        normalMinHPcount = hpCount[1];
        normalMaxHPcount = hpCount[2];
        maxNormalHP = SnpHP::GERMLINE_H2;
    }

    //the similarity of HP types
    //tumor variants
    tumHPsimilarity = (tumorMaxHPcount == 0) ? 0.0 : tumorMaxHPcount/(tumorMaxHPcount+tumorMinHPcount);
    norHPsimilarity = (normalMaxHPcount == 0) ? 0.0 : normalMaxHPcount/(normalMaxHPcount+normalMinHPcount);


    // determine the haplotype of the read
    int hpResult = ReadHP::unTag;

    // read have tumor variant
    if(tumorMaxHPcount != 0){
        //check the similarity of HP types in the tumor
        if(tumHPsimilarity >= percentageThreshold){
            //check the similarity of HP types in the normal
            if(norHPsimilarity >= percentageThreshold){
                switch(maxTumorHP){
                    case SnpHP::SOMATIC_H3:
                        if(maxNormalHP == SnpHP::GERMLINE_H1){
                            hpResult = ReadHP::H1_1;
                        }else if(maxNormalHP == SnpHP::GERMLINE_H2){
                            hpResult = ReadHP::H2_1;
                        }
                        break;
                    case SnpHP::SOMATIC_H4:
                        if(maxNormalHP == SnpHP::GERMLINE_H1){
                            hpResult = ReadHP::H1_2;
                        }else if(maxNormalHP == SnpHP::GERMLINE_H2){
                            hpResult = ReadHP::H2_2;
                        }
                        break;
                    default:
                        std::cerr << "ERROR: Unexpected haplotype : tumor Max HP= "<< maxTumorHP << std::endl;
                        exit(1);
                        break;
                }
            }
            // can't determine HP from germline
            else{
                //The read has a somatic SNP, but the germline haplotype(H1/H2) it came from is unknown
                switch(maxTumorHP){
                    case SnpHP::SOMATIC_H3:
                        hpResult = ReadHP::H3; break;
                    case SnpHP::SOMATIC_H4:
                        hpResult = ReadHP::H4; break;
                    default:
                        std::cerr << "ERROR: Unexpected haplotype : tumor Max HP= "<< maxTumorHP << std::endl;
                        exit(1);
                        break;
                }
            }
        }else{
            // no tag
            pqValue = 0;
            if(totalHighSimilarity != nullptr) (*totalHighSimilarity)++;
        }
    }
    // read haven't tumor variant
    else if(normalMaxHPcount != 0){
        if(norHPsimilarity >= percentageThreshold){
            hpResult = maxNormalHP;
        }else{
            // no tag
            pqValue = 0;
            if(totalHighSimilarity != nullptr) (*totalHighSimilarity)++;
        }
    }

    // cross two block
    // had at least one tumor-unique SNP in the current read
    if(hpCount[3] != 0){
        if(norCountPS.size() > 1){
            hpResult = ReadHP::unTag;
            if(totalCrossTwoBlock != nullptr) (*totalCrossTwoBlock)++;
        }
    // all SNPs are germline variants in the current read
    }else{
        if(norCountPS.size() > 1){
            hpResult = ReadHP::unTag;
            if(totalCrossTwoBlock != nullptr) (*totalCrossTwoBlock)++;
        }
    }

    // determine the quality of the read
    // There haven't been any variants in the current read
    if(normalMaxHPcount == 0 && tumorMaxHPcount == 0){
        //std::cerr<< "Read max = 0 : " << bam_get_qname(&aln) << "\n";
        if(totalWithOutVaraint != nullptr) (*totalWithOutVaraint)++;
        pqValue=0;
    }
    // read have tumor variant
    else if(tumorMaxHPcount != 0){
        if( tumorMaxHPcount == ( tumorMaxHPcount + tumorMinHPcount ) ){
            pqValue=40;
        }
        else{
            pqValue=-10*(std::log10((double)tumorMinHPcount/double(tumorMaxHPcount+tumorMinHPcount)));
        }
    }
    else if(normalMaxHPcount != 0){
        if( normalMaxHPcount == ( normalMaxHPcount + normalMinHPcount ) ){
            pqValue=40;
        }
        else{
            pqValue=-10*(std::log10((double)normalMinHPcount/double(normalMaxHPcount+normalMinHPcount)));
        }
    }

    return hpResult;
}


void ExtractSomaticDataStragtegy::judgeTumorOnlySnpHap(const std::string &chrName, int &curPos, MultiGenomeVar &curVar, std::string base, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *tumorAllelePosVec){
 //the tumor SNP GT is phased heterozygous
    //all bases of the same type at the current position in normal.bam
    if(tumorAllelePosVec == nullptr){
        std::cerr << "ERROR (SomaticDetectJudgeHP) => tumorAllelePosVec pointer cannot be nullptr"<< std::endl;
        exit(1);
    }

    std::string& TumorRefBase = curVar.Variant[TUMOR].allele.Ref;
    std::string& TumorAltBase = curVar.Variant[TUMOR].allele.Alt; 
    
    if(base == TumorAltBase){
        hpCount[3]++;
        if(variantsHP != nullptr) (*variantsHP)[curPos] = SnpHP::SOMATIC_H3;

        //record postions that tagged as HP3 for calculating the confidence of somatic positions
        (*tumorAllelePosVec).push_back(curPos);

    //base is not match to TumorRefBase & TumorAltBase (other HP)
    }else if(base != TumorRefBase && base != TumorAltBase){
        //hpCount[4]++;
        //if(variantsHP != nullptr) (*variantsHP)[curPos] = SnpHP::SOMATIC_H4;
        //(*tumorAllelePosVec).push_back(curPos);
    }

    if(tumCountPS != nullptr) (*tumCountPS)[curVar.Variant[TUMOR].PhasedSet]++;
}

void SomaticHaplotagStrategy::judgeTumorOnlySnpHap(const std::string &chrName, int &curPos, MultiGenomeVar &curVar, std::string base, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *tumorAllelePosVec){

    if(curVar.isSomaticVariant){
        std::string& TumorRefBase = curVar.Variant[TUMOR].allele.Ref;
        std::string& TumorAltBase = curVar.Variant[TUMOR].allele.Alt; 

        if(base == TumorAltBase){
            hpCount[3]++;
            if(variantsHP != nullptr) (*variantsHP)[curPos] = SnpHP::SOMATIC_H3;

        //base is not match to TumorRefBase & TumorAltBase (other HP)
        }else if(base != TumorRefBase && base != TumorAltBase){
            //hpCount[4]++;
            //if(variantsHP != nullptr) (*variantsHP)[curPos] = SnpHP::SOMATIC_H4;
            //std::cerr << "Somatic SNP: " << curPos << " " << base << " -> " << TumorRefBase << "|" << TumorAltBase << std::endl;
        }

        if(curVar.Variant[TUMOR].GT == GenomeType::PHASED_HETERO){
            if(tumCountPS != nullptr) (*tumCountPS)[curVar[TUMOR].PhasedSet]++;
        }
    }
}