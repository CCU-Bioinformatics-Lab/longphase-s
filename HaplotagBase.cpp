#include "HaplotagBase.h"

BamFileRAII::BamFileRAII(
    const std::string& BamFile
  , const std::string& fastaFile
  , htsThreadPool& threadPool
  , const HaplotagParameters& params
  , const bool writeOutputBam
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
        sam_hdr_add_pg(bamHdr, "longphase", "VN", params.version.c_str(), "CL", params.command.c_str(), NULL);

        // check bam file index
        idx = sam_index_load(in, BamFile.c_str());
        checkNullPointer(idx, "Cannot open index for bam file " + BamFile);

        // set thread
        if (hts_set_opt(in, HTS_OPT_THREAD_POOL, &threadPool) != 0) {
            throw std::runtime_error("Cannot set thread pool for input bam file " + BamFile);
        }

        if (writeOutputBam) {
            // output file mangement
            std::string writeBamFile = params.resultPrefix + "." + params.outputFormat;
            std::cerr << "set output bam file : " + writeBamFile << std::endl;
            // open output bam file
            out = hts_open(writeBamFile.c_str(), (params.outputFormat == "bam" ? "wb" : "wc" ));
            // load reference file
            hts_set_fai_filename(out, params.fastaFile.c_str());
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
    ParsingBamMode mode,
    bool writeOutputBam,
    bool mappingQualityFilter
): mode(mode), writeOutputBam(writeOutputBam), mappingQualityFilter(mappingQualityFilter)
{

}

HaplotagBamParser::~HaplotagBamParser(){

}

void HaplotagBamParser::parsingBam(
    const std::string &BamFile, 
    const HaplotagParameters &params, 
    const std::vector<std::string> &chrVec, 
    const std::map<std::string, int> &chrLength, 
    std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat, 
    std::map<Genome, VCF_Info> &vcfSet, 
    const Genome& genmoeType
){
    try{

        if(mode == ParsingBamMode::MULTI_THREAD && writeOutputBam) {
            throw std::runtime_error("Cannot set Multi Thread mode when writeOutputBam is true");
        }

        if(chrVec.size() == 0){
            throw std::runtime_error("chrVec is empty");
        }

        if(chrLength.size() == 0){
            throw std::runtime_error("chrLength is empty");
        }
        
        std::vector<int> last_pos;
        // get the last variant position of the reference
        getLastVarPos(last_pos, chrVec,mergedChrVarinat, genmoeType);
        // reference fasta parser
        FastaParser fastaParser(params.fastaFile, chrVec, last_pos, params.numThreads);

        // init data structure and get core n
        htsThreadPool threadPool = {NULL, 0};
        // creat thread pool
        if (!(threadPool.pool = hts_tpool_init(params.numThreads))) {
            throw std::runtime_error("Error creating thread pool");
        }

        switch(mode){
            case ParsingBamMode::SINGLE_THREAD:
                // Process the BAM file using a single thread
                std::cerr<< "[single thread]";
                processBamWithOutput(BamFile, params, chrVec, chrLength, fastaParser, threadPool, mergedChrVarinat, vcfSet, genmoeType);
                break;

            case ParsingBamMode::MULTI_THREAD:
                // Process the BAM file using multiple threads
                std::cerr<< "[multi thread]";
                processBamParallel(BamFile, params, chrVec, chrLength, fastaParser, threadPool, mergedChrVarinat, vcfSet, genmoeType);
                break;

            default:
                throw std::runtime_error("Unsupported parsing bam mode: " + std::to_string(mode));
        }

        hts_tpool_destroy(threadPool.pool);
    }catch(const std::exception& e){
        std::cerr << "[ERROR](HaplotagBamParser): " << e.what() << std::endl;
        exit(1);
    }

    return;
}

void HaplotagBamParser::processBamParallel(
    const std::string &BamFile, 
    const HaplotagParameters &params, 
    const std::vector<std::string> &chrVec, 
    const std::map<std::string, int> &chrLength, 
    const FastaParser &fastaParser,
    htsThreadPool &threadPool,
    std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat, 
    std::map<Genome, VCF_Info> &vcfSet, 
    const Genome& genmoeType
){
    // loop all chromosome
    #pragma omp parallel for schedule(dynamic) num_threads(params.numThreads) 
    for(auto chr : chrVec ){
        // bam file resource allocation
        BamFileRAII bamRAII(BamFile, params.fastaFile, threadPool, params);
        //create the chromosome processor
        auto chrProcessor = createProcessor(chr);
        chrProcessor->processSingleChromosome(chr, chrLength, params, fastaParser, mergedChrVarinat, bamRAII, genmoeType, vcfSet);
    }
}

void HaplotagBamParser::processBamWithOutput(
    const std::string &BamFile, 
    const HaplotagParameters &params, 
    const std::vector<std::string> &chrVec, 
    const std::map<std::string, int> &chrLength, 
    const FastaParser &fastaParser,
    htsThreadPool &threadPool,
    std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat, 
    std::map<Genome, VCF_Info> &vcfSet, 
    const Genome& genmoeType
){
    // bam file resource allocation
    BamFileRAII bamRAII(BamFile, params.fastaFile, threadPool, params, writeOutputBam);
    // loop all chromosome
    for(auto chr : chrVec ){
        //create the chromosome processor
        std::time_t begin = time(NULL);
        std::cerr<<"chr: " << chr << " ... " ;
        auto chrProcessor = createProcessor(chr);
        chrProcessor->processSingleChromosome(chr, chrLength, params, fastaParser, mergedChrVarinat, bamRAII, genmoeType, vcfSet);
        std::cerr<< difftime(time(NULL), begin) << "s\n";
    }
}

void HaplotagBamParser::getLastVarPos(
    std::vector<int>& last_pos, 
    const std::vector<std::string>& chrVec, 
    std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat,
    const Genome& genmoeType
){
    for( auto chr : chrVec ){
        bool existLastPos = false;

        for (auto lastVariantIter = mergedChrVarinat[chr].rbegin(); lastVariantIter != mergedChrVarinat[chr].rend(); ++lastVariantIter) {
            if (genmoeType == NORMAL) {
                if ((*lastVariantIter).second.isExists(NORMAL) && (*lastVariantIter).second.Variant[NORMAL].isExistPhasedSet()) {
                    last_pos.push_back((*lastVariantIter).first);
                        existLastPos = true;
                        break;
                    }
            }
            else if (genmoeType == TUMOR) {
                if ((*lastVariantIter).second.isExists(TUMOR) || 
                    ((*lastVariantIter).second.isExists(NORMAL) && (*lastVariantIter).second.Variant[NORMAL].isExistPhasedSet())) {
                    last_pos.push_back((*lastVariantIter).first);
                    existLastPos = true;
                    break;
                }
            }else{
                std::cerr << "ERROR (germlineGetRefLastVarPos) => unsupported genome type: " << genmoeType << std::endl;
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
    const std::string& chr,
    const std::map<std::string, int>& chrLength,
    const HaplotagParameters& params, 
    const FastaParser& fastaParser,
    std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat, 
    BamFileRAII& bamRAII,
    const Genome& genmoeType,
    std::map<Genome, VCF_Info> &vcfSet
){
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

    std::string region = !params.region.empty() ? params.region : std::string(chr + ":1-" + std::to_string(chrLength.at(chr)));
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
            processRead(*bamRAII.aln, *bamRAII.bamHdr, chr, params, genmoeType, currentVariants, firstVariantIter, vcfSet, ref_string);
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

CigarParser::CigarParser(int& ref_pos, int& query_pos)
: ref_pos(ref_pos), query_pos(query_pos){

    aln = nullptr;
    bamHdr = nullptr;
    chrName = nullptr;
    ref_string = nullptr;
    hpCount = nullptr;
    norCountPS = nullptr;
    variantsHP = nullptr;
}

CigarParser::~CigarParser(){

};


void CigarParser::parsingCigar(
    const bam1_t& aln,
    const bam_hdr_t& bamHdr,
    const std::string& chrName,
    const HaplotagParameters& params,
    std::map<int, MultiGenomeVar>::iterator& firstVariantIter,
    std::map<int, MultiGenomeVar>& currentVariants,
    const std::string& ref_string,
    std::map<int, int>& hpCount,
    std::map<int, int>& variantsHP,
    std::map<int, int>& norCountPS
){
    this->aln = &aln;
    this->bamHdr = &bamHdr;
    this->chrName = &chrName;
    this->params = &params;
    this->ref_string = &ref_string;
    this->hpCount = &hpCount;
    this->variantsHP = &variantsHP;
    this->norCountPS = &norCountPS;

    // Skip variants that are to the left of this read
    while (firstVariantIter != currentVariants.end() && (*firstVariantIter).first < aln.core.pos) {
        firstVariantIter++;
    }

    if (firstVariantIter == currentVariants.end()) {
        return;
    }

    currentVariantIter = firstVariantIter;

    // position relative to reference
    ref_pos = aln.core.pos;
    // position relative to read
    query_pos = 0;

    // reading cigar to detect snp on this read
    int aln_core_n_cigar = int(aln.core.n_cigar);
    for (int i = 0; i < aln_core_n_cigar; i++) {
        uint32_t* cigar = bam_get_cigar(&aln);
        int cigar_op = bam_cigar_op(cigar[i]);
        int length = bam_cigar_oplen(cigar[i]);

        // iterator next variant
        while (currentVariantIter != currentVariants.end() && (*currentVariantIter).first < ref_pos) {
            currentVariantIter++;
        }

        // CIGAR operators: MIDNSHP=X correspond 012345678
        // 0: alignment match (can be a sequence match or mismatch)
        // 7: sequence match
        // 8: sequence mismatch
        if( cigar_op == 0 || cigar_op == 7 || cigar_op == 8 ){ 
                while( currentVariantIter != currentVariants.end() && (*currentVariantIter).first < ref_pos + length){
                    int offset = (*currentVariantIter).first - ref_pos;
                    if( offset < 0){
                    }
                    else{
                        uint8_t *q = bam_get_seq(&aln);
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
                while( currentVariantIter != currentVariants.end() && (*currentVariantIter).first < ref_pos + length){
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
            std::cerr << "Alignment find unsupported CIGAR operation from read: " << bam_get_qname(&aln) << "\n";
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



void GermlineJudgeBase::germlineJudgeSnpHap(
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
                std::cerr << "ERROR (germlineJudgeSnpHap) => can't find the position:" 
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


void GermlineJudgeBase::germlineJudgeDeletionHap(
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

void GermlineJudgeBase::germlineJudgeSVHap(const bam1_t &aln, std::map<Genome, VCF_Info> &vcfSet, std::map<int, int>& hpCount, const int& genmoeType){
    auto readIter = vcfSet[Genome::NORMAL].readSVHapCount.find(bam_get_qname(&aln));
    if( readIter != vcfSet[Genome::NORMAL].readSVHapCount.end()){
        hpCount[SnpHP::GERMLINE_H1] += vcfSet[Genome::NORMAL].readSVHapCount[bam_get_qname(&aln)][0];
        hpCount[SnpHP::GERMLINE_H2] += vcfSet[Genome::NORMAL].readSVHapCount[bam_get_qname(&aln)][1];
    }
}

int GermlineJudgeBase::germlineDetermineReadHap(
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

void SomaticJudgeBase::SomaticJudgeSnpHP(std::map<int, MultiGenomeVar>::iterator &currentVariantIter, std::string chrName, std::string base, std::map<int, int> &hpCount, std::map<int, int> &norCountPS, std::map<int, int> &tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *tumorAllelePosVec){
    int curPos = (*currentVariantIter).first;
    auto& curVar = (*currentVariantIter).second;

    // normal & tumor SNP at the current position (base on normal phased SNPs)
    // both normal and tumor samples that do not exist in the high-confidence set
    if(curVar.isExists(NORMAL) && curVar.isExists(TUMOR)){

        // the tumor & normal SNP GT are phased heterozygous 
        if((curVar.Variant[NORMAL].is_phased_hetero) && (curVar.Variant[TUMOR].is_phased_hetero)){   
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

                if(norHP1 == base){
                    hpCount[1]++;
                    //Germline mutation
                    if(variantsHP != nullptr) (*variantsHP)[curPos] = SnpHP::GERMLINE_H1;
                }else if(norHP2 == base){                
                    hpCount[2]++;
                    //Germline mutation
                    if(variantsHP != nullptr) (*variantsHP)[curPos] = SnpHP::GERMLINE_H2;

                }else{
                    std::cerr<< "ERROR : (phased hetero)normal & (phased hetero)tumor not match base at position => chr:" << chrName << " pos: " << curPos << "\n";
                    exit(1);
                }
                norCountPS[curVar.Variant[NORMAL].PhasedSet]++;
            }
        //the normal SNP GT is phased heterozgous & the tumor SNP GT is unphased heterozgous 
        }else if((curVar.Variant[NORMAL].is_phased_hetero) && (curVar.Variant[TUMOR].is_unphased_hetero)){   
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

                if(norHP1 == base){
                    hpCount[1]++;
                    //Germline mutation
                    if(variantsHP != nullptr) (*variantsHP)[curPos] = SnpHP::GERMLINE_H1;
                }else if(norHP2 == base){                
                    hpCount[2]++;
                    //Germline mutation
                    if(variantsHP != nullptr) (*variantsHP)[curPos] = SnpHP::GERMLINE_H2;
                }else{
                    std::cerr<< "ERROR : (phased hetero)normal & (unphased hetero)tumor not match base at position => chr:" << chrName << " pos: " << curPos << "\n";
                    exit(1);
                }
                norCountPS[curVar.Variant[NORMAL].PhasedSet]++;
            }

        //the normal SNP GT is phased heterozgous & the tumor SNP GT is homozygous 
        }else if((curVar.Variant[NORMAL].is_phased_hetero) && (curVar.Variant[TUMOR].is_homozygous)){   
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

                if(norHP1 == base){
                    hpCount[1]++;
                    //Germline mutation
                    if(variantsHP != nullptr) (*variantsHP)[curPos] = SnpHP::GERMLINE_H1;
                }else if(norHP2 == base){                
                    hpCount[2]++;
                    //Germline mutation
                    if(variantsHP != nullptr) (*variantsHP)[curPos] = SnpHP::GERMLINE_H2;
                }else{
                    std::cerr<< "ERROR : (phased hetero)normal & (Homo)tumor not match base at position => chr:" << chrName << " pos: " << curPos << "\n";
                    exit(1);
                }
                norCountPS[curVar.Variant[NORMAL].PhasedSet]++;
            }
        }
    // only normal SNP at the current position
    }else if(curVar.isExists(NORMAL)){
        // the normal SNP GT is phased heterozgous SNP
        if((curVar.Variant[NORMAL].is_phased_hetero)){
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
    // only tumor SNP at the current position
    }else if(curVar.isExists(TUMOR)){
        //the tumor SNP GT is phased heterozygous
        if(curVar.Variant[TUMOR].is_phased_hetero == true){
            if(curVar.Variant[TUMOR].allele.Ref == base || curVar.Variant[TUMOR].allele.Alt == base){
                if(!curVar.Variant[TUMOR].isExistPhasedSet()){
                    std::cerr<< curPos << "\t"
                             << curVar.Variant[TUMOR].allele.Ref << "\t"
                             << curVar.Variant[TUMOR].allele.Alt << "\n";
                    exit(EXIT_SUCCESS);
                }else{
                    OnlyTumorSNPjudgeHP(chrName, curPos, curVar, base, hpCount, &tumCountPS, variantsHP, tumorAllelePosVec);
                }
            }
        //the tumor SNP GT is unphased heterozygous
        }else if(curVar.Variant[TUMOR].is_unphased_hetero == true){
            if(curVar.Variant[TUMOR].allele.Ref == base || curVar.Variant[TUMOR].allele.Alt == base){
                OnlyTumorSNPjudgeHP(chrName, curPos, curVar, base, hpCount, nullptr, variantsHP, tumorAllelePosVec);
            }           
        //the tumor SNP GT is homozygous
        }else if(curVar.Variant[TUMOR].is_homozygous == true){
            if(curVar.Variant[TUMOR].allele.Ref == base || curVar.Variant[TUMOR].allele.Alt == base){
                OnlyTumorSNPjudgeHP(chrName, curPos, curVar, base, hpCount, nullptr, variantsHP, tumorAllelePosVec);
            }
        }
    }
}

void SomaticJudgeBase::OnlyTumorSNPjudgeHP(const std::string &chrName, int &curPos, MultiGenomeVar &curVar, std::string base, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *tumorAllelePosVec){

}

int SomaticJudgeBase::determineReadHP(std::map<int, int> &hpCount, int &pqValue, std::map<int, int> &norCountPS, double &norHPsimilarity, double &tumHPsimilarity, double percentageThreshold, int *totalHighSimilarity, int *totalCrossTwoBlock, int *totalWithOutVaraint){
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

void chrReadHpResult::recordReadHp(int &pos, int &hpResult, int &BaseHP){
    posReadHpResult[pos].readHpCounter[hpResult]++;
    
    if(hpResult != ReadHP::unTag){
        if(BaseHP == SnpHP::SOMATIC_H3){
            if(hpResult != ReadHP::H1_1 && hpResult != ReadHP::H2_1 && hpResult != ReadHP::H3){
                std::cerr << "Error(recordReadHp) => error read hp : BaseHP: " <<BaseHP << " readHP: " << hpResult << " pos: " << pos+1 << std::endl; 
                exit(1);
            }
            posReadHpResult[pos].somaticSnpH3count++;
            posReadHpResult[pos].somaticBaseReadHpCounter[hpResult]++;
        }
    }
}

void chrReadHpResult::recordDeriveHp(int &pos, int &deriveHP, float deriveHPsimilarity){
    if(deriveHP != SnpHP::GERMLINE_H1 && deriveHP != SnpHP::GERMLINE_H2 && deriveHP != SnpHP::NONE_SNP){
        std::cerr << "Error(recordDeriveHp) => error derive hp : pos: " <<pos+1 << " deriveHP: " << deriveHP << std::endl; 
        exit(1);        
    }
    posReadHpResult[pos].deriveHP = deriveHP;
    if(deriveHPsimilarity != 0.0){
        posReadHpResult[pos].deriveHPsimilarVec.emplace_back(deriveHPsimilarity);
        if(deriveHPsimilarity != 1.0){
            // std::cout << "deriveHPsimilarity: " << deriveHPsimilarity << "\n";
            // std::cout << "deriveHPsimilarityVec: " << varReadHpResult[pos].deriveHPsimilarVec.back() << "\n";
        }
    }
}

void chrReadHpResult::recordAlignCoverRegion(int& curVarPos, int &startPos, int &endPos){
    if(posReadHpResult[curVarPos].coverRegionStartPos > startPos){
        posReadHpResult[curVarPos].coverRegionStartPos = startPos;
    }
    if(posReadHpResult[curVarPos].coverRegionEndPos < endPos){
        posReadHpResult[curVarPos].coverRegionEndPos = endPos;
    }
}


ReadHpDistriLog::ReadHpDistriLog(){

}

ReadHpDistriLog::~ReadHpDistriLog(){

}

void ReadHpDistriLog::loadChrKey(const std::string &chr){
    chrVarReadHpResult[chr] = chrReadHpResult();
}

chrReadHpResult* ReadHpDistriLog::getChrHpResultsPtr (const std::string &chr){
    return &(chrVarReadHpResult[chr]);
}

void ReadHpDistriLog::recordChrReadHp(const std::string &chr, int &pos, int &hpResult, int &BaseHP){
    //only use in single thread scenario
    chrVarReadHpResult[chr].recordReadHp(pos, hpResult, BaseHP);
}

void ReadHpDistriLog::recordChrDeriveHp(const std::string &chr, int &pos, int &deriveHP, float deriveHPsimilarity){
    //only use in single thread scenario
    chrVarReadHpResult[chr].recordDeriveHp(pos, deriveHP, deriveHPsimilarity);
}

void ReadHpDistriLog::recordChrAlignCoverRegion(const std::string &chr, int &pos, int &startPos, int &endPos){
    //only use in single thread scenario
    chrVarReadHpResult[chr].recordAlignCoverRegion(pos, startPos, endPos);
}


void ReadHpDistriLog::writeReadHpDistriLog(const HaplotagParameters &params, std::string logPosfix, const std::vector<std::string> &chrVec){
    std::ofstream *readHpDistriLog=NULL;
    readHpDistriLog=new std::ofstream(params.resultPrefix + logPosfix);

    int somaticSnpCount = 0;
    for(auto chr: chrVec){
        if(!chrVarReadHpResult[chr].posReadHpResult.empty()){
            somaticSnpCount += chrVarReadHpResult[chr].posReadHpResult.size();
        }
    }

    if(!readHpDistriLog->is_open()){
        std::cerr<< "Fail to open write file: " << params.resultPrefix + logPosfix << "\n";
        exit(1);
    }else{
        (*readHpDistriLog) << "####################################\n";
        (*readHpDistriLog) << "# Somatic SNP read HP distribution #\n";
        (*readHpDistriLog) << "####################################\n";
        (*readHpDistriLog) << "##MappingQualityThreshold:"        << params.qualityThreshold << "\n";
        (*readHpDistriLog) << "##SomaticSNP: " << somaticSnpCount << "\n";
        (*readHpDistriLog) << "#Chr\t"
                            << "Pos\t"
                            << "DeriveHP\t"
                            << "DeriveHPsimilarity\t\t"
                            << "AltCount\t"
                            << "somaticBase_H1-1\t"
                            << "somaticBase_H2-1\t"
                            << "somaticBase_H3\t\t"
                            << "HP1read\t"
                            << "HP2read\t"
                            << "HP1-1read\t"
                            << "HP2-1read\t"
                            << "HP3read\t"
                            << "untagRead\t"
                            << "HP1ratio\t"
                            << "HP2ratio\t"
                            << "HP1-1ratio\t"
                            << "HP2-1ratio\t"
                            << "HP3ratio\n";
    }

    for(auto chr: chrVec){
        std::map<int, ReadHpResult>::iterator curVarReadHpIter = chrVarReadHpResult[chr].posReadHpResult.begin();
        while(curVarReadHpIter != chrVarReadHpResult[chr].posReadHpResult.end()){
            int pos = (*curVarReadHpIter).first + 1;
            int HP1readCount = (*curVarReadHpIter).second.readHpCounter[ReadHP::H1];
            int HP1_1readCount = (*curVarReadHpIter).second.readHpCounter[ReadHP::H1_1];

            int HP2readCount = (*curVarReadHpIter).second.readHpCounter[ReadHP::H2];
            int HP2_1readCount = (*curVarReadHpIter).second.readHpCounter[ReadHP::H2_1];

            int HP3readCount = (*curVarReadHpIter).second.readHpCounter[ReadHP::H3];

            int totaltagRead = HP1readCount + HP2readCount + HP3readCount + HP1_1readCount + HP2_1readCount;

            float HP1readRatio = (float)HP1readCount / (float)totaltagRead; 
            float HP1_1readRatio = (float)HP1_1readCount / (float)totaltagRead; 

            float HP2readRatio = (float)HP2readCount / (float)totaltagRead; 
            float HP2_1readRatio = (float)HP2_1readCount / (float)totaltagRead; 

            float HP3readRatio = (float)HP3readCount / (float)totaltagRead; 

            float meanDeriveHPsimilarity = 0.0;
            if((*curVarReadHpIter).second.deriveHPsimilarVec.size() != 0){
                float size = (*curVarReadHpIter).second.deriveHPsimilarVec.size();
                for(float deriveHPsimilarity: (*curVarReadHpIter).second.deriveHPsimilarVec){
                    meanDeriveHPsimilarity += deriveHPsimilarity;
                    // std::cout << deriveHPsimilarity << "\t";
                }
                meanDeriveHPsimilarity /= size;
            }
            (*readHpDistriLog) << std::fixed << std::setprecision(3) 
                                << chr << "\t"
                                << pos << "\t"
                                << "H" << (*curVarReadHpIter).second.deriveHP << "\t"
                                << meanDeriveHPsimilarity << "\t\t"
                                << (*curVarReadHpIter).second.somaticSnpH3count << "\t"
                                << (*curVarReadHpIter).second.somaticBaseReadHpCounter[ReadHP::H1_1] << "\t"
                                << (*curVarReadHpIter).second.somaticBaseReadHpCounter[ReadHP::H2_1] << "\t"
                                << (*curVarReadHpIter).second.somaticBaseReadHpCounter[ReadHP::H3] << "\t\t"
                                << HP1readCount << "\t"
                                << HP2readCount << "\t\t"
                                << HP1_1readCount << "\t"
                                << HP2_1readCount << "\t"
                                << HP3readCount << "\t"
                                << (*curVarReadHpIter).second.readHpCounter[ReadHP::unTag] << "\t"
                                << HP1readRatio << "\t"
                                << HP2readRatio << "\t"
                                << HP1_1readRatio << "\t"
                                << HP2_1readRatio << "\t"
                                << HP3readRatio << "\n";
            curVarReadHpIter++;
        }
    }
    (*readHpDistriLog).close();
    delete readHpDistriLog;
    readHpDistriLog = nullptr;
}

void ReadHpDistriLog::writePosCoverRegionLog(const HaplotagParameters &params, std::string logPosfix, const std::vector<std::string> &chrVec){
    std::ofstream *posCoverRegionLog=NULL;
    posCoverRegionLog=new std::ofstream(params.resultPrefix + logPosfix);

    int somaticSnpCount = 0;
    for(auto chr: chrVec){
        if(!chrVarReadHpResult[chr].posReadHpResult.empty()){
            somaticSnpCount += chrVarReadHpResult[chr].posReadHpResult.size();
        }
    }

    if(!posCoverRegionLog->is_open()){
        std::cerr<< "Fail to open write file: " << params.resultPrefix + logPosfix << "\n";
        exit(1);
    }else{
        (*posCoverRegionLog) << "############################\n";
        (*posCoverRegionLog) << "# Somatic SNP cover region #\n";
        (*posCoverRegionLog) << "############################\n";
        (*posCoverRegionLog) << "##MappingQualityThreshold:"        << params.qualityThreshold << "\n";
        (*posCoverRegionLog) << "##SomaticSNP: " << somaticSnpCount << "\n";
        (*posCoverRegionLog) << "#Chr\t"
                              << "Pos\t"
                              << "Type\t"
                              << "StartPos\t"
                              << "EndPos\n";
    }

    for(auto chr: chrVec){
        std::map<int, ReadHpResult>::iterator curVarReadHpIter = chrVarReadHpResult[chr].posReadHpResult.begin();
        while(curVarReadHpIter != chrVarReadHpResult[chr].posReadHpResult.end()){
            int pos = (*curVarReadHpIter).first + 1;

            (*posCoverRegionLog) << std::fixed << std::setprecision(3) 
                                << chr << "\t"
                                << pos << "\t"
                                << "somatic" << "\t"
                                << (*curVarReadHpIter).second.coverRegionStartPos << "\t"
                                << (*curVarReadHpIter).second.coverRegionEndPos << "\n";
            curVarReadHpIter++;
        }
    }
    (*posCoverRegionLog).close();
    delete posCoverRegionLog;
    posCoverRegionLog = nullptr;
}

void ReadHpDistriLog::writeTagReadCoverRegionLog(const HaplotagParameters &params, std::string logPosfix, const std::vector<std::string> &chrVec, std::map<std::string, int> &chrLength){
    std::ofstream *tagReadCoverRegionLog=NULL;
    tagReadCoverRegionLog=new std::ofstream(params.resultPrefix + logPosfix);

    std::map<std::string, std::vector<coverRegionInfo>> coverRegion;

    // merge cover region from different SNPs
    for (const auto& chr : chrVec){
        auto curVarReadHpIter = chrVarReadHpResult[chr].posReadHpResult.begin();

        if (curVarReadHpIter == chrVarReadHpResult[chr].posReadHpResult.end()) {
            continue;
        }

        int curStartPos = curVarReadHpIter->second.coverRegionStartPos;
        int curEndPos = curVarReadHpIter->second.coverRegionEndPos;

        while (curVarReadHpIter != chrVarReadHpResult[chr].posReadHpResult.end()){

            auto nextVarReadHpIter = std::next(curVarReadHpIter);
            if(nextVarReadHpIter != chrVarReadHpResult[chr].posReadHpResult.end()){
                int nextStartPos = nextVarReadHpIter->second.coverRegionStartPos;
                int nextEndPos = nextVarReadHpIter->second.coverRegionEndPos;

                // region not overlap
                if (curEndPos < nextStartPos) {
                    coverRegionInfo regionInfo;
                    regionInfo.startPos = curStartPos;
                    regionInfo.endPos = curEndPos;
                    regionInfo.length = curEndPos - curStartPos + 1;
                    coverRegion[chr].emplace_back(regionInfo);
                    curStartPos = nextStartPos;
                    curEndPos = nextEndPos;
                }else{
                    // regions overlap
                    curStartPos = std::min(curStartPos, nextStartPos);
                    curEndPos = std::max(curEndPos, nextEndPos);
                }
            // last region
            }else{
                coverRegionInfo regionInfo;
                regionInfo.startPos = curStartPos;
                regionInfo.endPos = curEndPos;
                regionInfo.length = curEndPos - curStartPos + 1;
                coverRegion[chr].emplace_back(regionInfo);
                break;
            }

            curVarReadHpIter++;
        }
    }

    std::map<std::string, float> coverRegionRatio;
    long long totalChrLength = 0;
    long long totalChrCoverLength = 0;
    double totalChrCoverageRatio = 0.0;

    // calculate cover region ratio
    for(auto chr: chrVec){
        int totalCoverLength = 0;
        auto coverRegionIter = coverRegion[chr].begin();
        while(coverRegionIter != coverRegion[chr].end()){
            totalCoverLength += (*coverRegionIter).length;
            coverRegionIter++;
        }
        coverRegionRatio[chr] = (float)totalCoverLength / (float)chrLength[chr];
        totalChrLength += chrLength[chr];
        totalChrCoverLength += totalCoverLength;
    }
    totalChrCoverageRatio = (double)totalChrCoverLength / (double)totalChrLength;

    if(!tagReadCoverRegionLog->is_open()){
        std::cerr<< "Fail to open write file: " << params.resultPrefix + logPosfix << "\n";
        exit(1);
    }else{
        (*tagReadCoverRegionLog) << "##################################\n";
        (*tagReadCoverRegionLog) << "# Somatic reads cover region bed #\n";
        (*tagReadCoverRegionLog) << "##################################\n";
        (*tagReadCoverRegionLog) << "##MappingQualityThreshold: "   << params.qualityThreshold << "\n";
        (*tagReadCoverRegionLog) << "##----Chr coverage ratio----\n";
        (*tagReadCoverRegionLog) << "##Total chr coverage ratio: " << totalChrCoverageRatio << "\n";
        for(auto chr: chrVec){
            (*tagReadCoverRegionLog) <<"##" << chr << ":" << coverRegionRatio[chr] << "\n";
        }
        (*tagReadCoverRegionLog) << "#Chr\t"
                                 << "StartPos\t"
                                 << "EndPos\n";
                                // << "length\n";
    }

    for(auto chr: chrVec){
        auto coverRegionIter = coverRegion[chr].begin();
        while(coverRegionIter != coverRegion[chr].end()){

            (*tagReadCoverRegionLog) << std::fixed << std::setprecision(3) 
                                     << chr << "\t"
                                     << (*coverRegionIter).startPos << "\t"
                                     << (*coverRegionIter).endPos << "\n";
                                     //<< (*coverRegionIter).length << "\n";
            coverRegionIter++;
        }
    }
    (*tagReadCoverRegionLog).close();
    delete tagReadCoverRegionLog;
    tagReadCoverRegionLog = nullptr;
}


void ReadHpDistriLog::removeNotDeriveByH1andH2pos(const std::vector<std::string> &chrVec){
    for (const auto& chr : chrVec) {
        auto& chrResult = chrVarReadHpResult[chr];
        for (auto it = chrResult.posReadHpResult.begin(); it != chrResult.posReadHpResult.end(); ) {
            if (!it->second.existDeriveByH1andH2) {
                //std::cerr << "Removed position not derived by H1 and H2: " << chr << " " << it->first << std::endl;
                it = chrResult.posReadHpResult.erase(it);
            } else {
                ++it;
            }
        }
    }
}


VcfParser::VcfParser(bool tagTumorMode){
    this->tagTumorMode = tagTumorMode;
    reset();
}

VcfParser::VcfParser(){
    reset();
}

VcfParser::~VcfParser(){

}

void VcfParser::reset(){
    parseSnpFile = false;
    parseSVFile = false;
    parseMODFile = false;
    integerPS = false;
    psIndex.clear();
}

void VcfParser::variantParser(std::string &variantFile, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat){

    if( variantFile.find("gz") != std::string::npos ){
        // .vcf.gz 
        compressParser(variantFile, Info, mergedChrVarinat);
    }
    else if( variantFile.find("vcf") != std::string::npos ){
        // .vcf
        unCompressParser(variantFile, Info, mergedChrVarinat);
    }
    else{
        std::cerr<<"file: "<< variantFile << "\nnot vcf file. please check filename extension\n";
        exit(EXIT_FAILURE);
    }
    return;
}

void VcfParser::compressParser(std::string &variantFile, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat){
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
                parserProcess(input, Info, mergedChrVarinat);
            }
            // any trailing data in [eol, end) now is a partial line
            offset = std::copy(cur, end, buffer);
        }
        gzclose (file);
        free(buffer);
    }    
}

void VcfParser::unCompressParser(std::string &variantFile, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat){
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
            parserProcess(input, Info, mergedChrVarinat);
        }
    }
}

void VcfParser::setParseSnpFile(bool parseSnpFile){
    this->parseSnpFile = parseSnpFile;
}

void VcfParser::setParseSVFile(bool parseSVFile){
    this->parseSVFile = parseSVFile;
} 

void VcfParser::setParseMODFile(bool parseMODFile){
    this->parseMODFile = parseMODFile;
}

bool VcfParser::getParseSnpFile(){
    return this->parseSnpFile;
}

void VcfParser::parserProcess(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat){
    if( input.substr(0, 2) == "##" && parseSnpFile){
        if( input.find("contig=")!= std::string::npos ){
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
            if( input.find("Type=Integer")!= std::string::npos ){
                integerPS = true;
            }
            else if( input.find("Type=String")!= std::string::npos ){
                integerPS = false;
                std::cerr<< "PS type is String. Auto index to integer ... ";
            }
            else{
                std::cerr<< "ERROR: not found PS type (Type=Integer or Type=String).\n"; 
                exit(EXIT_SUCCESS);
            }
        }
    }
    else if ( input.substr(0, 1) == "#" ){
        
    }
    else{
        std::istringstream iss(input);
        std::vector<std::string> fields((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());

        if( fields.size() == 0 )
            return;

        // trans to 0-base
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
                varData.is_phased_hetero = true;
                varData.setVariantType();
                
                if(integerPS){
                    varData.PhasedSet = std::stoi(psValue);
                }
                else{
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

                if(Info.gene_type == NORMAL){
                    mergedChrVarinat[chr][pos].Variant[NORMAL] = varData;
                }else if(Info.gene_type == TUMOR){
                    mergedChrVarinat[chr][pos].Variant[TUMOR] = varData;
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
                // In which haplotype does SV occur
                if( fields[9][modifu_start] == '0' && fields[9][modifu_start+2] == '1' ){
                    modHaplotype = 1;
                }
                else if( fields[9][modifu_start] == '1' && fields[9][modifu_start+2] == '0' ){
                    modHaplotype = 0;
                }
                
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
        else if((tagTumorMode == true)){
            //homozygous SNPs
            if( fields[9][modifu_start] == '1' && fields[9][modifu_start+1] == '/' && fields[9][modifu_start+2] == '1' ){
                if(parseSnpFile){

                    VarData varData;
                    varData.allele.Ref = fields[3];
                    varData.allele.Alt = fields[4];
                    varData.is_homozygous = true;
                    varData.setVariantType();

                    if(Info.gene_type == NORMAL){
                        mergedChrVarinat[chr][pos].Variant[NORMAL] = varData;
                    }else if(Info.gene_type == TUMOR){
                        mergedChrVarinat[chr][pos].Variant[TUMOR] = varData;
                    }
                }
            //unphased heterozygous
            }else if( fields[9][modifu_start] == '0' && fields[9][modifu_start+1] == '/' && fields[9][modifu_start+2] == '1' ){
                if(parseSnpFile){

                    VarData varData;
                    varData.allele.Ref = fields[3];
                    varData.allele.Alt = fields[4];
                    
                    varData.is_unphased_hetero = true;
                    varData.setVariantType();

                    if(Info.gene_type == NORMAL){
                        mergedChrVarinat[chr][pos].Variant[NORMAL] = varData;
                    }else if(Info.gene_type == TUMOR){
                        mergedChrVarinat[chr][pos].Variant[TUMOR] = varData;
                    }
                }
            }
        }
    }
}
