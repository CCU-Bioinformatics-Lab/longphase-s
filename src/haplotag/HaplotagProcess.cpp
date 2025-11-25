#include "HaplotagProcess.h"

HaplotagProcess::HaplotagProcess(HaplotagParameters &params):
params(params),chrVec(nullptr),chrLength(nullptr),readStats(),processBegin(time(NULL))
{
    //initialize variable
    vcfSet[Genome::NORMAL] = VCF_Info{.sample = Genome::NORMAL};
    vcfSet[Genome::TUMOR] = VCF_Info{.sample = Genome::TUMOR};
    vcfSet[Genome::TRUTH_SOMATIC] = VCF_Info{.sample = Genome::TRUTH_SOMATIC};

    chrMultiVariants = new std::map<std::string, std::map<int, MultiGenomeVar>>();

}

HaplotagProcess::~HaplotagProcess(){
    delete chrMultiVariants;
};

void HaplotagProcess::printParamsMessage(){
    std::cerr<< "LongPhase-S v" << params.bamCfg.version << " - Haplotag\n";
    std::cerr<< "\n";
    std::cerr<< "phased SNP file:   " << params.snpFile             << "\n";
    std::cerr<< "phased SV file:    " << params.svFile              << "\n";
    std::cerr<< "phased MOD file:   " << params.modFile             << "\n";
    std::cerr<< "input bam file:    " << params.bamFile             << "\n";
    std::cerr<< "input ref file:    " << params.fastaFile           << "\n";
    std::cerr<< "output bam file:   " << params.bamCfg.resultPrefix + "." + params.bamCfg.outputFormat << "\n";
    std::cerr<< "number of threads: " << params.bamCfg.numThreads          << "\n";
    std::cerr<< "write log file:    " << (params.bamCfg.writeReadLog ? "true" : "false") << "\n";
    std::cerr<< "log file:          " << (params.bamCfg.writeReadLog ? (params.bamCfg.resultPrefix+".out") : "") << "\n";
    std::cerr<< "-------------------------------------------\n";
    std::cerr<< "tag region:                    " << (!params.bamCfg.region.empty() ? params.bamCfg.region : "all") << "\n";
    std::cerr<< "filter mapping quality below:  " << params.bamCfg.qualityThreshold    << "\n";
    std::cerr<< "percentage threshold:          " << params.bamCfg.percentageThreshold << "\n";
    std::cerr<< "tag supplementary:             " << (params.bamCfg.tagSupplementary ? "true" : "false") << "\n";
    std::cerr<< "-------------------------------------------\n";
}

void HaplotagProcess::pipelineProcess()
{
    printParamsMessage();
    // decide on the type of tagging for VCF and BAM files
    Genome tagSample = NORMAL;

    VcfParser vcfParser(tagSample);
    // load SNP, SV, MOD vcf file
    parseVariantFiles(vcfParser);
    //decide which genome sample chrVec and chrLength belong to
    setChrVecAndChrLength();
    // update chromosome processing based on region
    setProcessingChromRegion();
    // tag read
    tagRead(params,params.bamFile, tagSample);
    // postprocess after haplotag
    postprocessForHaplotag();

    printExecutionReport();
    
    return;
};

void HaplotagProcess::parseVariantFiles(VcfParser& vcfParser){
    // load SNP vcf file
    std::time_t begin = time(NULL);
    std::cerr<< getNormalSnpParsingMessage();

    vcfParser.setParseSnpFile(true);
    vcfParser.parsingVCF(params.snpFile, vcfSet[Genome::NORMAL], *chrMultiVariants);
    vcfParser.reset();
    std::cerr<< difftime(time(NULL), begin) << "s\n";

    // load SV vcf file
    if(params.svFile!=""){
        begin = time(NULL);
        std::cerr<< "parsing SV VCF ... ";
        vcfParser.setParseSVFile(true);
        vcfParser.parsingVCF(params.svFile, vcfSet[Genome::NORMAL], *chrMultiVariants);
        vcfParser.reset();
        std::cerr<< difftime(time(NULL), begin) << "s\n";    
    }
    
    // load MOD vcf file
    if(params.modFile!=""){
        begin = time(NULL);
        std::cerr<< "parsing MOD VCF ... ";
        vcfParser.setParseMODFile(true);
        vcfParser.parsingVCF(params.modFile, vcfSet[Genome::NORMAL], *chrMultiVariants);
        vcfParser.reset();
        std::cerr<< difftime(time(NULL), begin) << "s\n";    
    }
}

void HaplotagProcess::setChrVecAndChrLength(){
    chrVec = &(vcfSet[Genome::NORMAL].chrVec);
    chrLength = &(vcfSet[Genome::NORMAL].chrLength); 
}

void HaplotagProcess::setProcessingChromRegion(){
    if (!params.bamCfg.region.empty()) {//if the region is not empty, set the chromosome vector to the region
        auto colonPos = params.bamCfg.region.find(":");//find the colon position
        std::string regionChr;//the chromosome name
        if (colonPos != std::string::npos) {
            regionChr = params.bamCfg.region.substr(0, colonPos);
        }
        else {
            regionChr = params.bamCfg.region;
        }
        auto chrVecIter = std::find((*chrVec).begin(), (*chrVec).end(), regionChr);
        if (chrVecIter != (*chrVec).end()) {
            (*chrVec) = std::vector<std::string>{regionChr};
        } else {
            std::cerr << "[ERROR] Incorrect chromosome for input region: " << regionChr << std::endl;
            exit(1);
        }
    }

    // remove variant on chromosome not in chrVec
    std::map<std::string, std::map<int, MultiGenomeVar>>::iterator iter = chrMultiVariants->begin();
    while(iter != chrMultiVariants->end()){
        if(std::find((*chrVec).begin(), (*chrVec).end(), iter->first) == (*chrVec).end()){
            chrMultiVariants->erase(iter++);
        }else{
            iter++;
        }
    }
}


void HaplotagProcess::tagRead(HaplotagParameters &params, std::string& tagBamFile, const Genome& geneSample){

    // tag read
    std::time_t begin = time(NULL);
    std::cerr<< getTagReadStartMessage();

    ParsingBamConfig& config = params.bamCfg;

    ParsingBamControl control;
    control.mode = ParsingBamMode::SINGLE_THREAD;
    control.writeOutputBam = true;
    control.mappingQualityFilter = true;

    BamParserContext ctx(tagBamFile, params.fastaFile, *chrVec, *chrLength, *chrMultiVariants, vcfSet, geneSample);
    GermlineHaplotagBamParser* haplotagBamParser = createHaplotagBamParser(config, control, readStats);
    haplotagBamParser->createTagLog();
    haplotagBamParser->parsingBam(ctx);
    delete haplotagBamParser;

    std::cerr<< "tag read " << difftime(time(NULL), begin) << "s\n";
    return;
}

void HaplotagProcess::printExecutionReport(){
    std::cerr<< "-------------------------------------------\n";
    std::cerr<< "total process time        : " << difftime(time(NULL), processBegin) << "s\n";
    std::cerr<< "total alignment           : " << readStats.totalAlignment     << "\n";
    std::cerr<< "total supplementary       : " << readStats.totalSupplementary << "\n";
    std::cerr<< "total secondary           : " << readStats.totalSecondary     << "\n";
    std::cerr<< "total unmapped            : " << readStats.totalUnmapped      << "\n";
    std::cerr<< "total tagged alignments   : " << readStats.totalTagCount     << "\n";
    std::cerr<< "    L----total HP1        : " << readStats.totalHpCount[ReadHP::H1]     << "\n";   
    std::cerr<< "    L----total HP2        : " << readStats.totalHpCount[ReadHP::H2]     << "\n";   
    std::cerr<< "    L----total HP1-1      : " << readStats.totalHpCount[ReadHP::H1_1]   << "\n";   
    std::cerr<< "    L----total HP2-1      : " << readStats.totalHpCount[ReadHP::H2_1]   << "\n";   
    std::cerr<< "    L----total HP3        : " << readStats.totalHpCount[ReadHP::H3]     << "\n";   
    std::cerr<< "         L----only H3 SNP : " << readStats.totalreadOnlyH3Snp << "\n";  
    std::cerr<< "total untagged            : " << readStats.totalUnTagCount   << "\n";
    std::cerr<< "    L----lower mapping quality        : " << readStats.totalLowerQuality   << "\n"; 
    std::cerr<< "    L----no variant                   : " << readStats.totalEmptyVariant   << "\n";
    std::cerr<< "    L----start pos > last variant pos : " << readStats.totalOtherCase   << "\n";  
    std::cerr<< "    L----judge to untag               : " << readStats.totalHpCount[ReadHP::unTag] << "\n"; 
    std::cerr<< "         L----high similarity         : " << readStats.totalHighSimilarity   << "\n";   
    std::cerr<< "         L----cross two block         : " << readStats.totalCrossTwoBlock   << "\n";   
    std::cerr<< "         L----no variant judge HP     : " << readStats.totalWithOutVaraint   << "\n";   
    std::cerr<< "-------------------------------------------\n";
}

GermlineTagLog::GermlineTagLog(const HaplotagParameters& params) 
    : HaplotagReadLog<HaplotagParameters, TagReadLog>(params, params.bamCfg.resultPrefix+".out"){
}

GermlineTagLog::~GermlineTagLog(){}

void GermlineTagLog::addParamsMessage(){
    *tagReadLog << "##snpFile:" << params.snpFile << "\n"
                << "##svFile:" << params.svFile << "\n"
                << "##bamFile:" << params.bamFile << "\n"
                << "##resultPrefix:" << params.bamCfg.resultPrefix << "\n"
                << "##numThreads:" << params.bamCfg.numThreads << "\n"
                << "##region:" << params.bamCfg.region << "\n"
                << "##qualityThreshold:" << params.bamCfg.qualityThreshold << "\n"
                << "##percentageThreshold:" << params.bamCfg.percentageThreshold << "\n"
                << "##tagSupplementary:" << params.bamCfg.tagSupplementary << "\n";
}

void GermlineTagLog::writeBasicColumns(){
    *tagReadLog << "#ReadID\t"
                << "CHROM\t"
                << "ReadStart\t"
                << "Confidnet(%)\t"
                << "Haplotype\t"
                << "PhaseSet\t"
                << "TotalAllele\t"
                << "HP1Allele\t"
                << "HP2Allele\t"
                << "phasingQuality(PQ)\t"
                << "(Variant,HP)\t"
                << "(PhaseSet,Variantcount)\n";
}

void GermlineTagLog::writeTagReadLog(TagReadLog& data){

    // write tag log file
    (*tagReadLog) << bam_get_qname(&data.aln) << "\t"
                << data.bamHdr.target_name[data.aln.core.tid] << "\t"
                << data.aln.core.pos << "\t"
                << data.norHPsimilarity << "\t"
                << "H" << data.hpResultStr << "\t"
                << data.psResultStr << "\t"
                << data.hpCount[SnpHP::GERMLINE_H1] + data.hpCount[SnpHP::GERMLINE_H2] << "\t"
                << data.hpCount[SnpHP::GERMLINE_H1] << "\t"
                << data.hpCount[SnpHP::GERMLINE_H2] << "\t"
                << data.pqValue << "\t";

    // print position and HP
    for (const auto& v : data.variantsHP) {
        (*tagReadLog) << " " << v.first << "," << v.second;
    }

    (*tagReadLog) << "\t";

    // belong PS, number of variant
    for (const auto& v : data.norCountPS) {
        (*tagReadLog) << " " << v.first << "," << v.second;
    }

    (*tagReadLog) << "\n"; 
}


GermlineHaplotagBamParser::GermlineHaplotagBamParser(
    const ParsingBamConfig &config,
    const ParsingBamControl &control,
    ReadStatistics& readStats,
    const HaplotagParameters& params
):HaplotagBamParser(config, control),
    params(params), readStats(readStats), tagResult(nullptr)
{}

void GermlineHaplotagBamParser::createTagLog(){
    // Read log can only be written in single thread mode, which is enforced when writeOutputBam is true.
    // This is because log writing requires sequential processing of reads.
    if(params.bamCfg.writeReadLog && control.writeOutputBam){
        tagResult = createTagReadLog();
        tagResult->writeHeader();
    }
}

GermlineHaplotagBamParser::~GermlineHaplotagBamParser(){
    if(tagResult){
        delete tagResult;
    }
}



GermlineHaplotagChrProcessor::GermlineHaplotagChrProcessor(
    bool writeOutputBam,
    bool mappingQualityFilter,
    ReadStatistics& readStats,
    GermlineTagLog *tagResult
):ChromosomeProcessor(writeOutputBam, mappingQualityFilter),
    readStats(readStats),
    tagResult(tagResult) 
{
    localReadStats = ReadStatistics();
}

GermlineHaplotagChrProcessor::~GermlineHaplotagChrProcessor(){

}

void GermlineHaplotagChrProcessor::processLowMappingQuality(){
    localReadStats.totalLowerQuality++;
    localReadStats.totalUnTagCount++;
    localReadStats.totalAlignment++;
}

void GermlineHaplotagChrProcessor::processUnmappedRead(){
    localReadStats.totalUnmapped++;
    localReadStats.totalUnTagCount++;
    localReadStats.totalAlignment++;
}

void GermlineHaplotagChrProcessor::processSecondaryAlignment(){
    localReadStats.totalSecondary++;
    localReadStats.totalUnTagCount++;
    localReadStats.totalAlignment++;
}

void GermlineHaplotagChrProcessor::processSupplementaryAlignment(){
    localReadStats.totalSupplementary++;
    localReadStats.totalUnTagCount++;
    localReadStats.totalAlignment++;
}

void GermlineHaplotagChrProcessor::processEmptyVariants(){
    localReadStats.totalEmptyVariant++;
    localReadStats.totalUnTagCount++;
    localReadStats.totalAlignment++;
}

void GermlineHaplotagChrProcessor::processOtherCase(){
    localReadStats.totalOtherCase++;
    localReadStats.totalUnTagCount++;
    localReadStats.totalAlignment++;
}

void GermlineHaplotagChrProcessor::processRead(
    bam1_t &aln, 
    const bam_hdr_t &bamHdr,
    const std::string &ref_string,
    std::map<int, MultiGenomeVar> &currentVariants,
    std::map<int, MultiGenomeVar>::iterator &firstVariantIter,
    ChrProcContext& ctx
){

    if( (aln.core.flag & 0x800) != 0 ){
        localReadStats.totalSupplementary++;
    }

    int pqValue = 0;
    int psValue = 0; 
    int haplotype = ReadHP::unTag;

    haplotype = judgeHaplotype(bamHdr, aln, ctx.chrName, ctx.params.percentageThreshold, tagResult, pqValue, psValue, ctx.genomeSample, ref_string, ctx.params, firstVariantIter, currentVariants, ctx.vcfSet);

    initFlag(&aln, "HP");
    initFlag(&aln, "PS");
    initFlag(&aln, "PQ");


    if (haplotype != ReadHP::unTag){

        localReadStats.totalHpCount[haplotype]++;
        localReadStats.totalTagCount++;
        addAuxiliaryTags(&aln, haplotype, pqValue, psValue);
        
    }
    else{
        localReadStats.totalHpCount[ReadHP::unTag]++;
        localReadStats.totalUnTagCount++;
    }
    localReadStats.totalAlignment++;
}


void GermlineHaplotagChrProcessor::addAuxiliaryTags(bam1_t *aln, int& haplotype, int& pqValue, int& psValue){
    bam_aux_append(aln, "HP", 'i', sizeof(haplotype), (uint8_t*) &haplotype);
    bam_aux_append(aln, "PS", 'i', sizeof(psValue), (uint8_t*) &psValue);
    bam_aux_append(aln, "PQ", 'i', sizeof(pqValue), (uint8_t*) &pqValue);
}

int GermlineHaplotagChrProcessor::judgeHaplotype(
    const bam_hdr_t &bamHdr,
    const bam1_t &aln,
    std::string chrName,
    double percentageThreshold,
    GermlineTagLog *tagResult,
    int &pqValue,
    int &psValue,
    const int tagSample,
    const std::string &ref_string,
    const ParsingBamConfig &params,
    std::map<int, MultiGenomeVar>::iterator &firstVariantIter,
    std::map<int, MultiGenomeVar> &currentChrVariants,
    std::map<Genome, VCF_Info> &vcfSet
){

    std::map<int, int> hpCount;
    hpCount[SnpHP::GERMLINE_H1] = 0;
    hpCount[SnpHP::GERMLINE_H2] = 0;

    //record variants on this read
    std::map<int, int> variantsHP;

    std::map<int,int> countPS;

    // position relative to reference
    int ref_pos = aln.core.pos;
    // position relative to read
    int query_pos = 0;
    
    /// Create a CIGAR parser using polymorphism design
    // Use GermlineHaplotagCigarParser to process CIGAR strings for normal samples 
    CigarParserContext cigarCtx(aln, bamHdr, chrName, params, firstVariantIter, currentChrVariants, ref_string);   
    CigarParser* cigarParser = new GermlineHaplotagCigarParser(cigarCtx, ref_pos, query_pos);
    cigarParser->parsingCigar(hpCount, variantsHP, countPS);
    delete cigarParser;

    // get the number of SVs occurring on different haplotypes in a read
    judger.judgeSVHap(aln, vcfSet, hpCount, tagSample);

    double min,max;
    int hpResult = ReadHP::unTag;

    // determine the haplotype of the read
    hpResult = judger.judgeReadHap(hpCount, min, max, percentageThreshold, pqValue, psValue, countPS, &localReadStats.totalHighSimilarity, &localReadStats.totalWithOutVaraint);
     
    //write tag log file
    if(tagResult != nullptr){
        std::string hpResultStr = ((hpResult == ReadHP::unTag) ? "." : std::to_string(hpResult));
        std::string psResultStr = ".";

        if (hpResultStr != ".") {
            auto psIter = countPS.begin();
            psResultStr = std::to_string(psIter->first);
        }

        double norHPsimilarity = (max / (max + min));

        TagReadLog data{
            .aln = aln,
            .bamHdr = bamHdr,
            .norHPsimilarity = norHPsimilarity,
            .hpResultStr = hpResultStr,
            .psResultStr = psResultStr,
            .hpCount = hpCount,
            .pqValue = pqValue,
            .variantsHP = variantsHP,
            .norCountPS = countPS,
            .tumCountPS = nullptr,
            .deriveByHpSimilarity = 0.0
        };
        tagResult->writeTagReadLog(data);
    }

    return hpResult;
}

void GermlineHaplotagChrProcessor::initFlag(bam1_t *aln, std::string flag){

    uint8_t *hpTag = bam_aux_get(aln, flag.c_str() );

    if( hpTag != NULL )
        bam_aux_del(aln, hpTag);

    return;
}

void GermlineHaplotagChrProcessor::postProcess(const std::string &chr, std::map<int, MultiGenomeVar> &currentVariants){
    #pragma omp critical
    {
        readStats.totalAlignment += localReadStats.totalAlignment;
        readStats.totalHpCount[ReadHP::unTag] += localReadStats.totalHpCount[ReadHP::unTag];
        readStats.totalHpCount[ReadHP::H1] += localReadStats.totalHpCount[ReadHP::H1];
        readStats.totalHpCount[ReadHP::H2] += localReadStats.totalHpCount[ReadHP::H2];
        readStats.totalHpCount[ReadHP::H3] += localReadStats.totalHpCount[ReadHP::H3];
        readStats.totalHpCount[ReadHP::H1_1] += localReadStats.totalHpCount[ReadHP::H1_1];
        readStats.totalHpCount[ReadHP::H2_1] += localReadStats.totalHpCount[ReadHP::H2_1];
        readStats.totalSupplementary += localReadStats.totalSupplementary;
        readStats.totalSecondary += localReadStats.totalSecondary;
        readStats.totalUnmapped += localReadStats.totalUnmapped;
        readStats.totalTagCount += localReadStats.totalTagCount;
        readStats.totalUnTagCount += localReadStats.totalUnTagCount;
        readStats.totalLowerQuality += localReadStats.totalLowerQuality;
        readStats.totalOtherCase += localReadStats.totalOtherCase;
        readStats.totalunTag_HP0 += localReadStats.totalunTag_HP0;
        readStats.totalreadOnlyH3Snp += localReadStats.totalreadOnlyH3Snp;
        readStats.totalHighSimilarity += localReadStats.totalHighSimilarity;
        readStats.totalCrossTwoBlock += localReadStats.totalCrossTwoBlock;
        readStats.totalEmptyVariant += localReadStats.totalEmptyVariant;
        readStats.totalWithOutVaraint += localReadStats.totalWithOutVaraint;
    }
}

GermlineHaplotagCigarParser::GermlineHaplotagCigarParser(CigarParserContext& ctx, int& ref_pos, int& query_pos)
:CigarParser(ctx, ref_pos, query_pos)
{

}

GermlineHaplotagCigarParser::~GermlineHaplotagCigarParser(){

}

void GermlineHaplotagCigarParser::processMatchOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, std::string& base, bool& isAlt, int& offset){
    auto norVar = (*currentVariantIter).second.Variant[NORMAL];//get the current normal variant
    judger.judgeSnpHap(ctx.chrName, norVar, base, ref_pos, length, i, aln_core_n_cigar, cigar, currentVariantIter, *hpCount, *variantsHP, *norCountPS, isAlt);
    //judge the haplotype of the current normal variant
}

void GermlineHaplotagCigarParser::processDeletionOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, bool& alreadyJudgeDel){
    // only execute at the first phased normal snp
    if ((*currentVariantIter).second.isExists(NORMAL) && !alreadyJudgeDel){
        if((*currentVariantIter).second.Variant[NORMAL].GT == GenomeType::PHASED_HETERO){
            // longphase v1.73 only execute once
            alreadyJudgeDel = true;
            judger.judgeDeletionHap(ctx.chrName, ctx.ref_string, ref_pos, length, query_pos, currentVariantIter, &ctx.aln, *hpCount, *variantsHP, *norCountPS);
        }
    }
}



