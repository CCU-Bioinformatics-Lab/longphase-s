#include "SomaticHaplotagProcess.h"


SomaticHaplotagProcess::SomaticHaplotagProcess(SomaticHaplotagParameters &params)
    : HaplotagProcess(params.basic),sParams(params), somaticBenchmark(params.benchmarkVcf, params.benchmarkBedFile, params.basic.bamCfg.qualityThreshold)
{
    hpBeforeInheritance = new ReadHpDistriLog();
    hpAfterInheritance = new ReadHpDistriLog();
}

SomaticHaplotagProcess::~SomaticHaplotagProcess(){
    delete hpBeforeInheritance;
    delete hpAfterInheritance;
}

void SomaticHaplotagProcess::printParamsMessage(){
    std::cerr<< "LongPhase-S v" << sParams.basic.bamCfg.version << " - Somatic Haplotag\n";
    std::cerr<< "\n";
    std::cerr<< "[Input Files]\n";
    std::cerr<< "phased normal SNP file       : " << sParams.basic.snpFile << "\n";
    std::cerr<< "tumor SNP file               : " << sParams.tumorSnpFile << "\n";
    std::cerr<< "normal BAM file              : " << sParams.basic.bamFile << "\n";
    std::cerr<< "tumor BAM file               : " << sParams.tumorBamFile << "\n";
    std::cerr<< "reference file               : " << sParams.basic.fastaFile << "\n";
    std::cerr<< "\n";
    std::cerr<< "[Benchmark Files]\n";
    std::cerr<< "truth VCF file               : " << sParams.benchmarkVcf << "\n";
    std::cerr<< "truth BED file               : " << sParams.benchmarkBedFile << "\n";
    std::cerr<< "\n";
    std::cerr<< "[Output Files]\n";
    std::cerr<< "tagged tumor BAM file        : " << sParams.basic.bamCfg.resultPrefix + "." + sParams.basic.bamCfg.outputFormat << "\n";
    std::cerr<< "purity estimation file       : " << (sParams.callerCfg.estimateTumorPurity ? sParams.basic.bamCfg.resultPrefix + "_purity.out" : "") << "\n";
    std::cerr<< "somatic calling VCF file     : " << (sParams.writeSomaticCallingVcf ? sParams.basic.bamCfg.resultPrefix + "_sc.vcf" : "") << "\n";
    std::cerr<< "benchmark metrics file       : " << (sParams.benchmarkVcf != "" ? sParams.basic.bamCfg.resultPrefix + sParams.metricsSuffix : "") << "\n";
    std::cerr<< "log file                     : " << (sParams.basic.bamCfg.writeReadLog ? (sParams.basic.bamCfg.resultPrefix+".out") : "") << "\n";
    std::cerr<< "-------------------------------------------\n";
    std::cerr<< "[Somatic Haplotagging Params] " << "\n";
    std::cerr<< "number of threads            : " << sParams.basic.bamCfg.numThreads << "\n";
    std::cerr<< "tag region                   : " << (!sParams.basic.bamCfg.region.empty() ? sParams.basic.bamCfg.region : "all") << "\n";
    std::cerr<< "filter mapping quality below : " << sParams.basic.bamCfg.qualityThreshold << "\n";
    std::cerr<< "percentage threshold         : " << sParams.basic.bamCfg.percentageThreshold << "\n";
    std::cerr<< "write log file               : " << (sParams.basic.bamCfg.writeReadLog ? "enabled" : "disabled") << "\n";
    std::cerr<< "tag supplementary            : " << (sParams.basic.bamCfg.tagSupplementary ? "enabled" : "disabled") << "\n";
    std::cerr<< "\n";
    std::cerr<< "[Somatic Variant Calling Params] " << "\n";
    std::cerr<< "mapping quality              : " << sParams.basic.bamCfg.qualityThreshold << "\n";
    std::cerr<< "tumor purity value           : " << (sParams.callerCfg.estimateTumorPurity ? "automatic estimation" : std::to_string(sParams.callerCfg.tumorPurity)) << "\n";
    std::cerr<< "variant filtering            : " << (sParams.callerCfg.enableFilter ? "enabled" : "disabled") << "\n";
    std::cerr<< "write calling VCF            : " << (sParams.writeSomaticCallingVcf ? "enabled" : "disabled") << "\n";
    std::cerr<< "write calling log            : " << (sParams.callerCfg.writeCallingLog ? "enabled" : "disabled") << "\n";
    std::cerr<< "-------------------------------------------\n";
}

void SomaticHaplotagProcess::pipelineProcess()
{
    printParamsMessage();
    // decide on the type of tagging for VCF and BAM files
    Genome tagSample = TUMOR;

    VcfParser vcfParser(tagSample);
    // load SNP, SV, MOD vcf file
    parseVariantFiles(vcfParser);
    //decide which genome sample chrVec and chrLength belong to
    setChrVecAndChrLength();
    // [debug] calculate SNP counts
    displaySnpCounts();
    // update chromosome processing based on region
    setProcessingChromRegion();

    //somatic variant calling
    CallerContext ctx(sParams.basic.bamFile, sParams.tumorBamFile, 
                    sParams.basic.snpFile, sParams.tumorSnpFile, 
                    sParams.basic.fastaFile);

    SomaticVarCaller *somaticVarCaller = new SomaticVarCaller(sParams.callerCfg, sParams.basic.bamCfg, *chrVec);
    somaticVarCaller->variantCalling(ctx, *chrVec, *chrLength, (*chrMultiVariants), vcfSet);
    somaticVarCaller->getSomaticFlag(*chrVec, *chrMultiVariants);
    delete somaticVarCaller;

    // write somatic calling vcf
    if(sParams.writeSomaticCallingVcf){
        std::cerr<< "writing somatic variants to vcf file ... ";
        std::time_t begin = time(NULL);
        vcfParser.setCommandLine(sParams.basic.bamCfg.command);
        vcfParser.setVersion(sParams.basic.bamCfg.version);
        vcfParser.writingResultVCF(sParams.tumorSnpFile, vcfSet[Genome::TUMOR], *chrMultiVariants, sParams.basic.bamCfg.resultPrefix);
        std::cerr<< difftime(time(NULL), begin) << "s\n";
    }

    // remove tumor & benchmark variants outside bed regions
    if(somaticBenchmark.isLoadBedFile() && somaticBenchmark.isEnabled()){
        std::cerr << "[Benchmark] removing tumor & truth somatic variants outside bed regions ... ";
        std::time_t begin = time(NULL);
        somaticBenchmark.removeVariantsOutBedRegion(*chrMultiVariants);
        std::cerr<< difftime(time(NULL), begin) << "s\n";
        // [debug] display variants in bed region count
        // somaticBenchmark.displayBedRegionCount(*chrVec);
    }

    // return;
    // tag read
    tagRead(sParams.basic, sParams.tumorBamFile, tagSample);

    // postprocess after haplotag
    postprocessForHaplotag();

    printExecutionReport();

    return;
};

void SomaticHaplotagProcess::parseVariantFiles(VcfParser& vcfParser){
    // normal SNP, SV, MOD vcf file
    HaplotagProcess::parseVariantFiles(vcfParser);

    //load tumor snp vcf
    if(sParams.tumorSnpFile != ""){
        std::time_t begin = time(NULL);
        std::cerr<< "parsing tumor SNP VCF ... ";
        vcfParser.setParseSnpFile(true);
        vcfParser.parsingVCF(sParams.tumorSnpFile, vcfSet[Genome::TUMOR], *chrMultiVariants);
        vcfParser.reset();
        std::cerr<< difftime(time(NULL), begin) << "s\n";
    }

    //parse benchmark vcf for benchmarking
    if(sParams.benchmarkVcf != ""){
        std::time_t begin = time(NULL);
        std::cerr<< "[Benchmark] parsing truth VCF ... ";
        somaticBenchmark.setEnabled(true);
        somaticBenchmark.loadTruthSomaticVCF(sParams.benchmarkVcf, vcfSet[Genome::TRUTH_SOMATIC], *chrMultiVariants);
        std::cerr<< difftime(time(NULL), begin) << "s\n";

        // [debug] display somatic variant count
        // somaticBenchmark.displaySomaticVarCount(vcfSet[Genome::TRUTH_SOMATIC].chrVec, *chrMultiVariants);
    }

    // parse benchmark bed file
    if(sParams.benchmarkBedFile != "" && somaticBenchmark.isEnabled()){
        std::time_t begin = time(NULL);
        std::cerr<< "[Benchmark] parsing truth BED file ... ";
        somaticBenchmark.parseBedFile(sParams.benchmarkBedFile);
        std::cerr<< difftime(time(NULL), begin) << "s\n";

        // mark variants in bed regions
        somaticBenchmark.markVariantsInBedRegions(vcfSet[TUMOR].chrVec, *chrMultiVariants);
    }
}

void SomaticHaplotagProcess::setChrVecAndChrLength(){
    //check normal & tumor chr & length
    try{
        for(auto& chrIter : vcfSet[Genome::TUMOR].chrLength){
            if(vcfSet[Genome::NORMAL].chrLength.find(chrIter.first) != vcfSet[Genome::NORMAL].chrLength.end()){
                if(chrIter.second != vcfSet[Genome::NORMAL].chrLength[chrIter.first]){
                    throw std::runtime_error("tumor & normal VCFs chromosome length are not the same => chr: " + chrIter.first + " length: " + std::to_string(chrIter.second));
                }
            }else{
                throw std::runtime_error("tumor & normal VCFs chromosome count are not the same");
            }
        }

        // if tumor VCFs chromosome count are empty, use normal VCFs chromosome count
        if(vcfSet[Genome::TUMOR].chrVec.empty()){
            std::cerr<< "[WARNING] tumor VCF chromosome count is empty" << std::endl;
            if(vcfSet[Genome::NORMAL].chrVec.empty()){
                throw std::runtime_error("tumor & normal VCFs chromosome count are empty");
            }else{
                std::cerr<< "[INFO] use normal VCF chromosome count" << std::endl;
                chrVec = &(vcfSet[Genome::NORMAL].chrVec);
                chrLength = &(vcfSet[Genome::NORMAL].chrLength); 
            }
        }else{
            chrVec = &(vcfSet[Genome::TUMOR].chrVec);
            chrLength = &(vcfSet[Genome::TUMOR].chrLength); 
        }

        // if tumor VCFs chromosome length are empty, use normal VCFs chromosome length
        if(vcfSet[Genome::TUMOR].chrLength.empty()){
            std::cerr<< "[WARNING] tumor VCF chromosome length is empty" << std::endl;
            if(vcfSet[Genome::NORMAL].chrLength.empty()){
                throw std::runtime_error("tumor & normal VCFs chromosome length are empty");
            }else{
                std::cerr<< "[INFO] use normal VCF chromosome length" << std::endl;
                chrLength = &(vcfSet[Genome::NORMAL].chrLength); 
            }
        }else{
            chrLength = &(vcfSet[Genome::TUMOR].chrLength); 
        }

    }catch(const std::exception& e){
        std::cerr << "[ERROR] (setChrVecAndChrLength) :" << e.what() << std::endl;
        return;
    }
}

void SomaticHaplotagProcess::displaySnpCounts(){
    int tumor_snp_count = 0;
    int normal_snp_count = 0;
    int overlap_snp_count = 0;
    int tumor_insert_count = 0;
    int tumor_delete_count = 0;
    for(auto& chrIter : (*chrVec)){
        auto chrVarIter = (*chrMultiVariants)[chrIter].begin();
        while(chrVarIter != (*chrMultiVariants)[chrIter].end()){
            if((*chrVarIter).second.isExists(TUMOR)){
                if((*chrVarIter).second.Variant[TUMOR].variantType == VariantType::SNP){
                    tumor_snp_count++;
                }
                if((*chrVarIter).second.Variant[TUMOR].variantType == VariantType::INSERTION){
                    tumor_insert_count++;
                }
                if((*chrVarIter).second.Variant[TUMOR].variantType == VariantType::DELETION){
                    tumor_delete_count++;
                }
            }
            if((*chrVarIter).second.isExists(NORMAL)){
                normal_snp_count++;
            }
            if((*chrVarIter).second.isExists(TUMOR) && (*chrVarIter).second.isExists(NORMAL)){
                overlap_snp_count++;
            }
            chrVarIter++;
        }
    }
    std::cerr << "Normal SNP count: " << normal_snp_count << std::endl;
    std::cerr << "Tumor SNP count: " << tumor_snp_count << std::endl;
    std::cerr << "Overlap SNP count: " << overlap_snp_count << std::endl;
    std::cerr << "Tumor Insert count: " << tumor_insert_count << std::endl;
    std::cerr << "Tumor Delete count: " << tumor_delete_count << std::endl;
}

void SomaticHaplotagProcess::postprocessForHaplotag(){
    somaticBenchmark.writeTaggedSomaticReadReport(*chrVec, sParams.basic.bamCfg.resultPrefix + sParams.metricsSuffix);

    // other log files
    if(sParams.callerCfg.writeCallingLog){
        hpBeforeInheritance->writeReadHpDistriLog(sParams.basic.bamCfg.resultPrefix + "_read_distri_before_inheritance.out", *chrVec);
        hpAfterInheritance->writeReadHpDistriLog(sParams.basic.bamCfg.resultPrefix + "_read_distri_after_inheritance.out", *chrVec);
        //write snp cover region
        hpAfterInheritance->writePosCoverRegionLog(sParams.basic.bamCfg.resultPrefix + "_snp_cover_region.out", *chrVec);
        //write read cover region in whole genome
        hpAfterInheritance->writeTagReadCoverRegionLog(sParams.basic.bamCfg.resultPrefix + "_read_cover_region.bed", *chrVec, *chrLength);
    }

    if(sParams.writeBenchmarkLog){
        //write somatic read log
        somaticBenchmark.writeTotalTruthSomaticReadReport(*chrVec, sParams.basic.bamCfg.resultPrefix + "_total_truth_somatic_read.out");
        somaticBenchmark.writeTaggedReadReport(*chrVec, sParams.basic.bamCfg.resultPrefix + "_total_tagged_read.out");
        somaticBenchmark.writePosAlleleCountLog(*chrVec, sParams.basic.bamCfg.resultPrefix + "_allele_count.out", *chrMultiVariants);
        somaticBenchmark.writeBedRegionLog(*chrVec, *chrMultiVariants, sParams.basic.bamCfg.resultPrefix);
    }
}

SomaticHaplotagBamParser::SomaticHaplotagBamParser(
    const ParsingBamConfig &config,
    const ParsingBamControl &control,
    ReadStatistics& readStats,
    SomaticReadBenchmark& somaticBenchmark,
    ReadHpDistriLog& hpBeforeInheritance,
    ReadHpDistriLog& hpAfterInheritance,
    const SomaticHaplotagParameters& sParams
):GermlineHaplotagBamParser(config, control, readStats, sParams.basic),
    somaticBenchmark(somaticBenchmark),
    hpBeforeInheritance(hpBeforeInheritance),
    hpAfterInheritance(hpAfterInheritance),
    sParams(sParams)
{

}

SomaticHaplotagBamParser::~SomaticHaplotagBamParser(){

}

SomaticHaplotagChrProcessor::SomaticHaplotagChrProcessor(
    bool writeOutputBam,
    bool mappingQualityFilter,
    ReadStatistics& readStats,
    GermlineTagLog *tagResult,
    SomaticReadBenchmark& somaticBenchmark,
    ReadHpDistriLog& hpBeforeInheritance,
    ReadHpDistriLog& hpAfterInheritance,
    const std::string &chr
):GermlineHaplotagChrProcessor(writeOutputBam, mappingQualityFilter, readStats, tagResult),
    somaticReadCounter(nullptr)
{

    bool enableBenckmark = somaticBenchmark.isEnabled();
    #pragma omp critical
    {
        hpBeforeInheritance.loadChrKey(chr);
        hpAfterInheritance.loadChrKey(chr);
        somaticBenchmark.loadChrKey(chr);
        localHpBeforeInheritance = hpBeforeInheritance.getChrHpResultsPtr(chr);
        localHpAfterInheritance = hpAfterInheritance.getChrHpResultsPtr(chr);
        somaticReadCounter = new SomaticReadVerifier(enableBenckmark, somaticBenchmark.getMetricsPtr(chr));
    }
    this->chr = chr;
    begin = time(NULL);
}

SomaticHaplotagChrProcessor::~SomaticHaplotagChrProcessor(){
    delete somaticReadCounter;
    // #pragma omp critical
    // {
    //     std::cerr << chr << " : " << difftime(time(NULL), begin) << "s\n";
    // }
}

int SomaticHaplotagChrProcessor::judgeHaplotype(
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
    hpCount[1] = 0;
    hpCount[2] = 0;
    hpCount[3] = 0;
    hpCount[4] = 0;

    //record variants on this read
    std::map<int, int> variantsHP;

    // poition, <BaseHP,deriveByHp>
    std::map<int, std::pair<int , int>> somaticVarDeriveHP;

    //record PS count(vcf type, PS value, count)
    std::map<int, int> tumCountPS;
    std::map<int, int> norCountPS;

    // position relative to reference
    int ref_pos = aln.core.pos;
    // position relative to read
    int query_pos = 0;

    // Create a CIGAR parser using polymorphism design
    // Use GermlineHaplotagCigarParser to process CIGAR strings for normal samples
    CigarParserContext cigarCtx(aln, bamHdr, chrName, params, firstVariantIter, currentChrVariants, ref_string);   
    CigarParser* cigarParser = new SomaticHaplotagCigarParser(cigarCtx, ref_pos, query_pos, tumCountPS, somaticVarDeriveHP, *somaticReadCounter);
    cigarParser->parsingCigar(hpCount, variantsHP, norCountPS);
    delete cigarParser;
    //In the current version, only normal SVs are considered, without inclusion of tumor samples
    judger.judgeSVHap(aln, vcfSet, hpCount, Genome::NORMAL);

    int startPos = aln.core.pos + 1;
    int endPos = ref_pos;

    //the similarity of HP types
    //tumor variants
    double norHPsimilarity = 0.0;
    double tumHPsimilarity = 0.0;

    // determine the haplotype of the read
    int hpResult = ReadHP::unTag;
    hpResult = somaticJudger.judgeSomaticReadHap(hpCount, pqValue, norCountPS, norHPsimilarity, tumHPsimilarity, percentageThreshold, &localReadStats.totalHighSimilarity, &localReadStats.totalCrossTwoBlock, &localReadStats.totalWithOutVaraint);

    int hpResultBeforeInheritance = hpResult;

    float deriveByHpSimilarity = 0.0;

    if(hpResult == ReadHP::H3){
        // Inherit haplotype information for H3 reads based on somatic variant derived haplotype
        hpResult = inheritHaplotype(deriveByHpSimilarity, params.percentageThreshold, somaticVarDeriveHP, hpCount, hpResult, aln);
    }
 

    if(!somaticVarDeriveHP.empty()){
        for(auto somaticVarIter : somaticVarDeriveHP){
            int pos = somaticVarIter.first;
            int baseHP = somaticVarIter.second.first;
            int deriveHP = somaticVarIter.second.second;
            //Record read HP result before correction hp result
            localHpBeforeInheritance->recordReadHp(pos, hpResultBeforeInheritance, baseHP);
            localHpBeforeInheritance->recordDeriveHp(pos, deriveHP, 0.0);

            //Record read HP result after correction hp result
            localHpAfterInheritance->recordReadHp(pos, hpResult, baseHP);
            localHpAfterInheritance->recordDeriveHp(pos, deriveHP, deriveByHpSimilarity);

            if(hpResult != ReadHP::unTag){   

                // update cover region at somatic position
                localHpAfterInheritance->recordAlignCoverRegion(pos, startPos, endPos);
            }
        }
    }

    if(hpCount[1] == 0 && hpCount[2] == 0 && hpCount[3] != 0){
        if(hpResult == ReadHP::H3){
            localReadStats.totalreadOnlyH3Snp++;
        }
    }

    std::string hpResultStr = ".";
    std::string psResultStr = ".";

    hpResultStr = ReadHapUtil::readHapIntToString(hpResult);

    //determine the PS value of the read
    if( hpResult != ReadHP::unTag ){
        std::map<int, int>::iterator psIter;
        if(hpResult != ReadHP::H1 && hpResult != ReadHP::H2){

            if(norCountPS.size() != 0){
                psIter = norCountPS.begin();
                psResultStr = std::to_string((*psIter).first);
                psValue = (*psIter).first;
            }
            // the read only had unphased hetero tumor SNPs or homo SNPs
            else{
                psResultStr= "*";
                //for determining not write PS in the current read
                psValue = VarData::NONE_PHASED_SET;  
            }
        }else if(hpResult == ReadHP::H1 || hpResult == ReadHP::H2){
            psIter = norCountPS.begin();
            psResultStr = std::to_string((*psIter).first);
            psValue = (*psIter).first;
        }
    }

    std::string readID = bam_get_qname(&aln);
    
    // //record somatic read
    somaticReadCounter->recordTaggedRead(chrName, readID, hpResult, variantsHP, hpCount, norHPsimilarity, deriveByHpSimilarity, currentChrVariants);
    somaticReadCounter->recordCrossingTruthSomaticSnpRead(chrName, readID, hpResult, variantsHP, hpCount, norHPsimilarity, deriveByHpSimilarity, currentChrVariants);

    //write tag log file
    if(tagResult!=NULL){

        TagReadLog data{
            .aln = aln,
            .bamHdr = bamHdr,
            .norHPsimilarity = norHPsimilarity,
            .hpResultStr = hpResultStr,
            .psResultStr = psResultStr,
            .hpCount = hpCount,
            .pqValue = pqValue,
            .variantsHP = variantsHP,
            .norCountPS = norCountPS,
            .tumCountPS = &tumCountPS,
            .deriveByHpSimilarity = deriveByHpSimilarity
        };

        tagResult->writeTagReadLog(data);
    }

    return hpResult;
}

int SomaticHaplotagChrProcessor::inheritHaplotype(
    float &deriveByHpSimilarity,
    double percentageThreshold,
    std::map<int, std::pair<int , int>>& somaticVarDeriveHP,
    std::map<int, int>& hpCount,
    int &hpResult,
    const bam1_t &aln
){
    int deriveByH1 = 0;
    int deriveByH2 = 0;

    for(auto& somaticVarIter : somaticVarDeriveHP){
        int baseHP = somaticVarIter.second.first;
        int deriveHP = somaticVarIter.second.second;
        if(baseHP == SnpHP::SOMATIC_H3){
            if(deriveHP == SnpHP::GERMLINE_H1){
                deriveByH1++;
            }
            else if(deriveHP == SnpHP::GERMLINE_H2){
                deriveByH2++;
            }
        }
        // if(strcmp(bam_get_qname(&aln), "SRR25005626.4310873") == 0){
        //     std::cerr << "pos: " << somaticVarIter.first+1 << " BaseHP: " << BaseHP << " deriveHP: " << deriveHP << std::endl;
        // }
    }

    int max = 0;
    int min = 0;
    int maxHp = 0;

    if(deriveByH1 > deriveByH2){
        max = deriveByH1;
        min = deriveByH2;
        maxHp = SnpHP::GERMLINE_H1;
    }
    else{
        max = deriveByH2;
        min = deriveByH1;   
        maxHp = SnpHP::GERMLINE_H2;
    }

    deriveByHpSimilarity = (max == 0) ? 0.0 : ((float)max / ((float)max + (float)min));

    if(deriveByHpSimilarity != 1.0 && deriveByHpSimilarity != 0.0){
        // std::string readID = bam_get_qname(&aln);
        // for(auto& somaticVarIter : somaticVarDeriveHP){
        //     std::cerr << "pos: " << somaticVarIter.first+1 << " BaseHP: " << somaticVarIter.second.first << " deriveHP: " << somaticVarIter.second.second << std::endl;
        // }
        // std::cerr << "deriveByH1: " << deriveByH1 << " deriveByH2: " << deriveByH2 << std::endl;
        // std::cerr << "deriveByHpSimilarity: " << deriveByHpSimilarity << std::endl;
        // std::cerr << "readID: " << readID << std::endl;
    }

    if(deriveByHpSimilarity >= percentageThreshold){
        switch(maxHp){
            case SnpHP::GERMLINE_H1:
                hpResult = ReadHP::H1_1;
                break;
            case SnpHP::GERMLINE_H2:
                hpResult = ReadHP::H2_1;
                break; 
        }
    }

    return hpResult;
}

void SomaticHaplotagChrProcessor::addAuxiliaryTags(bam1_t *aln, int& haplotype, int& pqValue, int& psValue){
    std::string haplotype_str = "";
    haplotype_str = ReadHapUtil::readHapIntToString(haplotype);

    bam_aux_append(aln, "HP", 'Z', (haplotype_str.size() + 1), (const uint8_t*) haplotype_str.c_str());
    if(psValue != VarData::NONE_PHASED_SET) bam_aux_append(aln, "PS", 'i', sizeof(psValue), (uint8_t*) &psValue);
    bam_aux_append(aln, "PQ", 'i', sizeof(pqValue), (uint8_t*) &pqValue);
}

SomaticHaplotagCigarParser::SomaticHaplotagCigarParser(
    CigarParserContext& ctx,
    int& ref_pos,
    int& query_pos,
    std::map<int, int>& tumCountPS,
    std::map<int, std::pair<int, int>>& somaticVarDeriveHP,
    SomaticReadVerifier& somaticReadCounter
):CigarParser(ctx, ref_pos, query_pos),
    tumCountPS(tumCountPS),
    somaticVarDeriveHP(somaticVarDeriveHP),
    somaticReadCounter(somaticReadCounter)
{

}

SomaticHaplotagCigarParser::~SomaticHaplotagCigarParser(){

}

void SomaticHaplotagCigarParser::processMatchOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, std::string& base, bool& isAlt, int& offset){
    somaticJudger.judgeSomaticSnpHap(currentVariantIter, ctx.chrName, base, *hpCount, *norCountPS, tumCountPS, variantsHP, nullptr, isAlt);

    //record the somatic snp derive by which germline hp in this read
    if(currentVariantIter->second.isSomaticVariant){
        int& deriveByHp = currentVariantIter->second.somaticReadDeriveByHP;
        int BaseHp = SnpHP::NONE_SNP;

        if(variantsHP->find((*currentVariantIter).first) != variantsHP->end()){
            if(variantsHP->at((*currentVariantIter).first) == SnpHP::SOMATIC_H3){
                BaseHp = SnpHP::SOMATIC_H3;
            }
        }
        somaticVarDeriveHP[(*currentVariantIter).first] = std::make_pair(BaseHp, deriveByHp);

    } 
    
    somaticReadCounter.recordRefAltAlleleCount(ctx.chrName, base, currentVariantIter);
}

void SomaticHaplotagCigarParser::processDeletionOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, bool& alreadyJudgeDel){
    somaticReadCounter.recordDelReadCount(ctx.chrName, currentVariantIter);
}

void SomaticTagLog::addParamsMessage(){
    *tagReadLog << "##normalSnpFile:" << sParams.basic.snpFile << "\n"
                << "##tumorSnpFile:" << sParams.tumorSnpFile << "\n"
                << "##svFile:" << sParams.basic.svFile << "\n"
                << "##tumorBamFile:" << sParams.tumorBamFile << "\n"
                << "##bamFile:" << sParams.basic.bamFile << "\n"
                << "##resultPrefix:" << sParams.basic.bamCfg.resultPrefix << "\n"
                << "##numThreads:" << sParams.basic.bamCfg.numThreads << "\n"
                << "##region:" << sParams.basic.bamCfg.region << "\n"
                << "##qualityThreshold:" << sParams.basic.bamCfg.qualityThreshold << "\n"
                << "##somaticCallingThreshold:" << sParams.basic.bamCfg.qualityThreshold << "\n"
                << "##percentageThreshold:" << sParams.basic.bamCfg.percentageThreshold << "\n"
                << "##tagSupplementary:" << sParams.basic.bamCfg.tagSupplementary << "\n";
}

void SomaticTagLog::writeBasicColumns(){
    *tagReadLog << "#ReadID\t"
                << "CHROM\t"
                << "ReadStart\t"
                << "Confidnet(%)\t"
                << "deriveByHpSimilarity\t"
                << "Haplotype\t"
                << "PhaseSet\t"
                << "TotalAllele\t"
                << "HP1Allele\t"
                << "HP2Allele\t"
                << "HP3Allele\t"
                << "HP4Allele\t"
                << "phasingQuality(PQ)\t"
                << "(Variant,HP)\t"
                << "(PhaseSet,Variantcount)\n";
}

void SomaticTagLog::writeTagReadLog(TagReadLog& data) {

    if (data.tumCountPS == nullptr){
        std::cerr << "[ERROR](writeTagReadLog) => Tumor PS count is nullptr" << std::endl;
        exit(1);
    }

    (*tagReadLog)<< bam_get_qname(&data.aln)                                                    << "\t"
                << data.bamHdr.target_name[data.aln.core.tid]                                   << "\t"
                << data.aln.core.pos                                                            << "\t"
                << data.norHPsimilarity                                                         << "\t"
                << data.deriveByHpSimilarity                                                    << "\t"
                << "H" << data.hpResultStr                                                      << "\t"
                << data.psResultStr                                                             << "\t"
                << data.hpCount.at(1)+data.hpCount.at(2)+data.hpCount.at(3)+data.hpCount.at(4)  << "\t"
                << data.hpCount.at(1)                                                           << "\t"
                << data.hpCount.at(2)                                                           << "\t"
                << data.hpCount.at(3)                                                           << "\t"
                << data.hpCount.at(4)                                                           << "\t"
                << data.pqValue                                                                 << "\t\t";

    // print position and HP
    for(auto v : data.variantsHP ){
        (*tagReadLog)<< " " << (v.first + 1) << "," << v.second ;
    }

    (*tagReadLog) << "\t";

    // belong PS, number of variant
    (*tagReadLog)<< "NorPS:";
    for(auto v : data.norCountPS){
        (*tagReadLog)<< " " << v.first << "," << v.second ;
    }

    (*tagReadLog)<< " TumPS:";
    for(auto v : *data.tumCountPS){
        (*tagReadLog)<< " " << v.first << "," << v.second ;
    }

    (*tagReadLog)<< "\n";
}
