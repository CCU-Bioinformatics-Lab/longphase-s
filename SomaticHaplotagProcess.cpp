#include "SomaticHaplotagProcess.h"


SomaticHaplotagProcess::SomaticHaplotagProcess(HaplotagParameters &params)
    : HaplotagProcess(params), somaticParamsMessage(params)
{
    hpBeforeInheritance = new ReadHpDistriLog();
    hpAfterInheritance = new ReadHpDistriLog();
}

SomaticHaplotagProcess::~SomaticHaplotagProcess(){
    delete hpBeforeInheritance;
    delete hpAfterInheritance;
}

void SomaticHaplotagProcess::taggingProcess()
{
    somaticParamsMessage.addParamsMessage();
    somaticParamsMessage.printParamsMessage();

    // decide on the type of tagging for VCF and BAM files
    Genome tagGeneType = TUMOR;

    VcfParser vcfParser(tagGeneType);
    // load SNP, SV, MOD vcf file
    parseVariantFiles(vcfParser);
    //decide which genome type chrVec and chrLength belong to
    setChrVecAndChrLength();
    // calculate SNP counts
    calculateSnpCounts();
    // update chromosome processing based on region
    setProcessingChromRegion();

    //somatic variant calling
    SomaticVarCaller *somaticVarCaller = new SomaticVarCaller(*chrVec, params);
    somaticVarCaller->VariantCalling(params, *chrVec, *chrLength, (*mergedChrVarinat), vcfSet, Genome::TUMOR);
    somaticVarCaller->getSomaticFlag(*chrVec, *mergedChrVarinat);
    delete somaticVarCaller;
    // return;

    // tag read
    tagRead(params, tagGeneType);

    // postprocess after haplotag
    postprocessForHaplotag();

    return;
};

void SomaticHaplotagProcess::parseVariantFiles(VcfParser& vcfParser){
    // parse common variant files
    // normal SNP, SV, MOD vcf file
    parseCommonVariantFiles(vcfParser);

    //load tumor snp vcf
    if(params.tumorSnpFile != ""){
        std::time_t begin = time(NULL);
        std::cerr<< "parsing tumor SNP VCF ... ";
        vcfParser.setParseSnpFile(true);
        vcfParser.variantParser(params.tumorSnpFile, vcfSet[Genome::TUMOR], *mergedChrVarinat);
        vcfParser.reset();
        std::cerr<< difftime(time(NULL), begin) << "s\n";
    }

    //load seqc high con file for benchmarking
    if(params.benchmarkVcf != ""){
        std::time_t begin = time(NULL);
        std::cerr<< "loading high confidence SNP ... ";
        somaticBenchmark.setTestingFunc(true);
        somaticBenchmark.loadHighConSomatic(params.benchmarkVcf, vcfSet[Genome::HIGH_CON_SOMATIC], *mergedChrVarinat);
        std::cerr<< difftime(time(NULL), begin) << "s\n";
        somaticBenchmark.parseBedFile(params.benchmarkBedFile);
        somaticBenchmark.displaySomaticVarCount(vcfSet[HIGH_CON_SOMATIC].chrVec, *mergedChrVarinat);
    }

    if(params.benchmarkBedFile != ""){
        std::time_t begin = time(NULL);
        std::cerr<< "marking variants in bed regions ... \n";
        somaticBenchmark.markVariantsInBedRegions(vcfSet[TUMOR].chrVec, *mergedChrVarinat);
        std::cerr<< difftime(time(NULL), begin) << "s\n";

        somaticBenchmark.writeBedRegionLog(vcfSet[TUMOR].chrVec, *mergedChrVarinat, params.resultPrefix);

        // std::cerr<< "removing variants out bed regions ... \n";
        // begin = time(NULL);
        // somaticBenchmark.removeVariantsOutBedRegion(*mergedChrVarinat);
        // std::cerr<< difftime(time(NULL), begin) << "s\n";
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

        // chrVec = &(vcfSet[TUMOR].chrVec);
        // chrLength = &(vcfSet[TUMOR].chrLength); 
        chrVec = &(vcfSet[Genome::NORMAL].chrVec);
        chrLength = &(vcfSet[Genome::NORMAL].chrLength); 
    }catch(const std::exception& e){
        std::cerr << "[Error] (setChrVecAndChrLength) :" << e.what() << std::endl;
        return;
    }
}

void SomaticHaplotagProcess::calculateSnpCounts(){
    int tumor_snp_count = 0;
    int normal_snp_count = 0;
    int overlap_snp_count = 0;
    for(auto& chrIter : (*chrVec)){
        auto chrVarIter = (*mergedChrVarinat)[chrIter].begin();
        while(chrVarIter != (*mergedChrVarinat)[chrIter].end()){
            if((*chrVarIter).second.isExists(TUMOR)){
                tumor_snp_count++;
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
}

void SomaticHaplotagProcess::postprocessForHaplotag(){
    std::cerr<< "postprocess for haplotag ...\n";
    if(params.writeReadLog){
        hpBeforeInheritance->writeReadHpDistriLog(params, "_readDistri_beforeInheritance.out", *chrVec);
        hpAfterInheritance->writeReadHpDistriLog(params, "_readDistri_afterInheritance.out", *chrVec);
        //write snp cover region
        hpAfterInheritance->writePosCoverRegionLog(params, "_SnpCoverRegion.out", *chrVec);
        //write read cover region in whole genome
        hpAfterInheritance->writeTagReadCoverRegionLog(params, "_readCoverRegion.bed", *chrVec, *chrLength);
        //write somatic read log
        somaticBenchmark.writeTaggedReadLog(*chrVec, params, "_totalRead.out");
        somaticBenchmark.writeTaggedSomaticReadLog(*chrVec, params, "_somaticRead.out");
        somaticBenchmark.writeCrossHighConSnpReadLog(*chrVec, params, "_crossHighConSnpRead.out");
        somaticBenchmark.writePosAlleleCountLog(*chrVec, params, "_alleleCount.out", *mergedChrVarinat);
    }
}

SomaticHaplotagBamParser::SomaticHaplotagBamParser(
    ParsingBamMode mode,
    bool writeOutputBam,
    bool mappingQualityFilter,
    ReadStatistics& readStats,
    SomaticReadBenchmark& somaticBenchmark,
    ReadHpDistriLog& hpBeforeInheritance,
    ReadHpDistriLog& hpAfterInheritance
):GermlineHaplotagBamParser(mode, writeOutputBam, mappingQualityFilter, readStats),
    somaticBenchmark(somaticBenchmark),
    hpBeforeInheritance(hpBeforeInheritance),
    hpAfterInheritance(hpAfterInheritance)
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

    bool openTestingFunc = somaticBenchmark.getTestingFunc();
    #pragma omp critical
    {
        hpBeforeInheritance.loadChrKey(chr);
        hpAfterInheritance.loadChrKey(chr);
        somaticBenchmark.loadChrKey(chr);
        localHpBeforeInheritance = hpBeforeInheritance.getChrHpResultsPtr(chr);
        localHpAfterInheritance = hpAfterInheritance.getChrHpResultsPtr(chr);
        somaticReadCounter = new SomaticReadVerifier(openTestingFunc, somaticBenchmark.getMetricsPtr(chr));
    }
    this->chr = chr;
    begin = time(NULL);
}

SomaticHaplotagChrProcessor::~SomaticHaplotagChrProcessor(){
    delete somaticReadCounter;
    std::cerr << chr << " : " << difftime(time(NULL), begin) << "s\n";
}

int SomaticHaplotagChrProcessor::judgeHaplotype(
    const bam_hdr_t &bamHdr,
    const bam1_t &aln,
    std::string chrName,
    double percentageThreshold,
    GermlineTagLog *tagResult,
    int &pqValue,
    int &psValue,
    const int tagGeneType,
    const std::string &ref_string,
    const HaplotagParameters &params,
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
    CigarParser* cigarParser = new SomaticHaplotagCigarParser(ref_pos, query_pos, tumCountPS, somaticVarDeriveHP, *somaticReadCounter);
    cigarParser->parsingCigar(aln, bamHdr, chrName, params, firstVariantIter, currentChrVariants, ref_string, hpCount, variantsHP, norCountPS);
    delete cigarParser;
    //In the current version, only normal SVs are considered, without inclusion of tumor samples
    germlineJudgeSVHap(aln, vcfSet, hpCount, Genome::NORMAL);

    int startPos = aln.core.pos + 1;
    int endPos = ref_pos;

    //the similarity of HP types
    //tumor variants
    double norHPsimilarity = 0.0;
    double tumHPsimilarity = 0.0;

    // determine the haplotype of the read
    int hpResult = ReadHP::unTag;
    hpResult = determineReadHP(hpCount, pqValue, norCountPS, norHPsimilarity, tumHPsimilarity, percentageThreshold, &localReadStats.totalHighSimilarity, &localReadStats.totalCrossTwoBlock, &localReadStats.totalWithOutVaraint);

    int hpResultBeforeInheritance = hpResult;

    float deriveByHpSimilarity = 0.0;

    if(hpResult == ReadHP::H3){
        // Inherit haplotype information for H3 reads based on somatic variant derived haplotype
        hpResult = inheritHaplotype(deriveByHpSimilarity, params.percentageThreshold, somaticVarDeriveHP, hpCount, hpResult);
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

    //std::cout << "---------------------------------------" << std::endl;
    //std::cout << "HP1: "<< hpCount[1]<< " HP2: "<< hpCount[2]<< " HP3: "<< hpCount[3]<< " HP4: "<< hpCount[4] << std::endl;
    //std::cout << "TmaxC: "<< tumorMaxHPcount<< " TminC: "<< tumorMinHPcount<< " NmaxC: "<< normalMaxHPcount<< " NminC: "<< normalMinHPcount << std::endl;
    //std::cout << "TmaxHP: "<< maxTumorHP<< " NmaxHP: "<< maxNormalHP << std::endl;
    //std::cout << "Tsimilar: "<< tumorHPsimilarity<< " Nsimilar: "<< normalHPsimilarity << std::endl;
    //std::cout << "TpsCountSize: "<< tumCountPS.size()<< " NpsCountSize: "<< norCountPS.size() << std::endl;
    //std::cout << "hpResult: "<< hpResult << std::endl;

    std::string hpResultStr = ".";
    std::string psResultStr = ".";

    hpResultStr = convertHpResultToString(hpResult);

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
    somaticReadCounter->recordTaggedRead(chrName, readID, hpResultStr, variantsHP, hpCount, norHPsimilarity, deriveByHpSimilarity, currentChrVariants);
    somaticReadCounter->recordCrossingHighConSnpRead(chrName, readID, hpResultStr, variantsHP, hpCount, norHPsimilarity, deriveByHpSimilarity, currentChrVariants);

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
    int &hpResult
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

    //float deriveByHpSimilarity = (max == 0) ? 0.0 : ((float)max / ((float)max + (float)min));
    deriveByHpSimilarity = (max == 0) ? 0.0 : ((float)max / ((float)max + (float)min));

    if(deriveByHpSimilarity != 1.0 && deriveByHpSimilarity != 0.0){
        // std::cerr << "deriveByH1: " << deriveByH1 << " deriveByH2: " << deriveByH2 << std::endl;
        // std::cerr << "deriveByHpSimilarity: " << deriveByHpSimilarity << std::endl;
        // std::cerr << "readID: " << bam_get_qname(&aln) << std::endl;
        //exit(1);
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
    haplotype_str = convertHpResultToString(haplotype);

    bam_aux_append(aln, "HP", 'Z', (haplotype_str.size() + 1), (const uint8_t*) haplotype_str.c_str());
    if(psValue != VarData::NONE_PHASED_SET) bam_aux_append(aln, "PS", 'i', sizeof(psValue), (uint8_t*) &psValue);
    bam_aux_append(aln, "PQ", 'i', sizeof(pqValue), (uint8_t*) &pqValue);
}

std::string SomaticHaplotagChrProcessor::convertHpResultToString(int hpResult){
    switch(hpResult){
        case ReadHP::unTag:
            return ".";
        case ReadHP::H1:
            return "1";
        case ReadHP::H2:
            return "2";
        case ReadHP::H3:
            return "3";
        case ReadHP::H4:
            return "4";
        case ReadHP::H1_1:
            return "1-1";
        case ReadHP::H2_1:
            return "2-1";
        case ReadHP::H1_2:
            return "1-2";
        case ReadHP::H2_2:
            return "2-2";
        default:
            std::cerr << "ERROR (convertHpResultToString) => Unsupported HP result: " << hpResult << std::endl;
            exit(1);
    }
}

SomaticHaplotagCigarParser::SomaticHaplotagCigarParser(
    int& ref_pos,
    int& query_pos,
    std::map<int, int>& tumCountPS,
    std::map<int, std::pair<int, int>>& somaticVarDeriveHP,
    SomaticReadVerifier& somaticReadCounter
):CigarParser(ref_pos, query_pos),
    tumCountPS(tumCountPS),
    somaticVarDeriveHP(somaticVarDeriveHP),
    somaticReadCounter(somaticReadCounter)
{

}

SomaticHaplotagCigarParser::~SomaticHaplotagCigarParser(){

}

void SomaticHaplotagCigarParser::OnlyTumorSNPjudgeHP(const std::string &chrName, int &curPos, MultiGenomeVar &curVar, std::string base, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *tumorAllelePosVec){

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

        if(curVar.Variant[TUMOR].is_phased_hetero){
            if(tumCountPS != nullptr) (*tumCountPS)[curVar[TUMOR].PhasedSet]++;
        }
    }
}

void SomaticHaplotagCigarParser::processMatchOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, std::string& base){
    SomaticJudgeSnpHP(currentVariantIter, *chrName, base, *hpCount, *norCountPS, tumCountPS, variantsHP, nullptr);

    //record the somatic snp derive by which germline hp in this read
    if(currentVariantIter->second.isSomaticVariant){
        int& deriveByHp = currentVariantIter->second.somaticReadDeriveByHP;
        int BaseHp = SnpHP::NONE_SNP;
        // if(base == (*currentVariantIter).second.Variant[TUMOR].Alt){
        //     BaseHp = SnpHP::SOMATIC_H3;
        // }
        if(variantsHP->find((*currentVariantIter).first) != variantsHP->end()){
            if(variantsHP->at((*currentVariantIter).first) == SnpHP::SOMATIC_H3){
                BaseHp = SnpHP::SOMATIC_H3;
            }
        }
        somaticVarDeriveHP[(*currentVariantIter).first] = std::make_pair(BaseHp, deriveByHp);

    } 
    
    somaticReadCounter.recordRefAltAlleleCount(*chrName, base, currentVariantIter);
}

void SomaticHaplotagCigarParser::processDeletionOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, bool& alreadyJudgeDel){
    somaticReadCounter.recordDelReadCount(*chrName, currentVariantIter);
}

void SomaticTagLog::writeTagReadLog(TagReadLog& data) {

    if (data.tumCountPS == nullptr){
        std::cerr << "ERROR (writeTagReadLog) => Tumor PS count is nullptr" << std::endl;
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
