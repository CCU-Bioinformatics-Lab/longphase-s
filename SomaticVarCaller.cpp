#include "SomaticVarCaller.h"

/**
 * @brief Calculate common base information for tumor SNPs
 * 
 * Computes VAF, depth ratios, and haplotype imbalance metrics for tumor-normal analysis.
 * This function is used by both normal and tumor data extraction to ensure consistent
 * calculation of base statistics.
 * 
 * @param baseInfo Base information structure to update with calculated metrics
 * @param tumorAltBase Alternative base in tumor sample for VAF calculation
 */
void tumor_normal_analysis::calculateBaseCommonInfo(PosBase& baseInfo, std::string& tumorAltBase){
    int &depth = baseInfo.depth;
    int AltCount = baseInfo.getBaseCount(tumorAltBase);

    int filteredMpqDepth = baseInfo.filteredMpqDepth;
    int filteredMpqAltCount = baseInfo.getMpqBaseCount(tumorAltBase);

    baseInfo.VAF = base_analysis::calculateVAF(AltCount, depth);
    baseInfo.filteredMpqVAF = base_analysis::calculateVAF(filteredMpqAltCount, filteredMpqDepth);
    baseInfo.nonDelVAF = base_analysis::calculateVAF(AltCount, (depth - baseInfo.delCount));

    baseInfo.lowMpqReadRatio = base_analysis::calculateLowMpqReadRatio(depth, filteredMpqDepth);
    baseInfo.delRatio = base_analysis::calculateDelRatio(baseInfo.delCount, depth);

    //read hp count in the normal bam
    int H1readCount = baseInfo.ReadHpCount[ReadHP::H1];
    int H2readCount = baseInfo.ReadHpCount[ReadHP::H2];
    int germlineReadHpCount = H1readCount + H2readCount;

    baseInfo.germlineHaplotypeImbalanceRatio = base_analysis::calculateHaplotypeImbalanceRatio(H1readCount, H2readCount, germlineReadHpCount);
    baseInfo.percentageOfGermlineHp = base_analysis::calculatePercentageOfGermlineHp(germlineReadHpCount, depth);
}


double statisticsUtils::calculateMean(const std::map<int, double>& data){
    double size = data.size();
    if(size == 0) return 0.0;
    
    double sum = 0.0;
    for(auto meanIter : data){
        sum += meanIter.second;
    }
    return sum / size;
}

double statisticsUtils::calculateStandardDeviation(const std::map<int, double>& data, double mean){
    double variance = 0.0;
    for (auto meanIter : data) {
        double value = meanIter.second;
        variance += (value - mean) * (value - mean);
    }
    return std::sqrt(variance / data.size());
}

void statisticsUtils::calculateZScores(const std::map<int, double>& data, double mean, double stdDev, std::map<int, double> &zScores){
    for (auto meanIter : data) {
        double value = meanIter.second;
        if(stdDev == 0){
            zScores[meanIter.first] = 0.0;
        }else{
            zScores[meanIter.first] = ((value - mean) / stdDev); // calculate z-score
        }
    }
}


ExtractNorDataBamParser::ExtractNorDataBamParser(
    const ParsingBamConfig &config, 
    const ParsingBamControl &control,
    std::map<std::string, std::map<int, PosBase>>& chrPosNorBase
)
: HaplotagBamParser(config, control), chrPosNorBase(chrPosNorBase){}

ExtractNorDataBamParser::~ExtractNorDataBamParser(){
};

void ExtractNorDataBamParser::displayPosInfo(std::string chr, int pos){
    if(chrPosNorBase[chr].find(pos) == chrPosNorBase[chr].end()){
        std::cerr << "[ERROR](displayPosInfo) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }else{
        PosBase VarBase = chrPosNorBase[chr][pos];
        std::cout << " chr : " <<  chr <<" Pos :" << pos << std::endl;
        std::cout << " A_count : " << VarBase.A_count << std::endl;
        std::cout << " C_count : " << VarBase.C_count << std::endl;
        std::cout << " G_count : " << VarBase.G_count << std::endl;
        std::cout << " T_count : " << VarBase.T_count << std::endl;
        std::cout << " depth :   " << VarBase.depth << std::endl;
        std::cout << " unknow :  " << VarBase.unknow << std::endl;
        std::cout << " deletion count :  " << VarBase.delCount << std::endl;
        std::cout << " VAF :  " << VarBase.VAF << std::endl;
        std::cout << " filterd MPQ VAF :  " << VarBase.filteredMpqVAF << std::endl;
        std::cout << " Low MPQ read ratio :  " << VarBase.lowMpqReadRatio << std::endl;
    }
}

ExtractNorDataChrProcessor::ExtractNorDataChrProcessor(std::map<std::string, std::map<int, PosBase>> &chrPosNorBase, const std::string &chr){
    #pragma omp critical
    {
        variantBase = &((chrPosNorBase)[chr]);
    }
}

ExtractNorDataChrProcessor::~ExtractNorDataChrProcessor(){
    variantBase = nullptr;
}

void ExtractNorDataChrProcessor::processRead(
    bam1_t &aln, 
    const bam_hdr_t &bamHdr,
    const std::string &ref_string,
    std::map<int, MultiGenomeVar> &currentVariants,
    std::map<int, MultiGenomeVar>::iterator &firstVariantIter,
    ChrProcContext& ctx
){

    std::map<int, int> hpCount;

    hpCount[SnpHP::GERMLINE_H1] = 0; 
    hpCount[SnpHP::GERMLINE_H2] = 0;

    //record variants on this read
    std::map<int, int> variantsHP;

    // normal variant phased set, count
    std::map<int,int> norCountPS;

    std::vector<int> tumVarPosVec;

    // position relative to reference
    int ref_pos = aln.core.pos;
    // position relative to read
    int query_pos = 0;

    /// Create a CIGAR parser using polymorphism design
    // Use extractNorDataCigarParser to process CIGAR strings for normal samples    
    CigarParserContext cigarCtx(aln, bamHdr, ctx.chrName, ctx.params, firstVariantIter, currentVariants, ref_string);
    CigarParser* cigarParser = new ExtractNorDataCigarParser(cigarCtx, *variantBase, tumVarPosVec, ref_pos, query_pos, ctx.params.qualityThreshold);
    cigarParser->parsingCigar(hpCount, variantsHP, norCountPS);
    delete cigarParser;

    // get the number of SVs occurring on different haplotypes in a read
    if( aln.core.qual >= ctx.params.qualityThreshold ){
        judger.judgeSVHap(aln, ctx.vcfSet, hpCount, ctx.genomeSample);
    }

    double min = 0.0;
    double max = 0.0;
    int hpResult = ReadHP::unTag;
    int pqValue = 0;
    int psValue = 0;
    double percentageThreshold = ctx.params.percentageThreshold;
    // determine the haplotype of the read
    hpResult = judger.judgeReadHap(hpCount, min, max, percentageThreshold, pqValue, psValue, norCountPS, nullptr, nullptr);

    //record read hp to tumor SNP position
    for(auto pos : tumVarPosVec){
        (*variantBase)[pos].ReadHpCount[hpResult]++;
    }
}

void ExtractNorDataChrProcessor::postProcess(
    const std::string &chr,
    std::map<int, MultiGenomeVar> &currentVariants
){
    std::map<int, PosBase>::iterator currentPosIter = variantBase->begin();

    while (currentPosIter != variantBase->end()) {


        if(!currentVariants[(*currentPosIter).first].isExists(TUMOR)){
            std::cerr << "[ERROR](extractNorData:postProcess) => can't find the position : chr:" << chr << " pos: " << (*currentPosIter).first;
            exit(1);
        }

        PosBase &baseInfo = (currentPosIter->second);
        auto curVar = currentVariants[(*currentPosIter).first].Variant[TUMOR];

        // current variant is SNP
        if(curVar.variantType == VariantType::SNP){
            std::string& tumAltBase = curVar.allele.Alt;

            // calculate the base information of the tumor SNP
            tumor_normal_analysis::calculateBaseCommonInfo(baseInfo, tumAltBase);
        }
        
        currentPosIter++;
    } 
}

ExtractNorDataCigarParser::ExtractNorDataCigarParser(
    CigarParserContext& ctx,
    std::map<int, PosBase>& variantBase,
    std::vector<int>& tumVarPosVec,
    int& ref_pos,
    int& query_pos,
    const int& mappingQualityThr)
:CigarParser(ctx, ref_pos, query_pos), variantBase(variantBase), tumVarPosVec(tumVarPosVec), mappingQualityThr(mappingQualityThr){

}

ExtractNorDataCigarParser::~ExtractNorDataCigarParser(){

}

void ExtractNorDataCigarParser::processMatchOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, std::string& base){
    //statistically analyze SNP information exclusive to the tumor
    if((*currentVariantIter).second.isExists(TUMOR)){
        int curPos = (*currentVariantIter).first;

        auto curVar = (*currentVariantIter).second.Variant[TUMOR];
        
        // the variant is SNP
        if(curVar.variantType == VariantType::SNP){
            //record tumor SNP position
            tumVarPosVec.push_back(curPos);

            //counting current tumor SNP base and depth           
            countBaseNucleotide(variantBase[curPos], base, ctx.aln, mappingQualityThr);
        }
        // the indel(deletion) SNP position is the start position, and the deletion occurs at the next position
        else if(curVar.variantType == VariantType::DELETION){
            // the indel SNP start position is at the end of the deletion, and the next cigar operator is deletion
            if(curPos == (ref_pos + length - 1) && bam_cigar_op(cigar[i+1]) == 2 && i+1 < aln_core_n_cigar){

            }
        }
    }       

    if ( ctx.aln.core.qual >= mappingQualityThr && (*currentVariantIter).second.isExists(NORMAL)){
        // only judge the heterozygous SNP
        if((*currentVariantIter).second.Variant[NORMAL].GT == GenomeType::PHASED_HETERO){
            auto norVar = (*currentVariantIter).second.Variant[NORMAL];
            judger.judgeSnpHap(ctx.chrName, norVar, base, ref_pos, length, i, aln_core_n_cigar, cigar, currentVariantIter, *hpCount, *variantsHP, *norCountPS);
        }
    }
}

void ExtractNorDataCigarParser::processDeletionOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, bool& alreadyJudgeDel){
    //statistically analyze SNP information exclusive to the tumor
    if((*currentVariantIter).second.isExists(TUMOR)){

        auto curVar = (*currentVariantIter).second.Variant[TUMOR];

        int curPos = (*currentVariantIter).first;
        
        //record tumor SNP position
        tumVarPosVec.push_back(curPos);

        // the variant is SNP
        if(curVar.variantType == VariantType::SNP){
            countDeletionBase(variantBase[curPos]);
        }
        // the variant is deletion
        else if(curVar.variantType == VariantType::DELETION){

        }
    }

    // only execute at the first phased normal snp
    if ( ctx.aln.core.qual >= mappingQualityThr && (*currentVariantIter).second.isExists(NORMAL) && !alreadyJudgeDel){
        if((*currentVariantIter).second.Variant[NORMAL].GT == GenomeType::PHASED_HETERO){
            // longphase v1.73 only execute once
            alreadyJudgeDel = true;
            judger.judgeDeletionHap(ctx.chrName, ctx.ref_string, ref_pos, length, query_pos, currentVariantIter, &ctx.aln, *hpCount, *variantsHP, *norCountPS);
        }
    }
}


ExtractTumDataBamParser::ExtractTumDataBamParser(
    const ParsingBamConfig &config,
    const ParsingBamControl &control,
    std::map<std::string, std::map<int, SomaticData>>& chrPosSomaticInfo,
    std::map<std::string, std::map<std::string, ReadVarHpCount>>& chrReadHpResultSet,
    std::map<std::string, std::map<int, std::map<std::string, int>>>& chrTumorPosReadCorrBaseHP
):
  HaplotagBamParser(config, control),
  chrPosSomaticInfo(chrPosSomaticInfo),
  chrReadHpResultSet(chrReadHpResultSet),
  chrTumorPosReadCorrBaseHP(chrTumorPosReadCorrBaseHP){

}

ExtractTumDataBamParser::~ExtractTumDataBamParser(){

}



ExtractTumDataChrProcessor::ExtractTumDataChrProcessor(
    std::map<std::string, std::map<int, SomaticData>>& chrPosSomaticInfo,
    std::map<std::string, std::map<std::string, ReadVarHpCount>>& chrReadHpResultSet,
    std::map<std::string, std::map<int, std::map<std::string, int>>>& chrTumorPosReadCorrBaseHP,
    const std::string &chr
){
    #pragma omp critical
    {
        somaticPosInfo = &((chrPosSomaticInfo)[chr]);
        readHpResultSet = &((chrReadHpResultSet)[chr]);
        tumorPosReadCorrBaseHP = &((chrTumorPosReadCorrBaseHP)[chr]);
    }
}

ExtractTumDataChrProcessor::~ExtractTumDataChrProcessor(){

}

void ExtractTumDataChrProcessor::processRead(
    bam1_t &aln, 
    const bam_hdr_t &bamHdr,
    const std::string &ref_string,
    std::map<int, MultiGenomeVar> &currentVariants,
    std::map<int, MultiGenomeVar>::iterator &firstVariantIter,
    ChrProcContext& ctx
){
   
    std::map<int, int> hpCount;
    hpCount[SnpHP::GERMLINE_H1] = 0; 
    hpCount[SnpHP::GERMLINE_H2] = 0;
    hpCount[SnpHP::SOMATIC_H3] = 0;
    hpCount[SnpHP::SOMATIC_H4] = 0;
    //record variants on this read
    std::map<int, int> variantsHP;

    //record tumor-unique variants with low VAF in the normal.bam on this read
    std::vector<int> tumorAllelePosVec;

    //record tumor SNPs on this read
    std::vector<int> tumorSnpPosVec;

    //record PS count( PS value, count)
    std::map<int, int> tumCountPS;
    std::map<int, int> norCountPS;

    // // position relative to reference
    int ref_pos = aln.core.pos;
    // position relative to read
    int query_pos = 0;

    // Create CIGAR parser for tumor data extraction
    CigarParserContext cigarCtx(aln, bamHdr, ctx.chrName, ctx.params, firstVariantIter, currentVariants, ref_string);
    ExtractTumDataCigarParser cigarParser(cigarCtx, *somaticPosInfo, tumorAllelePosVec, tumorSnpPosVec, tumCountPS, ref_pos, query_pos, ctx.params.qualityThreshold);
    cigarParser.parsingCigar(hpCount, variantsHP, norCountPS);

    int pqValue = 0;   
    double normalHPsimilarity = 0.0;
    double tumorHPsimilarity = 0.0;

    //calculate germline read HP result for estimate tumor purity
    //calculate somatic read HP result for somatic read consistency filter
    int hpResult = somaticJudger.judgeSomaticReadHap(hpCount, pqValue, norCountPS, normalHPsimilarity, tumorHPsimilarity, ctx.params.percentageThreshold, nullptr, nullptr, nullptr);
   
    //classify read cases where tumor SNPs have low VAF in normal samples
    if(!tumorAllelePosVec.empty()){
        classifyReadsByCase(tumorAllelePosVec, norCountPS, hpCount, *somaticPosInfo);
        
        for(auto pos : tumorAllelePosVec){
            int baseHP = SnpHP::NONE_SNP;
            if(variantsHP.find(pos) != variantsHP.end()){
                baseHP = variantsHP[pos];
            }else{
                std::cerr << "[ERROR](SomaticStatisticSomaticPosInfo) => can't find the position" << std::endl;
                std::cerr << "chr:" << ctx.chrName << " pos: " << pos+1 << std::endl;
                std::cerr << "readID: " << bam_get_qname(&aln) << std::endl;
                exit(1);
            }
            if(baseHP != SnpHP::SOMATIC_H3){
                std::cerr << "[ERROR](SomaticStatisticSomaticPosInfo) => baseHP is not HP3 : chr:" << ctx.chrName << " pos: " << pos+1 << " baseHP: " << baseHP << std::endl;
                exit(1);
            }
            if(hpResult == ReadHP::H1_1 || hpResult == ReadHP::H2_1 || hpResult == ReadHP::H3 || hpResult == ReadHP::unTag){
                // record the somatic read HP for somatic read consistency filter
                (*somaticPosInfo)[pos].somaticReadHpCount[hpResult]++;
            }else if(hpResult == ReadHP::H1 || hpResult == ReadHP::H2){
                std::cerr << "[ERROR](SomaticStatisticSomaticPosInfo) => error somatic read HP : chr:" << ctx.chrName << " pos: " << pos+1 << " hpResult: " << hpResult << std::endl;
                exit(1);      
            }
        }
    }

    //record variants HP count for each read in tumor-only position
    if(!tumorSnpPosVec.empty()){

        std::string readID = bam_get_qname(&aln);
        // Prevent override of read ID by adding counter suffix
        if((*readHpResultSet).find(readID) != (*readHpResultSet).end()){
            (*readHpResultSet)[readID].readIDcount++;
            readID = readID + "-" + std::to_string((*readHpResultSet)[readID].readIDcount);
        }

        // Store haplotype counts and read metadata
        (*readHpResultSet)[readID].HP1= hpCount[1];
        (*readHpResultSet)[readID].HP2= hpCount[2];
        (*readHpResultSet)[readID].HP3= hpCount[3];
        (*readHpResultSet)[readID].HP4= hpCount[4];

        (*readHpResultSet)[readID].norCountPS = norCountPS;
        (*readHpResultSet)[readID].startPos = aln.core.pos + 1;
        (*readHpResultSet)[readID].endPos = ref_pos;
        (*readHpResultSet)[readID].readLength = query_pos;
        
        for(auto pos : tumorSnpPosVec){
            int tumorSnpBaseHP = SnpHP::NONE_SNP;

            if(variantsHP.find(pos) != variantsHP.end()){
                tumorSnpBaseHP = variantsHP[pos];
            }
            // record the base HP whatever it is germline or somatic for statistic read distribution at current position
            (*tumorPosReadCorrBaseHP)[pos][readID] = tumorSnpBaseHP;
            // record the base HP for estimate tumor purity
            (*somaticPosInfo)[pos].base.ReadHpCount[hpResult]++;
        } 
    }
}

void ExtractTumDataChrProcessor::classifyReadsByCase(std::vector<int> &tumorAllelePosVec, std::map<int, int> &norCountPS, std::map<int, int> &hpCount, std::map<int, SomaticData> &somaticPosInfo){
    
    //decide whether to tag the read or not
    bool recordRead = true;

    //only exist tumor homo SNP
    if( norCountPS.size() == 0 && hpCount[3] != 0){
        // std::cerr << "[ERROR](classifyReadsByCase) => read only had tumor homo SNP " << bam_get_qname(aln) << "\n";
    }

    // germline SNPs cross two block
    if( norCountPS.size() > 1){
        recordRead = false;
    }

    int zero_count = 0;
    float CleanHP3Threshold = 1.0;
    bool tagCleanHP3Read = false;

    if (hpCount[1] == 0) zero_count++;
    if (hpCount[2] == 0) zero_count++;

    if(hpCount[3] == 0 && hpCount[4] == 0){
        std::cerr << "[ERROR](classifyReadsByCase) => hp3 or hp4 count is 0" << std::endl;
        std::cerr << "hp3 count:" << hpCount[3] << " hp4 count:" << hpCount[4];
        exit(1);
    }
    
    // only HP3 or one type of HP SNP
    if((zero_count == 1 || zero_count == 2) && hpCount[3] != 0){
        tagCleanHP3Read = true;
    }else if(hpCount[1] + hpCount[2] != 0){
        float hp1Ratio = static_cast<float>(hpCount[1]) / (hpCount[1] + hpCount[2]);
        float hp2Ratio = static_cast<float>(hpCount[2]) / (hpCount[1] + hpCount[2]);

        if ((hp1Ratio >= CleanHP3Threshold) || (hp2Ratio >= CleanHP3Threshold)) {
            tagCleanHP3Read = true;
        }
    }

    for (int somaticPos : tumorAllelePosVec) {
        if(recordRead == false){
            somaticPosInfo[somaticPos].unTag++;
        }else if(tagCleanHP3Read){
            somaticPosInfo[somaticPos].totalCleanHP3Read++;
            if(hpCount[1] == 0 && hpCount[2] == 0 && hpCount[3] != 0){
                somaticPosInfo[somaticPos].pure_H3_read++;
            }else if(hpCount[1] !=0  && hpCount[2] == 0){ 
                somaticPosInfo[somaticPos].pure_H1_1_read++;      
            }else if(hpCount[1] == 0 && hpCount[2] != 0){ 
                somaticPosInfo[somaticPos].pure_H2_1_read++;                            
            }
        }else{
            somaticPosInfo[somaticPos].Mixed_HP_read++;
        }
    }
}

void ExtractTumDataChrProcessor::postProcess( const std::string &chr, std::map<int, MultiGenomeVar> &currentVariants){
    //calculate the information of the somatic positon 
    std::map<int, SomaticData>::iterator somaticVarIter = somaticPosInfo->begin();
    while( somaticVarIter != somaticPosInfo->end()){
        
        if(!currentVariants[(*somaticVarIter).first].isExists(TUMOR)){
            std::cerr << "[ERROR](extractTumData:postProcess) => can't find the position : chr:" << chr << " pos: " << (*somaticVarIter).first;
            exit(1);
        }

        auto curVar = currentVariants[(*somaticVarIter).first].Variant[TUMOR];

        // current variant is SNP
        if (curVar.variantType == VariantType::SNP) {

            std::string RefBase = curVar.allele.Ref;
            std::string AltBase = curVar.allele.Alt;

            //calculating Read Case Ratio
            int totalCleanHP3Read = (*somaticVarIter).second.totalCleanHP3Read;
            int totalHP1WithHP3Read = (*somaticVarIter).second.pure_H1_1_read;
            int totalHP2WithHP3Read = (*somaticVarIter).second.pure_H2_1_read;
            int totalOnlyHP3Read = (*somaticVarIter).second.pure_H3_read;
            int totalMessyHPRead = (*somaticVarIter).second.Mixed_HP_read;

            (*somaticVarIter).second.CaseReadCount = totalCleanHP3Read + totalMessyHPRead;

            // calculate ratios if there are case reads
            if((*somaticVarIter).second.CaseReadCount != 0){
                (*somaticVarIter).second.Mixed_HP_readRatio = (float)totalMessyHPRead / ((float)totalCleanHP3Read + (float)totalMessyHPRead);
                (*somaticVarIter).second.pure_H1_1_readRatio = (float)totalHP1WithHP3Read / ((float)totalCleanHP3Read + (float)totalMessyHPRead);
                (*somaticVarIter).second.pure_H2_1_readRatio = (float)totalHP2WithHP3Read / ((float)totalCleanHP3Read + (float)totalMessyHPRead);
                (*somaticVarIter).second.pure_H3_readRatio = (float)totalOnlyHP3Read / ((float)totalCleanHP3Read + (float)totalMessyHPRead);
            }else{
                (*somaticVarIter).second.Mixed_HP_readRatio = 0.0;
                (*somaticVarIter).second.pure_H1_1_readRatio = 0.0;
                (*somaticVarIter).second.pure_H2_1_readRatio = 0.0;
                (*somaticVarIter).second.pure_H3_readRatio = 0.0;
            }

            // calculate the base information of the tumor SNP
            tumor_normal_analysis::calculateBaseCommonInfo((*somaticVarIter).second.base, AltBase);

            // the distribution of reads HP at the current position
            int H1readCount = (*somaticVarIter).second.base.ReadHpCount[ReadHP::H1];
            int H2readCount = (*somaticVarIter).second.base.ReadHpCount[ReadHP::H2];

            int H1_1readCount = (*somaticVarIter).second.base.ReadHpCount[ReadHP::H1_1];
            int H2_1readCount = (*somaticVarIter).second.base.ReadHpCount[ReadHP::H2_1];

            // the count of reads based on germline HP
            int BaseOnH1ReadCount = H1readCount + H1_1readCount;
            int BaseOnH2ReadCount = H2readCount + H2_1readCount;
            int totalBaseOnGermlineReadHpCount = BaseOnH1ReadCount + BaseOnH2ReadCount;

            // the ratio of the read count of H1 and H2 based on germline HP
            (*somaticVarIter).second.allelicImbalanceRatio = base_analysis::calculateHaplotypeImbalanceRatio(BaseOnH1ReadCount, BaseOnH2ReadCount, totalBaseOnGermlineReadHpCount);

            // the ratio of the read count of H1 and H2 based on somatic HP
            int totalSomaticReadHpCount = H1_1readCount + H2_1readCount;
            (*somaticVarIter).second.somaticHaplotypeImbalanceRatio = base_analysis::calculateHaplotypeImbalanceRatio(H1_1readCount, H2_1readCount, totalSomaticReadHpCount);

            // current SNP GT type
            if(curVar.GT == GenomeType::UNPHASED_HOMO){
                (*somaticVarIter).second.GTtype = "Homo";
            }else if(curVar.GT == GenomeType::PHASED_HETERO){
                (*somaticVarIter).second.GTtype = "Hetero";
            }else if(curVar.GT == GenomeType::UNPHASED_HETERO){
                (*somaticVarIter).second.GTtype = "UnphasedHetero";
            }else{
                std::cerr << "[ERROR](GTtype) => can't find GTtype at chr: " << chr << " pos: " << ((*somaticVarIter).first) + 1;
                exit(1);  
            }    
        }
        somaticVarIter++;
    }
}

ExtractTumDataCigarParser::ExtractTumDataCigarParser(
    CigarParserContext& ctx,
    std::map<int, SomaticData> &somaticPosInfo, 
    std::vector<int>& tumorAllelePosVec, 
    std::vector<int>& tumorSnpPosVec, 
    std::map<int, int>& tumCountPS,
    int& ref_pos,
    int& query_pos,
    const int& mappingQualityThr
):CigarParser(ctx, ref_pos, query_pos),
  somaticPosInfo(somaticPosInfo),
  tumorAllelePosVec(tumorAllelePosVec),
  tumorSnpPosVec(tumorSnpPosVec),
  tumCountPS(tumCountPS),
  mappingQualityThr(mappingQualityThr){

}

ExtractTumDataCigarParser::~ExtractTumDataCigarParser(){

}

void ExtractTumDataCigarParser::processMatchOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, std::string& base){
    //waring : using ref length to split SNP and indel that will be effect case ratio result 
    if ( ctx.aln.core.qual >= mappingQualityThr ){
        somaticJudger.judgeSomaticSnpHap(currentVariantIter, ctx.chrName, base, *hpCount, *norCountPS, tumCountPS, variantsHP, &tumorAllelePosVec);
        if((*currentVariantIter).second.isExists(TUMOR)){
            tumorSnpPosVec.push_back((*currentVariantIter).first);
        }
    }

    //statistically analyze SNP information exclusive to the tumor
    if((*currentVariantIter).second.isExists(TUMOR)){

        auto curVar = (*currentVariantIter).second.Variant[TUMOR];

        if(curVar.variantType == VariantType::SNP){
            //counting current tumor SNP base and depth           
            countBaseNucleotide(somaticPosInfo[(*currentVariantIter).first].base, base, ctx.aln, mappingQualityThr);
        }
    }
}

void ExtractTumDataCigarParser::processDeletionOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, bool& alreadyJudgeDel){
    //statistically analyze SNP information exclusive to the tumor
    if((*currentVariantIter).second.isExists(TUMOR)){

        auto curVar = (*currentVariantIter).second.Variant[TUMOR];

        int curPos = (*currentVariantIter).first;

        if(curVar.variantType == VariantType::SNP){
            countDeletionBase(somaticPosInfo[curPos].base);
        }
        // the indel SNP start position isn't at the end of the deletion
        else if(curVar.variantType == VariantType::DELETION){

        }
    }
}

SomaticVarCaller::SomaticVarCaller(const CallerConfig &callerCfg, const ParsingBamConfig &bamCfg, const std::vector<std::string> &chrVec)
:callerCfg(callerCfg),
 bamCfg(bamCfg),
 chrVec(chrVec)
{
    chrPosSomaticInfo = new std::map<std::string, std::map<int, SomaticData>>();
    chrPosNorBase = new std::map<std::string, std::map<int, PosBase>>();
    callerReadHpDistri = new ReadHpDistriLog();
    denseTumorSnpInterval = new std::map<std::string, std::map<int, std::pair<int, DenseSnpInterval>>>();
    chrReadHpResultSet = new std::map<std::string, std::map<std::string, ReadVarHpCount>>();
    chrTumorPosReadCorrBaseHP = new std::map<std::string, std::map<int, std::map<std::string, int>>>();

    for(auto chr : chrVec){
        (*chrPosSomaticInfo)[chr] = std::map<int, SomaticData>();
        (*chrPosNorBase)[chr] = std::map<int, PosBase>();
        callerReadHpDistri->loadChrKey(chr);
        (*denseTumorSnpInterval)[chr] = std::map<int, std::pair<int, DenseSnpInterval>>();
        (*chrReadHpResultSet)[chr] = std::map<std::string, ReadVarHpCount>();
        (*chrTumorPosReadCorrBaseHP)[chr] = std::map<int, std::map<std::string, int>>();
    }
}

SomaticVarCaller::~SomaticVarCaller(){
    releaseMemory();
}

void SomaticVarCaller::releaseMemory(){
    delete chrPosSomaticInfo;
    delete chrPosNorBase;
    delete callerReadHpDistri;
    delete denseTumorSnpInterval;
    delete chrReadHpResultSet;
    delete chrTumorPosReadCorrBaseHP;
}

void SomaticVarCaller::variantCalling(
    const CallerContext &ctx,
    const std::vector<std::string> &chrVec,
    const std::map<std::string, int> &chrLength,
    std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants,
    std::map<Genome, VCF_Info> &vcfSet
){

    // extract somatic data from normal and tumor BAM
    extractSomaticData(ctx.normalBamFile, ctx.tumorBamFile, ctx.fastaFile, bamCfg, chrVec, chrLength, chrMultiVariants, vcfSet);

    double tumorPurity = 0.0;

    //estimate tumor purity or use the value from parameters
    if(callerCfg.estimateTumorPurity){
        tumorPurity = runTumorPurityEstimator(callerCfg.writeCallingLog, bamCfg.resultPrefix);
    }else{
        tumorPurity = callerCfg.tumorPurity;
    }

    // set filter params with tumor purity
    SetFilterParamsWithPurity(somaticParams, tumorPurity);
    
    std::cerr<< "calling somatic variants ... ";
    std::time_t begin = time(NULL);
    // loop all chromosome
    #pragma omp parallel for schedule(dynamic) num_threads(bamCfg.numThreads) 
    for(auto chr : chrVec ){
        // record the position that tagged as HP3
        // chr, variant position
        std::map<int, SomaticData> *somaticPosInfo = nullptr;
        // read ID, reads hpResult 
        chrReadHpResult *localCallerReadHpDistri = nullptr;
        
        // pos, (endPos, denseSnpInterval)
        std::map<int, std::pair<int, DenseSnpInterval>> *localDenseTumorSnpInterval = nullptr;

        // read ID, SNP HP count 
        std::map<std::string, ReadVarHpCount> *readHpResultSet = nullptr;
        // position, read ID, baseHP 
        std::map<int, std::map<std::string, int>> *tumorPosReadCorrBaseHP = nullptr;

        // records all variants within this chromosome.
        std::map<int, MultiGenomeVar> currentChrVariants;


        #pragma omp critical
        {
            somaticPosInfo = &((*chrPosSomaticInfo)[chr]);
            localCallerReadHpDistri = callerReadHpDistri->getChrHpResultsPtr(chr);
            localDenseTumorSnpInterval = &((*denseTumorSnpInterval)[chr]);
            readHpResultSet = &((*chrReadHpResultSet)[chr]);
            tumorPosReadCorrBaseHP = &((*chrTumorPosReadCorrBaseHP)[chr]);
            currentChrVariants = chrMultiVariants[chr];
        }

        //get close somatic SNP interval
        getDenseTumorSnpInterval(*somaticPosInfo, *readHpResultSet, *tumorPosReadCorrBaseHP, *localDenseTumorSnpInterval);

        //calculate information and filter somatic SNPs
        somaticFeatureFilter(somaticParams, currentChrVariants, chr, *somaticPosInfo, tumorPurity);

        //calibrate read HP to remove low confidence H3 SNP
        calibrateReadHP(chr, *somaticPosInfo, *readHpResultSet, *tumorPosReadCorrBaseHP);

        //calculate all read HP result
        calculateReadSetHP(chr, *readHpResultSet, *tumorPosReadCorrBaseHP, bamCfg.percentageThreshold);

        //statistic all read HP in somatic SNP position
        statisticSomaticPosReadHP(chr, *somaticPosInfo, *tumorPosReadCorrBaseHP, *readHpResultSet, *localCallerReadHpDistri);
        
        //find other somatic variants Haplotype (HP4/HP5)
        // findOtherSomaticSnpHP(chr, *somaticPosInfo, currentChrVariants);

    }
    std::cerr<< difftime(time(NULL), begin) << "s\n";

    // [debug] calling somatic SNPs count
    displayCallingSnpCount();
    
    if(callerCfg.writeCallingLog){
        //write the log file for variants with positions tagged as HP3
        writeSomaticVarCallingLog(ctx ,somaticParams, chrVec, chrMultiVariants);
        writeDenseTumorSnpIntervalLog(bamCfg.resultPrefix + "_dense_tumor_snp_interval.log", chrVec);

        callerReadHpDistri->writeReadHpDistriLog(bamCfg.resultPrefix + "_read_distri_scaller.out", chrVec);
        // remove the position that derived by HP1 or HP2
        callerReadHpDistri->removeNotDeriveByH1andH2pos(chrVec);
        // record the position that can't derived because exist H1-1 and H2-1 reads
        callerReadHpDistri->writeReadHpDistriLog(bamCfg.resultPrefix + "_read_distri_scaller_derive_by_H1_H2.out", chrVec);

        // writeOtherSomaticHpLog(bamCfg.resultPrefix + "_other_allele_somatic_var.log", chrVec, chrMultiVariants);
    }
    return;
}

void SomaticVarCaller::extractSomaticData(
    const std::string &normalBamFile,
    const std::string &tumorBamFile,
    const std::string &fastaFile,
    const ParsingBamConfig &config,
    const std::vector<std::string> &chrVec,
    const std::map<std::string, int> &chrLength,
    std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants,
    std::map<Genome, VCF_Info> &vcfSet
){
    //Count each base numbers at tumor SNP position in the Normal.bam
    std::cerr<< "extracting data from normal BAM ... ";
    std::time_t begin = time(NULL);


    ParsingBamControl control;

    BamParserContext norBamCtx(normalBamFile, fastaFile, chrVec, chrLength, chrMultiVariants, vcfSet, Genome::NORMAL);
    ExtractNorDataBamParser normalBamParser(config, control, *chrPosNorBase);
    normalBamParser.parsingBam(norBamCtx);
    std::cerr<< difftime(time(NULL), begin) << "s\n";

    std::cerr << "extracting data from tumor BAM ... ";
    begin = time(NULL);
    BamParserContext tumBamCtx(tumorBamFile, fastaFile, chrVec, chrLength, chrMultiVariants, vcfSet, Genome::TUMOR);
    ExtractTumDataBamParser tumorBamParser(config, control, *chrPosSomaticInfo, *chrReadHpResultSet, *chrTumorPosReadCorrBaseHP);
    tumorBamParser.parsingBam(tumBamCtx);
    std::cerr<< difftime(time(NULL), begin) << "s\n";
}

double SomaticVarCaller::runTumorPurityEstimator(bool writeReadLog, const std::string resultPrefix){
    double tumorPurity = 0.0;

    TumorPurityEstimator* purityEstimator = new TumorPurityEstimator(chrVec, *chrPosNorBase, *chrPosSomaticInfo, writeReadLog, resultPrefix);
    // estimate tumor purity
    tumorPurity = purityEstimator->estimateTumorPurity();
    //flag the position that is used for purity estimation
    purityEstimator->markStatisticFlag(*chrPosSomaticInfo);
    delete purityEstimator;

    return tumorPurity;
}


void SomaticVarCaller::SetFilterParamsWithPurity(SomaticVarFilterParams &somaticParams, double &tumorPurity){

    somaticParams.tumorPurity = tumorPurity;
    // display the filter tier
    FilterTier filterTier = TIER_0_2;

    if(tumorPurity >= 0.9 && tumorPurity <= 1.0){ filterTier = TIER_1_0;}
    else if(tumorPurity >= 0.7 && tumorPurity < 0.9){ filterTier = TIER_0_8;}
    else if(tumorPurity >= 0.5 && tumorPurity < 0.7){ filterTier = TIER_0_6;}
    else if(tumorPurity >= 0.3 && tumorPurity < 0.5){ filterTier = TIER_0_4;}
    else{ filterTier = TIER_0_2;}

    switch(filterTier){
        // Tier 1.0: High purity tumors (0.9-1.0) - Most strict filtering
        case TIER_1_0:
            somaticParams.norVAF_maxThr = 0.13;
            somaticParams.norDepth_minThr = 1;

            somaticParams.MessyReadRatioThreshold = 1.0;    
            somaticParams.ReadCount_minThr = 3.0;

            somaticParams.HapConsistency_ReadCount_maxThr = 12.0;
            somaticParams.HapConsistency_VAF_maxThr = 0.144;
            somaticParams.HapConsistency_somaticRead_minThr = 0.0;

            somaticParams.IntervalSnpCount_ReadCount_maxThr = 12.0;
            somaticParams.IntervalSnpCount_VAF_maxThr = 0.189;
            somaticParams.IntervalSnpCount_minThr = 4.0;
            somaticParams.zScore_maxThr = 5.233;
            break;
        // Tier 0.8: Medium-high purity tumors (0.7-0.9)
        case TIER_0_8:
            somaticParams.norVAF_maxThr = 0.13;
            somaticParams.norDepth_minThr = 1;

            somaticParams.MessyReadRatioThreshold = 1.0;    
            somaticParams.ReadCount_minThr = 3.0;

            somaticParams.HapConsistency_ReadCount_maxThr = 10.0;
            somaticParams.HapConsistency_VAF_maxThr = 0.130;
            somaticParams.HapConsistency_somaticRead_minThr = 1.0;

            somaticParams.IntervalSnpCount_ReadCount_maxThr = 10.0;
            somaticParams.IntervalSnpCount_VAF_maxThr = 0.133;
            somaticParams.IntervalSnpCount_minThr = 4.0;
            somaticParams.zScore_maxThr = 2.676;
            break;
        // Tier 0.6: Medium purity tumors (0.5-0.7)
        case TIER_0_6:
            somaticParams.norVAF_maxThr = 0.105;
            somaticParams.norDepth_minThr = 1;

            somaticParams.MessyReadRatioThreshold = 1.0;    
            somaticParams.ReadCount_minThr = 1.0;

            somaticParams.HapConsistency_ReadCount_maxThr = 10.0;
            somaticParams.HapConsistency_VAF_maxThr = 0.071;
            somaticParams.HapConsistency_somaticRead_minThr = 0.0;

            somaticParams.IntervalSnpCount_ReadCount_maxThr = 10.0;
            somaticParams.IntervalSnpCount_VAF_maxThr = 0.105;
            somaticParams.IntervalSnpCount_minThr = 4.0;
            somaticParams.zScore_maxThr = 5.683;
            break;
        // Tier 0.4: Medium-low purity tumors (0.3-0.5)
        case TIER_0_4:
            somaticParams.norVAF_maxThr = 0.117;
            somaticParams.norDepth_minThr = 1;

            somaticParams.MessyReadRatioThreshold = 1.0;    
            somaticParams.ReadCount_minThr = 1.0;

            somaticParams.HapConsistency_ReadCount_maxThr = 8.0;
            somaticParams.HapConsistency_VAF_maxThr = 0.035;
            somaticParams.HapConsistency_somaticRead_minThr = 1.0;

            somaticParams.IntervalSnpCount_ReadCount_maxThr = 8.0;
            somaticParams.IntervalSnpCount_VAF_maxThr = 0.049;
            somaticParams.IntervalSnpCount_minThr = 4.0;
            somaticParams.zScore_maxThr = 3.043;
            break;
        // Tier 0.2: Low purity tumors (0.0-0.3) - Most lenient filtering
        case TIER_0_2:
            somaticParams.norVAF_maxThr = 0.130;
            somaticParams.norDepth_minThr = 1;

            somaticParams.MessyReadRatioThreshold = 1.0;    
            somaticParams.ReadCount_minThr = 1.0;

            somaticParams.HapConsistency_ReadCount_maxThr = 8.0;
            somaticParams.HapConsistency_VAF_maxThr = 0.020;
            somaticParams.HapConsistency_somaticRead_minThr = 1.0;

            somaticParams.IntervalSnpCount_ReadCount_maxThr = 8.0;
            somaticParams.IntervalSnpCount_VAF_maxThr = 0.025;
            somaticParams.IntervalSnpCount_minThr = 8.0;
            somaticParams.zScore_maxThr = 1.953;
            break;
        default:
            std::cerr << "[ERROR] Unexpected filter tier: " << filterTier << std::endl;
            exit(1);
    }

    if (tumorPurity <= 0 || tumorPurity > 1.0) {
        std::cerr << "[WARNING] tumor purity is not in the range of 0.0 to 1.0: " << tumorPurity << std::endl;
        std::cerr << "[WARNING] setting default parameters (tier " <<  FilterTierUtils::getTierValue(filterTier) << ")" << std::endl;
    }else{
        std::cerr << "setting filter params " << "(tier " << FilterTierUtils::getTierValue(filterTier) << ")" << " with tumor purity: " << tumorPurity << std::endl;
    }
}

void SomaticVarCaller::somaticFeatureFilter(const SomaticVarFilterParams &somaticParams, std::map<int, MultiGenomeVar> &currentChrVariants,const std::string &chr, std::map<int, SomaticData> &somaticPosInfo, double& tumorPurity){
    //calculate the information of the somatic positon 
    std::map<int, SomaticData>::iterator somaticVarIter = somaticPosInfo.begin();

    while( somaticVarIter != somaticPosInfo.end()){
                
        if(!currentChrVariants[(*somaticVarIter).first].isExists(TUMOR)){
            std::cerr << "[ERROR](somaticFeatureFilter) => can't find the position : chr:" << chr << " pos: " << (*somaticVarIter).first;
            exit(1);
        }

        auto curVar = currentChrVariants[(*somaticVarIter).first].Variant[TUMOR];

        // current variant is SNP
        if (curVar.variantType == VariantType::SNP) {
                
            //normal VAF filter parameter
            float norVAF_maxThr = somaticParams.norVAF_maxThr;
            int norDepth_minThr = somaticParams.norDepth_minThr;

            //messy read filter parameter
            float messyReadRatioThreshold = somaticParams.MessyReadRatioThreshold;
            int readCountThreshold = somaticParams.ReadCount_minThr;

            //haplotype consistency filter parameter
            float HapConsistency_VAF_maxThr=somaticParams.HapConsistency_VAF_maxThr;
            int HapConsistency_ReadCount_maxThr=somaticParams.HapConsistency_ReadCount_maxThr;
            int HapConsistency_somaticRead_minThr=somaticParams.HapConsistency_somaticRead_minThr;

            //interval snp count filter parameter
            float IntervalSnpCount_VAF_maxThr=somaticParams.IntervalSnpCount_VAF_maxThr;
            int IntervalSnpCount_ReadCount_maxThr=somaticParams.IntervalSnpCount_ReadCount_maxThr;
            int IntervalSnpCount_minThr=somaticParams.IntervalSnpCount_minThr;
            float zScore_maxThr = somaticParams.zScore_maxThr;

            //filter out the somatic SNP
            (*somaticVarIter).second.isFilterOut = false;
            
            // Check all filters conditions
            bool stage1_filtered = false;
            bool messy_read_filtered = false;
            bool read_count_filtered = false;

            
            //stage 1 filter
            float norVAF = (*chrPosNorBase)[chr][(*somaticVarIter).first].VAF;
            float norDepth = (*chrPosNorBase)[chr][(*somaticVarIter).first].depth;

            if (!(norVAF <= norVAF_maxThr && norDepth > norDepth_minThr)) {
                stage1_filtered = true;
            }

            //stage 2 filter
            if((*somaticVarIter).second.Mixed_HP_readRatio >= messyReadRatioThreshold){
                messy_read_filtered = true;
            }
            if((*somaticVarIter).second.CaseReadCount <= readCountThreshold){
                read_count_filtered = true;
            }
            
            // Haplotype consistency filter check
            bool haplotype_filtered = false;

            int somaticReadH1_1 = (*somaticVarIter).second.somaticReadHpCount[ReadHP::H1_1];
            int somaticReadH2_1 = (*somaticVarIter).second.somaticReadHpCount[ReadHP::H2_1];

            if ((*somaticVarIter).second.CaseReadCount <= HapConsistency_ReadCount_maxThr && 
                (*somaticVarIter).second.base.VAF <= HapConsistency_VAF_maxThr) {
                if (somaticReadH1_1 > HapConsistency_somaticRead_minThr && 
                    somaticReadH2_1 > HapConsistency_somaticRead_minThr) {
                    haplotype_filtered = true;
                }
            }
            
            // interval snp count filter check
            bool zscore_filtered = false;
            if ((*somaticVarIter).second.CaseReadCount <= IntervalSnpCount_ReadCount_maxThr && 
                (*somaticVarIter).second.base.VAF <= IntervalSnpCount_VAF_maxThr) {

                int intervalSnpCount = (*somaticVarIter).second.intervalSnpCount;
                float zScore = (*somaticVarIter).second.zScore;
                if (intervalSnpCount > IntervalSnpCount_minThr && zScore <= zScore_maxThr && zScore >= 0.0) {
                    zscore_filtered = true;
                }
            }

            // Whether the tag is filtered by any filter
            if (stage1_filtered || messy_read_filtered || read_count_filtered || 
                haplotype_filtered || zscore_filtered) {
                (*somaticVarIter).second.isFilterOut = true;
            }

            // If a filter needs to be applied, skip the filtered variant
            if (callerCfg.enableFilter && (*somaticVarIter).second.isFilterOut) {
                somaticVarIter++;
                continue;
            }

            (*somaticVarIter).second.isHighConSomaticSNP = true;
        }

        somaticVarIter++;          
    }
}

void SomaticVarCaller::calculateIntervalData(bool &isStartPos, int &startPos, int &endPos, DenseSnpInterval &denseSnp, std::map<int, std::pair<int, DenseSnpInterval>> &localDenseTumorSnpInterval){
    double mean = statisticsUtils::calculateMean(denseSnp.snpAltMean);
    double stdDev = statisticsUtils::calculateStandardDeviation(denseSnp.snpAltMean, mean);
    statisticsUtils::calculateZScores(denseSnp.snpAltMean, mean, stdDev, denseSnp.snpZscore);

    denseSnp.totalAltMean = mean;
    denseSnp.StdDev = stdDev;   

    localDenseTumorSnpInterval[startPos] = std::make_pair(endPos, denseSnp); // save the interval
}

void SomaticVarCaller::getDenseTumorSnpInterval(std::map<int, SomaticData> &somaticPosInfo, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, std::map<std::string, int>> &tumorPosReadCorrBaseHP, std::map<int, std::pair<int, DenseSnpInterval>> &localDenseTumorSnpInterval){
    //calculate the mean of HP3 reads at each tumor SNP position
    for(auto somaticPosIter = tumorPosReadCorrBaseHP.begin(); somaticPosIter != tumorPosReadCorrBaseHP.end(); somaticPosIter++){
        float readCount = 0.0;
        float altMean = 0.0;
        int minStartPos = INT_MAX;
        int maxEndPos = INT_MIN;

        for(auto readHPIter = somaticPosIter->second.begin(); readHPIter != somaticPosIter->second.end(); readHPIter++){
            int baseHP = readHPIter->second;
            if(baseHP != SnpHP::SOMATIC_H3) continue;
            readCount++;

            if(readHpResultSet.find(readHPIter->first) != readHpResultSet.end()){
                altMean += readHpResultSet[readHPIter->first].HP3;
                minStartPos = std::min(minStartPos, readHpResultSet[readHPIter->first].startPos);
                maxEndPos = std::max(maxEndPos, readHpResultSet[readHPIter->first].endPos);
            }else{
                std::cerr << "[ERROR](getDenseTumorSnpInterval) => readID not found in readHpResultSet: " << readHPIter->first << std::endl;
                exit(1);
            }
        }
        if(altMean != 0) altMean /= readCount;

        if(somaticPosInfo.find(somaticPosIter->first) != somaticPosInfo.end()){
            somaticPosInfo[somaticPosIter->first].MeanAltCountPerVarRead = altMean;
        }else{
            std::cerr << "[ERROR](getDenseTumorSnpInterval) => somaticPosInfo not found: " << somaticPosIter->first << std::endl;
            exit(1);
        }
    }

    //find the interval of tumor SNPs
    auto somaticPosIter = somaticPosInfo.begin();
    bool isRecordStartPos = false;
    int startPos = 0;
    int dense_distance = INTERVAL_SNP_MAX_DISTANCE;

    DenseSnpInterval denseSnp;

    while (somaticPosIter != somaticPosInfo.end()){
        int curPos = somaticPosIter->first;

        // ensure not out of bounds
        auto nextIter = std::next(somaticPosIter);
        if (nextIter != somaticPosInfo.end()){
            int nextPos = nextIter->first;

            int snpDistance = nextPos - curPos;

            if (snpDistance <= dense_distance) {
                //record the start position of the dense tumor interval
                if (!isRecordStartPos) {
                    isRecordStartPos = true;
                    startPos = curPos;
                    denseSnp.snpAltMean[curPos] = somaticPosIter->second.MeanAltCountPerVarRead;
                    denseSnp.minDistance[curPos] = snpDistance;
                    denseSnp.snpCount++;
                }

                // compare the distance between the previous and next position
                if(snpDistance < denseSnp.minDistance[curPos]){
                    denseSnp.minDistance[curPos] = snpDistance;
                }

                //record the next position in the dense tumor interval
                denseSnp.snpAltMean[nextPos] = nextIter->second.MeanAltCountPerVarRead;
                denseSnp.minDistance[nextPos] = snpDistance;

                denseSnp.snpCount++;
               
            } else {
                if (isRecordStartPos) {
                    //record the end position of the dense tumor interval
                    calculateIntervalData(isRecordStartPos, startPos, curPos, denseSnp, localDenseTumorSnpInterval);
                    isRecordStartPos = false;
                    startPos = 0;
                    denseSnp.clear();
                }
            }
        }

        somaticPosIter = nextIter;
    }

    // check if there is an unfinished interval
    if (isRecordStartPos) {
        int endPos = somaticPosInfo.rbegin()->first;
        if(endPos - startPos <= dense_distance){
            calculateIntervalData(isRecordStartPos, startPos, endPos, denseSnp, localDenseTumorSnpInterval);
        }
    }

    // record the zScore of each SNP in the dense tumor interval
    for (const auto& entry : localDenseTumorSnpInterval) {
        if(entry.second.second.snpCount > 1){
            for(auto snpIter: entry.second.second.snpZscore){
                int pos = snpIter.first;
                somaticPosInfo[pos].inDenseTumorInterval = true;
                somaticPosInfo[pos].zScore = abs(snpIter.second);
                somaticPosInfo[pos].intervalSnpCount = entry.second.second.snpCount;
            }

            for(auto minDistanceIter: entry.second.second.minDistance){
                int pos = minDistanceIter.first;
                somaticPosInfo[pos].minDistance = minDistanceIter.second;
            }
        }
    }
}

/**
 * @brief Calibrate read haplotype assignments
 * 
 * Removes low confidence H3/H4 SNP assignments from reads to improve
 * the quality of somatic variant calling results.
 * 
 * @param chr Chromosome name being processed
 * @param somaticPosInfo Somatic position information
 * @param readHpResultSet Read haplotype results to calibrate
 * @param tumorPosReadCorrBaseHP Position-read mapping for validation
 */
void SomaticVarCaller::calibrateReadHP(const std::string &chr, std::map<int, SomaticData> &somaticPosInfo, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, std::map<std::string, int>> &tumorPosReadCorrBaseHP){
        // Calibrate read HP to remove low confidence H3/H4 SNP
        std::map<int, SomaticData>::iterator somaticVarIter = somaticPosInfo.begin();
        while( somaticVarIter != somaticPosInfo.end()){
            if(!(*somaticVarIter).second.isHighConSomaticSNP){
                int pos = (*somaticVarIter).first;

                // Only reads with HP3 SNPs will be calibrated
                if(tumorPosReadCorrBaseHP.find(pos) != tumorPosReadCorrBaseHP.end()){

                    for(auto readIdBaseHP : tumorPosReadCorrBaseHP[pos]){
                        std::string readID = readIdBaseHP.first;
                        int baseHP = readIdBaseHP.second;

                        switch (baseHP) {
                            case SnpHP::GERMLINE_H1:
                                break;
                            case SnpHP::GERMLINE_H2:
                                break;
                            case SnpHP::SOMATIC_H3:
                                readHpResultSet[readID].HP3--; break;
                            default:
            break;
    }

                        if(readHpResultSet[readID].HP3 < 0){
                            std::cerr << "[ERROR](calibrate read HP) => read HP3 or HP4 SNP count < 0 :" << std::endl;
                            std::cerr << "readID: "<< readID << " chr: "<< chr << " pos: " << pos+1<< std::endl;
                            std::cerr << "HP3: "<< readHpResultSet[readID].HP3 << std::endl;
                            exit(1);
                        }
                    }
                }else{
                    std::cerr << "[ERROR](calibrate read HP) => can't find pos in tumorPosReadCorrBaseHP : chr: "<< chr << " pos: " << pos+1 <<std::endl;
                    exit(1);
                }
            }
            somaticVarIter++;
        }
}

/**
 * @brief Calculate read set haplotype results
 * 
 * Determines final haplotype assignments for all reads based on
 * variant patterns and quality thresholds.
 * 
 * @param chr Chromosome name being processed
 * @param readHpResultSet Read haplotype results to update
 * @param tumorPosReadCorrBaseHP Position-read mapping for validation
 * @param percentageThreshold Percentage threshold for haplotype determination
 */
void SomaticVarCaller::calculateReadSetHP(const std::string &chr, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, std::map<std::string, int>> &tumorPosReadCorrBaseHP, const double& percentageThreshold){
        std::map<std::string, ReadVarHpCount>::iterator readTotalHPcountIter = readHpResultSet.begin();
        // Calculate all read HP result 
        while(readTotalHPcountIter != readHpResultSet.end()){
            std::string readID = (*readTotalHPcountIter).first;
            
            std::map<int, int> hpCount;
            hpCount[1] = (*readTotalHPcountIter).second.HP1;
            hpCount[2] = (*readTotalHPcountIter).second.HP2;
            hpCount[3] = (*readTotalHPcountIter).second.HP3;
            hpCount[4] = (*readTotalHPcountIter).second.HP4;
            int pqValue = 0;
            
            double normalHPsimilarity = 0.0;
            double tumorHPsimilarity = 0.0;
            std::map<int, int> norCountPS = (*readTotalHPcountIter).second.norCountPS;

            (*readTotalHPcountIter).second.hpResult = somaticJudger.judgeSomaticReadHap(hpCount, pqValue, norCountPS, normalHPsimilarity, tumorHPsimilarity, percentageThreshold, nullptr, nullptr, nullptr);
            
            readTotalHPcountIter++;
        }
}

void SomaticVarCaller::statisticSomaticPosReadHP(
    const std::string &chr, std::map<int, SomaticData> &somaticPosInfo,
    std::map<int, std::map<std::string, int>> &tumorPosReadCorrBaseHP,
    std::map<std::string, ReadVarHpCount> &readHpResultSet,
    chrReadHpResult &localReadHpDistri
){
    std::map<int, SomaticData>::iterator somaticVarIter = somaticPosInfo.begin();
    while(somaticVarIter != somaticPosInfo.end()){
        if((*somaticVarIter).second.isHighConSomaticSNP){
            int pos = (*somaticVarIter).first;

            if(tumorPosReadCorrBaseHP.find(pos) != tumorPosReadCorrBaseHP.end()){
                localReadHpDistri.posReadHpResult[pos] = ReadHpResult();

                // record the number of HP1-1 or HP2-1 derived from Base HP3
                std::map<int, int> deriveByHPfromBaseHp3;
                deriveByHPfromBaseHp3[ReadHP::H1_1] = 0;
                deriveByHPfromBaseHp3[ReadHP::H2_1] = 0;

                // record read hp in current position and calculate somatic read derive by HP1 or HP2
                for(std::pair<std::string, int> readIdBaseHP : tumorPosReadCorrBaseHP[pos]){
                    std::string readID = readIdBaseHP.first;
                    int baseHP = readIdBaseHP.second;
                    int hpResult = readHpResultSet[readID].hpResult;

                    //record read HP
                    localReadHpDistri.recordReadHp(pos, hpResult, baseHP);

                    if(baseHP == SnpHP::SOMATIC_H3){
                        switch (hpResult) {
                            case ReadHP::H1_1:
                                deriveByHPfromBaseHp3[ReadHP::H1_1]++; break;
                            case ReadHP::H2_1:
                                deriveByHPfromBaseHp3[ReadHP::H2_1]++; break;
                            default:
                                break;
                        }
                    }
                }

                //calculate HP1-1 and HP2-1 ratio
                int totalH3readHadGermlineSnp = deriveByHPfromBaseHp3[ReadHP::H1_1] + deriveByHPfromBaseHp3[ReadHP::H2_1];
                
                float HP1_1ratio = 0.0;
                float HP2_1ratio = 0.0;

                if(totalH3readHadGermlineSnp > 0){
                    if(deriveByHPfromBaseHp3[ReadHP::H1_1] > 0){
                        HP1_1ratio = (float)deriveByHPfromBaseHp3[ReadHP::H1_1] / (float)totalH3readHadGermlineSnp;
                    }
                    if(deriveByHPfromBaseHp3[ReadHP::H2_1] > 0){
                        HP2_1ratio = (float)deriveByHPfromBaseHp3[ReadHP::H2_1] / (float)totalH3readHadGermlineSnp;
                    }
                }

                localReadHpDistri.posReadHpResult[pos].existDeriveByH1andH2 = false;

                if(HP1_1ratio >= 1.0){
                    (*somaticVarIter).second.somaticReadDeriveByHP = SnpHP::GERMLINE_H1;
                }else if(HP2_1ratio >= 1.0){
                    (*somaticVarIter).second.somaticReadDeriveByHP = SnpHP::GERMLINE_H2;
                }else{
                    (*somaticVarIter).second.somaticReadDeriveByHP = SnpHP::NONE_SNP;

                    if((0 < HP1_1ratio && HP1_1ratio < 1.0)  || (0 < HP2_1ratio && HP2_1ratio < 1.0)){
                        localReadHpDistri.posReadHpResult[pos].existDeriveByH1andH2 = true;
                    }
                }
                //record derive HP of current snp
                localReadHpDistri.recordDeriveHp(pos, (*somaticVarIter).second.somaticReadDeriveByHP, 0.0);
            }else{
                std::cerr << "[ERROR](statistic all read HP) => can't find pos in tumorPosReadCorrBaseHP : chr: "<< chr << " pos: " << pos+1 <<std::endl;
                exit(1); 
            }

        }
        somaticVarIter++;
    }
}

void SomaticVarCaller::findOtherSomaticSnpHP(const std::string &chr, std::map<int, SomaticData> &somaticPosInfo, std::map<int, MultiGenomeVar> &currentChrVariants){
    std::map<int, SomaticData>::iterator somaticVarIter = somaticPosInfo.begin();
    while(somaticVarIter != somaticPosInfo.end()){
        if((*somaticVarIter).second.isHighConSomaticSNP){
            int pos = (*somaticVarIter).first;

            std::map<int, int> NucCount;
            NucCount[Nitrogenous::A] = (*somaticVarIter).second.base.MPQ_A_count;
            NucCount[Nitrogenous::C] = (*somaticVarIter).second.base.MPQ_C_count;
            NucCount[Nitrogenous::T] = (*somaticVarIter).second.base.MPQ_T_count;
            NucCount[Nitrogenous::G] = (*somaticVarIter).second.base.MPQ_G_count;

            int refAllele;
            int altAllele;

            if(currentChrVariants.find(pos) != currentChrVariants.end()){
                refAllele = NucUtil::convertStrNucToInt(currentChrVariants[pos].Variant[TUMOR].allele.Ref);
                altAllele = NucUtil::convertStrNucToInt(currentChrVariants[pos].Variant[TUMOR].allele.Alt);
            }else{
                std::cerr << "[ERROR](FindOtherSomaticSnpHP) => can't find position in currentChrVariants : chr:" << chr << " pos: " << pos + 1;
                exit(1);
            }

            NucCount.erase(refAllele);
            NucCount.erase(altAllele);

            int maxNuc = Nitrogenous::UNKOWN;
            int maxNucCount = 0;

            int minNuc = Nitrogenous::UNKOWN;
            int minNucCount = 0;
            for(auto nuc : NucCount){
                if(nuc.second > 1){
                    if(nuc.second > maxNucCount){
                        minNuc = maxNuc;
                        minNucCount = maxNucCount;
                        maxNuc = nuc.first;
                        maxNucCount = nuc.second;
                    }
                }
            }

            if(maxNuc != Nitrogenous::UNKOWN){
                (*somaticVarIter).second.somaticHp4Base = maxNuc;
                (*somaticVarIter).second.somaticHp4BaseCount = maxNucCount;
            }

            if(minNuc != Nitrogenous::UNKOWN){
                (*somaticVarIter).second.somaticHp5Base = minNuc;
                (*somaticVarIter).second.somaticHp5BaseCount = minNucCount;
            }
        }
        somaticVarIter++;
    }
}

void SomaticVarCaller::writeSomaticVarCallingLog(const CallerContext &ctx, const SomaticVarFilterParams &somaticParams, const std::vector<std::string> &chrVec, std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants){
    std::ofstream *tagHP3Log = new std::ofstream(bamCfg.resultPrefix+"_somatic_var.out");

    if(!tagHP3Log->is_open()){
        std::cerr<< "Fail to open write file: " << bamCfg.resultPrefix+"_somatic_var.out" << "\n";
        exit(1);
    }

    std::cerr << "writing somatic variants calling log ... ";
    std::time_t begin = time(NULL);

    int totalSomaticSNP = 0;
    for(auto chr : chrVec){
        std::map<int, SomaticData>::iterator somaticVarIter = (*chrPosSomaticInfo)[chr].begin();
        while( somaticVarIter != (*chrPosSomaticInfo)[chr].end()){
            if((*somaticVarIter).second.isHighConSomaticSNP){
                totalSomaticSNP++;
            }
            somaticVarIter++;
        }
    }

    //write header
    (*tagHP3Log) << "####################################\n"
                 << "#   Somatic Variants Calling Log   #\n"
                 << "####################################\n";
    (*tagHP3Log) << "##normalSnpFile:"       << ctx.normalSnpFile                    << "\n"
                 << "##tumorSnpFile:"        << ctx.tumorSnpFile                     << "\n" 
                 << "##bamFile:"             << ctx.normalBamFile                    << "\n"
                 << "##tumorBamFile:"        << ctx.tumorBamFile                     << "\n" 
                 << "##resultPrefix:"        << bamCfg.resultPrefix        << "\n"
                 << "##numThreads:"          << bamCfg.numThreads          << "\n"
                 << "##region:"              << bamCfg.region              << "\n"
                 << "##qualityThreshold:"    << bamCfg.qualityThreshold    << "\n"
                 << "##percentageThreshold:" << bamCfg.percentageThreshold << "\n"
                 << "##tagSupplementary:"    << bamCfg.tagSupplementary    << "\n"
                 << "##\n";

    //write filter parameter
    (*tagHP3Log) << "##======== Filter Parameters =========\n"
                 << "##Enable filter : " << callerCfg.enableFilter << "\n"
                 << "##Calling mapping quality :" << bamCfg.qualityThreshold << "\n"
                 << "##Tumor purity : " << somaticParams.tumorPurity << "\n"
                 << "##Normal VAF maximum threshold : " << somaticParams.norVAF_maxThr << "\n"
                 << "##Normal depth minimum threshold : " << somaticParams.norDepth_minThr << "\n"
                 << "##Messy read ratio threshold : " << somaticParams.MessyReadRatioThreshold << "\n"
                 << "##Somatic read count minimum threshold : " << somaticParams.ReadCount_minThr << "\n"
                 << "##Haplotag consistency filter VAF threshold : " << somaticParams.HapConsistency_VAF_maxThr << "\n"
                 << "##Haplotag consistency filter read count threshold : " << somaticParams.HapConsistency_ReadCount_maxThr << "\n"
                 << "##Haplotag consistency somatic read count minimum threshold : " << somaticParams.HapConsistency_somaticRead_minThr << "\n"
                 << "##Interval SNP count filter threshold : " << somaticParams.IntervalSnpCount_VAF_maxThr << "\n"
                 << "##Interval SNP count filter read count threshold : " << somaticParams.IntervalSnpCount_ReadCount_maxThr << "\n"
                 << "##Interval SNP count minimum threshold : " << somaticParams.IntervalSnpCount_minThr << "\n"
                 << "##Z-score maximum threshold : " << somaticParams.zScore_maxThr << "\n"
                 << "##==================================== \n"
                 << "##\n"
                 << "##Total Somatic SNPs: " << totalSomaticSNP << "\n"
                 << "##\n"; 
    (*tagHP3Log) << "#CHROM\t" 
                 << "POS\t"
                 << "ID\t"
                 << "REF\t" 
                 << "ALT\t"
                 << "AltCount\t"
                 << "ReadCount\t" 
                 << "NorAltCount\t"
                 << "PureH1-1\t"
                 << "PureH2-1\t"
                 << "PureH3\t"
                 << "MixedHpRead\t"
                 << "UnTag\t" 
                 << "PureH1-1ratio\t" 
                 << "PureH2-1ratio\t" 
                 << "PureH3ratio\t" 
                 << "MixedHpReadRatio\t" 
                 << "NorVAF\t"  
                 << "TumVAF\t"
                 << "NorMpqVAF\t"
                 << "TumMpqVAF\t"
                 << "NorVAF_substract\t"
                 << "TumVAF_substract\t"
                 << "NorDepth\t"  
                 << "TumDepth\t"
                 << "Subtract_Depth\t" 
                 << "NorDeletionCount\t"
                 << "TumDeletionCount\t"
                 << "NorDeletionRatio\t"
                 << "TumDeletionRatio\t"
                 << "NorMpqReadRatio\t"
                 << "TumMpqReadRatio\t"  
                 << "ShannonEntropy\t"
                 << "HomopolymerLength\t"
                 << "H1readCount\t"
                 << "H2readCount\t"
                 << "H1_1readCount\t"
                 << "H2_1readCount\t"
                 << "H3readCount\t"
                 << "GermlineReadHpCount\t"
                 << "GermlineReadHpImbalanceRatio\t"
                 << "SomaticReadHpImbalanceRatio\t"
                 << "BaseGermlineReadHpImbalanceRatio\t"
                 << "PercentageOfGermlineHp\t"
                 << "H1readCountInNorBam\t"
                 << "H2readCountInNorBam\t"
                 << "GermlineReadHpCountInNorBam\t"
                 << "GermlineReadHpImbalanceRatioInNorBam\t"
                 << "PercentageOfGermlineHpInNorBam\t"
                 << "GermlineReadHpImbalanceRatioDifference\t"
                 << "PercentageOfGermlineHpDifference\t"
                 << "SomaticRead_H1-1\t"
                 << "SomaticRead_H2-1\t"
                 << "SomaticRead_H3\t"
                 << "SomaticRead_unTag\t"
                 << "AltMeanCountPerVarRead\t"
                 << "zScore\t"
                 << "IntervalSnpCount\t"
                 << "IntervalMinDistance\t"
                 << "ExistNorSnp\t"
                 << "StatisticPurity\t"
                 << "isFilterOut\t"
                 << "NorNonDelAF\t"
                 << "TumNonDelAF\t"
                 << "GT\n";

    //write variants information
    for(auto chr : chrVec){
    
        std::map<int, SomaticData>::iterator somaticVarIter = (*chrPosSomaticInfo)[chr].begin();
        while( somaticVarIter != (*chrPosSomaticInfo)[chr].end()){
            
            if((*somaticVarIter).second.isHighConSomaticSNP == false){
                somaticVarIter++;
                continue;
            }

            //calculating Case Ratio
            int pure_H1_1_read = (*somaticVarIter).second.pure_H1_1_read;
            int pure_H2_1_read = (*somaticVarIter).second.pure_H2_1_read;
            int pure_H3_read = (*somaticVarIter).second.pure_H3_read;
            int Mixed_HP_read = (*somaticVarIter).second.Mixed_HP_read;
            int unTagCount = (*somaticVarIter).second.unTag;
            int readCount = (*somaticVarIter).second.CaseReadCount;

            float mixed_HP_readRatio = (*somaticVarIter).second.Mixed_HP_readRatio;

            float pure_H1_1_readRatio = (*somaticVarIter).second.pure_H1_1_readRatio;
            float pure_H2_1_readRatio = (*somaticVarIter).second.pure_H2_1_readRatio;
            float pure_H3_readRatio = (*somaticVarIter).second.pure_H3_readRatio;

            int norDepth = (*chrPosNorBase)[chr][(*somaticVarIter).first].depth;
            // int norMpqDepth = (*chrPosNorBase)[chr][(*somaticVarIter).first].filteredMpqDepth;

            int tumDepth = (*somaticVarIter).second.base.depth;
            // int tumMpqDepth = (*somaticVarIter).second.base.filteredMpqDepth;
           
            //calculate the Subtract in depth between normal and tumor
            int subtractDepth = tumDepth - norDepth;


            std::string RefBase;
            std::string AltBase;
                
            if(chrMultiVariants[chr][(*somaticVarIter).first].isExists(TUMOR) == true){
                RefBase = chrMultiVariants[chr][(*somaticVarIter).first].Variant[TUMOR].allele.Ref;
                AltBase = chrMultiVariants[chr][(*somaticVarIter).first].Variant[TUMOR].allele.Alt;
            }else{
                std::cerr << "[ERROR](write tag HP3 log file) => can't find the position : chr:" << chr << " pos: " << ((*somaticVarIter).first) + 1;
                exit(1);
            }

            if(RefBase == "" || AltBase == ""){
                std::cerr << "[ERROR](write tag HP3 log file) => can't find RefBase or AltBase : chr:" << chr << " pos: " << ((*somaticVarIter).first) + 1 << " RefBase:" << RefBase << " AltBase:" << AltBase;
                exit(1);
            }

            //tmp
            int tumMpqAltCount = (*somaticVarIter).second.base.getMpqBaseCount(AltBase);
            int norAltCount = (*chrPosNorBase)[chr][(*somaticVarIter).first].getBaseCount(AltBase);

            //calculating tumor VAF
            float tumVAF = (*somaticVarIter).second.base.VAF;
            float tumMpqVAF = (*somaticVarIter).second.base.filteredMpqVAF;

            float norVAF = (*chrPosNorBase)[chr][(*somaticVarIter).first].VAF;
            float norMpqVAF = (*chrPosNorBase)[chr][(*somaticVarIter).first].filteredMpqVAF;

            float norVAF_substract = (norMpqVAF - norVAF);

            float norMPQReadRatio = (*chrPosNorBase)[chr][(*somaticVarIter).first].lowMpqReadRatio; 

            // Calculating the difference in VAF between VAF and filtered low MPQ VAF
            float tumVAF_substract = (tumMpqVAF -tumVAF);
            float tumLowMpqReadRatio = (*somaticVarIter).second.base.lowMpqReadRatio;
            
            float norNonDelAF = (*chrPosNorBase)[chr][(*somaticVarIter).first].nonDelVAF;
            float tumNonDelAF = (*somaticVarIter).second.base.nonDelVAF;
                
            //current SNP GT type
            std::string GTtype = (*somaticVarIter).second.GTtype;

            if(GTtype == ""){
                std::cerr << "[ERROR](GTtype) => chr: " << chr << " pos: " << ((*somaticVarIter).first) + 1 << " GTtype: " << GTtype;
                exit(1);                
            }

            //The count of deletions occurring at the SNP position
            int norDeletionCount = (*chrPosNorBase)[chr][(*somaticVarIter).first].delCount;

            int tumDeletionCount = (*somaticVarIter).second.base.delCount;

            //the ratio of deletions occurring
            float norDeletionRatio = (*chrPosNorBase)[chr][(*somaticVarIter).first].delRatio;

            float tumDeletionRatio = (*somaticVarIter).second.base.delRatio;

            // the count of reads have somatic base at current position
            int somaticVarReadH1_1 = (*somaticVarIter).second.somaticReadHpCount[ReadHP::H1_1];
            int somaticVarReadH2_1 = (*somaticVarIter).second.somaticReadHpCount[ReadHP::H2_1];
            int somaticVarReadH3 = (*somaticVarIter).second.somaticReadHpCount[ReadHP::H3];
            int untaggedReadCount = (*somaticVarIter).second.somaticReadHpCount[ReadHP::unTag];

            // the distribution of reads HP at the current position
            int H1readCount = (*somaticVarIter).second.base.ReadHpCount[ReadHP::H1];
            int H2readCount = (*somaticVarIter).second.base.ReadHpCount[ReadHP::H2];
            int H1_1readCount = (*somaticVarIter).second.base.ReadHpCount[ReadHP::H1_1];
            int H2_1readCount = (*somaticVarIter).second.base.ReadHpCount[ReadHP::H2_1];
            int H3readCount = (*somaticVarIter).second.base.ReadHpCount[ReadHP::H3];
              
            // the ratio of the read count of H1 and H2
            int germlineReadHpCount = H1readCount + H2readCount;
            double germlineReadHpImbalanceRatio = (*somaticVarIter).second.base.germlineHaplotypeImbalanceRatio;

            // // the ratio of the germlineHP based on depth   
            double percentageOfGermlineHp = (*somaticVarIter).second.base.percentageOfGermlineHp;

            // // the ratio of the read count of H1 and H2 based on germline HP
            double baseOnGermlineReadHpImbalanceRatio = (*somaticVarIter).second.allelicImbalanceRatio;

            double somaticReadHpImbalanceRatio = (*somaticVarIter).second.somaticHaplotypeImbalanceRatio;

            //read hp count in the normal bam
            int H1readCountInNorBam = (*chrPosNorBase)[chr][(*somaticVarIter).first].ReadHpCount[ReadHP::H1];
            int H2readCountInNorBam = (*chrPosNorBase)[chr][(*somaticVarIter).first].ReadHpCount[ReadHP::H2];
            int germlineReadHpCountInNorBam = H1readCountInNorBam + H2readCountInNorBam;

            double germlineReadHpImbalanceRatioInNorBam = (*chrPosNorBase)[chr][(*somaticVarIter).first].germlineHaplotypeImbalanceRatio;

            double percentageOfGermlineHpInNorBam = (*chrPosNorBase)[chr][(*somaticVarIter).first].percentageOfGermlineHp;

            //difference of the germline consistency ratio in tumor and normal bam
            double germlineReadHpImbalanceRatioDifference = 0.0;
            germlineReadHpImbalanceRatioDifference = germlineReadHpImbalanceRatio - germlineReadHpImbalanceRatioInNorBam;
            
            //difference of the percentage of germline HP in tumor and normal bam
            double percentageOfGermlineHpDifference = 0.0;
            percentageOfGermlineHpDifference = percentageOfGermlineHp - percentageOfGermlineHpInNorBam;

            // z-score of the interval snp filter 
            double zScore = -1.0;
            if((*somaticVarIter).second.inDenseTumorInterval){
                if((*somaticVarIter).second.zScore < 0.0){
                    std::cerr << "[ERROR](zScore) => chr: " << chr << " pos: " << ((*somaticVarIter).first) + 1 << " zScore: " << (*somaticVarIter).second.zScore << std::endl;
                    exit(1);
                }else{
                    zScore = (*somaticVarIter).second.zScore;
                }
            }
            
            //HP3 SNP position (1-base)
            int HP3pos = (*somaticVarIter).first + 1;
            
            //write information
            (*tagHP3Log)<< chr << " \t"  //1 
                        << HP3pos << "\t"  //2
                        << "." << "\t"  //3
                        << RefBase << "\t" //4
                        << AltBase << "\t" //5
                        << tumMpqAltCount << "\t" //6
                        << readCount << "\t\t"  //7
                        << norAltCount << "\t" //8
                        << pure_H1_1_read << "\t" //9
                        << pure_H2_1_read << "\t" //10
                        << pure_H3_read << "\t" //11
                        << Mixed_HP_read << "\t" //12
                        << unTagCount << "\t\t"  //13
                        << pure_H1_1_readRatio <<"\t"  //14
                        << pure_H2_1_readRatio << "\t"  //15
                        << pure_H3_readRatio << "\t"  //16
                        << mixed_HP_readRatio << "\t\t"  //17
                        << norVAF << "\t"  //18
                        << tumVAF << "\t\t"  //19   
                        << norMpqVAF << "\t" //20
                        << tumMpqVAF << "\t\t" //21
                        << norVAF_substract << "\t" //22
                        << tumVAF_substract << "\t\t" //23
                        << norDepth << "\t"  //24
                        << tumDepth << "\t"  //25
                        << subtractDepth << "\t"  //26 
                        << norDeletionCount << "\t"  //27
                        << tumDeletionCount << "\t"  //28
                        << norDeletionRatio << "\t"  //29
                        << tumDeletionRatio << "\t"  //30
                        << norMPQReadRatio << "\t" //31
                        << tumLowMpqReadRatio << "\t" //32
                        << (*somaticVarIter).second.shannonEntropy << "\t" //33
                        << (*somaticVarIter).second.homopolymerLength << "\t\t" //34
                        << H1readCount << "\t" //35
                        << H2readCount << "\t" //36
                        << H1_1readCount << "\t" //37
                        << H2_1readCount << "\t" //38
                        << H3readCount << "\t" //39
                        << germlineReadHpCount << "\t" //40
                        << germlineReadHpImbalanceRatio << "\t" //41
                        << somaticReadHpImbalanceRatio << "\t" //42
                        << baseOnGermlineReadHpImbalanceRatio << "\t" //43
                        << percentageOfGermlineHp << "\t" //44
                        << H1readCountInNorBam << "\t" //45
                        << H2readCountInNorBam << "\t" //46
                        << germlineReadHpCountInNorBam << "\t" //47
                        << germlineReadHpImbalanceRatioInNorBam << "\t" //48
                        << percentageOfGermlineHpInNorBam << "\t" //49
                        << germlineReadHpImbalanceRatioDifference << "\t" //50
                        << percentageOfGermlineHpDifference << "\t" //51
                        << somaticVarReadH1_1 << "\t" //52
                        << somaticVarReadH2_1 << "\t" //53
                        << somaticVarReadH3 << "\t" //54
                        << untaggedReadCount << "\t" //55
                        << (*somaticVarIter).second.MeanAltCountPerVarRead << "\t" //56
                        << zScore << "\t" //57
                        << (*somaticVarIter).second.intervalSnpCount << "\t" //58
                        << (*somaticVarIter).second.minDistance << "\t" //59
                        << chrMultiVariants[chr][(*somaticVarIter).first].isExists(NORMAL) << "\t" //60
                        << (*somaticVarIter).second.statisticPurity << "\t" //61
                        << (*somaticVarIter).second.isFilterOut << "\t" //62
                        << norNonDelAF << "\t" //63
                        << tumNonDelAF << "\t" //64
                        << GTtype <<"\n";  //65
                        
            somaticVarIter++;          
        }
    }

    tagHP3Log->close();
    delete tagHP3Log;
    tagHP3Log = nullptr;
    std::cerr<< difftime(time(NULL), begin) << "s\n";  
}

void SomaticVarCaller::writeOtherSomaticHpLog(const std::string logFileName, const std::vector<std::string> &chrVec, std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants){
    std::ofstream *OtherHpSomaticVarLog=NULL;
    OtherHpSomaticVarLog=new std::ofstream(logFileName);

    int totalOtherSomaticHpVar = 0;
    for(auto chr: chrVec){
        for(auto somaticVar: (*chrPosSomaticInfo)[chr]){
            if(somaticVar.second.somaticHp4Base != Nitrogenous::UNKOWN){
                totalOtherSomaticHpVar++;
            }
        }
    }

    if(!OtherHpSomaticVarLog->is_open()){
        std::cerr<< "Fail to open write file: " << logFileName << "\n";
        exit(1);
    }else{
        (*OtherHpSomaticVarLog) << "################################\n";
        (*OtherHpSomaticVarLog) << "# Other Somatic HP Variant Log #\n";
        (*OtherHpSomaticVarLog) << "################################\n";
        (*OtherHpSomaticVarLog) << "##Tatal other somatic HP variants:"  << totalOtherSomaticHpVar << "\n";
        (*OtherHpSomaticVarLog) << "#CHROM\t"
                                << "POS\t"
                                << "REF\t"
                                << "ALT\t"
                                << "HP4Base\t"
                                << "HP4BaseCount\t"
                                << "HP5Base\t"
                                << "HP5BaseCount\n";
    }

    for(auto chr: chrVec){
        auto currentChrVar = chrMultiVariants[chr];
        for(auto somaticVar: (*chrPosSomaticInfo)[chr]){
            int pos = somaticVar.first;

            if(somaticVar.second.somaticHp4Base == Nitrogenous::UNKOWN){
                continue;
            }
            if(currentChrVar.find(pos) == currentChrVar.end()){
                std::cerr << "[ERROR](writeOtherSomaticHpLog) => can't find the position : chr:" << chr << " pos: " << pos + 1;
                exit(1);
            }

            (*OtherHpSomaticVarLog) << chr << "\t"
                                    << pos + 1 << "\t"
                                    << currentChrVar[pos].Variant[TUMOR].allele.Ref << "\t"
                                    << currentChrVar[pos].Variant[TUMOR].allele.Alt << "\t"
                                    << NucUtil::convertIntNucToStr(somaticVar.second.somaticHp4Base) << "\t"
                                    << somaticVar.second.somaticHp4BaseCount<< "\t"
                                    << NucUtil::convertIntNucToStr(somaticVar.second.somaticHp5Base) << "\t"
                                    << somaticVar.second.somaticHp5BaseCount << "\n";
        }
    
    
    }
    (*OtherHpSomaticVarLog).close();
    delete OtherHpSomaticVarLog;
    OtherHpSomaticVarLog = nullptr;
}

void SomaticVarCaller::writeDenseTumorSnpIntervalLog(const std::string logFileName, const std::vector<std::string> &chrVec){
    std::ofstream *closeSomaticSnpIntervalLog=NULL;
    closeSomaticSnpIntervalLog=new std::ofstream(logFileName);

    int totalIntervalCount = 0;
    for(auto chr: chrVec){
        totalIntervalCount += (*denseTumorSnpInterval)[chr].size();
    }

    if(!closeSomaticSnpIntervalLog->is_open()){
        std::cerr<< "Fail to open write file: " << logFileName << "\n";
        exit(1);
    }else{
        (*closeSomaticSnpIntervalLog) << "################################\n";
        (*closeSomaticSnpIntervalLog) << "# Dense Tumor SNP Interval Log #\n";
        (*closeSomaticSnpIntervalLog) << "################################\n";
        (*closeSomaticSnpIntervalLog) << "##Tatal intervals:"  << totalIntervalCount << "\n";
        (*closeSomaticSnpIntervalLog) << "#CHROM\t"
                                      << "startPos-endPos\t"
                                      << "snpCount\t"
                                      << "totalAltMean\t"
                                      << "stdDev\t"
                                      << "zScore\n"; //standard deviation
    }

    for(auto chr: chrVec){
        for(auto intervalIter: (*denseTumorSnpInterval)[chr]){
            int startPos = intervalIter.first;
            int endPos = intervalIter.second.first;
            int snpCount = intervalIter.second.second.snpCount;

            (*closeSomaticSnpIntervalLog) << chr << ":"
                                    << startPos + 1 << "-"
                                    << endPos + 1 << "\t"
                                    << snpCount << "\t"
                                    << intervalIter.second.second.totalAltMean << "\t"
                                    << intervalIter.second.second.StdDev << "\n";
            for(auto snpIter: intervalIter.second.second.snpAltMean){
                double zScore = intervalIter.second .second.snpZscore[snpIter.first];
                int minDistance = intervalIter.second.second.minDistance[snpIter.first];
                (*closeSomaticSnpIntervalLog) <<"#snp:altMean:zScore:minDistance=>  " << snpIter.first + 1 << " : " << snpIter.second << " : " << zScore << " : " << minDistance << "\n";
            }
            (*closeSomaticSnpIntervalLog) << "#\n";
            
        }
    }
    (*closeSomaticSnpIntervalLog).close();
    delete closeSomaticSnpIntervalLog;
    
    closeSomaticSnpIntervalLog = nullptr;
}

/**
 * @brief Get somatic variant flags
 * 
 * Marks variants as somatic based on calling results and sets appropriate flags
 * for downstream analysis and reporting.
 * 
 * @param chrVec Vector of chromosome names to process
 * @param chrMultiVariants Multi-genome chromosome variants to update
 */
void SomaticVarCaller::getSomaticFlag(const std::vector<std::string> &chrVec, std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants){
    int somaticVariantCount = 0;
    for(auto chr: chrVec){
        for(auto somaticVar: (*chrPosSomaticInfo)[chr]){
            // [debug] logic for somatic tagging validation
            // if(!somaticVar.second.isFilterOut && somaticVar.second.isHighConSomaticSNP){
            if(somaticVar.second.isHighConSomaticSNP){
                chrMultiVariants[chr][somaticVar.first].isSomaticVariant = true;
                chrMultiVariants[chr][somaticVar.first].somaticReadDeriveByHP = somaticVar.second.somaticReadDeriveByHP;
                somaticVariantCount++;
            }
        }
    }
    // [debug] somatic variant count
    std::cerr << "somatic variant count(Flag): " << somaticVariantCount << "\n";
}

/**
 * @brief Display calling SNP count
 * 
 * Prints the number of somatic SNPs called for debugging and validation purposes.
 * This function provides a summary of the calling results.
 */
void SomaticVarCaller::displayCallingSnpCount(){
    int somaticSNP = 0;

    for(auto chr : chrVec){
        auto somaticPosIter = (*chrPosSomaticInfo)[chr].begin();
        while(somaticPosIter != (*chrPosSomaticInfo)[chr].end()){
            if((*somaticPosIter).second.isHighConSomaticSNP){
                somaticSNP++;
            }
            somaticPosIter++;
        }
    }
    std::cerr << "calling somatic SNPs: " << somaticSNP << std::endl;
}