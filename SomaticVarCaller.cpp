#include "SomaticVarCaller.h"


ExtractNorDataBamParser::ExtractNorDataBamParser(std::map<std::string, std::map<int, PosBase>>& chrPosNorBase)
: chrPosNorBase(chrPosNorBase){}

ExtractNorDataBamParser::~ExtractNorDataBamParser(){
};

void ExtractNorDataBamParser::displayPosInfo(std::string chr, int pos){
    if(chrPosNorBase[chr].find(pos) == chrPosNorBase[chr].end()){
        std::cerr << "ERROR (displayPosInfo) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
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
    const std::string &chrName,
    const HaplotagParameters &params,
    const Genome& genmoeType,
    std::map<int, MultiGenomeVar> &currentVariants,
    std::map<int, MultiGenomeVar>::iterator &firstVariantIter,
    std::map<Genome, VCF_Info> &vcfSet,
    const std::string &ref_string
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
    CigarParser* cigarParser = new ExtractNorDataCigarParser(*variantBase, tumVarPosVec, ref_pos, query_pos);
    cigarParser->parsingCigar(aln, bamHdr, chrName, params, firstVariantIter, currentVariants, ref_string, hpCount, variantsHP, norCountPS);
    delete cigarParser;

    // get the number of SVs occurring on different haplotypes in a read
    if( aln.core.qual >= params.somaticCallingMpqThreshold ){
        germlineJudgeSVHap(aln, vcfSet, hpCount, genmoeType);
    }

    double min = 0.0;
    double max = 0.0;
    int hpResult = ReadHP::unTag;
    int pqValue = 0;
    int psValue = 0;
    double percentageThreshold = params.percentageThreshold;
    // determine the haplotype of the read
    hpResult = germlineDetermineReadHap(hpCount, min, max, percentageThreshold, pqValue, psValue, norCountPS, nullptr, nullptr);

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
            std::cerr << "Error(extractNorData:postProcess) => can't find the position : chr:" << chr << " pos: " << (*currentPosIter).first;
            exit(1);
        }

        PosBase &baseInfo = (currentPosIter->second);
        auto curVar = currentVariants[(*currentPosIter).first].Variant[TUMOR];

        // current variant is SNP
        if(curVar.variantType == VariantType::SNP){
            std::string& tumAltBase = curVar.allele.Alt;

            // calculate the base information of the tumor SNP
            calculateBaseCommonInfo(baseInfo, tumAltBase);
        }
        
        currentPosIter++;
    } 
}

ExtractNorDataCigarParser::ExtractNorDataCigarParser(std::map<int, PosBase>& variantBase, std::vector<int>& tumVarPosVec, int& ref_pos, int& query_pos)
:CigarParser(ref_pos, query_pos), variantBase(variantBase), tumVarPosVec(tumVarPosVec){

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
            countBaseNucleotide(variantBase[curPos], base, *aln, params->somaticCallingMpqThreshold);
        }
        // the indel(deletion) SNP position is the start position, and the deletion occurs at the next position
        else if(curVar.variantType == VariantType::DELETION){
            // the indel SNP start position is at the end of the deletion, and the next cigar operator is deletion
            if(curPos == (ref_pos + length - 1) && bam_cigar_op(cigar[i+1]) == 2 && i+1 < aln_core_n_cigar){

            }
        }
    }       

    if ( aln->core.qual >= params->somaticCallingMpqThreshold && (*currentVariantIter).second.isExists(NORMAL)){
        // only judge the heterozygous SNP
        if((*currentVariantIter).second.Variant[NORMAL].is_phased_hetero){
            auto norVar = (*currentVariantIter).second.Variant[NORMAL];
            germlineJudgeSnpHap(*chrName, norVar, base, ref_pos, length, i, aln_core_n_cigar, cigar, currentVariantIter, *hpCount, *variantsHP, *norCountPS);
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
    if ( aln->core.qual >= params->somaticCallingMpqThreshold && (*currentVariantIter).second.isExists(NORMAL) && !alreadyJudgeDel){
        if((*currentVariantIter).second.Variant[NORMAL].is_phased_hetero){
            // longphase v1.73 only execute once
            alreadyJudgeDel = true;
            germlineJudgeDeletionHap(*chrName, *ref_string, ref_pos, length, query_pos, currentVariantIter, aln, *hpCount, *variantsHP, *norCountPS);
        }
    }
    
    //std::cout << "read ID : "<<  bam_get_qname(&aln) <<" pos :" << curPos << " deletion count ++" << std::endl;
}


ExtractTumDataBamParser::ExtractTumDataBamParser(
    std::map<std::string, std::map<int, HP3_Info>>& chrPosSomaticInfo,
    std::map<std::string, std::map<std::string, ReadVarHpCount>>& chrReadHpResultSet,
    std::map<std::string, std::map<int, std::map<std::string, int>>>& chrTumorPosReadCorrBaseHP
):
  chrPosSomaticInfo(chrPosSomaticInfo),
  chrReadHpResultSet(chrReadHpResultSet),
  chrTumorPosReadCorrBaseHP(chrTumorPosReadCorrBaseHP){

}

ExtractTumDataBamParser::~ExtractTumDataBamParser(){

}



ExtractTumDataChrProcessor::ExtractTumDataChrProcessor(
    std::map<std::string, std::map<int, HP3_Info>>& chrPosSomaticInfo,
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
    const std::string &chrName, 
    const HaplotagParameters &params, 
    const Genome& genmoeType, 
    std::map<int, MultiGenomeVar> &currentVariants,
    std::map<int, MultiGenomeVar>::iterator &firstVariantIter, 
    std::map<Genome, VCF_Info> &vcfSet, 
    const std::string &ref_string
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

    ExtractTumDataCigarParser cigarParser(*somaticPosInfo, tumorAllelePosVec, tumorSnpPosVec, tumCountPS, ref_pos, query_pos);
    cigarParser.parsingCigar(aln, bamHdr, chrName, params, firstVariantIter, currentVariants, ref_string, hpCount, variantsHP, norCountPS);

    int pqValue = 0;   
    double normalHPsimilarity = 0.0;
    double tumorHPsimilarity = 0.0;
    //calculate germline read HP result for predict tumor purity
    //calculate somatic read HP result for somatic read consistency filter
    int hpResult =determineReadHP(hpCount, pqValue, norCountPS, normalHPsimilarity, tumorHPsimilarity, params.percentageThreshold, nullptr, nullptr, nullptr);
   
    //classify read cases where tumor SNPs have low VAF in normal samples
    if(!tumorAllelePosVec.empty()){
        ClassifyReadsByCase(tumorAllelePosVec, norCountPS, hpCount, params, *somaticPosInfo);
        
        for(auto pos : tumorAllelePosVec){
            int baseHP = SnpHP::NONE_SNP;
            if(variantsHP.find(pos) != variantsHP.end()){
                baseHP = variantsHP[pos];
            }else{
                std::cerr << "Error (SomaticStatisticSomaticPosInfo) => can't find the position" << std::endl;
                std::cerr << "chr:" << chrName << " pos: " << pos+1 << std::endl;
                std::cerr << "readID: " << bam_get_qname(&aln) << std::endl;
                exit(1);
            }
            if(baseHP != SnpHP::SOMATIC_H3){
                std::cerr << "Error (SomaticStatisticSomaticPosInfo) => baseHP is not HP3 : chr:" << chrName << " pos: " << pos+1 << " baseHP: " << baseHP << std::endl;
                exit(1);
            }
            if(hpResult == ReadHP::H1_1 || hpResult == ReadHP::H2_1 || hpResult == ReadHP::H3 || hpResult == ReadHP::unTag){
                // record the somatic read HP for somatic read consistency filter
                (*somaticPosInfo)[pos].somaticReadHpCount[hpResult]++;
            }else if(hpResult == ReadHP::H1 || hpResult == ReadHP::H2){
                std::cerr << "Error (SomaticStatisticSomaticPosInfo) => error somatic read HP : chr:" << chrName << " pos: " << pos+1 << " hpResult: " << hpResult << std::endl;
                exit(1);      
            }
        }
    }

    //record variants HP count for each read in tumor-only position
    if(!tumorSnpPosVec.empty()){

        std::string readID = bam_get_qname(&aln);
        //read ID overide
        if((*readHpResultSet).find(readID) != (*readHpResultSet).end()){
            (*readHpResultSet)[readID].readIDcount++;
            readID = readID + "-" + std::to_string((*readHpResultSet)[readID].readIDcount);
        }

        (*readHpResultSet)[readID].HP1= hpCount[1];
        (*readHpResultSet)[readID].HP2= hpCount[2];
        (*readHpResultSet)[readID].HP3= hpCount[3];
        (*readHpResultSet)[readID].HP4= hpCount[4];

        (*readHpResultSet)[readID].norCountPS = norCountPS;
        (*readHpResultSet)[readID].startPos = aln.core.pos + 1;
        (*readHpResultSet)[readID].endPos = ref_pos;
        (*readHpResultSet)[readID].readLength = query_pos;
        
        //if(readID == "SRR25005626.11816585"){
        //    std::cerr << "readID: "<< readID << " chr: "<< chr << std::endl;
        //    std::cout << "HP1: "<< readTotalHPcount[readID].HP1 << " HP2: "<< readTotalHPcount[readID].HP2 << " HP3: "<< readTotalHPcount[readID].HP3 << " HP4: "<< readTotalHPcount[readID].HP4  << std::endl;
        //}
        for(auto pos : tumorSnpPosVec){
            int tumorSnpBaseHP = SnpHP::NONE_SNP;

            if(variantsHP.find(pos) != variantsHP.end()){
                tumorSnpBaseHP = variantsHP[pos];
            }
            // record the base HP whatever it is germline or somatic for statistic read distribution at current position
            (*tumorPosReadCorrBaseHP)[pos][readID] = tumorSnpBaseHP;
            // record the base HP for predict tumor purity
            (*somaticPosInfo)[pos].base.ReadHpCount[hpResult]++;
        } 
    }
}

void ExtractTumDataChrProcessor::ClassifyReadsByCase(std::vector<int> &tumorAllelePosVec, std::map<int, int> &norCountPS, std::map<int, int> &hpCount, const HaplotagParameters &params, std::map<int, HP3_Info> &somaticPosInfo){
    
    //decide whether to tag the read or not
    bool recordRead = true;
    /*if( max/(max+min) < params.percentageThreshold){
        // no tag
        recordRead = false;
    }*/

    //only exist tumor homo SNP
    if( norCountPS.size() == 0 && hpCount[3] != 0){
        //std::cerr << "Error : read only had tumor homo SNP " << bam_get_qname(aln) << "\n";
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
        std::cerr << "Error : hp3 or hp4 count is 0" << std::endl;
        std::cerr << "hp3 count:" << hpCount[3] << " hp4 count:" << hpCount[4];
        exit(1);
    }
    
    // only HP3 or one type of HP SNP
    if((zero_count == 1 || zero_count == 2) && hpCount[3] != 0){
        tagCleanHP3Read = true;
        //std::cout <<"HP1: " <<hp1Count<<"  HP2: " << hpCount[2]<<"  HP3: " << hpCount[3] << std::endl;
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
    std::map<int, HP3_Info>::iterator somaticVarIter = somaticPosInfo->begin();
    while( somaticVarIter != somaticPosInfo->end()){
        
        if(!currentVariants[(*somaticVarIter).first].isExists(TUMOR)){
            std::cerr << "Error(extractTumData:postProcess) => can't find the position : chr:" << chr << " pos: " << (*somaticVarIter).first;
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
            calculateBaseCommonInfo((*somaticVarIter).second.base, AltBase);

            // // the distribution of reads HP at the current position
            int H1readCount = (*somaticVarIter).second.base.ReadHpCount[ReadHP::H1];
            int H2readCount = (*somaticVarIter).second.base.ReadHpCount[ReadHP::H2];

            int H1_1readCount = (*somaticVarIter).second.base.ReadHpCount[ReadHP::H1_1];
            int H2_1readCount = (*somaticVarIter).second.base.ReadHpCount[ReadHP::H2_1];

            // the count of reads based on germline HP
            int BaseOnH1ReadCount = H1readCount + H1_1readCount;
            int BaseOnH2ReadCount = H2readCount + H2_1readCount;
            int totalBaseOnGermlineReadHpCount = BaseOnH1ReadCount + BaseOnH2ReadCount;

            // the ratio of the read count of H1 and H2 based on germline HP
            (*somaticVarIter).second.allelicImbalanceRatio = calculateHaplotypeImbalanceRatio(BaseOnH1ReadCount, BaseOnH2ReadCount, totalBaseOnGermlineReadHpCount);

            // the ratio of the read count of H1 and H2 based on somatic HP
            int totalSomaticReadHpCount = H1_1readCount + H2_1readCount;
            (*somaticVarIter).second.somaticHaplotypeImbalanceRatio = calculateHaplotypeImbalanceRatio(H1_1readCount, H2_1readCount, totalSomaticReadHpCount);

            //current SNP GT type
            if(curVar.is_homozygous == true){
                (*somaticVarIter).second.GTtype = "Homo";
            }else if(curVar.is_phased_hetero == true){
                (*somaticVarIter).second.GTtype = "Hetero";
            }else if(curVar.is_unphased_hetero == true){
                (*somaticVarIter).second.GTtype = "UnphasedHetero";
            }else{
                std::cerr << "Error(GTtype) => can't find GTtype at chr: " << chr << " pos: " << ((*somaticVarIter).first) + 1;
                exit(1);  
            }    
        }
        somaticVarIter++;
    }
}

ExtractTumDataCigarParser::ExtractTumDataCigarParser(
    std::map<int, HP3_Info> &somaticPosInfo,
    std::vector<int>& tumorAllelePosVec,
    std::vector<int>& tumorSnpPosVec,
    std::map<int, int>& tumCountPS,
    int& ref_pos,
    int& query_pos
):CigarParser(ref_pos, query_pos),
  somaticPosInfo(somaticPosInfo), tumorAllelePosVec(tumorAllelePosVec), tumorSnpPosVec(tumorSnpPosVec), tumCountPS(tumCountPS){

}

ExtractTumDataCigarParser::~ExtractTumDataCigarParser(){

}

void ExtractTumDataCigarParser::processMatchOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, std::string& base){
    //waring : using ref length to split SNP and indel that will be effect case ratio result 
    if ( aln->core.qual >= params->somaticCallingMpqThreshold ){
        SomaticJudgeSnpHP(currentVariantIter, *chrName, base, *hpCount, *norCountPS, tumCountPS, variantsHP, &tumorAllelePosVec);
        if((*currentVariantIter).second.isExists(TUMOR)){
            tumorSnpPosVec.push_back((*currentVariantIter).first);
        }
    }

    //statistically analyze SNP information exclusive to the tumor
    if((*currentVariantIter).second.isExists(TUMOR)){

        auto curVar = (*currentVariantIter).second.Variant[TUMOR];

        if(curVar.variantType == VariantType::SNP){
            //counting current tumor SNP base and depth           
            countBaseNucleotide(somaticPosInfo[(*currentVariantIter).first].base, base, *aln, params->somaticCallingMpqThreshold);
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

void ExtractTumDataCigarParser::OnlyTumorSNPjudgeHP(const std::string &chrName, int &curPos, MultiGenomeVar &curVar, std::string base, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *tumorAllelePosVec){
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



SomaticVarCaller::SomaticVarCaller(const std::vector<std::string> &chrVec, const HaplotagParameters &params){
    chrPosSomaticInfo = new std::map<std::string, std::map<int, HP3_Info>>();
    chrPosNorBase = new std::map<std::string, std::map<int, PosBase>>();
    callerReadHpDistri = new ReadHpDistriLog();
    denseTumorSnpInterval = new std::map<std::string, std::map<int, std::pair<int, DenseSnpInterval>>>();
    chrReadHpResultSet = new std::map<std::string, std::map<std::string, ReadVarHpCount>>();
    chrTumorPosReadCorrBaseHP = new std::map<std::string, std::map<int, std::map<std::string, int>>>();

    for(auto chr : chrVec){
        (*chrPosSomaticInfo)[chr] = std::map<int, HP3_Info>();
        (*chrPosNorBase)[chr] = std::map<int, PosBase>();
        callerReadHpDistri->loadChrKey(chr);
        (*denseTumorSnpInterval)[chr] = std::map<int, std::pair<int, DenseSnpInterval>>();
        (*chrReadHpResultSet)[chr] = std::map<std::string, ReadVarHpCount>();
        (*chrTumorPosReadCorrBaseHP)[chr] = std::map<int, std::map<std::string, int>>();
    }

    // setting somatic calling filter params
    InitialSomaticFilterParams(params.enableFilter);  
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

void SomaticVarCaller::VariantCalling(
    const HaplotagParameters &params,
    const std::vector<std::string> &chrVec,
    const std::map<std::string, int> &chrLength,
    std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat,
    std::map<Genome, VCF_Info> &vcfSet,
    const Genome& genmoeType
){

    //Count each base numbers at tumor SNP position in the Normal.bam
    std::cerr<< "extracting data from normal BAM... ";
    ExtractNorDataBamParser normalBamParser(*chrPosNorBase);
    normalBamParser.parsingBam(params.bamFile, params, chrVec, chrLength, mergedChrVarinat, vcfSet, Genome::NORMAL);

    std::cerr << "extracting data from tumor BAM... ";
    ExtractTumDataBamParser tumorBamParser(*chrPosSomaticInfo, *chrReadHpResultSet, *chrTumorPosReadCorrBaseHP);
    tumorBamParser.parsingBam(params.tumorBamFile, params, chrVec, chrLength, mergedChrVarinat, vcfSet, Genome::TUMOR);

    TumorPurityPredictor* tumorPurityPredictor = new TumorPurityPredictor(params, chrVec, *chrPosNorBase, *chrPosSomaticInfo);
    // predict tumor purity
    double tumorPurity = tumorPurityPredictor->predictTumorPurity();
    //flag the position that is used for purity prediction
    tumorPurityPredictor->markStatisticFlag(*chrPosSomaticInfo);
    delete tumorPurityPredictor;

    // set filter params with tumor purity
    SetFilterParamsWithPurity(somaticParams, tumorPurity);
    
    std::cerr<< "calling somatic variants ... ";
    std::time_t begin = time(NULL);
    // loop all chromosome
    #pragma omp parallel for schedule(dynamic) num_threads(params.numThreads) 
    for(auto chr : chrVec ){
        // record the position that tagged as HP3
        // chr, variant position
        std::map<int, HP3_Info> *somaticPosInfo = nullptr;
        // read ID, reads hpResult 
        // std::map<int, ReadHpResult> *localCallerReadHpDistri = nullptr;
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
            currentChrVariants = mergedChrVarinat[chr];
        }

        //get close somatic SNP interval
        getDenseTumorSnpInterval(*somaticPosInfo, *readHpResultSet, *tumorPosReadCorrBaseHP, *localDenseTumorSnpInterval);

        //calculate information and filter somatic SNPs
        SomaticFeatureFilter(somaticParams, currentChrVariants, chr, *somaticPosInfo, tumorPurity);

        //find other somatic variants HP4/HP5 (temporary)
        // FindOtherSomaticSnpHP(chr, *somaticPosInfo, currentChrVariants);

        // Shannon entropy filter(temporary)
        // ShannonEntropyFilter(chr, *somaticPosInfo, currentChrVariants, ref_string);

        //calibrate read HP to remove low confidence H3/H4 SNP
        CalibrateReadHP(chr, somaticParams, *somaticPosInfo, *readHpResultSet, *tumorPosReadCorrBaseHP);

        //calculate all read HP result
        CalculateReadSetHP(params, chr, *readHpResultSet, *tumorPosReadCorrBaseHP);

        //statistic all read HP in somatic SNP position
        StatisticSomaticPosReadHP(chr, *somaticPosInfo, *tumorPosReadCorrBaseHP, *readHpResultSet, *localCallerReadHpDistri);

    }
    std::cerr<< difftime(time(NULL), begin) << "s\n";

    int tumorOnlySNP = 0;
    for(auto chr : chrVec){
        auto somaticPosIter = (*chrPosSomaticInfo)[chr].begin();
        while(somaticPosIter != (*chrPosSomaticInfo)[chr].end()){
            if((*somaticPosIter).second.isHighConSomaticSNP){
                tumorOnlySNP++;
            }
            somaticPosIter++;
        }
    }
    std::cerr << "calling somatic SNPs: " << tumorOnlySNP << std::endl;
    
    if(somaticParams.writeVarLog){
        //write the log file for variants with positions tagged as HP3
        WriteSomaticVarCallingLog(params ,somaticParams, chrVec, mergedChrVarinat);

        WriteOtherSomaticHpLog(params, chrVec, mergedChrVarinat);

        WriteDenseTumorSnpIntervalLog(params, chrVec);

        callerReadHpDistri->writeReadHpDistriLog(params, "_readDistri_Scaller.out", chrVec);

        // remove the position that derived by HP1 or HP2
        callerReadHpDistri->removeNotDeriveByH1andH2pos(chrVec);

        // record the position that can't derived because exist H1-1 and H2-1 reads
        callerReadHpDistri->writeReadHpDistriLog(params, "_readDistri_Scaller_derive_by_H1_H2.out", chrVec);

    }
    return;
}

void SomaticVarCaller::InitialSomaticFilterParams(bool enableFilter){
    
    // Determine whether to apply the filter
    somaticParams.applyFilter = enableFilter;
    somaticParams.writeVarLog = true;

    somaticParams.tumorPurity = 1.0;

    // Below the mapping quality read ratio threshold
    somaticParams.LowMpqRatioThreshold = 0.1;

    somaticParams.MessyReadRatioThreshold = 1.0;
    somaticParams.ReadCount_minThr = 3;

    somaticParams.HapConsistency_ReadCount_maxThr = 10;
    somaticParams.HapConsistency_VAF_maxThr = 0.2;

    somaticParams.IntervalSnpCount_ReadCount_maxThr = 10;
    somaticParams.IntervalSnpCount_VAF_maxThr = 0.15;
}


void SomaticVarCaller::SetFilterParamsWithPurity(SomaticFilterParaemter &somaticParams, double &tumorPurity){

    somaticParams.tumorPurity = tumorPurity;

    // tumor purity 1.0
    if (tumorPurity >= 0.9 && tumorPurity <= 1.0) 
    {
        somaticParams.norVAF_maxThr = 0.084;
        somaticParams.norDepth_minThr = 1;

        somaticParams.MessyReadRatioThreshold = 1.0;    
        somaticParams.ReadCount_minThr = 4.0;

        somaticParams.HapConsistency_ReadCount_maxThr = 12.0;
        somaticParams.HapConsistency_VAF_maxThr = 0.238;
        somaticParams.HapConsistency_somaticRead_minThr = 1.0;

        somaticParams.IntervalSnpCount_ReadCount_maxThr = 12.0;
        somaticParams.IntervalSnpCount_VAF_maxThr = 0.286;
        somaticParams.IntervalSnpCount_minThr = 4.0;
        somaticParams.zScore_maxThr = 4.621;
    } 
    // tumor purity 0.8
    else if (tumorPurity >= 0.7 && tumorPurity < 0.9) 
    {
        somaticParams.norVAF_maxThr = 0.095;
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
    } 
    // tumor purity 0.6
    else if (tumorPurity >= 0.5 && tumorPurity < 0.7) 
    {
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
    }
    // tumor purity 0.4 
    else if (tumorPurity >= 0.3 && tumorPurity < 0.5) 
    {
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
    }
    // tumor purity 0.2
    else if (tumorPurity > 0 && tumorPurity < 0.3) 
    {
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
    } 
    else 
    {
        std::cerr << "[Warning] tumor purity is not in the range of 0.0 to 1.0: " << tumorPurity << std::endl;
        std::cerr << "[Warning] setting default parameters" << std::endl;
        // default parameters (tumor purity = 0.2)
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
    }
}

void SomaticVarCaller::SomaticFeatureFilter(const SomaticFilterParaemter &somaticParams, std::map<int, MultiGenomeVar> &currentChrVariants,const std::string &chr, std::map<int, HP3_Info> &somaticPosInfo, double& tumorPurity){
    //calculate the information of the somatic positon 
    std::map<int, HP3_Info>::iterator somaticVarIter = somaticPosInfo.begin();

    while( somaticVarIter != somaticPosInfo.end()){
                
        if(!currentChrVariants[(*somaticVarIter).first].isExists(TUMOR)){
            std::cerr << "Error(SomaticFeatureFilter) => can't find the position : chr:" << chr << " pos: " << (*somaticVarIter).first;
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
            if (somaticParams.applyFilter && (*somaticVarIter).second.isFilterOut) {
                somaticVarIter++;
                continue;
            }

            (*somaticVarIter).second.isHighConSomaticSNP = true;
        }

        somaticVarIter++;          
    }
}

void SomaticVarCaller::ShannonEntropyFilter(const std::string &chr, std::map<int, HP3_Info> &somaticPosInfo, std::map<int, MultiGenomeVar> &currentChrVariants, std::string &ref_string){
    int window_size = 10;
    auto somaticVarIter = somaticPosInfo.begin();
    while(somaticVarIter != somaticPosInfo.end()){
        int read_count = (*somaticVarIter).second.CaseReadCount;
        if(read_count <= 10 && (*somaticVarIter).second.isHighConSomaticSNP){
            int snp_pos = (*somaticVarIter).first;
            int window_start = snp_pos - window_size/2;
            int window_end = snp_pos + window_size/2;

            int A_count=0;
            int C_count=0;
            int T_count=0;
            int G_count=0;

            char homopolymer_base = ' ';
            int homopolymer_length = 0;

            char max_homopolymer_base= ' ';
            int max_homopolymer_length=0;
            // std::cout << "------somatic SNP------: "<< snp_pos << std::endl;
            for(int pos = window_start; pos <= window_end; pos++){
                char base = ref_string.at(pos);
                A_count += base == 'A' ? 1 : 0;
                C_count += base == 'C' ? 1 : 0;
                T_count += base == 'T' ? 1 : 0;
                G_count += base == 'G' ? 1 : 0;


                // std::cout << "pos: " << pos << " base: " << base << std::endl;
                // std::cout << "homopolymer_base: " << homopolymer_base << " homoplerLen: " << homopolymer_length << std::endl;
                //count homopolymer length  
                if(homopolymer_base == base){
                    homopolymer_length++;
                }else{
                    if(max_homopolymer_length < homopolymer_length){
                        max_homopolymer_length = homopolymer_length;
                        max_homopolymer_base = homopolymer_base;
                        // std::cout << "---------------> max_base: " << max_homopolymer_base << " maxLen: " << max_homopolymer_length << std::endl;
                    }
                    homopolymer_base = base;
                    homopolymer_length = 1;

                }

            }

            // std::cout << "homopolymer_base: " << homopolymer_base << " homoplerLen: " << homopolymer_length << std::endl;
            if(max_homopolymer_length < homopolymer_length){
                max_homopolymer_length = homopolymer_length;
                max_homopolymer_base = homopolymer_base;
            }
            // std::cout << "---------------> max_base: " << max_homopolymer_base << " maxLen: " << max_homopolymer_length << std::endl;
            // std::cout << " " << std::endl;


            // if(homopolymer_length >= 5){
            // }
            // std::cout << "A_count: " << A_count << " C_count: " << C_count << " T_count: " << T_count << " G_count: " << G_count << std::endl;
            double entropy = calculateShannonEntropy(A_count, C_count, T_count, G_count);
            (*somaticVarIter).second.shannonEntropy = entropy;
            (*somaticVarIter).second.homopolymerLength = max_homopolymer_length;
            //filter out the low entropy somatic SNP
            if(entropy < 1.5){
                // (*somaticVarIter).second.isHighConSomaticSNP = false;
            }
        }
        somaticVarIter++;
    }
}

double SomaticVarCaller::entropyComponent(int count, int total){
    if (count == 0) return 0.0;
    double p = static_cast<double>(count) / total;
    return -p * std::log2(p);
}


double SomaticVarCaller::calculateShannonEntropy(int nA, int nC, int nT, int nG){
    int total = nA + nC + nT + nG;
    if (total == 0) return 0.0; 
    // calculate Shannon entropy
    double entropy = entropyComponent(nA, total) +
                     entropyComponent(nC, total) +
                     entropyComponent(nT, total) +
                     entropyComponent(nG, total);

    // std::cout << "total: " << total << std::endl;
    // std::cout << "entropy A: " << entropyComponent(nA, total) << std::endl;
    // std::cout << "entropy C: " << entropyComponent(nC, total) << std::endl;
    // std::cout << "entropy T: " << entropyComponent(nT, total) << std::endl;
    // std::cout << "entropy G: " << entropyComponent(nG, total) << std::endl;
    // std::cout << "entropy: " << entropy << std::endl;

    return entropy;
}

double SomaticVarCaller::calculateMean(const std::map<int, double>& data){
    double size = data.size();
    if(size == 0) return 0.0;
    
    double sum = 0.0;
    for(auto meanIter : data){
        sum += meanIter.second;
    }
    return sum / size;
}

double SomaticVarCaller::calculateStandardDeviation(const std::map<int, double>& data, double mean){
    double variance = 0.0;
    for (auto meanIter : data) {
        double value = meanIter.second;
        variance += (value - mean) * (value - mean);
    }
    return std::sqrt(variance / data.size());
}

void SomaticVarCaller::calculateZScores(const std::map<int, double>& data, double mean, double stdDev, std::map<int, double> &zScores){
    for (auto meanIter : data) {
        double value = meanIter.second;
        if(stdDev == 0){
            zScores[meanIter.first] = 0.0;
        }else{
            zScores[meanIter.first] = ((value - mean) / stdDev); // calculate z-score
        }
    }
}

void SomaticVarCaller::calculateIntervalZScore(bool &isStartPos, int &startPos, int &endPos, int &snpCount, DenseSnpInterval &denseSnp, std::map<int, std::pair<int, DenseSnpInterval>> &localDenseTumorSnpInterval){
    double mean = calculateMean(denseSnp.snpAltMean);
    double stdDev = calculateStandardDeviation(denseSnp.snpAltMean, mean);
    calculateZScores(denseSnp.snpAltMean, mean, stdDev, denseSnp.snpZscore);

    denseSnp.snpCount = snpCount;
    denseSnp.totalAltMean = mean;
    denseSnp.StdDev = stdDev;   

    localDenseTumorSnpInterval[startPos] = std::make_pair(endPos, denseSnp); // save the interval

    isStartPos = false;
    startPos = 0;
    snpCount = 0;
    denseSnp.snpAltMean.clear();
    denseSnp.snpZscore.clear();
}

void SomaticVarCaller::getDenseTumorSnpInterval(std::map<int, HP3_Info> &somaticPosInfo, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, std::map<std::string, int>> &tumorPosReadCorrBaseHP, std::map<int, std::pair<int, DenseSnpInterval>> &localDenseTumorSnpInterval){
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
                // if(somaticPosIter->first + 1 == 81435550 || somaticPosIter->first + 1 == 81435861){
                //     std::cout << "ReadID: " << readHPIter->first << " baseHP: " << baseHP << " HP3: " << readHpResultSet[readHPIter->first].HP3 << std::endl;
                // }
            }else{
                std::cerr << "ERROR: readID not found in readHpResultSet: " << readHPIter->first << std::endl;
                exit(1);
            }
        }
        if(altMean != 0) altMean /= readCount;

        if(somaticPosInfo.find(somaticPosIter->first) != somaticPosInfo.end()){
            somaticPosInfo[somaticPosIter->first].MeanAltCountPerVarRead = altMean;
        }else{
            std::cerr << "ERROR: somaticPosInfo not found: " << somaticPosIter->first << std::endl;
            exit(1);
        }
        // if(somaticPosIter->first + 1 == 81435550 || somaticPosIter->first + 1 == 81435861){
        //     std::cout << "pos: " << somaticPosIter->first + 1 << " readCount: " << readCount << " AltMeanCount: " << altMean << " minStartPos: " << minStartPos<< " maxEndPos: " << maxEndPos << std::endl;
        // }
    }

    //find the interval of tumor SNPs
    auto somaticPosIter = somaticPosInfo.begin();
    bool isStartPos = false;
    int startPos = 0;
    int snpCount = 0;
    int dense_distance = 5000;
    DenseSnpInterval denseSnp;

    while (somaticPosIter != somaticPosInfo.end()){
        int curPos = somaticPosIter->first;

        // ensure not out of bounds
        auto nextIter = std::next(somaticPosIter);
        if (nextIter != somaticPosInfo.end()){
            int nextPos = nextIter->first;

            if (nextPos - curPos <= dense_distance) {
                if (!isStartPos) {
                    isStartPos = true;
                    startPos = curPos;
                    snpCount++;
                    denseSnp.snpAltMean[curPos] = somaticPosIter->second.MeanAltCountPerVarRead;
                    // altMeanVec.push_back(somaticPosIter->second.altMeanCount);
                    // altPosVec.push_back(pos);
                    // std::cout << "------------->startPos: " << startPos << " AltMean: " << somaticPosIter->second.MeanAltCountPerVarRead << std::endl;
                }
                denseSnp.snpAltMean[nextPos] = nextIter->second.MeanAltCountPerVarRead;
                // altMeanVec.push_back(nextIter->second.altMeanCount);
                // altPosVec.push_back(nextPos);
                snpCount++;
                // std::cout << "pos: " << pos << " nextPos: " << nextPos << " snpCount: " << snpCount << " altVecSize: " << altMeanVec.size() << " nextAltMean: " << nextIter->second.altMeanCount << std::endl;
            } else {
                if (isStartPos) {
                    calculateIntervalZScore(isStartPos, startPos, curPos, snpCount, denseSnp, localDenseTumorSnpInterval);
                    // std::cout << "------------->endPos: " << pos << " stdDev: " << stdDev << "\n"<< std::endl;
                }else{
                    // std::cerr << "ERROR: isStartPos is false: pos: " << curPos << " nextPos: " << nextPos << " startPos: " << startPos << std::endl;
                }
            }
        }

        somaticPosIter = nextIter;
    }

    // check if there is an unfinished interval
    if (isStartPos) {
        int endPos = somaticPosInfo.rbegin()->first;
        if(endPos - startPos <= dense_distance){
            // std::cout << "Last interval ------------->" << "startPos: " << startPos << " endPos: " << endPos << std::endl;
            calculateIntervalZScore(isStartPos, startPos, endPos, snpCount, denseSnp, localDenseTumorSnpInterval);
        }
        // localDenseTumorSnpInterval[startPos] = std::make_pair(somaticPosInfo.rbegin()->first, denseSnpInterval()); // save the last interval
    }

    // record the zScore of each SNP in the dense tumor interval
    for (const auto& entry : localDenseTumorSnpInterval) {
        if(entry.second.second.snpCount > 1){
            for(auto snpIter: entry.second.second.snpZscore){
                somaticPosInfo[snpIter.first].inDenseTumorInterval = true;
                somaticPosInfo[snpIter.first].zScore = abs(snpIter.second);
                somaticPosInfo[snpIter.first].intervalSnpCount = entry.second.second.snpCount;
            }
        }
    }
}

void SomaticVarCaller::FindOtherSomaticSnpHP(const std::string &chr, std::map<int, HP3_Info> &somaticPosInfo, std::map<int, MultiGenomeVar> &currentChrVariants){
    std::map<int, HP3_Info>::iterator somaticVarIter = somaticPosInfo.begin();
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
                refAllele = convertStrNucToInt(currentChrVariants[pos].Variant[TUMOR].allele.Ref);
                altAllele = convertStrNucToInt(currentChrVariants[pos].Variant[TUMOR].allele.Alt);
            }else{
                std::cerr << "Error(FindOtherSomaticSnpHP) => can't find position in currentChrVariants : chr:" << chr << " pos: " << pos + 1;
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

            // if(maxNuc != Nitrogenous::UNKOWN && minNuc == Nitrogenous::UNKOWN){
            //     std::cout << "chr: " << chr << " pos: " << pos + 1 << "\n";
            //     std::cout << "maxNuc: " << convertIntNucToStr(maxNuc) << " minNuc: " << convertIntNucToStr(minNuc) << "\n";
            //     std::cout << "maxNucCount: " << maxNucCount << " minNucCount: " << minNucCount << "\n";
            // }
        }
        somaticVarIter++;
    }
}

int SomaticVarCaller::convertStrNucToInt(std::string &base){
    if(base == "A"){
        return Nitrogenous::A;
    }else if(base == "C"){
        return Nitrogenous::C;
    }else if(base == "G"){
        return Nitrogenous::G;
    }else if(base == "T"){
        return Nitrogenous::T;
    }else{
        std::cerr << "Error(convertNucleotideToInt) => can't find Allele : " << base << "\n";
        exit(1);
    }
}

std::string SomaticVarCaller::convertIntNucToStr(int base){
    if(base == Nitrogenous::A){
        return "A";
    }else if(base == Nitrogenous::C){
        return "C";
    }else if(base == Nitrogenous::G){
        return "G";
    }else if(base == Nitrogenous::T){
        return "T";
    }else if(base == Nitrogenous::UNKOWN){
        //std::cerr << "Error(convertIntNucToStr) => can't find Allele : " << base << "\n";
        return "UNKOWN";
    }else {
        std::cerr << "Error(convertIntNucToStr) => can't find Allele : " << base << "\n";
        exit(1);
    }
}

void SomaticVarCaller::CalibrateReadHP(const std::string &chr, const SomaticFilterParaemter &somaticParams, std::map<int, HP3_Info> &somaticPosInfo, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, std::map<std::string, int>> &tumorPosReadCorrBaseHP){
        //calibrate read HP to remove low confidence H3/H4 SNP
        std::map<int, HP3_Info>::iterator somaticVarIter = somaticPosInfo.begin();
        while( somaticVarIter != somaticPosInfo.end()){
            if(!(*somaticVarIter).second.isHighConSomaticSNP){
                int pos = (*somaticVarIter).first;

                //only reads with HP3 or HP4 SNPs will be calibrated
                if(tumorPosReadCorrBaseHP.find(pos) != tumorPosReadCorrBaseHP.end()){

                    for(auto readIdBaseHP : tumorPosReadCorrBaseHP[pos]){
                        std::string readID = readIdBaseHP.first;
                        int baseHP = readIdBaseHP.second;

                        switch (baseHP) {
                            case SnpHP::GERMLINE_H1:
                                break;
                            case SnpHP::GERMLINE_H2:
                                // std::cerr << "ERROR (calibrate read HP) => somatic Variant had HP1 or HP2 SNP :" << std::endl;
                                // std::cerr << "readID: " << readID << " chr: " << chr << " pos: " << pos + 1 << std::endl;
                                // std::cerr << "baseHP: " << baseHP << std::endl;
                                // exit(1);
                                break;
                            case SnpHP::SOMATIC_H3:
                                readHpResultSet[readID].HP3--; break;
                            default:
                                break;
                        }

                        if(readHpResultSet[readID].HP3 < 0){
                            std::cerr << "ERROR (calibrate read HP) => read HP3 or HP4 SNP count < 0 :" << std::endl;
                            std::cerr << "readID: "<< readID << " chr: "<< chr << " pos: " << pos+1<< std::endl;
                            std::cerr << "HP3: "<< readHpResultSet[readID].HP3 << std::endl;
                            exit(1);
                        }
                    }
                }else{
                    //std::cerr << "ERROR (calibrate read HP) => can't find pos in tumorPosReadCorrBaseHP : chr: "<< chr << " pos: " << pos+1 <<std::endl;
                    //exit(1);
                }
            }
            somaticVarIter++;
        }
}

void SomaticVarCaller::CalculateReadSetHP(const HaplotagParameters &params, const std::string &chr, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, std::map<std::string, int>> &tumorPosReadCorrBaseHP){
        std::map<std::string, ReadVarHpCount>::iterator readTotalHPcountIter = readHpResultSet.begin();
        //calculate all read HP result 
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

            (*readTotalHPcountIter).second.hpResult = determineReadHP(hpCount, pqValue, norCountPS, normalHPsimilarity, tumorHPsimilarity, params.percentageThreshold, nullptr, nullptr, nullptr);
            
            readTotalHPcountIter++;
        }
}

void SomaticVarCaller::StatisticSomaticPosReadHP(
    const std::string &chr, std::map<int, HP3_Info> &somaticPosInfo,
    std::map<int, std::map<std::string, int>> &tumorPosReadCorrBaseHP,
    std::map<std::string, ReadVarHpCount> &readHpResultSet,
    chrReadHpResult &localReadHpDistri
){
    std::map<int, HP3_Info>::iterator somaticVarIter = somaticPosInfo.begin();
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

                //error: haven't exist somatic HP read in current position 
                if( localReadHpDistri.posReadHpResult[pos].readHpCounter[ReadHP::H3] == 0 && 
                    localReadHpDistri.posReadHpResult[pos].readHpCounter[ReadHP::H1_1] == 0 &&
                    localReadHpDistri.posReadHpResult[pos].readHpCounter[ReadHP::H2_1] == 0){
                    // std::cerr << "ERROR (statistic all read HP) => hadn't exist somatic HP read in : chr: "<< chr << " pos: " << pos+1 <<std::endl;
                    // std::cerr << "HP1-1: "<< readHpDistributed[pos].hpResultCounter[ReadHP::H1_1] 
                    //           << " HP1-2: "<< readHpDistributed[pos].hpResultCounter[ReadHP::H1_2] 
                    //           << " HP2-1: "<< readHpDistributed[pos].hpResultCounter[ReadHP::H2_1] 
                    //           << " HP2-2: "<< readHpDistributed[pos].hpResultCounter[ReadHP::H2_2] 
                    //           << " HP3: "<< readHpDistributed[pos].hpResultCounter[ReadHP::H3] 
                    //           << " HP4: "<< readHpDistributed[pos].hpResultCounter[ReadHP::H4]
                    //           << std::endl;
                    // exit(1); 
                }

            }else{
                std::cerr << "ERROR (statistic all read HP) => can't find pos in tumorPosReadCorrBaseHP : chr: "<< chr << " pos: " << pos+1 <<std::endl;
                exit(1); 
            }

        }
        somaticVarIter++;
    }
}


void SomaticVarCaller::WriteSomaticVarCallingLog(const HaplotagParameters &params, const SomaticFilterParaemter &somaticParams, const std::vector<std::string> &chrVec, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat){
    std::ofstream *tagHP3Log = new std::ofstream(params.resultPrefix+"_HP3.out");

    if(!tagHP3Log->is_open()){
        std::cerr<< "Fail to open write file: " << params.resultPrefix+"_HP3.out" << "\n";
        exit(1);
    }

    std::cerr << "writing somatic varinats calling log ... ";
    std::time_t begin = time(NULL);

    int totalSomaticSNP = 0;
    for(auto chr : chrVec){
        std::map<int, HP3_Info>::iterator somaticVarIter = (*chrPosSomaticInfo)[chr].begin();
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
    (*tagHP3Log) << "##normalSnpFile:"       << params.snpFile             << "\n"
                 << "##tumorSnpFile:"        << params.tumorSnpFile        << "\n" 
                 << "##svFile:"              << params.svFile              << "\n"
                 << "##bamFile:"             << params.bamFile             << "\n"
                 << "##tumorBamFile:"        << params.tumorBamFile        << "\n" 
                 << "##resultPrefix:"        << params.resultPrefix        << "\n"
                 << "##numThreads:"          << params.numThreads          << "\n"
                 << "##region:"              << params.region              << "\n"
                 << "##qualityThreshold:"    << params.qualityThreshold    << "\n"
                 << "##percentageThreshold:" << params.percentageThreshold << "\n"
                 << "##tagSupplementary:"    << params.tagSupplementary    << "\n"
                 << "##tagTumor:"            << params.tagTumorSnp         << "\n"
                 << "##\n";

    //write filter parameter
    (*tagHP3Log) << "##======== Filter Parameters =========\n"
                 << "##Enable filter : " << somaticParams.applyFilter << "\n"
                 << "##Calling mapping quality :" << params.somaticCallingMpqThreshold << "\n"
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
                 //<< "##Low mapping quality read ratio threshold : " << somaticParams.LowMpqRatioThreshold << "\n"
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
                 << "GermlineReadHpConsistencyRatio\t"
                 << "SomaticReadHpConsistencyRatio\t"
                 << "BaseGermlineReadHpConsistencyRatio\t"
                 << "PercentageOfGermlineHp\t"
                 << "H1readCountInNorBam\t"
                 << "H2readCountInNorBam\t"
                 << "GermlineReadHpCountInNorBam\t"
                 << "GermlineReadHpConsistencyRatioInNorBam\t"
                 << "PercentageOfGermlineHpInNorBam\t"
                 << "GermlineReadHpConsistencyRatioDifference\t"
                 << "PercentageOfGermlineHpDifference\t"
                 << "SomaticRead_H1-1\t"
                 << "SomaticRead_H2-1\t"
                 << "SomaticRead_H3\t"
                 << "SomaticRead_unTag\t"
                 << "AltMeanCountPerVarRead\t"
                 << "zScore\t"
                 << "IntervalSnpCount\t"
                 << "ExistNorSnp\t"
                 << "StatisticPurity\t"
                 << "isFilterOut\t"
                 << "NorNonDelAF\t"
                 << "TumNonDelAF\t"
                 << "GT\n";

    //write variants information
    for(auto chr : chrVec){
    
        std::map<int, HP3_Info>::iterator somaticVarIter = (*chrPosSomaticInfo)[chr].begin();
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
                
            if(mergedChrVarinat[chr][(*somaticVarIter).first].isExists(TUMOR) == true){
                RefBase = mergedChrVarinat[chr][(*somaticVarIter).first].Variant[TUMOR].allele.Ref;
                AltBase = mergedChrVarinat[chr][(*somaticVarIter).first].Variant[TUMOR].allele.Alt;
            }else{
                std::cerr << "Error(write tag HP3 log file) => can't find the position : chr:" << chr << " pos: " << ((*somaticVarIter).first) + 1;
                exit(1);
            }

            if(RefBase == "" || AltBase == ""){
                std::cerr << "Error(write tag HP3 log file) => can't find RefBase or AltBase : chr:" << chr << " pos: " << ((*somaticVarIter).first) + 1 << " RefBase:" << RefBase << " AltBase:" << AltBase;
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
                std::cerr << "Error(GTtype) => chr: " << chr << " pos: " << ((*somaticVarIter).first) + 1 << " GTtype: " << GTtype;
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
            double baseOnGermlineReadHpConsistencyRatio = (*somaticVarIter).second.allelicImbalanceRatio;

            double somaticReadHpConsistencyRatio = (*somaticVarIter).second.somaticHaplotypeImbalanceRatio;

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
                    std::cerr << "Error(zScore) => chr: " << chr << " pos: " << ((*somaticVarIter).first) + 1 << " zScore: " << (*somaticVarIter).second.zScore << std::endl;
                    exit(1);
                }else{
                    zScore = (*somaticVarIter).second.zScore;
                }
            }
            
            //HP3 SNP position (1-base)
            int HP3pos = (*somaticVarIter).first + 1;
            
            //write information
            // (*tagHP3Log) << std::fixed << std::setprecision(3) 
                        // << chr << " \t"  
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
                        << somaticReadHpConsistencyRatio << "\t" //42
                        << baseOnGermlineReadHpConsistencyRatio << "\t" //43
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
                        << mergedChrVarinat[chr][(*somaticVarIter).first].isExists(NORMAL) << "\t" //59
                        << (*somaticVarIter).second.statisticPurity << "\t" //60
                        << (*somaticVarIter).second.isFilterOut << "\t" //61
                        << norNonDelAF << "\t" //62
                        << tumNonDelAF << "\t" //63
                        << GTtype <<"\n";  //64
                        
            somaticVarIter++;          
        }
    }

    tagHP3Log->close();
    delete tagHP3Log;
    tagHP3Log = nullptr;
    std::cerr<< difftime(time(NULL), begin) << "s\n";  
}

void SomaticVarCaller::WriteOtherSomaticHpLog(const HaplotagParameters &params, const std::vector<std::string> &chrVec, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat){
    std::ofstream *OtherHpSomaticVarLog=NULL;
    std::string logPosfix = "_otherHpSomaticVar.log";
    OtherHpSomaticVarLog=new std::ofstream(params.resultPrefix + logPosfix);

    int totalOtherSomaticHpVar = 0;
    for(auto chr: chrVec){
        for(auto somaticVar: (*chrPosSomaticInfo)[chr]){
            if(somaticVar.second.somaticHp4Base != Nitrogenous::UNKOWN){
                totalOtherSomaticHpVar++;
            }
        }
    }

    if(!OtherHpSomaticVarLog->is_open()){
        std::cerr<< "Fail to open write file: " << params.resultPrefix + logPosfix << "\n";
        exit(1);
    }else{
        (*OtherHpSomaticVarLog) << "################################\n";
        (*OtherHpSomaticVarLog) << "# Other Somatic HP Variant Log #\n";
        (*OtherHpSomaticVarLog) << "################################\n";
        (*OtherHpSomaticVarLog) << "##MappingQualityThreshold:"  << params.qualityThreshold << "\n";
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
        auto currentChrVar = mergedChrVarinat[chr];
        for(auto somaticVar: (*chrPosSomaticInfo)[chr]){
            int pos = somaticVar.first;

            if(somaticVar.second.somaticHp4Base == Nitrogenous::UNKOWN){
                continue;
            }
            if(currentChrVar.find(pos) == currentChrVar.end()){
                std::cerr << "Error(WriteOtherSomaticHpLog) => can't find the position : chr:" << chr << " pos: " << pos + 1;
                exit(1);
            }

            (*OtherHpSomaticVarLog) << chr << "\t"
                                    << pos + 1 << "\t"
                                    << currentChrVar[pos].Variant[TUMOR].allele.Ref << "\t"
                                    << currentChrVar[pos].Variant[TUMOR].allele.Alt << "\t"
                                    << convertIntNucToStr(somaticVar.second.somaticHp4Base) << "\t"
                                    << somaticVar.second.somaticHp4BaseCount<< "\t"
                                    << convertIntNucToStr(somaticVar.second.somaticHp5Base) << "\t"
                                    << somaticVar.second.somaticHp5BaseCount << "\n";
        }
    
    
    }
    (*OtherHpSomaticVarLog).close();
    delete OtherHpSomaticVarLog;
    OtherHpSomaticVarLog = nullptr;
}

void SomaticVarCaller::WriteDenseTumorSnpIntervalLog(const HaplotagParameters &params, const std::vector<std::string> &chrVec){
    std::ofstream *closeSomaticSnpIntervalLog=NULL;
    std::string logPosfix = "_denseTumorSnpInterval.log";
    closeSomaticSnpIntervalLog=new std::ofstream(params.resultPrefix + logPosfix);

    int totalIntervalCount = 0;
    for(auto chr: chrVec){
        totalIntervalCount += (*denseTumorSnpInterval)[chr].size();
    }

    if(!closeSomaticSnpIntervalLog->is_open()){
        std::cerr<< "Fail to open write file: " << params.resultPrefix + logPosfix << "\n";
        exit(1);
    }else{
        (*closeSomaticSnpIntervalLog) << "################################\n";
        (*closeSomaticSnpIntervalLog) << "# Dense Tumor SNP Interval Log #\n";
        (*closeSomaticSnpIntervalLog) << "################################\n";
        (*closeSomaticSnpIntervalLog) << "##MappingQualityThreshold:"  << params.qualityThreshold << "\n";
        (*closeSomaticSnpIntervalLog) << "##Tatal intervals:"  << totalIntervalCount << "\n";
        (*closeSomaticSnpIntervalLog) << "#CHROM\t"
                                      << "startPos-endPos\t"
                                      << "SnpCount\t"
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
                (*closeSomaticSnpIntervalLog) <<"#snp:altMean:zScore=>  " << snpIter.first + 1 << " : " << snpIter.second << " : " << zScore << "\n";
            }
            (*closeSomaticSnpIntervalLog) << "#\n";
            
        }
    }
    (*closeSomaticSnpIntervalLog).close();
    delete closeSomaticSnpIntervalLog;
    
    closeSomaticSnpIntervalLog = nullptr;
}

void SomaticVarCaller::getSomaticFlag(const std::vector<std::string> &chrVec, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat){
    for(auto chr: chrVec){
        for(auto somaticVar: (*chrPosSomaticInfo)[chr]){
            mergedChrVarinat[chr][somaticVar.first].isSomaticVariant = true;
            mergedChrVarinat[chr][somaticVar.first].somaticReadDeriveByHP = somaticVar.second.somaticReadDeriveByHP;
        }
    }
}