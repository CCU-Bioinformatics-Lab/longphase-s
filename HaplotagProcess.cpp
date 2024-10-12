#include "HaplotagProcess.h"



BamBaseCounter::BamBaseCounter(){
    ChrVariantBase = new std::map<std::string, std::map<int, PosBase>>();
};

BamBaseCounter::~BamBaseCounter(){
    delete ChrVariantBase;
};


void BamBaseCounter::CountingBamBase(const std::string &BamFile, const HaplotagParameters &params, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat, std::vector<std::string> &chrVec, std::map<std::string, int> &chrLength, int genmoeType){
    std::cerr<< "collecting data for the normal sample... ";
    std::time_t begin = time(NULL);

    if(chrVec.size() == 0){
        std::cerr<<"ERROR: chrVec is empty\n";
        exit(1);
    }

    if(chrLength.size() == 0){
        std::cerr<<"ERROR: chrLength is empty\n";
        exit(1);        
    }

    for(auto chr : chrVec){
        (*ChrVariantBase)[chr] = std::map<int, PosBase>();
    }

    // init data structure and get core n
    htsThreadPool threadPool = {NULL, 0};
    // creat thread pool
    if (!(threadPool.pool = hts_tpool_init(params.numThreads))) {
        fprintf(stderr, "Error creating thread pool\n");
        exit(1);
    }

    // loop all chromosome
    #pragma omp parallel for schedule(dynamic) num_threads(params.numThreads) 
    for(auto chr : chrVec ){

        // bam file resource allocation
        BamFileRAII bam(BamFile, params.fastaFile, threadPool);

        // store base information
        std::map<int, PosBase> *variantBase = nullptr;

        // variant position (0-base), allele haplotype set
        std::map<int, RefAltSet> currentVariants;

        #pragma omp critical
        {
            variantBase = &((*ChrVariantBase)[chr]);
            currentVariants = mergedChrVarinat[chr];
        }

        //inintial iterator
        std::map<int, RefAltSet>::iterator firstVariantIter = currentVariants.begin();

        std::map<int, RefAltSet>::reverse_iterator last = currentVariants.rbegin();

        std::string region = !params.region.empty() ? params.region : chr + ":1-" + std::to_string(chrLength[chr]);
        hts_itr_t* iter = sam_itr_querys(bam.idx, bam.bamHdr, region.c_str());

        while (sam_itr_multi_next(bam.in, iter, bam.aln) >= 0) {
            
            int flag = bam.aln->core.flag;

            //if ( aln->core.qual < params.qualityThreshold ){
            //    // mapping quality is lower than threshold
            //    continue;
            //}

            if( (flag & 0x4) != 0 ){
                // read unmapped
                continue;
            }
            if( (flag & 0x100) != 0 ){
                // secondary alignment. repeat.
                // A secondary alignment occurs when a given read could align reasonably well to more than one place.
                continue;
            }
            if( (flag & 0x800) != 0 && params.tagSupplementary == false ){
                // supplementary alignment
                // A chimeric alignment is represented as a set of linear alignments that do not have large overlaps.
                continue;
            }
            //currentVariants rbegin == rend
            if(last == currentVariants.rend()){ 
                //skip
            }
            else if(int(bam.aln->core.pos) <= (*last).first){
                StatisticBaseInfo(*bam.aln, chr, params, genmoeType, *variantBase, currentVariants, firstVariantIter);
            }
        }
        //std::cerr<<"calculate Base Information ...\n";
        CalculateBaseInfo(chr, *variantBase, currentVariants);
        
        variantBase = nullptr;
        hts_itr_destroy(iter);
    }
    hts_tpool_destroy(threadPool.pool);
    std::cerr<< difftime(time(NULL), begin) << "s\n";
    return;
}

void BamBaseCounter::StatisticBaseInfo(const bam1_t &aln, const std::string &chrName, const HaplotagParameters &params, int genmoeType, std::map<int, PosBase> &variantBase, std::map<int, RefAltSet> &currentVariants ,std::map<int, RefAltSet>::iterator &firstVariantIter){
    
    // Skip variants that are to the left of this read
    while( firstVariantIter != currentVariants.end() && (*(firstVariantIter)).first < aln.core.pos ){
        firstVariantIter++;
    }     

    if( firstVariantIter == currentVariants.end()){
        return ;
    }

    // position relative to reference
    int ref_pos = aln.core.pos;
    // position relative to read
    int query_pos = 0;
    // set variant start for current alignment
    std::map<int, RefAltSet>::iterator currentVariantIter = firstVariantIter;

    // reading cigar to detect snp on this read
    int aln_core_n_cigar = int(aln.core.n_cigar);
    for(int i = 0; i < aln_core_n_cigar ; i++ ){
        uint32_t *cigar = bam_get_cigar(&aln);
        int cigar_op = bam_cigar_op(cigar[i]);
        int length   = bam_cigar_oplen(cigar[i]);

        // iterator next variant
        while( currentVariantIter != currentVariants.end() && (*currentVariantIter).first < ref_pos ){
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

                    //statistically analyze SNP information exclusive to the tumor
                    if((*currentVariantIter).second.isExistTumor && !((*currentVariantIter).second.isExistNormal)){
                        int curPos = (*currentVariantIter).first;
                        //std::cout << "curPos :" << curPos << "is tumor "<<(*currentVariantIter).second.isExistTumor  << " isNormal: "<< (*currentVariantIter).second.isExistNormal<< std::endl;

                        //detect ref bsae length(temp :tumor SNP)
                        int tumRefLength = (*currentVariantIter).second.Variant[Genome::TUMOR].Ref.length();
                        int tumAltLength = (*currentVariantIter).second.Variant[Genome::TUMOR].Alt.length();
                        
                        // the variant is SNP
                        if(tumRefLength == 1 && tumAltLength == 1){

                            if(base == "A"){
                                variantBase[curPos].A_count++;
                                //std::cout << "base A read ID :" << bam_get_qname(aln) << std::endl;
                            }else if(base == "C"){
                                variantBase[curPos].C_count++;
                            }else if(base == "G"){
                                variantBase[curPos].G_count++;
                            }else if(base == "T"){
                                variantBase[curPos].T_count++;
                            }else{
                                variantBase[curPos].unknow++;
                            }
                            variantBase[curPos].depth++;  

                            // mapping quality is higher than threshold
                            if ( aln.core.qual >= params.somaticCallingMpqThreshold ){
                                if(base == "A"){
                                    variantBase[curPos].MPQ_A_count++;
                                }else if(base == "C"){
                                    variantBase[curPos].MPQ_C_count++;
                                }else if(base == "G"){
                                    variantBase[curPos].MPQ_G_count++;
                                }else if(base == "T"){
                                    variantBase[curPos].MPQ_T_count++;
                                }else{
                                    variantBase[curPos].MPQ_unknow++;
                                }
                                variantBase[curPos].filteredMpqDepth++;
                            }
                        }

                        // the indel(del) SNP position is the start position, and the deletion occurs at the next position
                        else if(tumRefLength > 1 && tumAltLength == 1){
                            // the indel SNP start position is at the end of the deletion, and the next cigar operator is deletion
                            if(curPos == (ref_pos + length - 1) && bam_cigar_op(cigar[i+1]) == 2 && i+1 < aln_core_n_cigar){

                            }
                        }
                    }               
                }
                currentVariantIter++;
            }
            query_pos += length;
            ref_pos += length;
        }
            // 1: insertion to the reference
        else if( cigar_op == 1 ){
            query_pos += length;
        }
        // 2: deletion from the reference
        else if( cigar_op == 2 ){
            while( currentVariantIter != currentVariants.end() && (*currentVariantIter).first < ref_pos + length){
                //statistically analyze SNP information exclusive to the tumor
                if((*currentVariantIter).second.isExistTumor && !((*currentVariantIter).second.isExistNormal)){

                    int curPos = (*currentVariantIter).first;

                    //detect ref bsae length(temp :tumor SNP)
                    int tumRefLength = (*currentVariantIter).second.Variant[Genome::TUMOR].Ref.length();
                    int tumAltLength = (*currentVariantIter).second.Variant[Genome::TUMOR].Alt.length();

                    // the variant is SNP
                    if(tumRefLength == 1 && tumAltLength == 1){
                        variantBase[curPos].delCount++;
                    }
                    // the variant is deletion
                    else if(tumRefLength > 1 && tumAltLength == 1){

                    }
                }
                
                //std::cout << "read ID : "<<  bam_get_qname(&aln) <<" pos :" << curPos << " deletion count ++" << std::endl;
                currentVariantIter++;
            }

            ref_pos += length;
        }
        // 3: skipped region from the reference
        else if( cigar_op == 3 ){
            ref_pos += length;
        }
        // 4: soft clipping (clipped sequences present in SEQ)
        else if( cigar_op == 4 ){
            query_pos += length;
        }
        // 5: hard clipping (clipped sequences NOT present in SEQ)
        // 6: padding (silent deletion from padded reference)
        else if( cigar_op == 5 || cigar_op == 6 ){
        // do nothing
        }
        else{
            std::cerr<< "alignment find unsupported CIGAR operation from read: " << bam_get_qname(&aln) << "\n";
            exit(1);
        }
    }
}

void BamBaseCounter::CalculateBaseInfo(const std::string &chr, std::map<int, PosBase> &VariantBase, std::map<int, RefAltSet> &currentVariants){
    std::map<int, PosBase>::iterator currentPosIter = VariantBase.begin();
    if(currentPosIter == VariantBase.end()){
        //exit(1);
        return;
    } 

    while (currentPosIter != VariantBase.end()) {

        PosBase *baseInfo = &(currentPosIter->second);

        // Determine whether this position is clean or not
        int zero_count = 0;
        if (baseInfo->A_count == 0) {
            zero_count++;
        }
        if (baseInfo->C_count == 0) {
            zero_count++;
        }
        if (baseInfo->G_count == 0) {
            zero_count++;
        }
        if (baseInfo->T_count == 0) {
            zero_count++;
        }

        // Which type of bases are the most and second most frequent
        int max_count = 0;
        int second_max_count = 0;
        std::string max_base = " ";
        std::string second_max_base = " ";

        // Determine max & second max base count
        if (baseInfo->A_count > max_count) {
            second_max_count = max_count;
            max_count = baseInfo->A_count;
            second_max_base = max_base;
            max_base = "A";
        } else if (baseInfo->A_count > second_max_count) {
            second_max_count = baseInfo->A_count;
            second_max_base = "A";
        }

        if (baseInfo->C_count > max_count) {
            second_max_count = max_count;
            max_count = baseInfo->C_count;
            second_max_base = max_base;
            max_base = "C";
        } else if (baseInfo->C_count > second_max_count) {
            second_max_count = baseInfo->C_count;
            second_max_base = "C";
        }

        if (baseInfo->G_count > max_count) {
            second_max_count = max_count;
            max_count = baseInfo->G_count;
            second_max_base = max_base;
            max_base = "G";
        } else if (baseInfo->G_count > second_max_count) {
            second_max_count = baseInfo->G_count;
            second_max_base = "G";
        }

        if (baseInfo->T_count > max_count) {
            second_max_count = max_count;
            max_count = baseInfo->T_count;
            second_max_base = max_base;
            max_base = "T";
        } else if (baseInfo->T_count > second_max_count) {
            second_max_count = baseInfo->T_count;
            second_max_base = "T";
        }

        int depth = baseInfo->depth;

        //record max base information
        baseInfo->max_count = max_count;
        baseInfo->max_base = max_base;
        baseInfo->max_ratio = (float)max_count / (float)depth;

        //record second max base information
        if(zero_count == 3){
            baseInfo->second_max_count = 0;
            baseInfo->second_max_base = ' ';
            baseInfo->second_max_ratio = 0.0;
        }else{
            baseInfo->second_max_count = second_max_count;
            baseInfo->second_max_base = second_max_base;
            baseInfo->second_max_ratio = (float)second_max_count / (float)depth;
        }

        // calculate VAF
        std::string tumRefBase = currentVariants[(*currentPosIter).first].Variant[Genome::TUMOR].Ref;
        std::string tumAltBase = currentVariants[(*currentPosIter).first].Variant[Genome::TUMOR].Alt;

        int AltCount = 0;
        int filteredMpqAltCount = 0;

        if(tumAltBase == "A"){
            AltCount = baseInfo->A_count;
            filteredMpqAltCount = baseInfo->MPQ_A_count;
        }else if(tumAltBase == "T"){
            AltCount = baseInfo->T_count;
            filteredMpqAltCount = baseInfo->MPQ_T_count;
        }else if(tumAltBase == "C"){
            AltCount = baseInfo->C_count;
            filteredMpqAltCount = baseInfo->MPQ_C_count;
        }else if(tumAltBase == "G"){
            AltCount = baseInfo->G_count;
            filteredMpqAltCount = baseInfo->MPQ_G_count;
        }


        if(AltCount != 0 && depth != 0){
            baseInfo->VAF = (float)AltCount / (float)depth;
        }

        int filteredMpqDepth = baseInfo->filteredMpqDepth;
        
        if(filteredMpqAltCount != 0 && filteredMpqDepth != 0){
            baseInfo->filteredMpqVAF = (float)filteredMpqAltCount / (float)filteredMpqDepth;
        }

        if(depth != 0){
            baseInfo->lowMpqReadRatio = (float)(depth - filteredMpqDepth) / (float)depth;
        }
        
        float refAllele_threshold = 0.90; //(VAF <= 0.1)      
        if(baseInfo->max_ratio >= refAllele_threshold && depth > 1 && baseInfo->lowMpqReadRatio <= 0.1){
            baseInfo->isHighRefAllelleFreq = true;
        } 
        
        currentPosIter++;
    } 
}

bool BamBaseCounter::isHighRefAllelleFreq(std::string chr, int pos){
    if((*ChrVariantBase)[chr].find(pos) == (*ChrVariantBase)[chr].end()){
        return false;
    }    
    return (*ChrVariantBase)[chr][pos].isHighRefAllelleFreq;
}

std::string BamBaseCounter::getMaxFreqBase(std::string chr, int pos){
    if((*ChrVariantBase)[chr].find(pos) == (*ChrVariantBase)[chr].end()){
        std::cerr << "ERROR (getMaxBase) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return (*ChrVariantBase)[chr][pos].max_base;
}

float BamBaseCounter::getMaxBaseRatio(std::string chr, int pos){
    if((*ChrVariantBase)[chr].find(pos) == (*ChrVariantBase)[chr].end()){
        std::cerr << "ERROR (getMaxBaseRatio) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return (*ChrVariantBase)[chr][pos].max_ratio;
}

float BamBaseCounter::getSecondMaxBaseRatio(std::string chr, int pos){
    if((*ChrVariantBase)[chr].find(pos) == (*ChrVariantBase)[chr].end()){
        std::cerr << "ERROR (getSecondMaxBaseRatio) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return (*ChrVariantBase)[chr][pos].second_max_ratio;
}

float BamBaseCounter::getLowMpqReadRatio(std::string chr, int pos){
    if((*ChrVariantBase)[chr].find(pos) == (*ChrVariantBase)[chr].end()){
        std::cerr << "ERROR (getLowMpqReadRatio) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return (*ChrVariantBase)[chr][pos].lowMpqReadRatio;
}

float BamBaseCounter::getVAF(std::string chr, int pos){
    if((*ChrVariantBase)[chr].find(pos) == (*ChrVariantBase)[chr].end()){
        std::cerr << "ERROR (getVAF) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return (*ChrVariantBase)[chr][pos].VAF;
}

float BamBaseCounter::getFilterdMpqVAF(std::string chr, int pos){
    if((*ChrVariantBase)[chr].find(pos) == (*ChrVariantBase)[chr].end()){
        std::cerr << "ERROR (getFilterdMpqVAF) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return (*ChrVariantBase)[chr][pos].filteredMpqVAF;
}

int BamBaseCounter::getBaseAcount(std::string chr, int pos){
    if((*ChrVariantBase)[chr].find(pos) == (*ChrVariantBase)[chr].end()){
        std::cerr << "ERROR (getBaseAcount) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return (*ChrVariantBase)[chr][pos].A_count;
}

int BamBaseCounter::getBaseCcount(std::string chr, int pos){
    if((*ChrVariantBase)[chr].find(pos) == (*ChrVariantBase)[chr].end()){
        std::cerr << "ERROR (getBaseCcount) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return (*ChrVariantBase)[chr][pos].C_count;
}

int BamBaseCounter::getBaseTcount(std::string chr, int pos){
    if((*ChrVariantBase)[chr].find(pos) == (*ChrVariantBase)[chr].end()){
        std::cerr << "ERROR (getBaseTcount) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return (*ChrVariantBase)[chr][pos].T_count;
}

int BamBaseCounter::getBaseGcount(std::string chr, int pos){
    if((*ChrVariantBase)[chr].find(pos) == (*ChrVariantBase)[chr].end()){
        std::cerr << "ERROR (getBaseGcount) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return (*ChrVariantBase)[chr][pos].G_count;
}

int BamBaseCounter::getDepth(std::string chr, int pos){
    if((*ChrVariantBase)[chr].find(pos) == (*ChrVariantBase)[chr].end()){
        std::cerr << "ERROR (getDepth) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return (*ChrVariantBase)[chr][pos].depth;
}

int BamBaseCounter::getMpqDepth(std::string chr, int pos){
    if((*ChrVariantBase)[chr].find(pos) == (*ChrVariantBase)[chr].end()){
        std::cerr << "ERROR (getMpqDepth) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return (*ChrVariantBase)[chr][pos].filteredMpqDepth;
}

int BamBaseCounter::getVarDeletionCount(std::string chr, int pos){
    if((*ChrVariantBase)[chr].find(pos) == (*ChrVariantBase)[chr].end()){
        std::cerr << "ERROR (getDelCount) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return (*ChrVariantBase)[chr][pos].delCount;
}

void BamBaseCounter::displayPosInfo(std::string chr, int pos){
    if((*ChrVariantBase)[chr].find(pos) == (*ChrVariantBase)[chr].end()){
        std::cerr << "ERROR (displayPosInfo) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }else{
        PosBase VarBase = (*ChrVariantBase)[chr][pos];
        std::cout << " chr : " <<  chr <<" Pos :" << pos << std::endl;
        std::cout << " A_count : " << VarBase.A_count << std::endl;
        std::cout << " C_count : " << VarBase.C_count << std::endl;
        std::cout << " G_count : " << VarBase.G_count << std::endl;
        std::cout << " T_count : " << VarBase.T_count << std::endl;
        std::cout << " max_base : " << VarBase.max_base << "  ratio : "<< VarBase.max_ratio << std::endl;
        std::cout << " second_max_base : " << VarBase.second_max_base << "  ratio : "<< VarBase.second_max_ratio << std::endl;
        std::cout << " depth :   " << VarBase.depth << std::endl;
        std::cout << " unknow :  " << VarBase.unknow << std::endl;
        std::cout << " deletion count :  " << VarBase.delCount << std::endl;
        std::cout << " VAF :  " << VarBase.VAF << std::endl;
        std::cout << " filterd MPQ VAF :  " << VarBase.filteredMpqVAF << std::endl;
        std::cout << " Low MPQ read ratio :  " << VarBase.lowMpqReadRatio << std::endl;
    }
}

void SomaticJudgeBase::SomaticJudgeSnpHP(std::map<int, RefAltSet>::iterator &currentVariantIter, VCF_Info *vcfSet, std::string chrName, std::string base, std::map<int, int> &hpCount, std::map<int, int> &norCountPS, std::map<int, int> &tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *readPosHP3, BamBaseCounter *NorBase, std::map<int, HP3_Info> *SomaticPos){
    int curPos = (*currentVariantIter).first;
    auto curVar = (*currentVariantIter).second;

    // normal & tumor SNP at the current position (base on normal phased SNPs)
    // both normal and tumor samples that do not exist in the high-confidence set
    if(curVar.isExistNormal == true && curVar.isExistTumor == true){

        // the tumor & normal SNP GT are phased heterozygous 
        if((curVar.Variant[Genome::NORMAL].is_phased_hetero) && (curVar.Variant[Genome::TUMOR].is_phased_hetero)){   
            if(curVar.Variant[Genome::NORMAL].Ref == base || curVar.Variant[Genome::NORMAL].Alt == base){
                //std::cerr<< "tag tumor normal SNP\n";
                std::map<int, int>::iterator NorPosPSiter = vcfSet[Genome::NORMAL].chrVariantPS[chrName].find(curPos);
                if( NorPosPSiter == vcfSet[Genome::NORMAL].chrVariantPS[chrName].end()){
                    std::cerr<< "Unable to locate the phase set of the current normal SNP\n"
                             << curPos << "\t"
                             << curVar.Variant[Genome::NORMAL].Ref << "\t"
                             << curVar.Variant[Genome::NORMAL].Alt  << "\n";
                    exit(EXIT_SUCCESS);
                }

                std::string norHP1 = vcfSet[Genome::NORMAL].chrVariantHP1[chrName][curPos];
                std::string norHP2 = vcfSet[Genome::NORMAL].chrVariantHP2[chrName][curPos];
                std::string tumHP1 = vcfSet[Genome::TUMOR].chrVariantHP1[chrName][curPos];
                std::string tumHP2 = vcfSet[Genome::TUMOR].chrVariantHP2[chrName][curPos];

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
                norCountPS[vcfSet[Genome::NORMAL].chrVariantPS[chrName][curPos]]++;
            }
        //the normal SNP GT is phased heterozgous & the tumor SNP GT is unphased heterozgous 
        }else if((curVar.Variant[Genome::NORMAL].is_phased_hetero) && (curVar.Variant[Genome::TUMOR].is_unphased_hetero)){   
            if(curVar.Variant[Genome::NORMAL].Ref == base || curVar.Variant[Genome::NORMAL].Alt == base){

                std::map<int, int>::iterator NorPosPSiter = vcfSet[Genome::NORMAL].chrVariantPS[chrName].find(curPos);
                if( NorPosPSiter == vcfSet[Genome::NORMAL].chrVariantPS[chrName].end()){
                    std::cerr<< "Unable to locate the phase set of the current normal SNP\n"
                             << curPos << "\t"
                             << curVar.Variant[Genome::NORMAL].Ref << "\t"
                             << curVar.Variant[Genome::NORMAL].Alt  << "\n";
                    exit(EXIT_SUCCESS);
                }

                std::string norHP1 = vcfSet[Genome::NORMAL].chrVariantHP1[chrName][curPos];
                std::string norHP2 = vcfSet[Genome::NORMAL].chrVariantHP2[chrName][curPos];

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
                norCountPS[vcfSet[Genome::NORMAL].chrVariantPS[chrName][curPos]]++;
            }
        //the normal SNP GT is phased heterozgous & the tumor SNP GT is homozygous 
        }else if((curVar.Variant[Genome::NORMAL].is_phased_hetero) && (curVar.Variant[Genome::TUMOR].is_homozygous)){   
            if(curVar.Variant[Genome::NORMAL].Ref == base || curVar.Variant[Genome::NORMAL].Alt == base){

                std::map<int, int>::iterator NorPosPSiter = vcfSet[Genome::NORMAL].chrVariantPS[chrName].find(curPos);
                if( NorPosPSiter == vcfSet[Genome::NORMAL].chrVariantPS[chrName].end()){
                    std::cerr<< "Unable to locate the phase set of the current normal SNP\n"
                             << curPos << "\t"
                             << curVar.Variant[Genome::NORMAL].Ref << "\t"
                             << curVar.Variant[Genome::NORMAL].Alt  << "\n";
                    exit(EXIT_SUCCESS);
                }

                std::string norHP1 = vcfSet[Genome::NORMAL].chrVariantHP1[chrName][curPos];
                std::string norHP2 = vcfSet[Genome::NORMAL].chrVariantHP2[chrName][curPos];

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
                norCountPS[vcfSet[Genome::NORMAL].chrVariantPS[chrName][curPos]]++;
            }
        }
    // only normal SNP at the current position
    }else if(curVar.isExistNormal == true){
        // the normal SNP GT is phased heterozgous SNP
        if((curVar.Variant[Genome::NORMAL].is_phased_hetero)){
            if(curVar.Variant[Genome::NORMAL].Ref == base || curVar.Variant[Genome::NORMAL].Alt == base){

                std::map<int, int>::iterator NorPosPSiter = vcfSet[Genome::NORMAL].chrVariantPS[chrName].find(curPos);
                if( NorPosPSiter == vcfSet[Genome::NORMAL].chrVariantPS[chrName].end()){
                    std::cerr<< "Unable to locate the phase set of the current normal SNP\n"
                             << curPos << "\t"
                             << curVar.Variant[Genome::NORMAL].Ref << "\t"
                             << curVar.Variant[Genome::NORMAL].Alt  << "\n";
                    exit(EXIT_SUCCESS);
                }

                std::string norHP1 = vcfSet[Genome::NORMAL].chrVariantHP1[chrName][curPos];
                std::string norHP2 = vcfSet[Genome::NORMAL].chrVariantHP2[chrName][curPos];

                if( base == norHP1){
                    hpCount[1]++;
                    if(variantsHP != nullptr) (*variantsHP)[curPos] = SnpHP::GERMLINE_H1;
                }
                if(base == norHP2){
                    hpCount[2]++;
                    if(variantsHP != nullptr) (*variantsHP)[curPos] = SnpHP::GERMLINE_H2;
                }
                norCountPS[vcfSet[Genome::NORMAL].chrVariantPS[chrName][curPos]]++;
            }
        }
    // only tumor SNP at the current position
    }else if(curVar.isExistTumor == true){

        //the tumor SNP GT is phased heterozygous
        if(curVar.Variant[Genome::TUMOR].is_phased_hetero == true){
            if(curVar.Variant[Genome::TUMOR].Ref == base || curVar.Variant[Genome::TUMOR].Alt == base){
                std::map<int, int>::iterator posPSiter = vcfSet[Genome::TUMOR].chrVariantPS[chrName].find(curPos);
                if( posPSiter == vcfSet[Genome::TUMOR].chrVariantPS[chrName].end() ){
                    std::cerr<< curPos << "\t"
                             << curVar.Variant[Genome::TUMOR].Ref << "\t"
                             << curVar.Variant[Genome::TUMOR].Alt << "\n";
                    exit(EXIT_SUCCESS);
                }else{
                    OnlyTumorSNPjudgeHP(chrName, curPos, curVar, base, vcfSet, hpCount, &tumCountPS, variantsHP, readPosHP3, NorBase, SomaticPos);
                }
            }
        //the tumor SNP GT is unphased heterozygous
        }else if(curVar.Variant[Genome::TUMOR].is_unphased_hetero == true){
            if(curVar.Variant[Genome::TUMOR].Ref == base || curVar.Variant[Genome::TUMOR].Alt == base){
                OnlyTumorSNPjudgeHP(chrName, curPos, curVar, base, vcfSet, hpCount, nullptr, variantsHP, readPosHP3, NorBase, SomaticPos);
            }           
        //the tumor SNP GT is homozygous
        }else if(curVar.Variant[Genome::TUMOR].is_homozygous == true){
            if(curVar.Variant[Genome::TUMOR].Ref == base || curVar.Variant[Genome::TUMOR].Alt == base){
                OnlyTumorSNPjudgeHP(chrName, curPos, curVar, base, vcfSet, hpCount, nullptr, variantsHP, readPosHP3, NorBase, SomaticPos);
            }
        }
    }
}

void SomaticJudgeBase::OnlyTumorSNPjudgeHP(const std::string &chrName, int &curPos, RefAltSet &curVar, std::string base, VCF_Info *vcfSet, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *readPosHP3, BamBaseCounter *NorBase, std::map<int, HP3_Info> *SomaticPos){

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
    if(hpCount[3] != 0 && hpCount[4] != 0 ){
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

int SomaticJudgeBase::convertStrNucToInt(std::string &base){
    if(base == "A"){
        return Nucleotide::A;
    }else if(base == "C"){
        return Nucleotide::C;
    }else if(base == "G"){
        return Nucleotide::G;
    }else if(base == "T"){
        return Nucleotide::T;
    }else{
        std::cerr << "Error(convertNucleotideToInt) => can't find Allele : " << base << "\n";
        exit(1);
    }
}

std::string SomaticJudgeBase::convertIntNucToStr(int base){
    if(base == Nucleotide::A){
        return "A";
    }else if(base == Nucleotide::C){
        return "C";
    }else if(base == Nucleotide::G){
        return "G";
    }else if(base == Nucleotide::T){
        return "T";
    }else if(base == Nucleotide::UNKOWN){
        //std::cerr << "Error(convertIntNucToStr) => can't find Allele : " << base << "\n";
        return "UNKOWN";
    }else {
        std::cerr << "Error(convertIntNucToStr) => can't find Allele : " << base << "\n";
        exit(1);
    }
}

void SomaticJudgeBase::recordReadHp(int &pos, int &hpResult, int &BaseHP, std::map<int, ReadHpResult> &varReadHpResult){
    varReadHpResult[pos].hpResultCounter[hpResult]++;
    
    if(hpResult != ReadHP::unTag){
        if(BaseHP == SnpHP::SOMATIC_H3){
            varReadHpResult[pos].somaticSnpH3count++;
        }
    }
}

readHpDistriLog::readHpDistriLog(){

}

readHpDistriLog::~readHpDistriLog(){

}

void readHpDistriLog::mergeLocalReadHp(const std::string &chr, std::map<int, ReadHpResult> &localReadHpResult){
    std::map<int, ReadHpResult>::iterator localReadHpIter = localReadHpResult.begin();
    while(localReadHpIter != localReadHpResult.end()){
        int pos = (*localReadHpIter).first;
        chrVarReadHpResult[chr][pos] = (*localReadHpIter).second;
        localReadHpIter++;
    }
}

void readHpDistriLog::writeReadHpDistriLog(HaplotagParameters &params, std::string logPosfix, const std::vector<std::string> &chrVec){
    std::ofstream *readHpDistriLog=NULL;
    readHpDistriLog=new std::ofstream(params.resultPrefix + logPosfix);

    int somaticSnpCount = 0;
    for(auto chr: chrVec){
        if(!chrVarReadHpResult[chr].empty()){
            somaticSnpCount += chrVarReadHpResult[chr].size();
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
                            << "Type\t"
                            << "AltCount\t\t"
                            << "HP1read\t"
                            << "HP1-1read\t"
                            << "HP1-2read\t"
                            << "HP2read\t"
                            << "HP2-1read\t"
                            << "HP2-2read\t"
                            << "HP3read\t"
                            << "HP4read\t"
                            << "untagRead\t"
                            << "HP1ratio\t"
                            << "HP1-1ratio\t"
                            << "HP1-2ratio\t"
                            << "HP2ratio\t"
                            << "HP2-1ratio\t"
                            << "HP2-2ratio\t"
                            << "HP3ratio\t"
                            << "HP4ratio\n";
    }

    for(auto chr: chrVec){
        std::map<int, ReadHpResult>::iterator curVarReadHpIter = chrVarReadHpResult[chr].begin();
        while(curVarReadHpIter != chrVarReadHpResult[chr].end()){
            int pos = (*curVarReadHpIter).first + 1;
            int HP1readCount = (*curVarReadHpIter).second.hpResultCounter[ReadHP::H1];
            int HP1_1readCount = (*curVarReadHpIter).second.hpResultCounter[ReadHP::H1_1];
            int HP1_2readCount = (*curVarReadHpIter).second.hpResultCounter[ReadHP::H1_2];

            int HP2readCount = (*curVarReadHpIter).second.hpResultCounter[ReadHP::H2];
            int HP2_1readCount = (*curVarReadHpIter).second.hpResultCounter[ReadHP::H2_1];
            int HP2_2readCount = (*curVarReadHpIter).second.hpResultCounter[ReadHP::H2_2];

            int HP3readCount = (*curVarReadHpIter).second.hpResultCounter[ReadHP::H3];
            int HP4readCount = (*curVarReadHpIter).second.hpResultCounter[ReadHP::H4];

            int totaltagRead = HP1readCount + HP2readCount + HP3readCount + HP4readCount + HP1_1readCount + HP1_2readCount + HP2_1readCount + HP2_2readCount;

            float HP1readRatio = (float)HP1readCount / (float)totaltagRead; 
            float HP1_1readRatio = (float)HP1_1readCount / (float)totaltagRead; 
            float HP1_2readRatio = (float)HP1_2readCount / (float)totaltagRead; 

            float HP2readRatio = (float)HP2readCount / (float)totaltagRead; 
            float HP2_1readRatio = (float)HP2_1readCount / (float)totaltagRead; 
            float HP2_2readRatio = (float)HP2_2readCount / (float)totaltagRead; 

            float HP3readRatio = (float)HP3readCount / (float)totaltagRead; 
            float HP4readRatio = (float)HP4readCount / (float)totaltagRead;

            (*readHpDistriLog) << std::fixed << std::setprecision(3) 
                                << chr << "\t"
                                << pos << "\t"
                                << "somatic" << "\t"
                                << (*curVarReadHpIter).second.somaticSnpH3count << "\t\t"
                                << HP1readCount << "\t"
                                << HP1_1readCount << "\t"
                                << HP1_2readCount << "\t\t"
                                << HP2readCount << "\t"
                                << HP2_1readCount << "\t"
                                << HP2_2readCount << "\t\t"
                                << HP3readCount << "\t"
                                << HP4readCount << "\t"
                                << (*curVarReadHpIter).second.hpResultCounter[ReadHP::unTag] << "\t"
                                << HP1readRatio << "\t"
                                << HP1_1readRatio << "\t"
                                << HP1_2readRatio << "\t\t"
                                << HP2readRatio << "\t"
                                << HP2_1readRatio << "\t"
                                << HP2_2readRatio << "\t\t"
                                << HP3readRatio << "\t"
                                << HP4readRatio << "\n";
            curVarReadHpIter++;
        }
    }
    (*readHpDistriLog).close();
    delete readHpDistriLog;
    readHpDistriLog = nullptr;
}

void readHpDistriLog::writePosCoverRegionLog(HaplotagParameters &params, std::string logPosfix, const std::vector<std::string> &chrVec){
    std::ofstream *posCoverRegionLog=NULL;
    posCoverRegionLog=new std::ofstream(params.resultPrefix + logPosfix);

    int somaticSnpCount = 0;
    for(auto chr: chrVec){
        if(!chrVarReadHpResult[chr].empty()){
            somaticSnpCount += chrVarReadHpResult[chr].size();
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
        std::map<int, ReadHpResult>::iterator curVarReadHpIter = chrVarReadHpResult[chr].begin();
        while(curVarReadHpIter != chrVarReadHpResult[chr].end()){
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

void readHpDistriLog::writeTagReadCoverRegionLog(HaplotagParameters &params, std::string logPosfix, const std::vector<std::string> &chrVec, std::map<std::string, int> &chrLength){
    std::ofstream *tagReadCoverRegionLog=NULL;
    tagReadCoverRegionLog=new std::ofstream(params.resultPrefix + logPosfix);

    std::map<std::string, std::vector<coverRegionInfo>> coverRegion;

    // merge cover region from different SNPs
    for (const auto& chr : chrVec){
        auto curVarReadHpIter = chrVarReadHpResult[chr].begin();

        if (curVarReadHpIter == chrVarReadHpResult[chr].end()) {
            continue;
        }

        int curStartPos = curVarReadHpIter->second.coverRegionStartPos;
        int curEndPos = curVarReadHpIter->second.coverRegionEndPos;

        while (curVarReadHpIter != chrVarReadHpResult[chr].end()){

            auto nextVarReadHpIter = std::next(curVarReadHpIter);
            if(nextVarReadHpIter != chrVarReadHpResult[chr].end()){
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
        (*tagReadCoverRegionLog) << "##MappingQualityThreshold:"        << params.qualityThreshold << "\n";
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


void readHpDistriLog::removeNotDeriveByH1andH2pos(const std::vector<std::string> &chrVec){
    for (const auto& chr : chrVec) {
        auto& chrResult = chrVarReadHpResult[chr];
        for (auto it = chrResult.begin(); it != chrResult.end(); ) {
            if (!it->second.existDeriveByH1andH2) {
                //std::cerr << "Removed position not derived by H1 and H2: " << chr << " " << it->first << std::endl;
                it = chrResult.erase(it);
            } else {
                ++it;
            }
        }
    }
}


BamFileRAII::BamFileRAII(const std::string& BamFile, const std::string& fastaFile, htsThreadPool &threadPool):
in(nullptr), bamHdr(nullptr), idx(nullptr), aln(nullptr)
{
    // open bam file
    in = hts_open(BamFile.c_str(), "r");
    if (in == nullptr) {
        std::cerr << "ERROR: Cannot open bam file " + BamFile << std::endl;
        exit(1);
    }

    // load reference file
    if (hts_set_fai_filename(in, fastaFile.c_str()) != 0) {
        std::cerr << "ERROR: Cannot set FASTA index file for " + fastaFile << std::endl;
        exit(1);
    }

    // input reader
    bamHdr = sam_hdr_read(in);
    if (bamHdr == nullptr) {
        std::cerr << "ERROR: Cannot read header from bam file " + BamFile << std::endl;
        exit(1);
    }

    // check bam file index
    idx = sam_index_load(in, BamFile.c_str());
    if (idx == nullptr) {
        std::cerr << "ERROR: Cannot open index for bam file " + BamFile << std::endl;
        exit(1);
    }

    // set thread
    if (hts_set_opt(in, HTS_OPT_THREAD_POOL, &threadPool) != 0) {
        std::cerr << "ERROR: Cannot set thread pool for bam file " + BamFile << std::endl;
        exit(1);
    }

    // initialize an alignment
    aln = bam_init1();
    if (aln == nullptr) {
        std::cerr << "ERROR: Cannot initialize alignment for bam file " + BamFile << std::endl;
        exit(1);
    }
}

BamFileRAII::~BamFileRAII(){
    if (aln) bam_destroy1(aln);
    if (idx) hts_idx_destroy(idx);
    if (bamHdr) bam_hdr_destroy(bamHdr);
    if (in) sam_close(in);
}


SomaticVarCaller::SomaticVarCaller(){
    SomaticChrPosInfo = new std::map<std::string, std::map<int, HP3_Info>>();
    chrVarReadHpResult = new std::map<std::string, std::map<int, ReadHpResult>>();
    callerReadHpDistri = new readHpDistriLog();
}

SomaticVarCaller::~SomaticVarCaller(){
    releaseMemory();
}

void SomaticVarCaller::releaseMemory(){
    delete SomaticChrPosInfo;
    delete chrVarReadHpResult;
    delete callerReadHpDistri;
}

void SomaticVarCaller::VariantCalling(const std::string BamFile, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat,const std::vector<std::string> &chrVec, std::map<std::string, int> &chrLength, HaplotagParameters &params, VCF_Info *vcfSet, BamBaseCounter &NorBase){
    std::cerr << "calling somatic variants ... ";
    std::time_t begin = time(NULL);

    if(chrVec.size() == 0){
        std::cerr<<"ERROR: chrVec is empty\n";
        exit(1);
    }

    if(chrLength.size() == 0){
        std::cerr<<"ERROR: chrLength is empty\n";
        exit(1);        
    }

    for(auto chr : chrVec){
        (*SomaticChrPosInfo)[chr] = std::map<int, HP3_Info>();
        (*chrVarReadHpResult)[chr] = std::map<int, ReadHpResult>();
    }

    // somatic calling filter params
    SomaticFilterParaemter somaticParams;

    // setting somatic calling filter params
    InitialSomaticFilterParams(somaticParams);  

    // init data structure and get core n
    htsThreadPool threadPool = {NULL, 0};
    // creat thread pool
    if (!(threadPool.pool = hts_tpool_init(params.numThreads))) {
        fprintf(stderr, "Error creating thread pool\n");
        exit(1);
    }

    // loop all chromosome
    #pragma omp parallel for schedule(dynamic) num_threads(params.numThreads) 
    for(auto chr : chrVec ){
        
        // bam file resource allocation
        BamFileRAII bam(BamFile, params.fastaFile, threadPool);

        // record the position that tagged as HP3
        // chr, variant position
        std::map<int, HP3_Info> *somaticPosInfo = nullptr;
        // read ID, reads hpResult 
        std::map<int, ReadHpResult> *readHpDistributed = nullptr;
        // read ID, SNP HP count 
        std::map<std::string, ReadVarHpCount> readHpResultSet;
        // position, read ID, baseHP 
        std::map<int, std::map<std::string, int>> somaticPosReadHPCount;

        // records all variants within this chromosome.
        std::map<int, RefAltSet> currentChrVariants;

        #pragma omp critical
        {
            somaticPosInfo = &((*SomaticChrPosInfo)[chr]);
            readHpDistributed = &((*chrVarReadHpResult)[chr]);
            currentChrVariants = mergedChrVarinat[chr];
        }

        // since each read is sorted based on the start coordinates, to save time, 
        // firstVariantIter keeps track of the first variant that each read needs to check.
        std::map<int, RefAltSet>::iterator firstVariantIter = currentChrVariants.begin();
        // get the coordinates of the last variant
        // the tagging process will not be perform if the read's start coordinate are over than last variant.
        std::map<int, RefAltSet>::reverse_iterator lastVariant = currentChrVariants.rbegin();
        
        // tagging will be attempted for reads within the specified coordinate range.
        const std::string region = !params.region.empty() ? params.region : chr + ":1-" + std::to_string(chrLength[chr]);
        hts_itr_t* Iter = sam_itr_querys(bam.idx, bam.bamHdr, region.c_str());
        
        // iter all reads
        while (sam_itr_multi_next(bam.in, Iter, bam.aln) >= 0) {

            int flag = bam.aln->core.flag;

            //if ( aln->core.qual < params.qualityThreshold ){
            //    // mapping quality is lower than threshold
            //    continue;
            //}            

            if( (flag & 0x4) != 0 ){
                // read unmapped
                // write this alignment to result bam file
                continue;
            }
            if( (flag & 0x100) != 0 ){
                // secondary alignment. repeat.
                // A secondary alignment occurs when a given read could align reasonably well to more than one place.
                continue;
            }
            if( (flag & 0x800) != 0 && params.tagSupplementary == false){
                // supplementary alignment
                // A chimeric alignment is represented as a set of linear alignments that do not have large overlaps.
                continue;
            }
            //currentVariants rbegin == rend
            if(lastVariant == currentChrVariants.rend()){
                //skip
                continue;
            }
            else if(int(bam.aln->core.pos) <= (*lastVariant).first){
                //statistics for tumor SNP base count and depth, and classify cases
                StatisticSomaticPosInfo(*bam.bamHdr, *bam.aln, chr, params, &NorBase, vcfSet, *somaticPosInfo, currentChrVariants, firstVariantIter, readHpResultSet, somaticPosReadHPCount);
            }
        }
        //calculate information and filter somatic SNPs
        SomaticFeatureFilter(somaticParams, currentChrVariants, chr, *somaticPosInfo);

        //find other somatic variants HP4/HP5
        FindOtherSomaticSnpHP(chr, *somaticPosInfo, currentChrVariants);

        //calibrate read HP to remove low confidence H3/H4 SNP
        CalibrateReadHP(chr, somaticParams, *somaticPosInfo, readHpResultSet, somaticPosReadHPCount);

        //calculate all read HP result
        CalculateChrReadHP(params, chr, readHpResultSet, somaticPosReadHPCount);

        //statistic all read HP in somatic SNP position
        StatisticSomaticPosReadHP(chr, *somaticPosInfo, somaticPosReadHPCount, readHpResultSet, *readHpDistributed);

        //merge local read hp result to global read hp result
        #pragma omp critical
        {
            callerReadHpDistri->mergeLocalReadHp(chr, *readHpDistributed);
        }

        somaticPosInfo = nullptr;
        readHpDistributed = nullptr;
        hts_itr_destroy(Iter);
    }
    hts_tpool_destroy(threadPool.pool);
    std::cerr<< difftime(time(NULL), begin) << "s\n";
    
    if(somaticParams.writeVarLog){
        //write the log file for variants with positions tagged as HP3
        WriteSomaticVarCallingLog(params ,somaticParams, chrVec, NorBase, mergedChrVarinat);

        WriteOtherSomaticHpLog(params, chrVec, mergedChrVarinat);

        callerReadHpDistri->writeReadHpDistriLog(params, "_readDistri_Scaller.out", chrVec);

        // remove the position that derived by HP1 or HP2
        callerReadHpDistri->removeNotDeriveByH1andH2pos(chrVec);

        // record the position that can't derived because exist H1-1 and H2-1 reads
        callerReadHpDistri->writeReadHpDistriLog(params, "_readDistri_Scaller_derive_by_H1_H2.out", chrVec);
    }
    return;
}

void SomaticVarCaller::InitialSomaticFilterParams(SomaticFilterParaemter &somaticParams){
    
    // Determine whether to apply the filter
    somaticParams.applyFilter = true;
    somaticParams.writeVarLog = true;

    // Below the mapping quality read ratio threshold
    somaticParams.LowMpqRatioThreshold = 0.1;

    // Phased heterozygous SNPs filter threshold
    somaticParams.Hetero_OnlyHP3ReadRatioThreshold = 1.0;
    somaticParams.Hetero_MessyReadRatioThreshold = 1.0;
    somaticParams.Hetero_readCountThreshold = 1;
    somaticParams.Hetero_VAF_upper_threshold = 1.0;
    somaticParams.Hetero_VAF_lower_threshold = 0.1;
    somaticParams.Hetero_tumDeletionRatio = 0.6;

    // Unphased Heterozygous SNPs filter threshold
    somaticParams.Unphased_Hetero_OnlyHP3ReadRatioThreshold = 1.0;
    somaticParams.Unphased_Hetero_MessyReadRatioThreshold = 1.0;
    somaticParams.Unphased_Hetero_readCountThreshold = 1;
    somaticParams.Unphased_Hetero_VAF_upper_threshold = 1.0;
    somaticParams.Unphased_Hetero_VAF_lower_threshold = 0.1;
    somaticParams.Unphased_Hetero_tumDeletionRatio = 0.6;

    // Homozygous SNPs filter threshold
    somaticParams.Homo_OnlyHP3ReadRatioThreshold = 1.0;
    somaticParams.Homo_MessyReadRatioThreshold = 1.0;
    somaticParams.Homo_readCountThreshold = 1;
    somaticParams.Homo_VAF_upper_threshold = 1.1; //not used
    somaticParams.Homo_VAF_lower_threshold = 0.4;
    somaticParams.Homo_tumDeletionRatio = 0.4;
}

void SomaticVarCaller::SetSomaticFilterParams(const SomaticFilterParaemter &somaticParams, std::string GTtype
                            , float &OnlyHP3ReadRatioThreshold, float &messyReadRatioThreshold,int &readCountThreshold
                            , float &VAFthreshold_upper, float &VAFthreshold_lower, float &tumDeletionRatioThreshold, float &tumLowMpqRatioThreshold){

    // Below the mapping quality read ratio threshold
    tumLowMpqRatioThreshold = somaticParams.LowMpqRatioThreshold;

    //setting threshold by Genotype
    if(GTtype == "Hetero"){
        OnlyHP3ReadRatioThreshold = somaticParams.Hetero_OnlyHP3ReadRatioThreshold;
        messyReadRatioThreshold = somaticParams.Hetero_MessyReadRatioThreshold;
        readCountThreshold = somaticParams.Hetero_readCountThreshold;
        VAFthreshold_upper = somaticParams.Hetero_VAF_upper_threshold;
        VAFthreshold_lower = somaticParams.Hetero_VAF_lower_threshold;
        tumDeletionRatioThreshold = somaticParams.Hetero_tumDeletionRatio;
    }else if(GTtype == "UnphasedHetero"){
        OnlyHP3ReadRatioThreshold = somaticParams.Unphased_Hetero_OnlyHP3ReadRatioThreshold;
        messyReadRatioThreshold = somaticParams.Unphased_Hetero_MessyReadRatioThreshold;
        readCountThreshold = somaticParams.Unphased_Hetero_readCountThreshold;
        VAFthreshold_upper = somaticParams.Unphased_Hetero_VAF_upper_threshold;
        VAFthreshold_lower = somaticParams.Unphased_Hetero_VAF_lower_threshold;
        tumDeletionRatioThreshold = somaticParams.Unphased_Hetero_tumDeletionRatio;                       
    }else if(GTtype == "Homo"){
        OnlyHP3ReadRatioThreshold = somaticParams.Homo_OnlyHP3ReadRatioThreshold;
        messyReadRatioThreshold  = somaticParams.Homo_MessyReadRatioThreshold;
        readCountThreshold = somaticParams.Homo_readCountThreshold;
        VAFthreshold_upper = somaticParams.Homo_VAF_upper_threshold;
        VAFthreshold_lower = somaticParams.Homo_VAF_lower_threshold;
        tumDeletionRatioThreshold = somaticParams.Homo_tumDeletionRatio;
    }
}

void SomaticVarCaller::StatisticSomaticPosInfo(const  bam_hdr_t &bamHdr,const bam1_t &aln, const std::string &chr, HaplotagParameters &params, BamBaseCounter *NorBase, VCF_Info *vcfSet, std::map<int, HP3_Info> &somaticPosInfo, std::map<int, RefAltSet> &currentChrVariants, std::map<int, RefAltSet>::iterator &firstVariantIter, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, std::map<std::string, int>> &somaticPosReadHPCount){
    
    std::map<int, int> hpCount;
    hpCount[1] = 0; 
    hpCount[2] = 0;
    hpCount[3] = 0;
    hpCount[4] = 0;

    //record variants on this read
    std::map<int, int> variantsHP;

    //record tumor-unique variants with low VAF in the normal.bam on this read
    std::vector<int> readPosHP3;

    std::vector<int> readTumorOnlyPos;

    //record PS count( PS value, count)
    std::map<int, int> TumCountPS;
    std::map<int, int> NorCountPS;

    // Skip variants that are to the left of this read
    while( firstVariantIter != currentChrVariants.end() && (*(firstVariantIter)).first < aln.core.pos ){
        firstVariantIter++;
    }     

    if(firstVariantIter == currentChrVariants.end()){
        return;
    }

    // position relative to reference
    int ref_pos = aln.core.pos;
    // position relative to read
    int query_pos = 0;
    // set variant start for current alignment
    std::map<int, RefAltSet>::iterator currentVariantIter = firstVariantIter;

    // reading cigar to detect snp on this read
    int aln_core_n_cigar = int(aln.core.n_cigar);
    for(int i = 0; i < aln_core_n_cigar ; i++ ){

        uint32_t *cigar = bam_get_cigar(&aln);
        int cigar_op = bam_cigar_op(cigar[i]);
        int length   = bam_cigar_oplen(cigar[i]);

        // iterator next variant
        while( currentVariantIter != currentChrVariants.end() && (*currentVariantIter).first < ref_pos ){
            currentVariantIter++;
        }

        // CIGAR operators: MIDNSHP=X correspond 012345678
        // 0: alignment match (can be a sequence match or mismatch)
        // 7: sequence match
        // 8: sequence mismatch
        if( cigar_op == 0 || cigar_op == 7 || cigar_op == 8 ){    
            while( currentVariantIter != currentChrVariants.end() && (*currentVariantIter).first < ref_pos + length){

                int offset = (*currentVariantIter).first - ref_pos;

                if( offset < 0){
                }
                else{
                    uint8_t *q = bam_get_seq(&aln);
                    char base_chr = seq_nt16_str[bam_seqi(q,query_pos + offset)];
                    std::string base(1, base_chr);


                    //waring : using ref length to split SNP and indel that will be effect case ratio result 
                    if ( aln.core.qual >= params.somaticCallingMpqThreshold ){
                        SomaticJudgeSnpHP(currentVariantIter, vcfSet , chr, base, hpCount, NorCountPS, TumCountPS, &variantsHP, &readPosHP3, NorBase, &somaticPosInfo);
                        if((*currentVariantIter).second.isExistTumor && !((*currentVariantIter).second.isExistNormal)){
                            readTumorOnlyPos.push_back((*currentVariantIter).first);
                        }
                    }

                    //statistically analyze SNP information exclusive to the tumor
                    if((*currentVariantIter).second.isExistTumor && !((*currentVariantIter).second.isExistNormal)){
                        int tumRefLength = (*currentVariantIter).second.Variant[Genome::TUMOR].Ref.length();
                        int tumAltLength = (*currentVariantIter).second.Variant[Genome::TUMOR].Alt.length();

                        if(tumRefLength == 1 && tumAltLength == 1){
                            if ( aln.core.qual >= params.somaticCallingMpqThreshold ){
                                //counting current tumor SNP basd and depth with filtered MPQ
                                if(base == "A"){
                                    somaticPosInfo[(*currentVariantIter).first].base.MPQ_A_count++;
                                }else if(base == "C"){
                                    somaticPosInfo[(*currentVariantIter).first].base.MPQ_C_count++;
                                }else if(base == "G"){
                                    somaticPosInfo[(*currentVariantIter).first].base.MPQ_G_count++;
                                }else if(base == "T"){
                                    somaticPosInfo[(*currentVariantIter).first].base.MPQ_T_count++;
                                }else{
                                    somaticPosInfo[(*currentVariantIter).first].base.MPQ_unknow++;
                                }
                                somaticPosInfo[(*currentVariantIter).first].base.filteredMpqDepth++;
                            } 
                            //counting current tumor SNP basd and depth without filtered MPQ
                            if(base == "A"){
                                somaticPosInfo[(*currentVariantIter).first].base.A_count++;
                            }else if(base == "C"){
                                somaticPosInfo[(*currentVariantIter).first].base.C_count++;
                            }else if(base == "G"){
                                somaticPosInfo[(*currentVariantIter).first].base.G_count++;
                            }else if(base == "T"){
                                somaticPosInfo[(*currentVariantIter).first].base.T_count++;
                            }else{
                                somaticPosInfo[(*currentVariantIter).first].base.unknow++;
                            }
                        }
                        somaticPosInfo[(*currentVariantIter).first].base.depth++;
                    }
                }
                currentVariantIter++;
            }
            query_pos += length;
            ref_pos += length;
        }
        // 1: insertion to the reference
        else if( cigar_op == 1 ){
            query_pos += length;
        }
        // 2: deletion from the reference
        else if( cigar_op == 2 ){
            while( currentVariantIter != currentChrVariants.end() && (*currentVariantIter).first < ref_pos + length){
                
                //statistically analyze SNP information exclusive to the tumor
                if((*currentVariantIter).second.isExistTumor && !((*currentVariantIter).second.isExistNormal)){

                    int curPos = (*currentVariantIter).first;

                    //detect ref bsae length(temp :tumor SNP)
                    int tumRefLength = (*currentVariantIter).second.Variant[Genome::TUMOR].Ref.length();
                    int tumAltLength = (*currentVariantIter).second.Variant[Genome::TUMOR].Alt.length();

                    if(tumRefLength == 1 && tumAltLength == 1){
                        somaticPosInfo[curPos].base.delCount++;
                    }
                    // the indel SNP start position isn't at the end of the deletion
                    else if(tumRefLength > 1 && tumAltLength == 1){

                    }
                }
                currentVariantIter++;
            }
            ref_pos += length;
        // 3: skipped region from the reference
        }else if( cigar_op == 3 ){
            ref_pos += length;
        }
        // 4: soft clipping (clipped sequences present in SEQ)
        else if( cigar_op == 4 ){
            query_pos += length;
        }
        // 5: hard clipping (clipped sequences NOT present in SEQ)
        // 6: padding (silent deletion from padded reference)
        else if( cigar_op == 5 || cigar_op == 6 ){
        // do nothing
        }
        else{
            std::cerr<< "alignment find unsupported CIGAR operation from read: " << bam_get_qname(&aln) << "\n";
            exit(1);
        }
    }


    //classify read cases where tumor SNPs have low VAF in normal samples
    if(!readPosHP3.empty()){
        ClassifyReadsByCase(readPosHP3, NorCountPS, hpCount, params, somaticPosInfo);
    }

    //record variants HP count for each read in tumor-only position
    if(!readTumorOnlyPos.empty()){

        std::string readID = bam_get_qname(&aln);
        //read ID overide
        if(readHpResultSet.find(readID) != readHpResultSet.end()){
            readHpResultSet[readID].readIDcount++;
            readID = readID + "-" + std::to_string(readHpResultSet[readID].readIDcount++);
        }

        readHpResultSet[readID].HP1= hpCount[1];
        readHpResultSet[readID].HP2= hpCount[2];
        readHpResultSet[readID].HP3= hpCount[3];
        readHpResultSet[readID].HP4= hpCount[4];

        readHpResultSet[readID].NorCountPS = NorCountPS;
        readHpResultSet[readID].startPos = aln.core.pos + 1;
        readHpResultSet[readID].endPos = ref_pos;
        readHpResultSet[readID].readLength = query_pos;
        
        //if(readID == "SRR25005626.11816585"){
        //    std::cerr << "readID: "<< readID << " chr: "<< chr << std::endl;
        //    std::cout << "HP1: "<< readTotalHPcount[readID].HP1 << " HP2: "<< readTotalHPcount[readID].HP2 << " HP3: "<< readTotalHPcount[readID].HP3 << " HP4: "<< readTotalHPcount[readID].HP4  << std::endl;
        //}

        for(auto pos : readTumorOnlyPos){
            
            int somaticSnpBaseHP = SnpHP::NONE_SNP;

            if(variantsHP.find(pos) != variantsHP.end()){
                somaticSnpBaseHP = variantsHP[pos];
            }

            somaticPosReadHPCount[pos][readID] = somaticSnpBaseHP;
        } 
    }
}

void SomaticVarCaller::OnlyTumorSNPjudgeHP(const std::string &chrName, int &curPos, RefAltSet &curVar, std::string base, VCF_Info *vcfSet, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *readPosHP3, BamBaseCounter *NorBase, std::map<int, HP3_Info> *SomaticPos){
    //the tumor SNP GT is phased heterozygous
    //all bases of the same type at the current position in normal.bam

    if(readPosHP3 == nullptr){
        std::cerr << "ERROR (SomaticDetectJudgeHP) => readPosHP3 pointer cannot be nullptr"<< std::endl;
        exit(1);
    }
    if(SomaticPos == nullptr){
        std::cerr << "ERROR (SomaticDetectJudgeHP) => SomaticPos pointer cannot be nullptr"<< std::endl;
        exit(1);
    }
    
    if((*NorBase).isHighRefAllelleFreq(chrName, curPos) == true){
        std::string TumorRefBase = curVar.Variant[Genome::TUMOR].Ref;
        std::string TumorAltBase = curVar.Variant[Genome::TUMOR].Alt;

        //max count base match to refBase in normal.bam
        if((*NorBase).getMaxFreqBase(chrName, curPos) == TumorRefBase){
            if(base == TumorAltBase){
                hpCount[3]++;
                if(variantsHP != nullptr) (*variantsHP)[curPos] = SnpHP::SOMATIC_H3;

                //record postions that tagged as HP3 for calculating the confidence of somatic positions
                (*readPosHP3).push_back(curPos);
                (*SomaticPos)[curPos].isNormalPosLowVAF = true;                                

            //base is not match to TumorRefBase & TumorAltBase (other HP)
            }else if(base != TumorRefBase && base != TumorAltBase){
                //hpCount[4]++;
                //if(variantsHP != nullptr) (*variantsHP)[curPos] = SnpHP::SOMATIC_H4;
                //(*readPosHP3).push_back(curPos);
                //(*SomaticPos)[curPos].isNormalPosLowVAF = true;  
            }

            if(tumCountPS != nullptr) (*tumCountPS)[vcfSet[Genome::TUMOR].chrVariantPS[chrName][curPos]]++;

        //max count base not match to tumorRefBase in normal.bam
        }else if((*NorBase).getMaxFreqBase(chrName, curPos) != TumorRefBase){
            //temp 
        }
    //exist more than one type of base at the current position in normal.bam 
    }else{
        //temp
    }
}

void SomaticVarCaller::ClassifyReadsByCase(std::vector<int> &readPosHP3, std::map<int, int> &NorCountPS, std::map<int, int> &hpCount, const HaplotagParameters &params, std::map<int, HP3_Info> &somaticPosInfo){
    
    /*double max, min;
    if(hpCount[3] >= hpCount[4]){
        max = hpCount[3];
        min = hpCount[4];
    }else{
        max = hpCount[4];
        min = hpCount[3];
    }*/

    //decide whether to tag the read or not
    bool recordRead = true;
    /*if( max/(max+min) < params.percentageThreshold){
        // no tag
        recordRead = false;
    }*/

    //only exist tumor homo SNP
    if( NorCountPS.size() == 0 && hpCount[3] != 0){
        //std::cerr << "Error : read only had tumor homo SNP " << bam_get_qname(aln) << "\n";
    }

    // germline SNPs cross two block
    if( NorCountPS.size() > 1){
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

    for (int somaticPos : readPosHP3) {
        if(recordRead == false){
            somaticPosInfo[somaticPos].unTag++;
        }else if(tagCleanHP3Read){
            somaticPosInfo[somaticPos].totalCleanHP3Read++;
            if(hpCount[1] == 0 && hpCount[2] == 0 && hpCount[3] != 0){
                somaticPosInfo[somaticPos].OnlyHP3Read++;
            }else if(hpCount[1] !=0  && hpCount[2] == 0){ 
                somaticPosInfo[somaticPos].HP1withHP3Read++;      
            }else if(hpCount[1] == 0 && hpCount[2] != 0){ 
                somaticPosInfo[somaticPos].HP2withHP3Read++;                            
            }
        }else{
            somaticPosInfo[somaticPos].MessyHPRead++;
        }
    }
}

void SomaticVarCaller::SomaticFeatureFilter(const SomaticFilterParaemter &somaticParams, std::map<int, RefAltSet> &currentChrVariants,const std::string &chr, std::map<int, HP3_Info> &somaticPosInfo){
//calculate the information of the somatic positon 
    std::map<int, HP3_Info>::iterator somaticVarIter = somaticPosInfo.begin();
    while( somaticVarIter != somaticPosInfo.end()){

        //Stage 1 filter
        if((*somaticVarIter).second.isNormalPosLowVAF == false){
            somaticVarIter++;
            continue;
        }

        //calculating Read Case Ratio
        int totalCleanHP3Read = (*somaticVarIter).second.totalCleanHP3Read;
        int totalHP1WithHP3Read = (*somaticVarIter).second.HP1withHP3Read;
        int totalHP2WithHP3Read = (*somaticVarIter).second.HP2withHP3Read;
        int totalOnlyHP3Read = (*somaticVarIter).second.OnlyHP3Read;
        int totalMessyHPRead = (*somaticVarIter).second.MessyHPRead;

        (*somaticVarIter).second.CaseReadCount = totalCleanHP3Read + totalMessyHPRead;

        (*somaticVarIter).second.MessyHPReadRatio = (float)totalMessyHPRead / ((float)totalCleanHP3Read + (float)totalMessyHPRead);

        (*somaticVarIter).second.HP1withHP3ReadRatio = (float)totalHP1WithHP3Read / ((float)totalCleanHP3Read + (float)totalMessyHPRead);
        (*somaticVarIter).second.HP2WithHP3ReadRatio = (float)totalHP2WithHP3Read / ((float)totalCleanHP3Read + (float)totalMessyHPRead);
        (*somaticVarIter).second.OnlyHP3ReadRatio = (float)totalOnlyHP3Read / ((float)totalCleanHP3Read + (float)totalMessyHPRead);


        //get tumor depth & filtered MPQ tumor depth
        int tumDepth = (*somaticVarIter).second.base.depth;
        int tumMpqDepth = (*somaticVarIter).second.base.filteredMpqDepth;

        //calculate VAF
        int tumAltCount = 0;
        std::string RefBase;
        std::string AltBase;
                
        if(currentChrVariants[(*somaticVarIter).first].isExistTumor == true){
            RefBase = currentChrVariants[(*somaticVarIter).first].Variant[Genome::TUMOR].Ref;
            AltBase = currentChrVariants[(*somaticVarIter).first].Variant[Genome::TUMOR].Alt;
        }else{
            std::cerr << "Error(calculate tumor VAF) => can't find the position : chr:" << chr << " pos: " << (*somaticVarIter).first;
            exit(1);
        }

        if(RefBase == "" || AltBase == ""){
            std::cerr << "Error(calculate tumor VAF) => can't find RefBase or AltBase : chr:" << chr << " pos: " << (*somaticVarIter).first << " RefBase:" << RefBase << " AltBase:" << AltBase;
            exit(1);
        }

        if(AltBase == "A"){
            //AltCount = TumBase.getBaseAcount(chr, (*somaticVarIter).first);
            tumAltCount = somaticPosInfo[(*somaticVarIter).first].base.A_count;
        }else if(AltBase == "C"){
            //AltCount = TumBase.getBaseCcount(chr, (*somaticVarIter).first);
            tumAltCount = somaticPosInfo[(*somaticVarIter).first].base.C_count;
        }else if(AltBase == "T"){
            //AltCount = TumBase.getBaseTcount(chr, (*somaticVarIter).first);
            tumAltCount = somaticPosInfo[(*somaticVarIter).first].base.T_count;
        }else if(AltBase == "G"){
            //AltCount = TumBase.getBaseGcount(chr, (*somaticVarIter).first);
            tumAltCount = somaticPosInfo[(*somaticVarIter).first].base.G_count;
        }else{
            std::cerr << "Error(calculate tumor VAF) => can't match RefBase or AltBase : chr:" << chr << " pos: " << ((*somaticVarIter).first)+1 << " RefBase:" << RefBase << " AltBase:" << AltBase;
            exit(1);
        }

        (*somaticVarIter).second.base.VAF =  (float)tumAltCount / (float)tumDepth;

        // Calculating the difference in VAF between VAF and filtered low MPQ VAF
        (*somaticVarIter).second.base.lowMpqReadRatio = (float)(tumDepth - tumMpqDepth) / (float)tumDepth;

        if((*somaticVarIter).second.base.lowMpqReadRatio < 0){
            std::cerr << "Error(tumMPQReadRatio) => chr: " << chr << " pos: " << ((*somaticVarIter).first) + 1;
            std::cerr << "tumMPQReadRatio " << (*somaticVarIter).second.base.lowMpqReadRatio << " norDepth: " << tumDepth << " norMPQDepth: " << tumMpqDepth;
        }
                
        //current SNP GT type
        if(currentChrVariants[(*somaticVarIter).first].Variant[Genome::TUMOR].is_homozygous == true){
            (*somaticVarIter).second.GTtype = "Homo";
        }else if(currentChrVariants[(*somaticVarIter).first].Variant[Genome::TUMOR].is_phased_hetero == true){
            (*somaticVarIter).second.GTtype = "Hetero";
        }else if(currentChrVariants[(*somaticVarIter).first].Variant[Genome::TUMOR].is_unphased_hetero == true){
            (*somaticVarIter).second.GTtype = "UnphasedHetero";
        }else{
            std::cerr << "Error(GTtype) => can't find GTtype at chr: " << chr << " pos: " << ((*somaticVarIter).first) + 1;
            exit(1);  
        }

        std::string GTtype = (*somaticVarIter).second.GTtype;
            
        int tumDeletionCount = somaticPosInfo[(*somaticVarIter).first].base.delCount;

        //the ratio of deletions occurring
        if(tumDeletionCount == 0){
            (*somaticVarIter).second.tumDelRatio = 0.0;
        }else{
            (*somaticVarIter).second.tumDelRatio = (float)tumDeletionCount / ((float)tumDeletionCount + (float)tumDepth);
        }

        //filter threshold 
        float OnlyHP3ReadRatioThreshold;
        float messyReadRatioThreshold;
        int readCountThreshold;
        float VAF_upper_threshold;
        float VAF_lower_threshold;
        float tumDeletionRatioThreshold;
        float tumLowMpqRatioThreshold;

        //setting somatic calling filter parameter
        SetSomaticFilterParams(somaticParams, GTtype, OnlyHP3ReadRatioThreshold
                             , messyReadRatioThreshold, readCountThreshold, VAF_upper_threshold
                             , VAF_lower_threshold, tumDeletionRatioThreshold, tumLowMpqRatioThreshold);

        //stage 2 filter
        if( (!somaticParams.applyFilter) ||
            (
            //((*somaticVarIter).second.OnlyHP3ReadRatio < OnlyHP3ReadRatioThreshold) && 
            ((*somaticVarIter).second.MessyHPReadRatio < messyReadRatioThreshold) && 
            ((*somaticVarIter).second.CaseReadCount > readCountThreshold) &&
            ((*somaticVarIter).second.base.VAF > VAF_lower_threshold) && 
            ((*somaticVarIter).second.base.VAF < VAF_upper_threshold) &&
            ((*somaticVarIter).second.tumDelRatio < tumDeletionRatioThreshold) &&
            ((*somaticVarIter).second.base.lowMpqReadRatio <= tumLowMpqRatioThreshold)
            ) ){
              
            (*somaticVarIter).second.isHighConSomaticSNP = true;
        }

        somaticVarIter++;          
    }
}

void SomaticVarCaller::FindOtherSomaticSnpHP(const std::string &chr, std::map<int, HP3_Info> &somaticPosInfo, std::map<int, RefAltSet> &currentChrVariants){
    std::map<int, HP3_Info>::iterator somaticVarIter = somaticPosInfo.begin();
    while(somaticVarIter != somaticPosInfo.end()){
        if((*somaticVarIter).second.isHighConSomaticSNP){
            int pos = (*somaticVarIter).first;

            std::map<int, int> NucCount;
            NucCount[Nucleotide::A] = (*somaticVarIter).second.base.MPQ_A_count;
            NucCount[Nucleotide::C] = (*somaticVarIter).second.base.MPQ_C_count;
            NucCount[Nucleotide::T] = (*somaticVarIter).second.base.MPQ_T_count;
            NucCount[Nucleotide::G] = (*somaticVarIter).second.base.MPQ_G_count;

            int refAllele;
            int altAllele;

            if(currentChrVariants.find(pos) != currentChrVariants.end()){
                refAllele = convertStrNucToInt(currentChrVariants[pos].Variant[Genome::TUMOR].Ref);
                altAllele = convertStrNucToInt(currentChrVariants[pos].Variant[Genome::TUMOR].Alt);
            }else{
                std::cerr << "Error(FindOtherSomaticSnpHP) => can't find position in currentChrVariants : chr:" << chr << " pos: " << pos + 1;
                exit(1);
            }

            NucCount.erase(refAllele);
            NucCount.erase(altAllele);

            int maxNuc = Nucleotide::UNKOWN;
            int maxNucCount = 0;

            int minNuc = Nucleotide::UNKOWN;
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

            if(maxNuc != Nucleotide::UNKOWN){
                (*somaticVarIter).second.somaticHp4Base = maxNuc;
                (*somaticVarIter).second.somaticHp4BaseCount = maxNucCount;
            }

            if(minNuc != Nucleotide::UNKOWN){
                (*somaticVarIter).second.somaticHp5Base = minNuc;
                (*somaticVarIter).second.somaticHp5BaseCount = minNucCount;
            }

            // if(maxNuc != Nucleotide::UNKOWN && minNuc == Nucleotide::UNKOWN){
            //     std::cout << "chr: " << chr << " pos: " << pos + 1 << "\n";
            //     std::cout << "maxNuc: " << convertIntNucToStr(maxNuc) << " minNuc: " << convertIntNucToStr(minNuc) << "\n";
            //     std::cout << "maxNucCount: " << maxNucCount << " minNucCount: " << minNucCount << "\n";
            // }
        }
        somaticVarIter++;
    }
}

void SomaticVarCaller::CalibrateReadHP(const std::string &chr, const SomaticFilterParaemter &somaticParams, std::map<int, HP3_Info> &somaticPosInfo, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, std::map<std::string, int>> &somaticPosReadHPCount){
        //calibrate read HP to remove low confidence H3/H4 SNP
        std::map<int, HP3_Info>::iterator somaticVarIter = somaticPosInfo.begin();
        while( somaticVarIter != somaticPosInfo.end()){
            if(!(*somaticVarIter).second.isHighConSomaticSNP){
                int pos = (*somaticVarIter).first;

                //only reads with HP3 or HP4 SNPs will be calibrated
                if(somaticPosReadHPCount.find(pos) != somaticPosReadHPCount.end()){

                    for(auto readIdBaseHP : somaticPosReadHPCount[pos]){
                        std::string readID = readIdBaseHP.first;
                        int baseHP = readIdBaseHP.second;

                        switch (baseHP) {
                            case SnpHP::GERMLINE_H1:
                            case SnpHP::GERMLINE_H2:
                                std::cerr << "ERROR (calibrate read HP) => somatic Variant had HP1 or HP2 SNP :" << std::endl;
                                std::cerr << "readID: " << readID << " chr: " << chr << " pos: " << pos + 1 << std::endl;
                                std::cerr << "baseHP: " << baseHP << std::endl;
                                exit(1);
                            case SnpHP::SOMATIC_H3:
                                readHpResultSet[readID].HP3--; break;
                            case SnpHP::SOMATIC_H4:
                                readHpResultSet[readID].HP4--; break;
                            default:
                                break;
                        }

                        if(readHpResultSet[readID].HP3 < 0 || readHpResultSet[readID].HP4 < 0){
                            std::cerr << "ERROR (calibrate read HP) => read HP3 or HP4 SNP count < 0 :" << std::endl;
                            std::cerr << "readID: "<< readID << " chr: "<< chr << " pos: " << pos+1<< std::endl;
                            std::cerr << "HP3: "<< readHpResultSet[readID].HP3 << " HP4: "<< readHpResultSet[readID].HP3  << std::endl;
                            exit(1);
                        }
                    }
                }else{
                    //std::cerr << "ERROR (calibrate read HP) => can't find pos in somaticPosReadHPCount : chr: "<< chr << " pos: " << pos+1 <<std::endl;
                    //exit(1);
                }
            }
            somaticVarIter++;
        }
}

void SomaticVarCaller::CalculateChrReadHP(const HaplotagParameters &params, const std::string &chr, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, std::map<std::string, int>> &somaticPosReadHPCount){
        std::map<std::string, ReadVarHpCount>::iterator readTotalHPcountIter = readHpResultSet.begin();
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
            std::map<int, int> NorCountPS = (*readTotalHPcountIter).second.NorCountPS;

            (*readTotalHPcountIter).second.hpResult = determineReadHP(hpCount, pqValue, NorCountPS, normalHPsimilarity, tumorHPsimilarity, params.percentageThreshold, nullptr, nullptr, nullptr);
            
            readTotalHPcountIter++;
        }
}

void SomaticVarCaller::StatisticSomaticPosReadHP(const std::string &chr, std::map<int, HP3_Info> &somaticPosInfo, std::map<int, std::map<std::string, int>> &somaticPosReadHPCount, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, ReadHpResult> &readHpDistributed){
    std::map<int, HP3_Info>::iterator somaticVarIter = somaticPosInfo.begin();
    while(somaticVarIter != somaticPosInfo.end()){
        if((*somaticVarIter).second.isHighConSomaticSNP){
            int pos = (*somaticVarIter).first;

            if(somaticPosReadHPCount.find(pos) != somaticPosReadHPCount.end()){
                readHpDistributed[pos] = ReadHpResult();

                // record the number of HP1-1 or HP2-1 derived from Base HP3
                std::map<int, int> deriveByHPfromBaseHp3;
                deriveByHPfromBaseHp3[ReadHP::H1_1] = 0;
                deriveByHPfromBaseHp3[ReadHP::H2_1] = 0;

                for(std::pair<std::string, int> readIdBaseHP : somaticPosReadHPCount[pos]){
                    std::string readID = readIdBaseHP.first;
                    int baseHP = readIdBaseHP.second;
                    int hpResult = readHpResultSet[readID].hpResult;

                    recordReadHp(pos, hpResult, baseHP, readHpDistributed);

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

                readHpDistributed[pos].existDeriveByH1andH2 = false;

                if(HP1_1ratio >= 1.0){
                    (*somaticVarIter).second.somaticReadDeriveByHP = SnpHP::GERMLINE_H1;
                }else if(HP2_1ratio >= 1.0){
                    (*somaticVarIter).second.somaticReadDeriveByHP = SnpHP::GERMLINE_H2;
                }else{
                    (*somaticVarIter).second.somaticReadDeriveByHP = SnpHP::NONE_SNP;

                    if((0 < HP1_1ratio && HP1_1ratio < 1.0)  || (0 < HP2_1ratio && HP2_1ratio < 1.0)){
                        readHpDistributed[pos].existDeriveByH1andH2 = true;
                        // std::cerr << "ERROR (statistic all read HP) => somatic Variant had HP1-1 and HP2-1 reads :"<< std::endl;
                        // std::cerr << "chr: "<< chr << " pos: " << pos+1 <<std::endl;
                        // std::cerr << "HP1-1_ratio: "<< HP1_1ratio << " HP2-1_ratio: "<< HP2_1ratio << std::endl;
                        // std::cerr << "HP1-1: "<< deriveByHPfromBaseHp3["1-1"] << " HP2-1: "<< deriveByHPfromBaseHp3["2-1"] << " HP3: "<< readHpDistributed[pos].HP3read << " HP4: "<< readHpDistributed[pos].HP4read<< std::endl;
                    }
                }

                if( readHpDistributed[pos].hpResultCounter[ReadHP::H3] == 0 && 
                    readHpDistributed[pos].hpResultCounter[ReadHP::H4] == 0 && 
                    readHpDistributed[pos].hpResultCounter[ReadHP::H1_1] == 0 &&
                    readHpDistributed[pos].hpResultCounter[ReadHP::H1_2] == 0 &&
                    readHpDistributed[pos].hpResultCounter[ReadHP::H2_1] == 0 &&  
                    readHpDistributed[pos].hpResultCounter[ReadHP::H2_2] == 0){
                    std::cerr << "ERROR (statistic all read HP) => hadn't exist somatic HP read in : chr: "<< chr << " pos: " << pos+1 <<std::endl;
                    std::cerr << "HP1-1: "<< readHpDistributed[pos].hpResultCounter[ReadHP::H1_1] 
                              << " HP1-2: "<< readHpDistributed[pos].hpResultCounter[ReadHP::H1_2] 
                              << " HP2-1: "<< readHpDistributed[pos].hpResultCounter[ReadHP::H2_1] 
                              << " HP2-2: "<< readHpDistributed[pos].hpResultCounter[ReadHP::H2_2] 
                              << " HP3: "<< readHpDistributed[pos].hpResultCounter[ReadHP::H3] 
                              << " HP4: "<< readHpDistributed[pos].hpResultCounter[ReadHP::H4]
                              << std::endl;
                    exit(1); 
                }

            }else{
                std::cerr << "ERROR (statistic all read HP) => can't find pos in somaticPosReadHPCount : chr: "<< chr << " pos: " << pos+1 <<std::endl;
                exit(1); 
            }

        }
        somaticVarIter++;
    }
}


void SomaticVarCaller::WriteSomaticVarCallingLog(const HaplotagParameters &params, const SomaticFilterParaemter &somaticParams, const std::vector<std::string> &chrVec, BamBaseCounter &NorBase, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat){
    std::ofstream *tagHP3Log = new std::ofstream(params.resultPrefix+"_HP3.out");

    if(!tagHP3Log->is_open()){
        std::cerr<< "Fail to open write file: " << params.resultPrefix+"_HP3.out" << "\n";
        exit(1);
    }

    std::cout << "writing somatic varinats calling log ... ";
    std::time_t begin = time(NULL);

    int totalSomaticSNP = 0;
    for(auto chr : chrVec){
        std::map<int, HP3_Info>::iterator somaticVarIter = (*SomaticChrPosInfo)[chr].begin();
        while( somaticVarIter != (*SomaticChrPosInfo)[chr].end()){
            if((*somaticVarIter).second.isHighConSomaticSNP){
                totalSomaticSNP++;
            }
            somaticVarIter++;
        }
    }

    //write header
    (*tagHP3Log) << "#####################################\n"
                 << "#   Somatic Varirants Calling Log   #\n"
                 << "#####################################\n";
    (*tagHP3Log) << "##NormalSnpFile:"       << params.snpFile             << "\n"
                 << "##TumorSnpFile:"        << params.tumorSnpFile        << "\n" 
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
    (*tagHP3Log) << "##----------- Filter Parameters -----------\n"
                 << "##\n"
                 << "##  Filtered : " << somaticParams.applyFilter << "\n"
                 << "##  Calling mapping quality :" << params.somaticCallingMpqThreshold << "\n"
                 << "##  Low mapping quality read ratio threshold : " << somaticParams.LowMpqRatioThreshold << "\n"
                 << "##\n"; 
    (*tagHP3Log) << "##---Phased Heterozygous SNP---\n"
                 //<< "##  Only HP3 Read Ratio threshold : " << somaticParams.Unphased_Hetero_OnlyHP3ReadRatioThreshold << "\n"
                 << "##  Messy Read Ratio threshold : " << somaticParams.Unphased_Hetero_MessyReadRatioThreshold << "\n"
                 << "##  Read Count threshold : " << somaticParams.Unphased_Hetero_readCountThreshold << "\n"
                 << "##  VAF upper bound threshold : " << somaticParams.Unphased_Hetero_VAF_upper_threshold << "\n"
                 << "##  VAF lower bound threshold : " << somaticParams.Unphased_Hetero_VAF_lower_threshold << "\n"
                 << "##  Tumor deletion ratio threshold : " << somaticParams.Unphased_Hetero_tumDeletionRatio << "\n"
                 << "##\n"; 
    (*tagHP3Log) << "##---Unphased Heterozygous SNP---\n"
                 //<< "##  Only HP3 Read Ratio threshold : " << somaticParams.Hetero_OnlyHP3ReadRatioThreshold << "\n"
                 << "##  Messy Read Ratio threshold : " << somaticParams.Hetero_MessyReadRatioThreshold << "\n"
                 << "##  Read Count threshold : " << somaticParams.Hetero_readCountThreshold << "\n"
                 << "##  VAF upper bound threshold : " << somaticParams.Hetero_VAF_upper_threshold << "\n"
                 << "##  VAF lower bound threshold : " << somaticParams.Hetero_VAF_lower_threshold << "\n"
                 << "##  Tumor deletion ratio threshold : " << somaticParams.Hetero_tumDeletionRatio << "\n"
                 << "##\n"; 
    (*tagHP3Log) << "##---Homozygous SNP---\n"
                 //<< "##  Only HP3 Read Ratio threshold : " << somaticParams.Homo_OnlyHP3ReadRatioThreshold << "\n"
                 << "##  Messy Read Ratio threshold : " << somaticParams.Homo_MessyReadRatioThreshold << "\n"
                 << "##  Read Count threshold : " << somaticParams.Homo_readCountThreshold << "\n"
                 << "##  VAF upper bound threshold : " << somaticParams.Homo_VAF_upper_threshold << "\n"
                 << "##  VAF lower bound threshold : " << somaticParams.Homo_VAF_lower_threshold << "\n"
                 << "##  Tumor deletion ratio threshold : " << somaticParams.Homo_tumDeletionRatio << "\n"
                 << "##---------------------------------------- \n"
                 << "##\n"
                 << "##Total Somatic SNPs: " << totalSomaticSNP << "\n"
                 << "##\n"; 
    (*tagHP3Log) << "#CHROM\t" 
                 << "POS\t" 
                 << "ReadCount\t" 
                 << "(HP1withHP3,HP2withHP3,OnlyHP3,MessyHP)\t"
                 << "unTag\t" 
                 << "HP1withHP3Ratio\t" 
                 << "HP2withHP3Ratio\t" 
                 << "OnlyHP3Ratio\t" 
                 << "MessyReadRatio\t" 
                 << "norVAF\t"  
                 << "tumVAF\t"
                 << "subtract Depth\t" 
                 << "NorDepth\t"  
                 << "TumDepth\t"
                 << "NorDeletionCount\t"
                 << "TumDeletionCount\t"
                 << "norDeletionRatio\t"
                 << "tumDeletionRatio\t"
                 << "norMpqVAF\t"
                 << "tumMpqVAF\t"
                 << "norVAF_substract\t"
                 << "tumVAF_substract\t"
                 << "norMpqReadRatio\t"
                 << "tumMpqReadRatio\t"  
                 << "GT\n";

    //write variants information
    for(auto chr : chrVec){
    
        std::map<int, HP3_Info>::iterator somaticVarIter = (*SomaticChrPosInfo)[chr].begin();
        while( somaticVarIter != (*SomaticChrPosInfo)[chr].end()){
            
            if((*somaticVarIter).second.isHighConSomaticSNP == false){
                somaticVarIter++;
                continue;
            }

            //calculating Case Ratio
            int totalHP1WithHP3Read = (*somaticVarIter).second.HP1withHP3Read;
            int totalHP2WithHP3Read = (*somaticVarIter).second.HP2withHP3Read;
            int totalOnlyHP3Read = (*somaticVarIter).second.OnlyHP3Read;
            int totalMessyHPRead = (*somaticVarIter).second.MessyHPRead;
            int unTagCount = (*somaticVarIter).second.unTag;
            int readCount = (*somaticVarIter).second.CaseReadCount;

            float messyHPReadRatio = (*somaticVarIter).second.MessyHPReadRatio;

            float HP1withHP3ReadRatio = (*somaticVarIter).second.HP1withHP3ReadRatio;
            float HP2WithHP3ReadRatio = (*somaticVarIter).second.HP2WithHP3ReadRatio;
            float OnlyHP3ReadRatio = (*somaticVarIter).second.OnlyHP3ReadRatio;

            int norDepth = NorBase.getDepth(chr, (*somaticVarIter).first);
            int norMpqDepth = NorBase.getMpqDepth(chr, (*somaticVarIter).first);

            int tumDepth = (*somaticVarIter).second.base.depth;
            int tumMpqDepth = (*somaticVarIter).second.base.filteredMpqDepth;
           
            //calculate the Subtract in depth between normal and tumor
            int subtractDepth = tumDepth - norDepth;

            //calculate VAF
            int tumMpqAltCount = 0;
            int norAltCount = 0;
            std::string RefBase;
            std::string AltBase;
                
            if(mergedChrVarinat[chr][(*somaticVarIter).first].isExistTumor == true){
                RefBase = mergedChrVarinat[chr][(*somaticVarIter).first].Variant[Genome::TUMOR].Ref;
                AltBase = mergedChrVarinat[chr][(*somaticVarIter).first].Variant[Genome::TUMOR].Alt;
            }else{
                std::cerr << "Error(write tag HP3 log file) => can't find the position : chr:" << chr << " pos: " << ((*somaticVarIter).first) + 1;
                exit(1);
            }

            if(RefBase == "" || AltBase == ""){
                std::cerr << "Error(write tag HP3 log file) => can't find RefBase or AltBase : chr:" << chr << " pos: " << ((*somaticVarIter).first) + 1 << " RefBase:" << RefBase << " AltBase:" << AltBase;
                exit(1);
            }

            if(AltBase == "A"){
                tumMpqAltCount = (*SomaticChrPosInfo)[chr][(*somaticVarIter).first].base.MPQ_A_count;
                norAltCount = NorBase.getBaseAcount(chr, (*somaticVarIter).first);
            }else if(AltBase == "C"){
                tumMpqAltCount = (*SomaticChrPosInfo)[chr][(*somaticVarIter).first].base.MPQ_C_count;
                norAltCount = NorBase.getBaseCcount(chr, (*somaticVarIter).first);
            }else if(AltBase == "T"){
                tumMpqAltCount = (*SomaticChrPosInfo)[chr][(*somaticVarIter).first].base.MPQ_T_count;
                norAltCount = NorBase.getBaseTcount(chr, (*somaticVarIter).first);
            }else if(AltBase == "G"){
                tumMpqAltCount = (*SomaticChrPosInfo)[chr][(*somaticVarIter).first].base.MPQ_G_count;
                norAltCount = NorBase.getBaseGcount(chr, (*somaticVarIter).first);
            }else{
                std::cerr << "Error(write tag HP3 log file) => can't match RefBase or AltBase : chr:" << chr << " pos: " << ((*somaticVarIter).first) + 1 << " RefBase:" << RefBase << " AltBase:" << AltBase;
                exit(1);
            }

            //calculating tumor VAF
            float tumVAF = (*somaticVarIter).second.base.VAF;
            float tumMpqVAF = (float)tumMpqAltCount / (float)tumMpqDepth;

            float norVAF =  (float)norAltCount / (float)norDepth;
            float norCountVAF = NorBase.getVAF(chr,(*somaticVarIter).first);
            float norMpqVAF = NorBase.getFilterdMpqVAF(chr,(*somaticVarIter).first);

            if(norVAF != norCountVAF){
                std::cerr << "Error(diff nor VAF) => can't find GTtype at chr: " << chr << " pos: " << ((*somaticVarIter).first) + 1;
                exit(1); 
            }
            float norVAF_substract = (norMpqVAF - norCountVAF);
            float norMPQReadRatio = (float)(norDepth - norMpqDepth) / (float)norDepth;
            float TmpNorMPQReadRatio = NorBase.getLowMpqReadRatio(chr,(*somaticVarIter).first);

            if(norMPQReadRatio != TmpNorMPQReadRatio){
                std::cerr << "Error(diff low MPQ read ratio) => can't find GTtype at chr: " << chr << " pos: " << ((*somaticVarIter).first) + 1 << "\n";
                std::cerr << "norMPQReadRatio: " << norMPQReadRatio << " TmpNorMPQReadRatio: " << TmpNorMPQReadRatio;
                exit(1); 
            }

            // Calculating the difference in VAF between VAF and filtered low MPQ VAF
            float tumVAF_substract = (tumMpqVAF -tumVAF);
            float tumLowMpqReadRatio = (*somaticVarIter).second.base.lowMpqReadRatio;

            if(norMPQReadRatio < 0){
                std::cerr << "Error(norMPQReadRatio) => chr: " << chr << " pos: " << ((*somaticVarIter).first) + 1;
                std::cerr << "norMPQReadRatio " << norMPQReadRatio << " norDepth: " << norDepth << " norMPQDepth: " << norMpqDepth;
                exit(1);
            }

            if(tumLowMpqReadRatio < 0){
                std::cerr << "Error(tumMPQReadRatio) => chr: " << chr << " pos: " << ((*somaticVarIter).first) + 1;
                std::cerr << "tumMPQReadRatio " << tumLowMpqReadRatio << " norDepth: " << tumDepth << " norMPQDepth: " << tumMpqDepth;
                exit(1);
            }

                
            //current SNP GT type
            std::string GTtype = (*somaticVarIter).second.GTtype;

            if(GTtype == ""){
                std::cerr << "Error(GTtype) => chr: " << chr << " pos: " << ((*somaticVarIter).first) + 1 << " GTtype: " << GTtype;
                exit(1);                
            }

            //The count of deletions occurring at the SNP position
            int norDeletionCount = NorBase.getVarDeletionCount(chr, (*somaticVarIter).first);

            //int tumDeletionCount = TumBase.getVarDeletionCount(chr, (*somaticVarIter).first);
            int tumDeletionCount = (*SomaticChrPosInfo)[chr][(*somaticVarIter).first].base.delCount;

            //the ratio of deletions occurring
            float norDeletionRatio;
            float tumDeletionRatio = (*somaticVarIter).second.tumDelRatio;

            if(norDeletionCount == 0){
                norDeletionRatio = 0.0;
            }else{
                norDeletionRatio = (float)norDeletionCount / ((float)norDeletionCount + (float)norDepth);
            }

            
            //HP3 SNP position (1-base)
            int HP3pos = (*somaticVarIter).first + 1;

            //write information
            (*tagHP3Log) << std::fixed << std::setprecision(3) 
                        << chr << " \t"  //1 
                        << HP3pos << "\t"  //2
                        << readCount << "\t\t"  //3
                        << "(" << totalHP1WithHP3Read << "," << totalHP2WithHP3Read << "," << totalOnlyHP3Read << "," << totalMessyHPRead << ")\t" //4 
                        << unTagCount << "\t\t"  //5
                        << HP1withHP3ReadRatio <<"\t"  //6
                        << HP2WithHP3ReadRatio << "\t"  //7
                        << OnlyHP3ReadRatio << "\t"  //8
                        << messyHPReadRatio << "\t\t"  //9
                        << norVAF << "\t"  //10
                        << tumVAF << "\t\t"  //11
                        << subtractDepth << "\t"  //12 
                        << norDepth << "\t"  //13
                        << tumDepth << "\t"  //14
                        << norDeletionCount << "\t"  //15
                        << tumDeletionCount << "\t"  //16
                        << norDeletionRatio << "\t"  //17
                        << tumDeletionRatio << "\t"  //18
                        << norMpqVAF << "\t" //19
                        << tumMpqVAF << "\t\t" //20
                        << norVAF_substract << "\t" //21
                        << tumVAF_substract << "\t\t" //22
                        << norMPQReadRatio << "\t" //23
                        << tumLowMpqReadRatio << "\t" //24
                        << GTtype <<"\n";  //25
                        
            somaticVarIter++;          
        }
    }

    tagHP3Log->close();
    delete tagHP3Log;
    tagHP3Log = nullptr;
    std::cerr<< difftime(time(NULL), begin) << "s\n";  
}

void SomaticVarCaller::WriteOtherSomaticHpLog(const HaplotagParameters &params, const std::vector<std::string> &chrVec, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat){
    std::ofstream *OtherHpSomaticVarLog=NULL;
    std::string logPosfix = "_otherHpSomaticVar.log";
    OtherHpSomaticVarLog=new std::ofstream(params.resultPrefix + logPosfix);

    int totalOtherSomaticHpVar = 0;
    for(auto chr: chrVec){
        for(auto somaticVar: (*SomaticChrPosInfo)[chr]){
            if(somaticVar.second.somaticHp4Base != Nucleotide::UNKOWN){
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
        for(auto somaticVar: (*SomaticChrPosInfo)[chr]){
            int pos = somaticVar.first;

            if(somaticVar.second.somaticHp4Base == Nucleotide::UNKOWN){
                continue;
            }
            if(currentChrVar.find(pos) == currentChrVar.end()){
                std::cerr << "Error(WriteOtherSomaticHpLog) => can't find the position : chr:" << chr << " pos: " << pos + 1;
                exit(1);
            }

            (*OtherHpSomaticVarLog) << chr << "\t"
                                    << pos + 1 << "\t"
                                    << currentChrVar[pos].Variant[Genome::TUMOR].Ref << "\t"
                                    << currentChrVar[pos].Variant[Genome::TUMOR].Alt << "\t"
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

std::map<std::string, std::map<int, HP3_Info>> SomaticVarCaller::getSomaticChrPosInfo(){
    return (*SomaticChrPosInfo);
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
}

void VcfParser::variantParser(std::string &variantFile, VCF_Info &Info, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat){

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

void VcfParser::compressParser(std::string &variantFile, VCF_Info &Info, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat){
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

void VcfParser::unCompressParser(std::string &variantFile, VCF_Info &Info, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat){
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

void VcfParser::parserProcess(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat){
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
                RefAlt tmp;
                tmp.Ref = fields[3]; 
                tmp.Alt = fields[4];

                tmp.is_phased_hetero = true;
                tmp.is_unphased_hetero = false;
                tmp.is_homozygous= false;

                Info.chrVariant[chr][pos]=tmp;

                if(Info.gene_type == Genome::NORMAL){
                    mergedChrVarinat[chr][pos].Variant[Genome::NORMAL] = tmp;
                    mergedChrVarinat[chr][pos].isExistNormal = true;    
                }else if(Info.gene_type == Genome::TUMOR){
                    mergedChrVarinat[chr][pos].Variant[Genome::TUMOR] = tmp;
                    mergedChrVarinat[chr][pos].isExistTumor = true; 
                }
                
                if(integerPS){
                    Info.chrVariantPS[chr][pos]=std::stoi(psValue);
                }
                else{
                    std::map<std::string, int>::iterator psIter = Info.psIndex.find(psValue);
                    
                    if( psIter == Info.psIndex.end() ){
                        Info.psIndex[psValue] = Info.psIndex.size();
                    }
                    Info.chrVariantPS[chr][pos]=Info.psIndex[psValue];
                }
                
                // record haplotype allele
                if( fields[9][modifu_start] == '0' && fields[9][modifu_start+2] == '1' ){
                    Info.chrVariantHP1[chr][pos]=fields[3];
                    Info.chrVariantHP2[chr][pos]=fields[4];
                }
                else if( fields[9][modifu_start] == '1' && fields[9][modifu_start+2] == '0' ){
                    Info.chrVariantHP1[chr][pos]=fields[4];
                    Info.chrVariantHP2[chr][pos]=fields[3];
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
                    RefAlt tmp;
                    tmp.Ref = fields[3]; 
                    tmp.Alt = fields[4];
                    tmp.is_phased_hetero = false;
                    tmp.is_unphased_hetero = false;
                    tmp.is_homozygous = true;

                    Info.chrVariant[chr][pos] = tmp;

                    if(Info.gene_type == Genome::NORMAL){
                        mergedChrVarinat[chr][pos].Variant[Genome::NORMAL] = tmp;
                        mergedChrVarinat[chr][pos].isExistNormal = true;    
                    }else if(Info.gene_type == Genome::TUMOR){
                        mergedChrVarinat[chr][pos].Variant[Genome::TUMOR] = tmp;
                        mergedChrVarinat[chr][pos].isExistTumor = true; 
                    }
                }
            //unphased heterozygous
            }else if( fields[9][modifu_start] == '0' && fields[9][modifu_start+1] == '/' && fields[9][modifu_start+2] == '1' ){
                if(parseSnpFile){
                    RefAlt tmp;
                    tmp.Ref = fields[3]; 
                    tmp.Alt = fields[4];

                    tmp.is_phased_hetero = false;
                    tmp.is_unphased_hetero = true;
                    tmp.is_homozygous = false;

                    Info.chrVariant[chr][pos] = tmp;

                    if(Info.gene_type == Genome::NORMAL){
                        mergedChrVarinat[chr][pos].Variant[Genome::NORMAL] = tmp;
                        mergedChrVarinat[chr][pos].isExistNormal = true;    
                    }else if(Info.gene_type == Genome::TUMOR){
                        mergedChrVarinat[chr][pos].Variant[Genome::TUMOR] = tmp;
                        mergedChrVarinat[chr][pos].isExistTumor = true; 
                    }
                }
            }
        }
    }
}

highConBenchmark::highConBenchmark(){
    setParseSnpFile(true);
    openTestingFunc = false;
}
highConBenchmark::~highConBenchmark(){

}

void highConBenchmark::setTestingFunc(bool openTestingFunc){
    this->openTestingFunc = openTestingFunc;
}

void highConBenchmark::loadHighConSomatic(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat){
    if(!openTestingFunc) return;
    variantParser(input, Info, mergedChrVarinat);
}

void highConBenchmark::parserProcess(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat){
    if(!openTestingFunc) return;
    
    if( input.substr(0, 2) == "##" && getParseSnpFile()){
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

        RefAlt tmp;
        tmp.Ref = fields[3]; 
        tmp.Alt = fields[4];
        mergedChrVarinat[chr][pos].Variant[Genome::SEQC_HIGH_CON] = tmp;
        mergedChrVarinat[chr][pos].isExistSeqcHighCon = true;
    }
}

void highConBenchmark::recordDelReadCount(const std::string &chr, std::map<int, RefAltSet>::iterator &currentVariantIter){
    if(!openTestingFunc) return;
    
    if(currentVariantIter->second.isExistSeqcHighCon){
        int pos = (*currentVariantIter).first;
        posAltRefDelCount[chr][pos].delCount++;

        //record somatic position for record crossing high con snp read
        highConSomaticPos.push_back(std::make_pair(pos, SnpHP::NONE_SNP));
    }        
}

void highConBenchmark::recordRefAltAlleleCount(const std::string &chr, std::string &base, std::map<int, RefAltSet>::iterator &currentVariantIter){
    if(!openTestingFunc) return;

    if(currentVariantIter->second.isExistSeqcHighCon){
        int pos = currentVariantIter->first;
        std::string refAllele = currentVariantIter->second.Variant[Genome::SEQC_HIGH_CON].Ref;
        std::string altAllele = currentVariantIter->second.Variant[Genome::SEQC_HIGH_CON].Alt;

        int baseHP = SnpHP::NONE_SNP;

        if(base == refAllele){
            posAltRefDelCount[chr][pos].refCount++;
        }else if(base == altAllele){
            posAltRefDelCount[chr][pos].altCount++;
            baseHP = SnpHP::SOMATIC_H3;
        }
        //record somatic position for record crossing high con snp read
        highConSomaticPos.push_back(std::make_pair(pos, baseHP));
    }
}

void highConBenchmark::recordCrossingHighConSnpRead(const std::string &chr, std::string &readID, std::string &hpResult, std::map<int, int> &variantsHP, std::map<int, RefAltSet> &currentChrVariants){
    // if not open testing function, return
    if(!openTestingFunc) return;

    somaticReadLog tmp;
    tmp.chr = chr;
    tmp.readID = readID;
    tmp.hpResult = hpResult;

    bool isCrossHighConSomatic = false;
    bool existHighConVaraints = false;

    for(auto varIter : highConSomaticPos){
        int pos = varIter.first;
        int baseHP = varIter.second;

        //transfer hpResult to H3
        if(baseHP == SnpHP::SOMATIC_H3){
            existHighConVaraints = true;
        }

        tmp.somaticSnpHp[pos] = baseHP;

        isCrossHighConSomatic = true;
    }

    if(isCrossHighConSomatic){
        //exist high con variants alt allele
        if(existHighConVaraints){
            //correction hpResult that exist high con variants
            if(hpResult == "1"){
                tmp.hpResult = "1-1";
            }else if(hpResult == "2"){
                tmp.hpResult = "2-1";
            }else if(hpResult == "."){
                tmp.hpResult = "3";
            }
        }else{
            //correction hpResult that not exist high con variants
            if(hpResult == "2-1"){
                tmp.hpResult = "2";
            }else if(hpResult == "1-1"){
                tmp.hpResult = "1";
            }else if(hpResult == "3"){
                tmp.hpResult = ".";
            }
        }
    }

    if(isCrossHighConSomatic){
        readsCrossingHighConSnpVec.push_back(tmp);
    }
    
    //clear high con somatic position in current read for next read
    if(!highConSomaticPos.empty()){
        highConSomaticPos.clear();
    }
}

void highConBenchmark::recordTaggedSomaticRead(const std::string &chr, std::string &readID, std::string &hpResult, std::map<int, int> &variantsHP, std::map<int, RefAltSet> &currentChrVariants){
    // if not open testing function, return
    if(!openTestingFunc) return;

    somaticReadLog tmp;
    tmp.chr = chr;
    tmp.readID = readID;
    tmp.hpResult = hpResult;

    bool readExistHighConSomatic = false;

    auto varIter = variantsHP.begin();
    while(varIter != variantsHP.end()){
        int pos = varIter->first;
        int snpHP = varIter->second;
        if(currentChrVariants.find(pos) != currentChrVariants.end()){
            if(currentChrVariants[pos].isExistSeqcHighCon && (snpHP == SnpHP::SOMATIC_H3 || snpHP == SnpHP::SOMATIC_H4)){
                tmp.somaticSnpHp[pos] = snpHP;
                readExistHighConSomatic = true;
            }
        }
        varIter++;
    }

    if(readExistHighConSomatic){
        taggedSomaticReadVec.push_back(tmp);
    }
}

void highConBenchmark::writePosAlleleCountLog(std::vector<std::string> &chrVec, HaplotagParameters &params, std::string logPosfix, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat){
    // if not open testing function, return
    if(!openTestingFunc) return;

    std::ofstream *refAltCountLog=NULL;
    refAltCountLog=new std::ofstream(params.resultPrefix + logPosfix);
    int totalVariantCount = 0;

    for(auto &chr : chrVec){
        totalVariantCount += posAltRefDelCount[chr].size();
    }

    if(!refAltCountLog->is_open()){
        std::cerr<< "Fail to open write file: " << params.resultPrefix + logPosfix << "\n";
        exit(1);
    }else{
        (*refAltCountLog) << "#############################\n";
        (*refAltCountLog) << "# Somatic SNP allele count #\n";
        (*refAltCountLog) << "#############################\n";
        (*refAltCountLog) << "##High confidence VCF:"  << params.seqcHighCon << "\n";
        (*refAltCountLog) << "##MappingQualityThreshold:"  << params.qualityThreshold << "\n";
        (*refAltCountLog) << "##Tatal variants:"  << totalVariantCount << "\n";
        (*refAltCountLog) << "#CHROM\t"
                          << "POS\t"
                          << "REF\t"
                          << "ALT\t"
                          << "REF_COUNT\t"
                          << "ALT_COUNT\t"
                          << "DEL_COUNT\n";
    }

    for(auto &chr : chrVec){
        for(auto &posIter : posAltRefDelCount[chr]){
            (*refAltCountLog) << chr << "\t"
                              << posIter.first << "\t"
                              << mergedChrVarinat[chr][posIter.first].Variant[Genome::SEQC_HIGH_CON].Ref << "\t"
                              << mergedChrVarinat[chr][posIter.first].Variant[Genome::SEQC_HIGH_CON].Alt << "\t"
                              << posIter.second.refCount << "\t"
                              << posIter.second.altCount << "\t"
                              << posIter.second.delCount << "\n";
        }
    }

    refAltCountLog->close();
    delete refAltCountLog;
    refAltCountLog = nullptr;
}

void highConBenchmark::writeTaggedSomaticReadLog(HaplotagParameters &params, std::string logPosfix){
    // if not open testing function, return
    if(!openTestingFunc) return;
    writeReadLog(params, logPosfix, taggedSomaticReadVec);
}

void highConBenchmark::writeCrossHighConSnpReadLog(HaplotagParameters &params, std::string logPosfix){
    // if not open testing function, return
    if(!openTestingFunc) return;
    writeReadLog(params, logPosfix, readsCrossingHighConSnpVec);
}

void highConBenchmark::writeReadLog(HaplotagParameters &params, std::string logPosfix, std::vector<somaticReadLog> &somaticReadVec){
    // if not open testing function, return
    if(!openTestingFunc) return;

    std::ofstream *somaticReadLog=NULL;
    somaticReadLog=new std::ofstream(params.resultPrefix + logPosfix);

    int totalTruthSomaticReads = 0;
    for(auto readIter :readsCrossingHighConSnpVec){
        if(readIter.hpResult == "1-1" || readIter.hpResult == "2-1" || readIter.hpResult == "3"){
            totalTruthSomaticReads++;
        }
    }

    int totalTaggedSomaticReads = 0;
    for(auto readIter :taggedSomaticReadVec){
        if(readIter.hpResult == "1-1" || readIter.hpResult == "2-1" || readIter.hpResult == "3"){
            totalTaggedSomaticReads++;
        }
    }

    if(!somaticReadLog->is_open()){
        std::cerr<< "Fail to open write file: " << params.resultPrefix + logPosfix << "\n";
        exit(1);
    }else{
        (*somaticReadLog) << "#####################\n";
        (*somaticReadLog) << "# Somatic Reads Log #\n";
        (*somaticReadLog) << "#####################\n";
        (*somaticReadLog) << "##High confidence VCF: "  << params.seqcHighCon << "\n";
        (*somaticReadLog) << "##MappingQualityThreshold: "  << params.qualityThreshold << "\n";
        (*somaticReadLog) << "##Tatal tagged somatic reads: "  << totalTaggedSomaticReads << "\n";
        (*somaticReadLog) << "##Tatal truth somatic reads: "  << totalTruthSomaticReads << "\n";
        (*somaticReadLog) << "##Truth somatic read ratio: "  << (float)totalTaggedSomaticReads / (float)totalTruthSomaticReads << "\n";
        (*somaticReadLog) << "#CHROM\t"
                          << "ReadID\t"
                          << "Haplotype\t"
                          << "somaticVariant,HP\n";
    }

    for(auto somaticRead: somaticReadVec){

        (*somaticReadLog) << somaticRead.chr << "\t"
                          << somaticRead.readID << "\t"
                          << "H" << somaticRead.hpResult << "\t";

        for(auto SnpHp: somaticRead.somaticSnpHp){
            (*somaticReadLog) << SnpHp.first+1 << "," << SnpHp.second << "\t";
        }
        (*somaticReadLog) << "\n";
    }
    (*somaticReadLog).close();
    delete somaticReadLog;
    somaticReadLog = nullptr;
}

void highConBenchmark::displaySomaticVarCount(std::vector<std::string> &chrVec, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat){
    if(!openTestingFunc) return;
    int totalVariantCount = 0;
    
    for(auto &chr : chrVec){
        std::map<int, RefAltSet>::iterator chrVariantIter = mergedChrVarinat[chr].begin();
        while(chrVariantIter != mergedChrVarinat[chr].end()){
            if(chrVariantIter->second.isExistSeqcHighCon){
                totalVariantCount++;
            }
            chrVariantIter++;
        }
    }
    std::cout << "Total somatic variants: " << totalVariantCount << "\n";
}


void HaplotagProcess::tagRead(HaplotagParameters &params, const int geneType){

    // input file management
    std::string openBamFile = params.bamFile;

    if(geneType == Genome::TUMOR){
        openBamFile = params.tumorBamFile;
    }else if(geneType == Genome::NORMAL){
        openBamFile = params.bamFile;
    }

    // open bam file
    samFile *in = hts_open(openBamFile.c_str(), "r");
    // load reference file
    hts_set_fai_filename(in, params.fastaFile.c_str() );
    // input reader
    bam_hdr_t *bamHdr = sam_hdr_read(in);
    // header add pg tag
    sam_hdr_add_pg(bamHdr, "longphase", "VN", params.version.c_str(), "CL", params.command.c_str(), NULL);
    // bam file index
    hts_idx_t *idx = NULL;
    // check input bam file
    if (in == NULL) {
        std::cerr<<"ERROR: Cannot open bam file " << openBamFile.c_str() << "\n";
    }
    // check bam file index
    if ((idx = sam_index_load(in, openBamFile.c_str())) == 0) {
        std::cerr<<"ERROR: Cannot open index for bam file\n";
        exit(1);
    }

    // output file mangement
    std::string writeBamFile = params.resultPrefix + "." + params.outputFormat;
    // open output bam file
    samFile *out = hts_open(writeBamFile.c_str(), (params.outputFormat == "bam" ? "wb" : "wc" ));
    // load reference file
    hts_set_fai_filename(out, params.fastaFile.c_str() );
    // output writer
    int result = sam_hdr_write(out, bamHdr);

    // record reference last variant pos
    std::vector<int> last_pos;
    for( auto chr : *chrVec ){
        auto lastVariantIter = vcfSet[geneType].chrVariantPS[chr].rbegin();
        if( lastVariantIter != vcfSet[geneType].chrVariantPS[chr].rend() ){
            last_pos.push_back(lastVariantIter->second);
        }
        else{
            last_pos.push_back(0);
        }
    }

    // reference fasta parser
    FastaParser fastaParser(params.fastaFile, *chrVec, last_pos, params.numThreads);

    // write tag read detail information
    std::ofstream *tagResult=NULL;
    if(params.writeReadLog){
        tagResult=new std::ofstream(params.resultPrefix+".out");
        if(!tagResult->is_open()){
            std::cerr<< "Fail to open write file: " << params.resultPrefix+".out" << "\n";
            exit(1);
        }
        else{
            (*tagResult) << "##snpFile:"                 << params.snpFile                    << "\n";
            if(tagTumorMode)
            (*tagResult) << "##TumorSnpFile:"            << params.tumorSnpFile               << "\n"; //new
            (*tagResult) << "##svFile:"                  << params.svFile                     << "\n";
            (*tagResult) << "##bamFile:"                 << params.bamFile                    << "\n";
            if(tagTumorMode)
            (*tagResult) << "##tumorBamFile:"            << params.tumorBamFile               << "\n"; //new
            (*tagResult) << "##resultPrefix:"            << params.resultPrefix               << "\n";
            (*tagResult) << "##numThreads:"              << params.numThreads                 << "\n";
            (*tagResult) << "##region:"                  << params.region                     << "\n";
            if(tagTumorMode)
            (*tagResult) << "##tagTumor:"                << params.tagTumorSnp                << "\n";  //new
            if(tagTumorMode)
            (*tagResult) << "##somaticCallingThreshold:" << params.somaticCallingMpqThreshold << "\n";  //new
            (*tagResult) << "##qualityThreshold:"        << params.qualityThreshold           << "\n";
            (*tagResult) << "##percentageThreshold:"     << params.percentageThreshold        << "\n";
            (*tagResult) << "##tagSupplementary:"        << params.tagSupplementary           << "\n";
            (*tagResult) << "#Read\t"
                         << "Chr\t"
                         << "ReadStart\t"
                         << "Confidnet(%)\t"
                         << "Haplotype\t"
                         << "PhaseSet\t"
                         << "TotalAllele\t"
                         << "HP1Allele\t"
                         << "HP2Allele\t";
            if(tagTumorMode){
                (*tagResult) << "HP3Allele\t"
                             << "HP4Allele\t";
            }
            (*tagResult) << "phasingQuality(PQ)\t"
                         << "(Variant,HP)\t"
                         << "(PhaseSet,Variantcount)\n";

        }
    }
    // init data structure and get core n
    htsThreadPool threadPool = {NULL, 0};
    // creat thread pool
    if (!(threadPool.pool = hts_tpool_init(params.numThreads))) {
        fprintf(stderr, "Error creating thread pool\n");
    }
    // set thread
    hts_set_opt(in, HTS_OPT_THREAD_POOL, &threadPool);
    hts_set_opt(out, HTS_OPT_THREAD_POOL, &threadPool);
    // initialize an alignment
    bam1_t *aln = bam_init1();

    // loop all chromosome
    for(auto chr : *chrVec ){
        std::time_t begin = time(NULL);
        std::cerr<<"chr: " << chr << " ... " ;
        // records all variants within this chromosome.
        currentChrVariants = (*mergedChrVarinat)[chr];
        // since each read is sorted based on the start coordinates, to save time, 
        // firstVariantIter keeps track of the first variant that each read needs to check.
        firstVariantIter = currentChrVariants.begin();
        // get the coordinates of the last variant
        // the tagging process will not be perform if the read's start coordinate are over than last variant.
        std::map<int, RefAltSet>::reverse_iterator last = currentChrVariants.rbegin();
        
        // fetch chromosome string
        std::string chr_reference = fastaParser.chrString.at(chr);
        
        // tagging will be attempted for reads within the specified coordinate region.
        std::string region = !params.region.empty() ? params.region : chr + ":1-" + std::to_string((*chrLength)[chr]);
        hts_itr_t *iter = sam_itr_querys(idx, bamHdr, region.c_str());
        // iter all reads
        while ((result = sam_itr_multi_next(in, iter, aln)) >= 0) {
            totalAlignment++;
            int flag = aln->core.flag;

            if ( aln->core.qual < params.qualityThreshold ){
                // mapping quality is lower than threshold
                totalUnTagCount++;
                totalLowerQuality++;
            }
            else if( (flag & 0x4) != 0 ){
                // read unmapped
                totalUnmapped++;
                totalUnTagCount++;
            }
            else if( (flag & 0x100) != 0 ){
                // secondary alignment. repeat.
                // A secondary alignment occurs when a given read could align reasonably well to more than one place.
                totalSecondary++;
                totalUnTagCount++;
            }
            else if( (flag & 0x800) != 0 && params.tagSupplementary == false ){
                // supplementary alignment
                // A chimeric alignment is represented as a set of linear alignments that do not have large overlaps.
                totalSupplementary++;
                totalUnTagCount++;
            }
            else if(last == currentChrVariants.rend()){
                // skip 
                totalUnTagCount++;
                totalEmptyVariant++;
            }
            else if(int(aln->core.pos) <= (*last).first){
                
                if( (flag & 0x800) != 0 ){
                    totalSupplementary++;
                }

                int pqValue = 0;
                int psValue = 0; 
                int haplotype = ReadHP::unTag;

                if(tagTumorMode){
                    haplotype = SomaticJudgeHaplotype(*bamHdr, *aln, chr, params.percentageThreshold, tagResult, pqValue, psValue, geneType, chr_reference);
                }else{
                    haplotype = judgeHaplotype(*bamHdr, *aln, chr, params.percentageThreshold, tagResult, pqValue, psValue, geneType, chr_reference);
                }

                initFlag(aln, "HP");
                initFlag(aln, "PS");
                initFlag(aln, "PQ");

                if (haplotype != ReadHP::unTag){

                    totalHpCount[haplotype]++;
                    totalTagCount++;

                    if(tagTumorMode){
                        std::string haplotype_str = "";
                        haplotype_str = convertHpResultToString(haplotype);

                        bam_aux_append(aln, "HP", 'Z', (haplotype_str.size() + 1), (const uint8_t*) haplotype_str.c_str());
                        if(psValue != -1) bam_aux_append(aln, "PS", 'i', sizeof(psValue), (uint8_t*) &psValue);
                        bam_aux_append(aln, "PQ", 'i', sizeof(pqValue), (uint8_t*) &pqValue);
                    }else{
                        bam_aux_append(aln, "HP", 'i', sizeof(haplotype), (uint8_t*) &haplotype);
                        bam_aux_append(aln, "PS", 'i', sizeof(psValue), (uint8_t*) &psValue);
                        bam_aux_append(aln, "PQ", 'i', sizeof(pqValue), (uint8_t*) &pqValue);
                    }
                }
                else{
                    totalunTag_HP0++;
                    totalUnTagCount++;
                }
            }
            else{
                totalOtherCase++;
                totalUnTagCount++;
            }

            // write this alignment to result bam file
            result = sam_write1(out, bamHdr, aln);
        }

        if (tagTumorMode && params.writeReadLog){
            hpBeforeInheritance->mergeLocalReadHp(chr, (*beforeCorrReadHpResult)[chr]);
            hpAfterInheritance->mergeLocalReadHp(chr, (*afterCorrReadHpResult)[chr]);
        }
        hts_itr_destroy(iter);  
        std::cerr<< difftime(time(NULL), begin) << "s\n";
    }
    if(tagResult!=NULL){
        (*tagResult).close();
        delete tagResult;
        tagResult = nullptr;
    }
    if(tagTumorMode && params.writeReadLog){
        hpBeforeInheritance->writeReadHpDistriLog(params, "_readDistri_beforeInheritance.out", *chrVec);
        hpAfterInheritance->writeReadHpDistriLog(params, "_readDistri_afterInheritance.out", *chrVec);
        //write snp cover region
        hpAfterInheritance->writePosCoverRegionLog(params, "_SnpCoverRegion.out", *chrVec);
        //write read cover region in whole genome
        hpAfterInheritance->writeTagReadCoverRegionLog(params, "_readCoverRegion.bed", *chrVec, *chrLength);
        //write somatic read log
        highConSomaticData.writeTaggedSomaticReadLog(params, "_somaticRead.out");
        highConSomaticData.writeCrossHighConSnpReadLog(params, "_crossHighConSnpRead.out");
        highConSomaticData.writePosAlleleCountLog(*chrVec, params, "_alleleCount.out", *mergedChrVarinat);
    }

    hts_idx_destroy(idx);
    bam_hdr_destroy(bamHdr);
    bam_destroy1(aln);
    sam_close(in);
    sam_close(out);
    hts_tpool_destroy(threadPool.pool);

    return;
}

void HaplotagProcess::initFlag(bam1_t *aln, std::string flag){

    uint8_t *hpTag = bam_aux_get(aln, flag.c_str() );

    if( hpTag != NULL )
        bam_aux_del(aln, hpTag);

    return;
}

int HaplotagProcess::judgeHaplotype(const bam_hdr_t &bamHdr,const bam1_t &aln, std::string chrName, double percentageThreshold, std::ofstream *tagResult, int &pqValue, int &psValue, const int tagGeneType, std::string &ref_string){

    int hp1Count = 0;
    int hp2Count = 0;
    //record variants on this read
    std::map<int, int> variantsHP;

    std::map<int,int> countPS;

    // Skip variants that are to the left of this read
    while( firstVariantIter != currentChrVariants.end() && (*firstVariantIter).first < aln.core.pos ){
        firstVariantIter++;
    }
    
    if( firstVariantIter == currentChrVariants.end() ){
        return 0;
    }

    // position relative to reference
    int ref_pos = aln.core.pos;
    // position relative to read
    int query_pos = 0;
    // set variant start for current alignment
    std::map<int, RefAltSet>::iterator currentVariantIter = firstVariantIter;

    // reading cigar to detect snp on this read
    int aln_core_n_cigar = int(aln.core.n_cigar);
    for(int i = 0; i < aln_core_n_cigar ; i++ ){
        uint32_t *cigar = bam_get_cigar(&aln);
        int cigar_op = bam_cigar_op(cigar[i]);
        int length   = bam_cigar_oplen(cigar[i]);

        // iterator next variant
        while( currentVariantIter != currentChrVariants.end() && (*currentVariantIter).first < ref_pos ){
            currentVariantIter++;
        }

        // CIGAR operators: MIDNSHP=X correspond 012345678
        // 0: alignment match (can be a sequence match or mismatch)
        // 7: sequence match
        // 8: sequence mismatch
        if( cigar_op == 0 || cigar_op == 7 || cigar_op == 8 ){

            while( currentVariantIter != currentChrVariants.end() && (*currentVariantIter).first < ref_pos + length){

                int curPos = (*currentVariantIter).first;
                auto norVar = (*currentVariantIter).second.Variant[Genome::NORMAL];
                int refAlleleLen = norVar.Ref.length();
                int altAlleleLen = norVar.Alt.length();

                int offset = (*currentVariantIter).first - ref_pos;

                if( offset < 0){
                }
                else{
                    uint8_t *q = bam_get_seq(&aln);
                    char base_chr = seq_nt16_str[bam_seqi(q,query_pos + offset)];
                    std::string base(1, base_chr);
                    
                    // currentVariant is SNP
                    if( refAlleleLen == 1 && altAlleleLen == 1 ){
                        // Detected that the base of the read is either REF or ALT. 
                        if( (base == norVar.Ref) || (base == norVar.Alt) ){

                            std::map<int, int>::iterator posPSiter = vcfSet[Genome::NORMAL].chrVariantPS[chrName].find((*currentVariantIter).first);

                            if( posPSiter == vcfSet[Genome::NORMAL].chrVariantPS[chrName].end() ){
                                std::cerr<< curPos << "\t"
                                         << norVar.Ref << "\t"
                                         << norVar.Alt << "\n";
                                exit(EXIT_SUCCESS);
                            }
                            else{
                                if( base == vcfSet[Genome::NORMAL].chrVariantHP1[chrName][curPos]){
                                    hp1Count++;
                                    variantsHP[curPos]=0;
                                }
                                if( base == vcfSet[Genome::NORMAL].chrVariantHP2[chrName][curPos]){
                                    hp2Count++;
                                    variantsHP[curPos]=1;
                                }
                                countPS[vcfSet[Genome::NORMAL].chrVariantPS[chrName][curPos]]++;
                            }
                            
                        }
                    }
                    // currentVariant is insertion
                    else if( refAlleleLen == 1 && altAlleleLen != 1 && i+1 < aln_core_n_cigar){
                        
                        int hp1Length = vcfSet[Genome::NORMAL].chrVariantHP1[chrName][curPos].length();
                        int hp2Length = vcfSet[Genome::NORMAL].chrVariantHP2[chrName][curPos].length();
                        
                        if ( ref_pos + length - 1 == (*currentVariantIter).first && bam_cigar_op(cigar[i+1]) == 1 ) {
                            // hp1 occur insertion
                            if( hp1Length != 1 && hp2Length == 1 ){
                                hp1Count++;
                                variantsHP[curPos]=0;
                            }
                            // hp2 occur insertion
                            else if( hp1Length == 1 && hp2Length != 1 ){
                                hp2Count++;
                                variantsHP[curPos]=1;
                            }
                        }
                        else {
                            // hp1 occur insertion
                            if( hp1Length != 1 && hp2Length == 1 ){
                                hp2Count++;
                                variantsHP[curPos]=1;
                            }
                            // hp2 occur insertion
                            else if( hp1Length == 1 && hp2Length != 1 ){
                                hp1Count++;
                                variantsHP[curPos]=0;
                            }
                        }
                        countPS[vcfSet[Genome::NORMAL].chrVariantPS[chrName][curPos]]++;
                    } 
                    // currentVariant is deletion
                    else if( refAlleleLen != 1 && altAlleleLen == 1 && i+1 < aln_core_n_cigar) {
                        
                        int hp1Length = vcfSet[Genome::NORMAL].chrVariantHP1[chrName][curPos].length();
                        int hp2Length = vcfSet[Genome::NORMAL].chrVariantHP2[chrName][curPos].length();
                        
                        if ( ref_pos + length - 1 == (*currentVariantIter).first && bam_cigar_op(cigar[i+1]) == 2 ) {
                            // hp1 occur deletion
                            if( hp1Length != 1 && hp2Length == 1 ){
                                hp1Count++;
                                variantsHP[curPos]=0;
                            }
                            // hp2 occur deletion
                            else if( hp1Length == 1 && hp2Length != 1 ){
                                hp2Count++;
                                variantsHP[curPos]=1;
                            }
                        }
                        else {
                            // hp2 occur deletion
                            if( hp1Length != 1 && hp2Length == 1 ){
                                hp2Count++;
                                variantsHP[curPos]=1;
                            }
                            // hp1 occur deletion
                            else if( hp1Length == 1 && hp2Length != 1 ){
                                hp1Count++;
                                variantsHP[curPos]=0;
                            }
                        }
                        countPS[vcfSet[Genome::NORMAL].chrVariantPS[chrName][curPos]]++;
                    } 

                }
                currentVariantIter++;
            }
            query_pos += length;
            ref_pos += length;
        }
            // 1: insertion to the reference
        else if( cigar_op == 1 ){
            query_pos += length;
        }
            // 2: deletion from the reference
        else if( cigar_op == 2 ){
            
            if(ref_string != ""){
                int del_len = length;
                if ( ref_pos + del_len + 1 == (*currentVariantIter).first ){
                    //if( homopolymerLength((*currentVariantIter).first , ref_string) >=3 ){
                        // special case
                    //}
                }
                else if( (*currentVariantIter).first >= ref_pos  && (*currentVariantIter).first < ref_pos + del_len ){
                    // check variant in homopolymer
                    if( homopolymerLength((*currentVariantIter).first , ref_string) >=3 ){
                        
                        int curPos = (*currentVariantIter).first;
                        auto norVar = (*currentVariantIter).second.Variant[Genome::NORMAL];
                        int refAlleleLen = norVar.Ref.length();
                        int altAlleleLen = norVar.Alt.length();
                        
                        // SNP
                        if( refAlleleLen == 1 && altAlleleLen == 1){
                            // get the next match
                            char base_chr = seq_nt16_str[bam_seqi(bam_get_seq(&aln), query_pos)];
                            std::string base(1, base_chr);

                            if( base == vcfSet[Genome::NORMAL].chrVariantHP1[chrName][curPos]){
                                hp1Count++;
                                variantsHP[curPos]=0;
                            }
                            if( base == vcfSet[Genome::NORMAL].chrVariantHP2[chrName][curPos]){
                                hp2Count++;
                                variantsHP[curPos]=1;
                            }
                            countPS[vcfSet[Genome::NORMAL].chrVariantPS[chrName][curPos]]++;
                        }
                        
                        // the read deletion contain VCF's deletion
                        else if( refAlleleLen != 1 && altAlleleLen == 1){

                            int hp1Length = vcfSet[Genome::NORMAL].chrVariantHP1[chrName][curPos].length();
                            int hp2Length = vcfSet[Genome::NORMAL].chrVariantHP2[chrName][curPos].length();
                            // hp1 occur deletion
                            if( hp1Length != 1 && hp2Length == 1 ){
                                hp1Count++;
                                variantsHP[curPos]=0;
                            }
                            // hp2 occur deletion
                            else if( hp1Length == 1 && hp2Length != 1 ){
                                hp2Count++;
                                variantsHP[curPos]=1;
                            }
                            countPS[vcfSet[Genome::NORMAL].chrVariantPS[chrName][curPos]]++;
                        }
                    }
                }
            }
            ref_pos += length;
        }
            // 3: skipped region from the reference
        else if( cigar_op == 3 ){
            ref_pos += length;
        }
            // 4: soft clipping (clipped sequences present in SEQ)
        else if( cigar_op == 4 ){
            query_pos += length;
        }
            // 5: hard clipping (clipped sequences NOT present in SEQ)
            // 6: padding (silent deletion from padded reference)
        else if( cigar_op == 5 || cigar_op == 6 ){
            // do nothing
        }
        else{
            std::cerr<< "alignment find unsupported CIGAR operation from read: " << bam_get_qname(&aln) << "\n";
            exit(1);
        }
    }

    double min,max;

    auto readIter = vcfSet[tagGeneType].readSVHapCount.find(bam_get_qname(&aln));
    if( readIter != vcfSet[tagGeneType].readSVHapCount.end() ){
        hp1Count += vcfSet[tagGeneType].readSVHapCount[bam_get_qname(&aln)][0];
        hp2Count += vcfSet[tagGeneType].readSVHapCount[bam_get_qname(&aln)][1];
    }

    if(hp1Count > hp2Count){
        min = hp2Count;
        max = hp1Count;
    }
    else{
        min = hp1Count;
        max = hp2Count;
    }

    int hpResult = ReadHP::unTag;
    if( max/(max+min) < percentageThreshold){
        // no tag
        pqValue = 0;
        totalHighSimilarity++;
    }
    else{
        if(hp1Count > hp2Count){
            hpResult = ReadHP::H1;
        }
        if(hp1Count < hp2Count){
            hpResult = ReadHP::H2;
        }
    }

    if( max == 0 ){
        pqValue=0;
        totalWithOutVaraint++;
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

    if(tagResult!=NULL){
        //write tag log file
        std::string hpResultStr = ((hpResult == ReadHP::unTag )? "." : std::to_string(hpResult) );
        std::string psResultStr = ".";

        if( hpResultStr != "." ){
            auto psIter = countPS.begin();
            psResultStr = std::to_string((*psIter).first);
        }

        (*tagResult)<< bam_get_qname(&aln)              << "\t"
                    << bamHdr.target_name[aln.core.tid] << "\t"
                    << aln.core.pos                     << "\t"
                    << max/(max+min)                    << "\t"
                    << hpResultStr                      << "\t"
                    << psResultStr                      << "\t"
                    << hp1Count+hp2Count                << "\t"
                    << hp1Count                         << "\t"
                    << hp2Count                         << "\t"
                    << pqValue                          << "\t";


        // print position and HP
        for(auto v : variantsHP ){
            (*tagResult)<< " " << v.first << "," << v.second ;
        }

        (*tagResult) << "\t";

        // belong PS, number of variant
        for(auto v : countPS ){
            (*tagResult)<< " " << v.first << "," << v.second ;
        }

        (*tagResult)<< "\n";
    }

    return hpResult;
}

int HaplotagProcess::SomaticJudgeHaplotype(const bam_hdr_t &bamHdr,const bam1_t &aln, const std::string &chrName, double percentageThreshold, std::ofstream *tagResult, int &pqValue, int &psValue, const int tagGeneType, std::string &ref_string){

    std::map<int, int> hpCount;
    hpCount[1] = 0;
    hpCount[2] = 0;
    hpCount[3] = 0;
    hpCount[4] = 0;

    //record variants on this read
    std::map<int, int> variantsHP;

    // poition, deriveByHp
    std::map<int, int> somaticVarDeriveHP;

    //record PS count(vcf type, PS value, count)
    std::map<int, int> tumCountPS;
    std::map<int, int> norCountPS;


    // Skip variants that are to the left of this read
    while( firstVariantIter != currentChrVariants.end() && (*(firstVariantIter)).first < aln.core.pos ){
        firstVariantIter++;
    }     

    if( firstVariantIter == currentChrVariants.end() )
        return 0;

    // position relative to reference
    int ref_pos = aln.core.pos;
    // position relative to read
    int query_pos = 0;
    // set variant start for current alignment
    std::map<int, RefAltSet>::iterator currentVariantIter = firstVariantIter;

    // reading cigar to detect snp on this read
    int aln_core_n_cigar = int(aln.core.n_cigar);
    for(int i = 0; i < aln_core_n_cigar ; i++ ){
        uint32_t *cigar = bam_get_cigar(&aln);
        int cigar_op = bam_cigar_op(cigar[i]);
        int length   = bam_cigar_oplen(cigar[i]);

        // iterator next variant
        while( currentVariantIter != currentChrVariants.end() && (*currentVariantIter).first < ref_pos ){
            currentVariantIter++;
        }

        // CIGAR operators: MIDNSHP=X correspond 012345678
        // 0: alignment match (can be a sequence match or mismatch)
        // 7: sequence match
        // 8: sequence mismatch
        if( cigar_op == 0 || cigar_op == 7 || cigar_op == 8 ){

            while( currentVariantIter != currentChrVariants.end() && (*currentVariantIter).first < ref_pos + length){
                
                int offset = (*currentVariantIter).first - ref_pos;

                if( offset < 0){
                }
                else{
                    uint8_t *q = bam_get_seq(&aln);
                    char base_chr = seq_nt16_str[bam_seqi(q,query_pos + offset)];
                    std::string base(1, base_chr);

                    //std::cout << "flag 1" << std::endl;
                    SomaticJudgeSnpHP(currentVariantIter, vcfSet , chrName, base, hpCount, norCountPS, tumCountPS, &variantsHP, nullptr, nullptr, &((*chrPosReadCase)[chrName]));
                    if((*chrPosReadCase)[chrName].find((*currentVariantIter).first) != (*chrPosReadCase)[chrName].end()){

                        //record the somatic snp derive by which germline hp in this read
                        if((*chrPosReadCase)[chrName][(*currentVariantIter).first].isHighConSomaticSNP){
                            int deriveByHp = (*chrPosReadCase)[chrName][(*currentVariantIter).first].somaticReadDeriveByHP;
                            somaticVarDeriveHP[(*currentVariantIter).first] = deriveByHp;
                        }
                    }   
                    
                    highConSomaticData.recordRefAltAlleleCount(chrName, base, currentVariantIter);
                }
                currentVariantIter++;
            }

            query_pos += length;
            ref_pos += length;
        }
            // 1: insertion to the reference
        else if( cigar_op == 1 ){
            query_pos += length;
        }
            // 2: deletion from the reference
        else if( cigar_op == 2 ){
            while( currentVariantIter != currentChrVariants.end() && (*currentVariantIter).first < ref_pos + length){
                highConSomaticData.recordDelReadCount(chrName, currentVariantIter);
                currentVariantIter++;
            }
            ref_pos += length;
        }
            // 3: skipped region from the reference
        else if( cigar_op == 3 ){
            ref_pos += length;
        }
            // 4: soft clipping (clipped sequences present in SEQ)
        else if( cigar_op == 4 ){
            query_pos += length;
        }
            // 5: hard clipping (clipped sequences NOT present in SEQ)
            // 6: padding (silent deletion from padded reference)
        else if( cigar_op == 5 || cigar_op == 6 ){
            // do nothing
        }
        else{
            std::cerr<< "alignment find unsupported CIGAR operation from read: " << bam_get_qname(&aln) << "\n";
            exit(1);
        }
    }

    //In the current version, only normal SVs are considered, without inclusion of tumor samples
    auto readIter = vcfSet[tagGeneType].readSVHapCount.find(bam_get_qname(&aln));
    if( readIter != vcfSet[tagGeneType].readSVHapCount.end() ){
        hpCount[1] += vcfSet[tagGeneType].readSVHapCount[bam_get_qname(&aln)][0];
        hpCount[2] += vcfSet[tagGeneType].readSVHapCount[bam_get_qname(&aln)][1];
    }

    int startPos = aln.core.pos + 1;
    int endPos = ref_pos;

    //the similarity of HP types
    //tumor variants
    double norHPsimilarity = 0.0;
    double tumHPsimilarity = 0.0;

    // determine the haplotype of the read
    int hpResult = ReadHP::unTag;
    hpResult = determineReadHP(hpCount, pqValue, norCountPS, norHPsimilarity, tumHPsimilarity, percentageThreshold, &totalHighSimilarity, &totalCrossTwoBlock, &totalWithOutVaraint);

    //Record read HP result before correction hp result
    if(!somaticVarDeriveHP.empty()){
        for(auto somaticVarIter : somaticVarDeriveHP){
            int pos = somaticVarIter.first;
            int baseHP = SnpHP::NONE_SNP;
            //if current read have somatic SNP then record the somatic SNP
            if(variantsHP.find(pos) != variantsHP.end()){
                baseHP = variantsHP[pos];
            }
            recordReadHp(pos, hpResult, baseHP, (*beforeCorrReadHpResult)[chrName]);
        }
    }

    //correction read HP result by deriveByHp 
    if(hpResult == ReadHP::H3){
        int deriveByH1 = 0;
        int deriveByH2 = 0;
        for(auto& somaticVarIter : somaticVarDeriveHP){
            if(somaticVarIter.second == SnpHP::GERMLINE_H1){
                deriveByH1++;
            }
            else if(somaticVarIter.second == SnpHP::GERMLINE_H2){
                deriveByH2++;
            }
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

        float deriveByHpSimilarity = (max == 0) ? 0.0 : (max / (max + min));

        if(deriveByHpSimilarity != 1.0 && deriveByHpSimilarity != 0.0){
            std::cerr << "deriveByHpSimilarity: " << deriveByHpSimilarity << std::endl;
            std::cerr << "deriveByH1: " << deriveByH1 << " deriveByH2: " << deriveByH2 << std::endl;
            std::cerr << "readID: " << bam_get_qname(&aln) << std::endl;
            exit(1);
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
    }

    if(hpCount[1] == 0 && hpCount[2] == 0){
        if(hpResult == ReadHP::H3 && hpCount[3] != 0 && hpCount[4] == 0){
            totalreadOnlyH3Snp++;
        }
    }

    //Record read HP result after correction hp result
    if(!somaticVarDeriveHP.empty()){
        for(auto somaticVarIter : somaticVarDeriveHP){
            int pos = somaticVarIter.first;
            int baseHP = SnpHP::NONE_SNP;
            //if current read have somatic SNP then record the somatic SNP
            if(variantsHP.find(pos) != variantsHP.end()){
                baseHP = variantsHP[pos];
            }
            recordReadHp(pos, hpResult, baseHP, (*afterCorrReadHpResult)[chrName]);

            // update cover region at somatic position
            if(hpResult != ReadHP::unTag){      
                
                if((*afterCorrReadHpResult)[chrName][pos].coverRegionStartPos > startPos){
                    (*afterCorrReadHpResult)[chrName][pos].coverRegionStartPos = startPos;
                }
                if((*afterCorrReadHpResult)[chrName][pos].coverRegionEndPos < endPos){
                    (*afterCorrReadHpResult)[chrName][pos].coverRegionEndPos = endPos;
                }
            }
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
            if(tumCountPS.size() != 0){
                psIter = tumCountPS.begin();
                psResultStr = std::to_string((*psIter).first);
                psValue = (*psIter).first;    
            }
            // the read only had unphased hetero tumor SNPs or homo SNPs
            else{
                psResultStr= "*";
                psValue = -1; //for determining not write PS in the current read 
            }
        }else if(hpResult == ReadHP::H1 || hpResult == ReadHP::H2){
            psIter = norCountPS.begin();
            psResultStr = std::to_string((*psIter).first);
            psValue = (*psIter).first;
        }
    }

    std::string readID = bam_get_qname(&aln);
    
    //record somatic read
    highConSomaticData.recordCrossingHighConSnpRead(chrName, readID, hpResultStr, variantsHP, currentChrVariants);
    highConSomaticData.recordTaggedSomaticRead(chrName, readID , hpResultStr, variantsHP, currentChrVariants);

    //write tag log file
    if(tagResult!=NULL){

        (*tagResult)<< readID                                       << "\t"
                    << bamHdr.target_name[aln.core.tid]             << "\t"
                    << aln.core.pos                                 << "\t"
                    << norHPsimilarity                              << "\t"
                    << "H" << hpResultStr                           << "\t"
                    << psResultStr                                  << "\t"
                    << hpCount[1]+hpCount[2]+hpCount[3]+hpCount[4]  << "\t" //modify
                    << hpCount[1]                                   << "\t"
                    << hpCount[2]                                   << "\t"
                    << hpCount[3]                                   << "\t" //new
                    << hpCount[4]                                   << "\t" //new
                    << pqValue                                      << "\t\t";
                    // << "(" << startPos << "," << endPos << ")" << "\t"
                    // << "(" << 0 << "," << query_pos << ")" << "\t\t";

        // print position and HP
        for(auto v : variantsHP ){
            (*tagResult)<< " " << v.first << "," << v.second ;
        }

        (*tagResult) << "\t";

        // belong PS, number of variant
        (*tagResult)<< "NorPS:";
        for(auto v : norCountPS){
            (*tagResult)<< " " << v.first << "," << v.second ;
        }

        (*tagResult)<< " TumPS:";
        for(auto v : tumCountPS){
            (*tagResult)<< " " << v.first << "," << v.second ;
        }

        (*tagResult)<< "\n";

    }

    return hpResult;
}

std::string HaplotagProcess::convertHpResultToString(int hpResult){
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

void HaplotagProcess::OnlyTumorSNPjudgeHP(const std::string &chrName, int &curPos, RefAltSet &curVar, std::string base, VCF_Info *vcfSet, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *readPosHP3, BamBaseCounter *NorBase, std::map<int, HP3_Info> *SomaticPos){

    if(SomaticPos == nullptr){
        std::cerr << "ERROR (SomaticTaggingJudgeHP) => SomaticPos pointer cannot be nullptr"<< std::endl;
        exit(1);
    }

    if((*SomaticPos).find(curPos) != (*SomaticPos).end()){
        if((*SomaticPos)[curPos].isHighConSomaticSNP){
            std::string TumorRefBase = curVar.Variant[Genome::TUMOR].Ref;
            std::string TumorAltBase = curVar.Variant[Genome::TUMOR].Alt; 

            if(base == TumorAltBase){
                hpCount[3]++;
                if(variantsHP != nullptr) (*variantsHP)[curPos] = SnpHP::SOMATIC_H3;

            //base is not match to TumorRefBase & TumorAltBase (other HP)
            }else if(base != TumorRefBase && base != TumorAltBase){
                //hpCount[4]++;
                //if(variantsHP != nullptr) (*variantsHP)[curPos] = SnpHP::SOMATIC_H4;
                //std::cerr << "Somatic SNP: " << curPos << " " << base << " -> " << TumorRefBase << "|" << TumorAltBase << std::endl;
            }

            if(curVar.Variant[Genome::TUMOR].is_phased_hetero){
                if(tumCountPS != nullptr) (*tumCountPS)[vcfSet[Genome::TUMOR].chrVariantPS[chrName][curPos]]++;
            }
        }
    }
}


HaplotagProcess::HaplotagProcess():
totalAlignment(0),totalSupplementary(0),totalSecondary(0),totalUnmapped(0),totalTagCount(0),totalUnTagCount(0),processBegin(time(NULL))
{
    //initialize variable
    chrVec = nullptr;
    chrLength = nullptr;

    vcfSet[Genome::NORMAL].gene_type = Genome::NORMAL;
    vcfSet[Genome::TUMOR].gene_type = Genome::TUMOR;

    mergedChrVarinat = new std::map<std::string, std::map<int, RefAltSet>>();;
    chrPosReadCase = new std::map<std::string, std::map<int, HP3_Info>>();

    beforeCorrReadHpResult = new std::map<std::string, std::map<int, ReadHpResult>>();
    afterCorrReadHpResult = new std::map<std::string, std::map<int, ReadHpResult>>();

    hpBeforeInheritance = new readHpDistriLog();
    hpAfterInheritance = new readHpDistriLog();

    //verification variable
    totalLowerQuality = 0;
    totalOtherCase = 0;
    totalunTag_HP0 = 0;
    totalreadOnlyH3Snp = 0;
    totalHighSimilarity = 0;
    totalCrossTwoBlock = 0;
    totalEmptyVariant = 0;
    totalWithOutVaraint = 0;
}

HaplotagProcess::~HaplotagProcess(){
    std::cerr<< "-------------------------------------------\n";
    std::cerr<< "total process time:    " << difftime(time(NULL), processBegin) << "s\n";
    std::cerr<< "total alignment:       " << totalAlignment     << "\n";
    std::cerr<< "total supplementary:   " << totalSupplementary << "\n";
    std::cerr<< "total secondary:       " << totalSecondary     << "\n";
    std::cerr<< "total unmapped:        " << totalUnmapped      << "\n";
    std::cerr<< "total tag alignment:   " << totalTagCount     << "\n";
    std::cerr<< "    L----total HP1   : " << totalHpCount[ReadHP::H1]     << "\n";   //new
    std::cerr<< "    L----total HP2   : " << totalHpCount[ReadHP::H2]     << "\n";   //new
    std::cerr<< "    L----total HP1-1 : " << totalHpCount[ReadHP::H1_1]   << "\n";   //new
    //std::cerr<< "    L----total HP1-2 : " << totalHpCount[ReadHP::H1_2]   << "\n";   //new
    std::cerr<< "    L----total HP2-1 : " << totalHpCount[ReadHP::H2_1]   << "\n";   //new
    //std::cerr<< "    L----total HP2-2 : " << totalHpCount[ReadHP::H2_2]   << "\n";   //new
    std::cerr<< "    L----total HP3   : " << totalHpCount[ReadHP::H3]     << "\n";   //new
    std::cerr<< "         L----total read only H3 Snp : " << totalreadOnlyH3Snp << "\n";   //new
    std::cerr<< "total untagged:        " << totalUnTagCount   << "\n";
    std::cerr<< "    L----total lower mapping quality:    " << totalLowerQuality   << "\n";   //new
    std::cerr<< "    L----total EmptyVariant:             " << totalEmptyVariant   << "\n";   //new
    std::cerr<< "    L----total start > last variant pos: " << totalOtherCase   << "\n";   //new
    std::cerr<< "    L----total judge to untag:           " << totalunTag_HP0   << "\n";   //new
    std::cerr<< "         L----total HighSimilarity:      " << totalHighSimilarity   << "\n";   //new
    std::cerr<< "         L----total CrossTwoBlock:       " << totalCrossTwoBlock   << "\n";   //new
    std::cerr<< "         L----total WithOut Variant:     " << totalWithOutVaraint   << "\n";   //new
    std::cerr<< "-------------------------------------------\n";

    delete mergedChrVarinat;
    delete chrPosReadCase;

    delete beforeCorrReadHpResult;
    delete afterCorrReadHpResult;

    delete hpBeforeInheritance;
    delete hpAfterInheritance;
};


void HaplotagProcess::TaggingProcess(HaplotagParameters &params)
{
    std::cerr<< "phased SNP file:       " << params.snpFile             << "\n";
    if(params.tagTumorSnp) 
    std::cerr<< "phased tumor SNP file: " << params.tumorSnpFile        << "\n";  
    std::cerr<< "phased SV file:        " << params.svFile              << "\n";
    std::cerr<< "phased MOD file:       " << params.modFile             << "\n";
    std::cerr<< "input bam file:        " << params.bamFile             << "\n";
    if(params.tagTumorSnp) 
    std::cerr<< "input tumor bam file:  " << params.tumorBamFile        << "\n"; 
    std::cerr<< "input ref file:        " << params.fastaFile           << "\n";
    std::cerr<< "output bam file:       " << params.resultPrefix + "." + params.outputFormat << "\n";
    std::cerr<< "number of threads:     " << params.numThreads          << "\n";
    std::cerr<< "write log file:        " << (params.writeReadLog ? "true" : "false") << "\n";
    std::cerr<< "log file:              " << (params.writeReadLog ? (params.resultPrefix+".out") : "") << "\n";
    std::cerr<< "-------------------------------------------\n";
    std::cerr<< "somatic mode:                    " << (params.tagTumorSnp ? "true" : "false") << "\n"; 
    std::cerr<< "tag region:                      " << (!params.region.empty() ? params.region : "all") << "\n";
    if(params.tagTumorSnp)
    std::cerr<< "somatic calling mapping quality: " << params.somaticCallingMpqThreshold    << "\n"; 
    std::cerr<< "filter mapping quality below:    " << params.qualityThreshold    << "\n";
    std::cerr<< "percentage threshold:            " << params.percentageThreshold << "\n";
    std::cerr<< "tag supplementary:               " << (params.tagSupplementary ? "true" : "false") << "\n";
    std::cerr<< "-------------------------------------------\n";
    
    tagTumorMode=params.tagTumorSnp;
    // decide on the type of tagging for VCF and BAM files
    int tagGeneType;

    if(tagTumorMode){
        tagGeneType = Genome::TUMOR;
    }else{
        tagGeneType = Genome::NORMAL;
    }

    VcfParser vcfParser(tagTumorMode);

    if(tagTumorMode){
        //load seqc high con file for benchmarking
        if(params.seqcHighCon != ""){
            std::time_t begin = time(NULL);
            std::cerr<< "loading high confidence SNP ... ";
            highConSomaticData.setTestingFunc(true);
            highConSomaticData.loadHighConSomatic(params.seqcHighCon, vcfSet[Genome::SEQC_HIGH_CON], *mergedChrVarinat);
            std::cerr<< difftime(time(NULL), begin) << "s\n";
            highConSomaticData.displaySomaticVarCount(vcfSet[Genome::SEQC_HIGH_CON].chrVec, *mergedChrVarinat);
        }
        //load tumor snp vcf
        if(params.tumorSnpFile != ""){
            std::time_t begin = time(NULL);
            std::cerr<< "parsing tumor SNP VCF ... ";
            vcfParser.setParseSnpFile(true);
            vcfParser.variantParser(params.tumorSnpFile, vcfSet[Genome::TUMOR], *mergedChrVarinat);
            vcfParser.reset();
            std::cerr<< difftime(time(NULL), begin) << "s\n";
        }
    }

    // load SNP vcf file
    std::time_t begin = time(NULL);

    if(tagTumorMode){
        std::cerr<< "parsing normal SNP VCF ... ";
    }else{
        std::cerr<< "parsing SNP VCF ... ";        
    }

    vcfParser.setParseSnpFile(true);
    vcfParser.variantParser(params.snpFile, vcfSet[Genome::NORMAL], *mergedChrVarinat);
    vcfParser.reset();
    std::cerr<< difftime(time(NULL), begin) << "s\n";

    // int tumor_snp_count = 0;
    // int normal_snp_count = 0;
    // for(auto& chrIter : vcfSet[Genome::NORMAL].chrVec){
    //     auto chrVarIter = (*mergedChrVarinat)[chrIter].begin();
    //     while(chrVarIter != (*mergedChrVarinat)[chrIter].end()){
    //         if((*chrVarIter).second.isExistTumor){
    //             tumor_snp_count++;
    //         }
    //         if((*chrVarIter).second.isExistNormal){
    //             normal_snp_count++;
    //         }
    //         chrVarIter++;
    //     }
    // }
    // std::cerr << "Normal SNP count: " << normal_snp_count << std::endl;
    // std::cerr << "Tumor SNP count: " << tumor_snp_count << std::endl;

    // load SV vcf file
    if(params.svFile!=""){
        begin = time(NULL);
        std::cerr<< "parsing SV VCF ... ";
        vcfParser.setParseSVFile(true);
        vcfParser.variantParser(params.svFile, vcfSet[tagGeneType], *mergedChrVarinat);
        vcfParser.reset();
        std::cerr<< difftime(time(NULL), begin) << "s\n";    
    }
    
    // load MOD vcf file
    if(params.modFile!=""){
        begin = time(NULL);
        std::cerr<< "parsing MOD VCF ... ";
        vcfParser.setParseMODFile(true);
        vcfParser.variantParser(params.modFile, vcfSet[tagGeneType], *mergedChrVarinat);
        vcfParser.reset();
        std::cerr<< difftime(time(NULL), begin) << "s\n";    
    }

    //decide which genome type chrVec and chrLength belong to
    if(tagGeneType == Genome::NORMAL){
        chrVec = &(vcfSet[Genome::NORMAL].chrVec);
        chrLength = &(vcfSet[Genome::NORMAL].chrLength); 
    }else{

        //check normal & tumor chr & length
        for(auto& chrIter : vcfSet[Genome::TUMOR].chrLength){
            if(vcfSet[Genome::NORMAL].chrLength.find(chrIter.first) != vcfSet[Genome::NORMAL].chrLength.end()){
                if(chrIter.second != vcfSet[Genome::NORMAL].chrLength[chrIter.first]){
                    std::cerr << "tumor & normal VCFs chromosome length are not the same" << std::endl;
                    std::cerr << "chr: " << chrIter.first << " length: " << chrIter.second << std::endl;
                    return ;
                }
            }else{
                std::cerr << "tumor & normal VCFs chromosome count are not the same" << std::endl;
                return ;
            }
        }

        chrVec = &(vcfSet[Genome::TUMOR].chrVec);
        chrLength = &(vcfSet[Genome::TUMOR].chrLength); 
    }

    // update chromosome processing based on region
    if (!params.region.empty()) {
        auto colonPos = params.region.find(":");
        std::string regionChr;
        if (colonPos != std::string::npos) {
            regionChr = params.region.substr(0, colonPos);
        }
        else {
            regionChr = params.region;
        }
        auto chrVecIter = std::find((*chrVec).begin(), (*chrVec).end(), regionChr);
        if (chrVecIter != (*chrVec).end()) {
            (*chrVec) = std::vector<std::string>{regionChr};
        } else {
            std::cerr << "ERROR: Incorrect chromosome for input region\n" << std::endl;
            exit(1);
        }
    }

    //somatic SNPs calling
    if(tagTumorMode){

        //Count each base numbers at tumor SNP position in the Normal.bam
        BamBaseCounter *NorBase = new BamBaseCounter();
        NorBase->CountingBamBase(params.bamFile, params, (*mergedChrVarinat), *chrVec, *chrLength, Genome::TUMOR);

        //record the HP3 confidence of each read
        SomaticVarCaller *SomaticVar = new SomaticVarCaller();
        SomaticVar->VariantCalling(params.tumorBamFile, (*mergedChrVarinat), *chrVec, *chrLength, params, vcfSet, *NorBase);
        (*chrPosReadCase) = SomaticVar->getSomaticChrPosInfo();

        delete NorBase;
        delete SomaticVar;
        NorBase = nullptr;
        SomaticVar = nullptr;
        //return;
    }

    // tag read
    begin = time(NULL);

    if(tagGeneType == Genome::TUMOR){
        std::cerr<< "somatic tagging start ...\n";
    }else{
        std::cerr<< "tag read start ...\n";
    }
    tagRead(params, tagGeneType);

    std::cerr<< "tag read " << difftime(time(NULL), begin) << "s\n";

    return;
};



