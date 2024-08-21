#include "HaplotagProcess.h"


BamBaseCounter::BamBaseCounter(std::string BamFile, const HaplotagParameters &params, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat, std::vector<std::string> &chrVec, std::map<std::string, int> &chrLength, int genmoeType){
    std::cerr<< "calculating base information in the normal genome bam ... ";
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
        ChrVariantBase[chr] = std::map<int, PosBase>();
    }

    // init data structure and get core n
    htsThreadPool threadPool = {NULL, 0};
    // creat thread pool
    if (!(threadPool.pool = hts_tpool_init(params.numThreads))) {
        fprintf(stderr, "Error creating thread pool\n");
    }

    // loop all chromosome
    #pragma omp parallel for schedule(dynamic) num_threads(params.numThreads) firstprivate(mergedChrVarinat)    
    for(auto chr : chrVec ){

        //int numThreads = omp_get_num_threads();
        //int threadId = omp_get_thread_num();
        //#pragma omp critical
        //{
        //    std::cout << chr << "  threads: " << numThreads << "  ID: "<< threadId <<std::endl;
        //}

        // open bam file
        samFile *in = hts_open(BamFile.c_str(), "r");
        // load reference file
        hts_set_fai_filename(in, params.fastaFile.c_str() );
        // input reader
        bam_hdr_t *bamHdr = sam_hdr_read(in);
        // bam file index
        hts_idx_t *idx = NULL;
        // check input bam file
        if (in == NULL) {
            std::cerr<<"ERROR: Cannot open bam file " << BamFile.c_str() << "\n";
            exit(1);
        }
        // check bam file index
        if ((idx = sam_index_load(in, BamFile.c_str())) == 0) {
            std::cerr<<"ERROR: Cannot open index for bam file " << BamFile.c_str() << "\n";
            exit(1);
        }

        // set thread
        hts_set_opt(in, HTS_OPT_THREAD_POOL, &threadPool);
        // initialize an alignment
        bam1_t *aln = bam_init1();

        //store base information
        std::map<int, PosBase> *variantBase = nullptr;

        variantBase = &(ChrVariantBase[chr]);

        //inintial iterator
        // variant position (0-base), allele haplotype set
        std::map<int, RefAltSet> currentVariants = mergedChrVarinat[chr];

        std::map<int, RefAltSet>::iterator firstVariantIter = currentVariants.begin();

        std::map<int, RefAltSet>::reverse_iterator last = currentVariants.rbegin();

        std::string region = !params.region.empty() ? params.region : chr + ":1-" + std::to_string(chrLength[chr]);
        hts_itr_t* iter = sam_itr_querys(idx, bamHdr, region.c_str());

        while (sam_itr_multi_next(in, iter, aln) >= 0) {
            
            int flag = aln->core.flag;

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
            else if(int(aln->core.pos) <= (*last).first){
                StatisticBaseInfo(*aln, chr, params, genmoeType, *variantBase, currentVariants, firstVariantIter);
            }
        }
        //std::cerr<<"calculate Base Information ...\n";
        CalculateBaseInfo(chr, *variantBase, currentVariants);
        
        variantBase = nullptr;
        hts_itr_destroy(iter);
        hts_idx_destroy(idx);
        bam_hdr_destroy(bamHdr);
        bam_destroy1(aln);
        sam_close(in);
    }
    hts_tpool_destroy(threadPool.pool);
    std::cerr<< difftime(time(NULL), begin) << "s\n";
    return;
};

BamBaseCounter::~BamBaseCounter(){

};


void BamBaseCounter::StatisticBaseInfo(const bam1_t &aln, std::string chrName, const HaplotagParameters &params, int genmoeType, std::map<int, PosBase> &variantBase, std::map<int, RefAltSet> &currentVariants ,std::map<int, RefAltSet>::iterator &firstVariantIter){
    
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
                        if(tumRefLength == 1 && tumRefLength == 1){

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

void BamBaseCounter::CalculateBaseInfo(std::string chr, std::map<int, PosBase> &VariantBase, std::map<int, RefAltSet> &currentVariants){
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
    if(ChrVariantBase[chr].find(pos) == ChrVariantBase[chr].end()){
        return false;
    }    
    return ChrVariantBase[chr][pos].isHighRefAllelleFreq;
}

std::string BamBaseCounter::getMaxFreqBase(std::string chr, int pos){
    if(ChrVariantBase[chr].find(pos) == ChrVariantBase[chr].end()){
        std::cout << "ERROR (getMaxBase) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return ChrVariantBase[chr][pos].max_base;
}

float BamBaseCounter::getMaxBaseRatio(std::string chr, int pos){
    if(ChrVariantBase[chr].find(pos) == ChrVariantBase[chr].end()){
        std::cout << "ERROR (getMaxBaseRatio) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return ChrVariantBase[chr][pos].max_ratio;
}

float BamBaseCounter::getSecondMaxBaseRatio(std::string chr, int pos){
    if(ChrVariantBase[chr].find(pos) == ChrVariantBase[chr].end()){
        std::cout << "ERROR (getSecondMaxBaseRatio) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return ChrVariantBase[chr][pos].second_max_ratio;
}

float BamBaseCounter::getLowMpqReadRatio(std::string chr, int pos){
    if(ChrVariantBase[chr].find(pos) == ChrVariantBase[chr].end()){
        std::cout << "ERROR (getLowMpqReadRatio) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return ChrVariantBase[chr][pos].lowMpqReadRatio;
}

float BamBaseCounter::getVAF(std::string chr, int pos){
    if(ChrVariantBase[chr].find(pos) == ChrVariantBase[chr].end()){
        std::cout << "ERROR (getVAF) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return ChrVariantBase[chr][pos].VAF;
}

float BamBaseCounter::getFilterdMpqVAF(std::string chr, int pos){
    if(ChrVariantBase[chr].find(pos) == ChrVariantBase[chr].end()){
        std::cout << "ERROR (getFilterdMpqVAF) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return ChrVariantBase[chr][pos].filteredMpqVAF;
}

int BamBaseCounter::getBaseAcount(std::string chr, int pos){
    if(ChrVariantBase[chr].find(pos) == ChrVariantBase[chr].end()){
        std::cout << "ERROR (getBaseAcount) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return ChrVariantBase[chr][pos].A_count;
}

int BamBaseCounter::getBaseCcount(std::string chr, int pos){
    if(ChrVariantBase[chr].find(pos) == ChrVariantBase[chr].end()){
        std::cout << "ERROR (getBaseCcount) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return ChrVariantBase[chr][pos].C_count;
}

int BamBaseCounter::getBaseTcount(std::string chr, int pos){
    if(ChrVariantBase[chr].find(pos) == ChrVariantBase[chr].end()){
        std::cout << "ERROR (getBaseTcount) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return ChrVariantBase[chr][pos].T_count;
}

int BamBaseCounter::getBaseGcount(std::string chr, int pos){
    if(ChrVariantBase[chr].find(pos) == ChrVariantBase[chr].end()){
        std::cout << "ERROR (getBaseGcount) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return ChrVariantBase[chr][pos].G_count;
}

int BamBaseCounter::getDepth(std::string chr, int pos){
    if(ChrVariantBase[chr].find(pos) == ChrVariantBase[chr].end()){
        std::cout << "ERROR (getDepth) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return ChrVariantBase[chr][pos].depth;
}

int BamBaseCounter::getMpqDepth(std::string chr, int pos){
    if(ChrVariantBase[chr].find(pos) == ChrVariantBase[chr].end()){
        std::cout << "ERROR (getMpqDepth) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return ChrVariantBase[chr][pos].filteredMpqDepth;
}

int BamBaseCounter::getVarDeletionCount(std::string chr, int pos){
    if(ChrVariantBase[chr].find(pos) == ChrVariantBase[chr].end()){
        std::cout << "ERROR (getDelCount) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return ChrVariantBase[chr][pos].delCount;
}

void BamBaseCounter::displayPosInfo(std::string chr, int pos){
    if(ChrVariantBase[chr].find(pos) == ChrVariantBase[chr].end()){
        std::cout << "ERROR (displayPosInfo) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }else{
        PosBase VarBase = ChrVariantBase[chr][pos];
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

void SomaticJudgeBase::SomaticJudgeSnpHP(std::map<int, RefAltSet>::iterator &currentVariantIter, VCF_Info *vcfSet, std::string chrName, std::string base, std::map<int, int> &hpCount, std::map<int, int> &norCountPS, std::map<int, int> &tumCountPS, std::map<int,std::string> *variantsHP, std::vector<int> *readPosHP3, BamBaseCounter *NorBase, std::map<int, HP3_Info> *SomaticPos){
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
                    if(tumHP1 == base){
                        if(variantsHP != nullptr) (*variantsHP)[curPos] = "(1,1)";
                    }else if(tumHP2 == base){
                        if(variantsHP != nullptr) (*variantsHP)[curPos] = "(1,2)";
                    }
                }else if(norHP2 == base){                
                    hpCount[2]++;
                    //Germline mutation
                    if(tumHP1 == base){
                        if(variantsHP != nullptr) (*variantsHP)[curPos] = "(2,1)";
                    }else if(tumHP2 == base){
                        if(variantsHP != nullptr) (*variantsHP)[curPos] = "(2,2)";
                    }
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
                    if(variantsHP != nullptr) (*variantsHP)[curPos] = "(1,*)";
                }else if(norHP2 == base){                
                    hpCount[2]++;
                    //Germline mutation
                    if(variantsHP != nullptr) (*variantsHP)[curPos] = "(2,*)";
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
                    if(variantsHP != nullptr) (*variantsHP)[curPos] = "(1,*)";
                }else if(norHP2 == base){                
                    hpCount[2]++;
                    //Germline mutation
                    if(variantsHP != nullptr) (*variantsHP)[curPos] = "(2,*)";
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
                    if(variantsHP != nullptr) (*variantsHP)[curPos] = '1';
                }
                if(base == norHP2){
                    hpCount[2]++;
                    if(variantsHP != nullptr) (*variantsHP)[curPos] = '2';
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
                    OnlyTumorSNPjudgeHP(curVar, curPos, NorBase, vcfSet, chrName, base, hpCount, &tumCountPS, variantsHP, readPosHP3, SomaticPos);
                }
            }
        //the tumor SNP GT is unphased heterozygous
        }else if(curVar.Variant[Genome::TUMOR].is_unphased_hetero == true){
            if(curVar.Variant[Genome::TUMOR].Ref == base || curVar.Variant[Genome::TUMOR].Alt == base){
                OnlyTumorSNPjudgeHP(curVar, curPos, NorBase, vcfSet, chrName, base, hpCount, nullptr, variantsHP, readPosHP3, SomaticPos);
            }           
        //the tumor SNP GT is homozygous
        }else if(curVar.Variant[Genome::TUMOR].is_homozygous == true){
            if(curVar.Variant[Genome::TUMOR].Ref == base || curVar.Variant[Genome::TUMOR].Alt == base){
                OnlyTumorSNPjudgeHP(curVar, curPos, NorBase, vcfSet, chrName, base, hpCount, nullptr, variantsHP, readPosHP3, SomaticPos);
            }
        }
    }
}

void SomaticJudgeBase::OnlyTumorSNPjudgeHP(RefAltSet &curVar, int &curPos, BamBaseCounter *NorBase, VCF_Info vcfSet[2], std::string chrName, std::string base, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int,std::string> *variantsHP, std::vector<int> *readPosHP3, std::map<int, HP3_Info> *SomaticPos){

}

std::string SomaticJudgeBase::determineReadHP(std::map<int, int> &hpCount, int &pqValue, std::map<int, int> &norCountPS, double percentageThreshold, int *totalHighSimilarity, int *totalCrossTwoBlock, int *totalWithOutVaraint){
    double normalMinHPcount = 0;
    double normalMaxHPcount = 0;
    std::string maxNormalHP = "0";

    double tumorMinHPcount = 0;
    double tumorMaxHPcount = 0;
    int maxTumorHP = 0;

    // determine max and min
    //tumor HP count
    if(hpCount[3] > hpCount[4]){
        tumorMinHPcount = hpCount[4];
        tumorMaxHPcount = hpCount[3];
        maxTumorHP = 3;
    }
    else{
        tumorMinHPcount = hpCount[3];
        tumorMaxHPcount = hpCount[4];
        maxTumorHP = 4;
    }

    //normal HP count
    if(hpCount[1] > hpCount[2]){
        normalMinHPcount = hpCount[2];
        normalMaxHPcount = hpCount[1];
        maxNormalHP = "1";
    }
    else{
        normalMinHPcount = hpCount[1];
        normalMaxHPcount = hpCount[2];
        maxNormalHP = "2";
    }

    //the similarity of HP types
    //tumor variants
    double tumorHPsimilarity = (tumorMaxHPcount == 0) ? 0.0 : tumorMaxHPcount/(tumorMaxHPcount+tumorMinHPcount);
    double normalHPsimilarity = (normalMaxHPcount == 0) ? 0.0 : normalMaxHPcount/(normalMaxHPcount+normalMinHPcount);


    // determine the haplotype of the read
    std::string hpResult = "0";

    // read have tumor variant
    if(tumorMaxHPcount != 0){
        //check the similarity of HP types in the tumor
        if(tumorHPsimilarity >= percentageThreshold){
            //check the similarity of HP types in the normal
            if(normalHPsimilarity >= percentageThreshold){
                switch(maxTumorHP){
                    case 3:
                        hpResult = maxNormalHP + "-1"; break;
                    case 4:
                        hpResult = maxNormalHP + "-2"; break;
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
                    case 3:
                        hpResult = "3"; break;
                    case 4:
                        hpResult = "4"; break;
                    default:
                        std::cerr << "ERROR: Unexpected haplotype : tumor Max HP= "<< maxTumorHP << std::endl;
                        exit(1);
                        break;
                }
            }
        }else{
            // no tag
            pqValue = 0;
            (*totalHighSimilarity)++;
        }
    }
    // read haven't tumor variant
    else if(normalMaxHPcount != 0){
        if(normalHPsimilarity >= percentageThreshold){
            hpResult = maxNormalHP;
        }else{
            // no tag
            pqValue = 0;
            (*totalHighSimilarity)++;
        }
    }

    // cross two block
    // had at least one tumor SNP in the current read
    if(hpCount[3] != 0 && hpCount[4] != 0 ){

    // all SNPs are germline variants in the current read
    }else{
        if(norCountPS.size() > 1){
            hpResult = "0";
            (*totalCrossTwoBlock)++;
        }
    }

    // determine the quality of the read
    // There haven't been any variants in the current read
    if(normalMaxHPcount == 0 && tumorMaxHPcount == 0){
        //std::cerr<< "Read max = 0 : " << bam_get_qname(&aln) << "\n";
        (*totalWithOutVaraint)++;
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


SomaticVarCaller::SomaticVarCaller(std::string BamFile, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat, std::vector<std::string> &chrVec, std::map<std::string, int> &chrLength, HaplotagParameters &params, VCF_Info *vcfSet, BamBaseCounter &NorBase){
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
        SomaticChrPosInfo[chr] = std::map<int, HP3_Info>();
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
    }

    // loop all chromosome
    #pragma omp parallel for schedule(dynamic) num_threads(params.numThreads) firstprivate(mergedChrVarinat)
    for(auto chr : chrVec ){

        // open bam file
        samFile *in = hts_open(BamFile.c_str(), "r");
        // load reference file
        hts_set_fai_filename(in, params.fastaFile.c_str() );
        // input reader
        bam_hdr_t *bamHdr = sam_hdr_read(in);
        // bam file index
        hts_idx_t *idx = NULL;
        // check input bam file
        if (in == NULL) {
            std::cerr<<"ERROR: Cannot open bam file " << BamFile.c_str() << "\n";
        }
        // check bam file index
        if ((idx = sam_index_load(in, BamFile.c_str())) == 0) {
            std::cerr<<"ERROR: Cannot open index for bam file " << BamFile.c_str() << "\n";
        }

        // set thread
        hts_set_opt(in, HTS_OPT_THREAD_POOL, &threadPool);
        // initialize an alignment
        bam1_t *aln = bam_init1();

        // record the position that tagged as HP3
        // chr, variant position
        std::map<int, HP3_Info> *somaticPosInfo = nullptr;

        somaticPosInfo = &(SomaticChrPosInfo[chr]);

        // records all variants within this chromosome.
        std::map<int, RefAltSet> currentChrVariants = mergedChrVarinat[chr];
        // since each read is sorted based on the start coordinates, to save time, 
        // firstVariantIter keeps track of the first variant that each read needs to check.
        std::map<int, RefAltSet>::iterator firstVariantIter = currentChrVariants.begin();
        // get the coordinates of the last variant
        // the tagging process will not be perform if the read's start coordinate are over than last variant.
        std::map<int, RefAltSet>::reverse_iterator lastVariant = currentChrVariants.rbegin();
        
        // tagging will be attempted for reads within the specified coordinate range.
        std::string region = !params.region.empty() ? params.region : chr + ":1-" + std::to_string(chrLength[chr]);
        hts_itr_t* Iter = sam_itr_querys(idx, bamHdr, region.c_str());
        
        // iter all reads
        while (sam_itr_multi_next(in, Iter, aln) >= 0) {

            int flag = aln->core.flag;

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
            else if(int(aln->core.pos) <= (*lastVariant).first){
                //statistics for tumor SNP base count and depth, and classify cases
                StatisticSomaticPosInfo(*bamHdr, *aln, chr, params, &NorBase, vcfSet, *somaticPosInfo, currentChrVariants, firstVariantIter);
            }
        }
        //calculate information and filter somatic SNPs
        SomaticFeatureFilter(somaticParams, currentChrVariants, chr, *somaticPosInfo);

        somaticPosInfo = nullptr;
        hts_itr_destroy(Iter);
        hts_idx_destroy(idx);
        bam_hdr_destroy(bamHdr);
        bam_destroy1(aln);
        sam_close(in);
    }
    hts_tpool_destroy(threadPool.pool);
    std::cerr<< difftime(time(NULL), begin) << "s\n";
    
    if(somaticParams.writeVarLog){
        std::ofstream *tagHP3Log = new std::ofstream(params.resultPrefix+"_HP3.out");

        if(!tagHP3Log->is_open()){
            std::cerr<< "Fail to open write file: " << params.resultPrefix+"_HP3.out" << "\n";
            exit(1);
        }

        //write the log file for variants with positions tagged as HP3
        WriteSomaticVarCallingLog(params ,somaticParams, tagHP3Log, chrVec, NorBase, mergedChrVarinat);
        tagHP3Log->close();

        delete tagHP3Log;
        tagHP3Log = nullptr;
    }

    return;
}

SomaticVarCaller::~SomaticVarCaller(){

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
    somaticParams.Hetero_VAF_lower_threshold = 0.2;
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

void SomaticVarCaller::StatisticSomaticPosInfo(const  bam_hdr_t &bamHdr,const bam1_t &aln, std::string chr, HaplotagParameters &params, BamBaseCounter *NorBase, VCF_Info *vcfSet, std::map<int, HP3_Info> &somaticPosInfo, std::map<int, RefAltSet> &currentChrVariants, std::map<int, RefAltSet>::iterator &firstVariantIter){
    
    std::map<int, int> hpCount;
    hpCount[1] = 0; 
    hpCount[2] = 0;
    hpCount[3] = 0;
    hpCount[4] = 0;

    //record variants on this read
    std::vector<int> readPosHP3;

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

                    //int tumRefLength = (*currentVariantIter).second.Variant[Gene::TUMOR].Ref.length();
                    //int tumAltLength = (*currentVariantIter).second.Variant[Gene::TUMOR].Alt.length();

                    //waring : using ref length to split SNP and indel that will be effect case ratio result 
                    if ( aln.core.qual >= params.somaticCallingMpqThreshold ){

                        SomaticJudgeSnpHP(currentVariantIter, vcfSet , chr, base, hpCount, NorCountPS, TumCountPS, nullptr, &readPosHP3, NorBase, &somaticPosInfo);
                        
                        //statistically analyze SNP information exclusive to the tumor
                        if((*currentVariantIter).second.isExistTumor && !((*currentVariantIter).second.isExistNormal)){
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
                    }  

                    //statistically analyze SNP information exclusive to the tumor
                    if((*currentVariantIter).second.isExistTumor && !((*currentVariantIter).second.isExistNormal)){
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
        ClassifyReadsByCase(readPosHP3, TumCountPS, hpCount, params, somaticPosInfo);
    }         
}

void SomaticVarCaller::OnlyTumorSNPjudgeHP(RefAltSet & curVar, int & curPos, BamBaseCounter * NorBase, VCF_Info * vcfSet, std::string chrName, std::string base, std::map<int,int>& hpCount, std::map<int,int>* tumCountPS, std::map<int,std::string>* variantsHP, std::vector<int>* readPosHP3, std::map<int,HP3_Info>* SomaticPos){
    //the tumor SNP GT is phased heterozygous
    //all bases of the same type at the current position in normal.bam

    if(readPosHP3 == nullptr){
        std::cout << "ERROR (SomaticDetectJudgeHP) => readPosHP3 pointer cannot be nullptr"<< std::endl;
        exit(1);
    }
    if(SomaticPos == nullptr){
        std::cout << "ERROR (SomaticDetectJudgeHP) => SomaticPos pointer cannot be nullptr"<< std::endl;
        exit(1);
    }
    
    if((*NorBase).isHighRefAllelleFreq(chrName, curPos) == true){
        std::string TumorRefBase = curVar.Variant[Genome::TUMOR].Ref;
        std::string TumorAltBase = curVar.Variant[Genome::TUMOR].Alt; 

        //max count base match to refBase in normal.bam
        if((*NorBase).getMaxFreqBase(chrName, curPos) == TumorRefBase){
            if(base == TumorAltBase){
                hpCount[3]++;
                if(variantsHP != nullptr) (*variantsHP)[curPos] = '3';

                //record postions that tagged as HP3 for calculating the confidence of somatic positions
                (*readPosHP3).push_back(curPos);
                (*SomaticPos)[curPos].isNormalPosLowVAF = true;                                

            //base is not match to TumorRefBase & TumorAltBase (other HP)
            }else if(base != TumorRefBase && base != TumorAltBase){
                hpCount[4]++;
                if(variantsHP != nullptr) (*variantsHP)[curPos] = '4';
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

void SomaticVarCaller::ClassifyReadsByCase(std::vector<int> &readPosHP3, std::map<int, int> &countPS, std::map<int, int> &hpCount, const HaplotagParameters &params, std::map<int, HP3_Info> &somaticPosInfo){
    
    double max, min;
    if(hpCount[3] >= hpCount[4]){
        max = hpCount[3];
        min = hpCount[4];
    }else{
        max = hpCount[4];
        min = hpCount[3];
    }

    //decide whether to tag the read or not
    bool tagRead = true;
    if( max/(max+min) < params.percentageThreshold){
        // no tag
        tagRead = false;
    }

    //only exist tumor homo SNP
    if( countPS.size() == 0 && hpCount[3] != 0){
        //std::cerr << "Error : read only had tumor homo SNP " << bam_get_qname(aln) << "\n";
    }

    // cross two block
    //if( countPS.size() > 1){
    //    tagRead = false;
    //}

    int zero_count = 0;
    float CleanHP3Threshold = 1.0;
    bool tagCleanHP3Read = false;

    if (hpCount[1] == 0) zero_count++;
    if (hpCount[2] == 0) zero_count++;

    if(hpCount[3] == 0){
        std::cerr << "Error : hp3 count = 0";
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
        if(tagRead == false){
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

void SomaticVarCaller::SomaticFeatureFilter(const SomaticFilterParaemter &somaticParams, std::map<int, RefAltSet> &currentChrVariants, std::string chr, std::map<int, HP3_Info> &somaticPosInfo){
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
            std::cerr << "Error(calculate tumor VAF) => can't match RefBase or AltBase : chr:" << chr << " pos: " << (*somaticVarIter).first << " RefBase:" << RefBase << " AltBase:" << AltBase;
            exit(1);
        }

        (*somaticVarIter).second.base.VAF =  (float)tumAltCount / (float)tumDepth;

        // Calculating the difference in VAF between VAF and filtered low MPQ VAF
        (*somaticVarIter).second.base.lowMpqReadRatio = (float)(tumDepth - tumMpqDepth) / (float)tumDepth;

        if((*somaticVarIter).second.base.lowMpqReadRatio < 0){
            std::cerr << "Error(tumMPQReadRatio) => chr: " << chr << " pos: " << (*somaticVarIter).first;
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
            std::cerr << "Error(GTtype) => can't find GTtype at chr: " << chr << " pos: " << (*somaticVarIter).first;
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

        //filter SNP
        if( (!somaticParams.applyFilter) ||
            (
            //((*somaticVarIter).second.OnlyHP3ReadRatio < OnlyHP3ReadRatioThreshold) && 
            ((*somaticVarIter).second.MessyHPReadRatio < messyReadRatioThreshold) && 
            ((*somaticVarIter).second.CaseReadCount > readCountThreshold) &&
            ((*somaticVarIter).second.base.VAF > VAF_lower_threshold) && 
            ((*somaticVarIter).second.base.VAF < VAF_upper_threshold) &&
            ((*somaticVarIter).second.tumDelRatio < tumDeletionRatioThreshold) &&
            ((*somaticVarIter).second.base.lowMpqReadRatio <= tumLowMpqRatioThreshold)) ){
                
            (*somaticVarIter).second.isHighConSomaticSNP = true;
        }
        somaticVarIter++;          
    }
}

void SomaticVarCaller::WriteSomaticVarCallingLog(const HaplotagParameters &params, const SomaticFilterParaemter &somaticParams, std::ofstream *tagHP3Log, const std::vector<std::string> &chrVec, BamBaseCounter &NorBase, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat){
    if(tagHP3Log == NULL){
        return;
    }

    std::cout << "writing somatic varinats calling log ... ";
    std::time_t begin = time(NULL);

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
    
        std::map<int, HP3_Info>::iterator somaticVarIter = SomaticChrPosInfo[chr].begin();
        while( somaticVarIter != SomaticChrPosInfo[chr].end()){
            
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
                std::cerr << "Error(write tag HP3 log file) => can't find the position : chr:" << chr << " pos: " << (*somaticVarIter).first;
                exit(1);
            }

            if(RefBase == "" || AltBase == ""){
                std::cerr << "Error(write tag HP3 log file) => can't find RefBase or AltBase : chr:" << chr << " pos: " << (*somaticVarIter).first << " RefBase:" << RefBase << " AltBase:" << AltBase;
                exit(1);
            }

            if(AltBase == "A"){
                tumMpqAltCount = SomaticChrPosInfo[chr][(*somaticVarIter).first].base.MPQ_A_count;
                norAltCount = NorBase.getBaseAcount(chr, (*somaticVarIter).first);
            }else if(AltBase == "C"){
                tumMpqAltCount = SomaticChrPosInfo[chr][(*somaticVarIter).first].base.MPQ_C_count;
                norAltCount = NorBase.getBaseCcount(chr, (*somaticVarIter).first);
            }else if(AltBase == "T"){
                tumMpqAltCount = SomaticChrPosInfo[chr][(*somaticVarIter).first].base.MPQ_T_count;
                norAltCount = NorBase.getBaseTcount(chr, (*somaticVarIter).first);
            }else if(AltBase == "G"){
                tumMpqAltCount = SomaticChrPosInfo[chr][(*somaticVarIter).first].base.MPQ_G_count;
                norAltCount = NorBase.getBaseGcount(chr, (*somaticVarIter).first);
            }else{
                std::cerr << "Error(write tag HP3 log file) => can't match RefBase or AltBase : chr:" << chr << " pos: " << (*somaticVarIter).first << " RefBase:" << RefBase << " AltBase:" << AltBase;
                exit(1);
            }

            //calculating tumor VAF
            float tumVAF = (*somaticVarIter).second.base.VAF;
            float tumMpqVAF = (float)tumMpqAltCount / (float)tumMpqDepth;

            float norVAF =  (float)norAltCount / (float)norDepth;
            float norCountVAF = NorBase.getVAF(chr,(*somaticVarIter).first);
            float norMpqVAF = NorBase.getFilterdMpqVAF(chr,(*somaticVarIter).first);

            if(norVAF != norCountVAF){
                std::cerr << "Error(diff nor VAF) => can't find GTtype at chr: " << chr << " pos: " << (*somaticVarIter).first;
                exit(1); 
            }
            float norVAF_substract = (norMpqVAF - norCountVAF);
            float norMPQReadRatio = (float)(norDepth - norMpqDepth) / (float)norDepth;
            float TmpNorMPQReadRatio = NorBase.getLowMpqReadRatio(chr,(*somaticVarIter).first);

            if(norMPQReadRatio != TmpNorMPQReadRatio){
                std::cerr << "Error(diff low MPQ read ratio) => can't find GTtype at chr: " << chr << " pos: " << (*somaticVarIter).first << "\n";
                std::cerr << "norMPQReadRatio: " << norMPQReadRatio << " TmpNorMPQReadRatio: " << TmpNorMPQReadRatio;
                exit(1); 
            }

            // Calculating the difference in VAF between VAF and filtered low MPQ VAF
            float tumVAF_substract = (tumMpqVAF -tumVAF);
            float tumLowMpqReadRatio = (*somaticVarIter).second.base.lowMpqReadRatio;

            if(norMPQReadRatio < 0){
                std::cerr << "Error(norMPQReadRatio) => chr: " << chr << " pos: " << (*somaticVarIter).first;
                std::cerr << "norMPQReadRatio " << norMPQReadRatio << " norDepth: " << norDepth << " norMPQDepth: " << norMpqDepth;
            }

            if(tumLowMpqReadRatio < 0){
                std::cerr << "Error(tumMPQReadRatio) => chr: " << chr << " pos: " << (*somaticVarIter).first;
                std::cerr << "tumMPQReadRatio " << tumLowMpqReadRatio << " norDepth: " << tumDepth << " norMPQDepth: " << tumMpqDepth;
            }

                
            //current SNP GT type
            std::string GTtype = (*somaticVarIter).second.GTtype;

            if(GTtype == ""){
                std::cerr << "Error(GTtype) => chr: " << chr << " pos: " << (*somaticVarIter).first << " GTtype: " << GTtype;
                exit(1);                
            }

            //The count of deletions occurring at the SNP position
            int norDeletionCount = NorBase.getVarDeletionCount(chr, (*somaticVarIter).first);

            //int tumDeletionCount = TumBase.getVarDeletionCount(chr, (*somaticVarIter).first);
            int tumDeletionCount = SomaticChrPosInfo[chr][(*somaticVarIter).first].base.delCount;

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
    std::cerr<< difftime(time(NULL), begin) << "s\n";  
}

bool SomaticVarCaller::isHighConfidneceSomaticSNP(std::string chr, int pos){
    if(SomaticChrPosInfo[chr].find(pos) == SomaticChrPosInfo[chr].end()){
        return false;
    }    
    return SomaticChrPosInfo[chr][pos].isHighConSomaticSNP;
}

std::map<int, HP3_Info>* SomaticVarCaller::getSomaticPosInfo(std::string chr){
    return &(SomaticChrPosInfo[chr]);
}

std::map<std::string, std::map<int, HP3_Info>> SomaticVarCaller::getSomaticChrPosInfo(){
    return SomaticChrPosInfo;
}


void HaplotagProcess::variantParser(std::string variantFile, VCF_Info &Info){

    if( variantFile.find("gz") != std::string::npos ){
        // .vcf.gz 
        compressParser(variantFile, Info);
    }
    else if( variantFile.find("vcf") != std::string::npos ){
        // .vcf
        unCompressParser(variantFile, Info);
    }
    else{
        std::cerr<<"file: "<< variantFile << "\nnot vcf file. please check filename extension\n";
        exit(EXIT_FAILURE);
    }
    return;
}

void HaplotagProcess::compressParser(std::string &variantFile, VCF_Info &Info){
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
                parserProcess(input, Info);
            }
            // any trailing data in [eol, end) now is a partial line
            offset = std::copy(cur, end, buffer);
        }
        gzclose (file);
        free(buffer);
    }    
}

void HaplotagProcess::unCompressParser(std::string &variantFile, VCF_Info &Info){
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
            parserProcess(input, Info);
        }
    }
}

void HaplotagProcess::parserProcess(std::string &input, VCF_Info &Info){
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
        std::cerr<<"ERROR: Cannot open index for bam file " << openBamFile.c_str() << "\n";
    }

    // output file mangement
    std::string writeBamFile = params.resultPrefix + "." + params.outputFormat;
    // open output bam file
    samFile *out = hts_open(writeBamFile.c_str(), (params.outputFormat == "bam" ? "wb" : "wc" ));
    // load reference file
    hts_set_fai_filename(out, params.fastaFile.c_str() );
    // output writer
    int result = sam_hdr_write(out, bamHdr);
    // check index file
    if ((idx = sam_index_load(in, openBamFile.c_str())) == 0) {
        std::cerr<<"ERROR: Cannot open index for bam file\n";
        exit(1);
    }

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
                         << "HP2Allele\t"
                         << "HP3Allele\t"   //new
                         << "HP4Allele\t"   //new
                         << "phasingQuality(PQ)\t"
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
        currentChrVariants = mergedChrVarinat[chr];
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
                std::string haplotype;

                if(tagTumorMode){
                    haplotype = SomaticJudgeHaplotype(*bamHdr, *aln, chr, params.percentageThreshold, tagResult, pqValue, psValue, geneType, chr_reference);
                }else{
                    haplotype = judgeHaplotype(*bamHdr, *aln, chr, params.percentageThreshold, tagResult, pqValue, psValue, geneType, chr_reference);
                }

                if(haplotype.empty()){
                    std::cerr<< "haplotag string is empty: " << bam_get_qname(aln) << "\n";
                    exit(1);
                }

                initFlag(aln, "HP");
                initFlag(aln, "PS");
                initFlag(aln, "PQ");

                if (haplotype != "0"){

                    bam_aux_append(aln, "HP", 'Z', (haplotype.size() + 1), (const uint8_t*) haplotype.c_str());
                    if(psValue != -1) bam_aux_append(aln, "PS", 'i', sizeof(psValue), (uint8_t*) &psValue);
                    bam_aux_append(aln, "PQ", 'i', sizeof(pqValue), (uint8_t*) &pqValue);
                    totalTagCount++;

                    /*switch (haplotype){
                        case 1 : totalHP1++; break;
                        case 2 : totalHP2++; break;
                        case 3 : totalHP3++; break;
                        case 4 : totalHP4++; break;
                        default: break;
                    }*/
                    if(haplotype == "1"){
                        totalHP1++;
                    }else if (haplotype == "2"){
                        totalHP2++;
                    }else if (haplotype == "3"){
                        totalHP3++;
                    }else if (haplotype == "4"){
                        totalHP4++;
                    }else if(haplotype == "1-1"){
                        totalHP1_1++;
                    }else if(haplotype == "1-2"){
                        totalHP1_2++;
                    }else if(haplotype == "2-1"){
                        totalHP2_1++;
                    }else if(haplotype == "2-2"){
                        totalHP2_2++;
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
        hts_itr_destroy(iter);  
        std::cerr<< difftime(time(NULL), begin) << "s\n";
    }
    if(tagResult!=NULL){
        (*tagResult).close();
        delete tagResult;
        tagResult = nullptr;
    }
    if(tagTumorMode && params.writeReadLog){
        std::ofstream *readHpDistriLog=NULL;
        readHpDistriLog=new std::ofstream(params.resultPrefix+"_readDistri.out");

        int somaticSnpCount = 0;
        for(auto chr: *chrVec){
            if(!chrVarReadHpResult[chr].empty()){
                somaticSnpCount += chrVarReadHpResult[chr].size();
            }
        }

        if(!readHpDistriLog->is_open()){
            std::cerr<< "Fail to open write file: " << params.resultPrefix+"_readDistri.out" << "\n";
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
        for(auto chr: *chrVec){
            std::map<int, ReadHpResult>::iterator curVarReadHpIter = chrVarReadHpResult[chr].begin();
            while(curVarReadHpIter != chrVarReadHpResult[chr].end()){
                int pos = (*curVarReadHpIter).first + 1;
                int HP1readCount = (*curVarReadHpIter).second.HP1readCount;
                int HP1_1readCount = (*curVarReadHpIter).second.HP1_1readCount;
                int HP1_2readCount = (*curVarReadHpIter).second.HP1_2readCount;

                int HP2readCount = (*curVarReadHpIter).second.HP2readCount;
                int HP2_1readCount = (*curVarReadHpIter).second.HP2_1readCount;
                int HP2_2readCount = (*curVarReadHpIter).second.HP2_2readCount;

                int HP3readCount = (*curVarReadHpIter).second.HP3readCount;
                int HP4readCount = (*curVarReadHpIter).second.HP4readCount;

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
                                   << "somatic" << "\t\t"
                                   << HP1readCount << "\t"
                                   << HP1_1readCount << "\t"
                                   << HP1_2readCount << "\t\t"
                                   << HP2readCount << "\t"
                                   << HP2_1readCount << "\t"
                                   << HP2_2readCount << "\t\t"
                                   << HP3readCount << "\t"
                                   << HP4readCount << "\t"
                                   << (*curVarReadHpIter).second.unTagReadCount << "\t"
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

std::string HaplotagProcess::judgeHaplotype(const bam_hdr_t &bamHdr,const bam1_t &aln, std::string chrName, double percentageThreshold, std::ofstream *tagResult, int &pqValue, int &psValue, const int tagGeneType, std::string &ref_string){

    int hp1Count = 0;
    int hp2Count = 0;
    //record variants on this read
    std::map<int,std::string> variantsHP;

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
                                    variantsHP[curPos]="0";
                                }
                                if( base == vcfSet[Genome::NORMAL].chrVariantHP2[chrName][curPos]){
                                    hp2Count++;
                                    variantsHP[curPos]="1";
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
                                variantsHP[curPos]="0";
                            }
                            // hp2 occur insertion
                            else if( hp1Length == 1 && hp2Length != 1 ){
                                hp2Count++;
                                variantsHP[curPos]="1";
                            }
                        }
                        else {
                            // hp1 occur insertion
                            if( hp1Length != 1 && hp2Length == 1 ){
                                hp2Count++;
                                variantsHP[curPos]="1";
                            }
                            // hp2 occur insertion
                            else if( hp1Length == 1 && hp2Length != 1 ){
                                hp1Count++;
                                variantsHP[curPos]="0";
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
                                variantsHP[curPos]="0";
                            }
                            // hp2 occur deletion
                            else if( hp1Length == 1 && hp2Length != 1 ){
                                hp2Count++;
                                variantsHP[curPos]="1";
                            }
                        }
                        else {
                            // hp2 occur deletion
                            if( hp1Length != 1 && hp2Length == 1 ){
                                hp2Count++;
                                variantsHP[curPos]="1";
                            }
                            // hp1 occur deletion
                            else if( hp1Length == 1 && hp2Length != 1 ){
                                hp1Count++;
                                variantsHP[curPos]="0";
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
                                variantsHP[curPos]="0";
                            }
                            if( base == vcfSet[Genome::NORMAL].chrVariantHP2[chrName][curPos]){
                                hp2Count++;
                                variantsHP[curPos]="1";
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
                                variantsHP[curPos]="0";
                            }
                            // hp2 occur deletion
                            else if( hp1Length == 1 && hp2Length != 1 ){
                                hp2Count++;
                                variantsHP[curPos]="1";
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

    std::string hpResult = "0";
    if( max/(max+min) < percentageThreshold){
        // no tag
        pqValue = 0;
        totalHighSimilarity++;
    }
    else{
        if(hp1Count > hp2Count){
            hpResult = "1";
        }
        if(hp1Count < hp2Count){
            hpResult = "2";
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
        hpResult = "0";
    }

    if(tagResult!=NULL){
        //write tag log file
        std::string hpResultStr = ((hpResult == "0" )? "." : hpResult );
        std::string psResultStr = ".";

        if( hpResultStr != "." ){
            auto psIter = countPS.begin();
            psResultStr = std::to_string((*psIter).first);
            psValue = (*psIter).first;
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

std::string HaplotagProcess::SomaticJudgeHaplotype(const  bam_hdr_t &bamHdr,const bam1_t &aln, std::string chrName, double percentageThreshold, std::ofstream *tagResult, int &pqValue, int &psValue, const int tagGeneType, std::string &ref_string){

    std::map<int, int> hpCount;
    hpCount[0] = 0;
    hpCount[1] = 0;
    hpCount[3] = 0;
    hpCount[4] = 0;

    //record variants on this read
    std::map<int,std::string> variantsHP;

    std::vector<int> somaticVar;

    //record PS count(vcf type, PS value, count)
    std::map<int, int> tumCountPS;
    std::map<int, int> norCountPS;

    //std::map<int, std::map<int, int>> mergedCountPS;

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
                    SomaticJudgeSnpHP(currentVariantIter, vcfSet , chrName, base, hpCount, norCountPS, tumCountPS, &variantsHP, nullptr, nullptr, &(chrPosReadCase[chrName]));
                    if(chrPosReadCase[chrName].find((*currentVariantIter).first) != chrPosReadCase[chrName].end()){
                        if(chrPosReadCase[chrName][(*currentVariantIter).first].isHighConSomaticSNP){
                            somaticVar.push_back((*currentVariantIter).first);
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

    // determine the haplotype of the read
    std::string hpResult = "0";
    hpResult = determineReadHP(hpCount, pqValue, norCountPS, percentageThreshold, &totalHighSimilarity, &totalCrossTwoBlock, &totalWithOutVaraint);

    /*double normalMinHPcount = 0;
    double normalMaxHPcount = 0;
    std::string maxNormalHP = "0";

    double tumorMinHPcount = 0;
    double tumorMaxHPcount = 0;
    int maxTumorHP = 0;

    // determine max and min
    //tumor HP count
    if(hpCount[3] > hpCount[4]){
        tumorMinHPcount = hpCount[4];
        tumorMaxHPcount = hpCount[3];
        maxTumorHP = 3;
    }
    else{
        tumorMinHPcount = hpCount[3];
        tumorMaxHPcount = hpCount[4];
        maxTumorHP = 4;
    }

    //normal HP count
    if(hpCount[1] > hpCount[2]){
        normalMinHPcount = hpCount[2];
        normalMaxHPcount = hpCount[1];
        maxNormalHP = "1";
    }
    else{
        normalMinHPcount = hpCount[1];
        normalMaxHPcount = hpCount[2];
        maxNormalHP = "2";
    }

    //the similarity of HP types
    //tumor variants
    double tumorHPsimilarity = (tumorMaxHPcount == 0) ? 0.0 : tumorMaxHPcount/(tumorMaxHPcount+tumorMinHPcount);
    double normalHPsimilarity = (normalMaxHPcount == 0) ? 0.0 : normalMaxHPcount/(normalMaxHPcount+normalMinHPcount);


    // determine the haplotype of the read
    std::string hpResult = "0";

    // read have tumor variant
    if(tumorMaxHPcount != 0){
        //check the similarity of HP types in the tumor
        if(tumorHPsimilarity >= percentageThreshold){
            //check the similarity of HP types in the normal
            if(normalHPsimilarity >= percentageThreshold){
                switch(maxTumorHP){
                    case 3:
                        hpResult = maxNormalHP + "-1"; break;
                    case 4:
                        hpResult = maxNormalHP + "-2"; break;
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
                    case 3:
                        hpResult = "3"; break;
                    case 4:
                        hpResult = "4"; break;
                    default:
                        std::cerr << "ERROR: Unexpected haplotype : tumor Max HP= "<< maxTumorHP << std::endl;
                        exit(1);
                        break;
                }
            }
        }else{
            // no tag
            pqValue = 0;
            totalHighSimilarity++;
        }
    }
    // read haven't tumor variant
    else if(normalMaxHPcount != 0){
        if(normalHPsimilarity >= percentageThreshold){
            hpResult = maxNormalHP;
        }else{
            // no tag
            pqValue = 0;
            totalHighSimilarity++;
        }
    }

    // cross two block
    // had at least one tumor SNP in the current read
    if(hpCount[3] != 0 && hpCount[4] != 0 ){

    // all SNPs are germline variants in the current read
    }else{
        if(norCountPS.size() > 1){
            hpResult = "0";
            totalCrossTwoBlock++;
        }
    }


    // determine the quality of the read
    // There haven't been any variants in the current read
    if(normalMaxHPcount == 0 && tumorMaxHPcount == 0){
        //std::cerr<< "Read max = 0 : " << bam_get_qname(&aln) << "\n";
        totalWithOutVaraint++;
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
    }*/

    //std::cout << "---------------------------------------" << std::endl;
    //std::cout << "HP1: "<< hpCount[1]<< " HP2: "<< hpCount[2]<< " HP3: "<< hpCount[3]<< " HP4: "<< hpCount[4] << std::endl;
    //std::cout << "TmaxC: "<< tumorMaxHPcount<< " TminC: "<< tumorMinHPcount<< " NmaxC: "<< normalMaxHPcount<< " NminC: "<< normalMinHPcount << std::endl;
    //std::cout << "TmaxHP: "<< maxTumorHP<< " NmaxHP: "<< maxNormalHP << std::endl;
    //std::cout << "Tsimilar: "<< tumorHPsimilarity<< " Nsimilar: "<< normalHPsimilarity << std::endl;
    //std::cout << "TpsCountSize: "<< tumCountPS.size()<< " NpsCountSize: "<< norCountPS.size() << std::endl;
    //std::cout << "hpResult: "<< hpResult << std::endl;


    std::string hpResultStr = ((hpResult == "0") ? "." : hpResult );
    std::string psResultStr = ".";
    
    std::pair<std::string, int> psLog = std::make_pair(".", 0);

    //determine the PS value of the read
    if( hpResultStr != "." ){
        std::map<int, int>::iterator psIter;
        if(hpResult != "1" && hpResult != "2"){
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
        }else if(hpResult == "1" || hpResult == "2"){
            psIter = norCountPS.begin();
            psResultStr = std::to_string((*psIter).first);
            psValue = (*psIter).first;
        }
    }

    //write tag log file
    if(tagResult!=NULL){

        (*tagResult)<< bam_get_qname(&aln)                          << "\t"
                    << bamHdr.target_name[aln.core.tid]             << "\t"
                    << aln.core.pos                                 << "\t"
                    //<< normalMaxHPcount/(normalMaxHPcount+normalMinHPcount)<< "\t"
                    << hpResultStr                                  << "\t"
                    << psResultStr                                  << "\t"
                    << hpCount[1]+hpCount[2]+hpCount[3]+hpCount[4]  << "\t" //modify
                    << hpCount[1]                                   << "\t"
                    << hpCount[2]                                   << "\t"
                    << hpCount[3]                                   << "\t" //new
                    << hpCount[4]                                   << "\t" //new
                    << pqValue                                      << "\t";


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

        (*tagResult)<< " ";

        (*tagResult)<< "TumPS:";
        for(auto v : tumCountPS){
            (*tagResult)<< " " << v.first << "," << v.second ;
        }

        (*tagResult)<< "\n";
    }

    //Record read HP result
    if(!somaticVar.empty()){
        //std::cout << "Size of recordPos: " << somaticVar.size() << std::endl;
        for(int pos : somaticVar){
            /*switch (hpResult){
                case 0:
                    chrVarReadHpResult[chrName][pos].unTagReadCount++;
                    break;
                case 1:
                    chrVarReadHpResult[chrName][pos].HP1readCount++;
                    break;
                case 2:
                    chrVarReadHpResult[chrName][pos].HP2readCount++;
                    break;
                case 3:
                    chrVarReadHpResult[chrName][pos].HP3readCount++;
                    break;
                case 4:
                    chrVarReadHpResult[chrName][pos].HP4readCount++;
                    break;
                default:
                    std::cerr<< "Error: Unexpected read HP result" << hpResult << "\n";
                    exit(1);
                    break;
            }*/
            if(hpResult == "0"){
                chrVarReadHpResult[chrName][pos].unTagReadCount++;
            }else if (hpResult == "1"){
                chrVarReadHpResult[chrName][pos].HP1readCount++;
            }else if (hpResult == "2"){
                chrVarReadHpResult[chrName][pos].HP2readCount++;
            }else if (hpResult == "3"){
                chrVarReadHpResult[chrName][pos].HP3readCount++;
            }else if (hpResult == "4"){
                chrVarReadHpResult[chrName][pos].HP4readCount++;
            }else if (hpResult == "1-1"){
                chrVarReadHpResult[chrName][pos].HP1_1readCount++;
            }else if (hpResult == "1-2"){
                chrVarReadHpResult[chrName][pos].HP1_2readCount++;
            }else if (hpResult == "2-1"){
                chrVarReadHpResult[chrName][pos].HP2_1readCount++;
            }else if (hpResult == "2-2"){
                chrVarReadHpResult[chrName][pos].HP2_2readCount++;
            }else{
                std::cerr << "Error: Unexpected read HP result " << hpResult << "\n";
                exit(1);
            }
        }
    }
    return hpResult;
}

void HaplotagProcess::OnlyTumorSNPjudgeHP(RefAltSet &curVar, int &curPos, BamBaseCounter *NorBase, VCF_Info *vcfSet, std::string chrName, std::string base, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int,std::string> *variantsHP, std::vector<int> *readPosHP3, std::map<int, HP3_Info> *SomaticPos){

    if(SomaticPos == nullptr){
        std::cout << "ERROR (SomaticTaggingJudgeHP) => SomaticPos pointer cannot be nullptr"<< std::endl;
        exit(1);
    }

    if((*SomaticPos).find(curPos) != (*SomaticPos).end()){
        if((*SomaticPos)[curPos].isHighConSomaticSNP){
            std::string TumorRefBase = curVar.Variant[Genome::TUMOR].Ref;
            std::string TumorAltBase = curVar.Variant[Genome::TUMOR].Alt; 

            if(base == TumorAltBase){
                hpCount[3]++;
                if(variantsHP != nullptr) (*variantsHP)[curPos] = '3';

            //base is not match to TumorRefBase & TumorAltBase (other HP)
            }else if(base != TumorRefBase && base != TumorAltBase){
                hpCount[4]++;
                if(variantsHP != nullptr) (*variantsHP)[curPos] = '4';
            }

            if(curVar.Variant[Genome::TUMOR].is_phased_hetero){
                if(tumCountPS != nullptr) (*tumCountPS)[vcfSet[Genome::TUMOR].chrVariantPS[chrName][curPos]]++;
            }
        }
    }
}

HaplotagProcess::HaplotagProcess(HaplotagParameters params):
totalAlignment(0),totalSupplementary(0),totalSecondary(0),totalUnmapped(0),totalTagCount(0),totalUnTagCount(0),processBegin(time(NULL)),integerPS(false)
{
    std::cerr<< "phased SNP file:       " << params.snpFile             << "\n";
    if(params.tagTumorSnp) 
    std::cerr<< "phased tumor SNP file: " << params.tumorSnpFile        << "\n";  //new
    std::cerr<< "phased SV file:        " << params.svFile              << "\n";
    std::cerr<< "phased MOD file:       " << params.modFile             << "\n";
    std::cerr<< "input bam file:        " << params.bamFile             << "\n";
    if(params.tagTumorSnp) 
    std::cerr<< "input tumor bam file:  " << params.tumorBamFile        << "\n";  //new
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

    //------------------------verification parameter----------------------------
    tagTumorMode=params.tagTumorSnp;
    vcfSet[Genome::NORMAL].gene_type = Genome::NORMAL;
    vcfSet[Genome::TUMOR].gene_type = Genome::TUMOR;

    parseSVFile = false;
    parseSnpFile = false;
    parseMODFile = false;

    totalLowerQuality = 0;
    totalOtherCase = 0;
    totalHP1 = 0;
    totalHP2 = 0;
    totalHP3 = 0;
    totalHP4 = 0;
    totalHP1_1 = 0;
    totalHP1_2 = 0;
    totalHP2_1 = 0;
    totalHP2_2 = 0;
    totalunTag_HP0 = 0;
    totalHighSimilarity = 0;
    totalCrossTwoBlock = 0;
    totalEmptyVariant = 0;
    totalWithOutVaraint = 0;
    
    // decide on the type of tagging for VCF and BAM files
    int tagGeneType;

    if(tagTumorMode){
        tagGeneType = Genome::TUMOR;
    }else{
        tagGeneType = Genome::NORMAL;
    }

    //---------------------------------------------------------------

    if(tagTumorMode){
        //load tumor snp vcf
        if(params.tumorSnpFile != ""){
            std::time_t begin = time(NULL);
            std::cerr<< "parsing tumor SNP VCF ... ";
            parseSnpFile = true;
            variantParser(params.tumorSnpFile, vcfSet[Genome::TUMOR]);
            parseSnpFile = false;
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

    parseSnpFile = true;
    variantParser(params.snpFile, vcfSet[Genome::NORMAL]);
    parseSnpFile = false;
    std::cerr<< difftime(time(NULL), begin) << "s\n";

    // load SV vcf file
    if(params.svFile!=""){
        begin = time(NULL);
        std::cerr<< "parsing SV VCF ... ";
        parseSVFile = true;
        variantParser(params.svFile, vcfSet[tagGeneType]);
        parseSVFile = false;
        std::cerr<< difftime(time(NULL), begin) << "s\n";    
    }
    
    // load MOD vcf file
    if(params.modFile!=""){
        begin = time(NULL);
        std::cerr<< "parsing MOD VCF ... ";
        parseMODFile = true;
        variantParser(params.modFile, vcfSet[tagGeneType]);
        parseMODFile = false;
        std::cerr<< difftime(time(NULL), begin) << "s\n";    
    }

    //decide which genome type chrVec and chrLength belong to
    if(tagGeneType == Genome::NORMAL){
        chrVec = &(vcfSet[Genome::NORMAL].chrVec);
        chrLength = &(vcfSet[Genome::NORMAL].chrLength); 
    }else{
        //check normal & tumor chr & length
        if(vcfSet[Genome::NORMAL].chrVec.size() != vcfSet[Genome::TUMOR].chrVec.size()){
            std::cout << "tumor & normal VCFs chromosome count are not the same" << std::endl;
            return ;
        }

        if(vcfSet[Genome::NORMAL].chrLength.size() != vcfSet[Genome::TUMOR].chrLength.size()){
            std::cout << "tumor & normal VCFs  chrLength count not same" << std::endl;
            return ;
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


    if(tagTumorMode){

        //Count each base numbers at tumor SNP position in the Normal.bam
        BamBaseCounter *NorBase = new BamBaseCounter(params.bamFile, params, mergedChrVarinat, *chrVec, *chrLength, Genome::TUMOR);

        //record the HP3 confidence of each read
        SomaticVarCaller *SomaticVar = new SomaticVarCaller(params.tumorBamFile, mergedChrVarinat, *chrVec, *chrLength, params, vcfSet, *NorBase);
        chrPosReadCase = SomaticVar->getSomaticChrPosInfo();

        delete NorBase;
        delete SomaticVar;

        NorBase = nullptr;
        SomaticVar = nullptr;
        //exit(1);
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

HaplotagProcess::~HaplotagProcess(){
    std::cerr<< "-------------------------------------------\n";
    std::cerr<< "total process time:    " << difftime(time(NULL), processBegin) << "s\n";
    std::cerr<< "total alignment:       " << totalAlignment     << "\n";
    std::cerr<< "total supplementary:   " << totalSupplementary << "\n";
    std::cerr<< "total secondary:       " << totalSecondary     << "\n";
    std::cerr<< "total unmapped:        " << totalUnmapped      << "\n";
    std::cerr<< "total tag alignment:   " << totalTagCount     << "\n";
    std::cerr<< "    L----total HP1   : " << totalHP1     << "\n";   //new
    std::cerr<< "    L----total HP1-1 : " << totalHP1_1   << "\n";   //new
    std::cerr<< "    L----total HP1-2 : " << totalHP1_2   << "\n";   //new
    std::cerr<< "    L----total HP2   : " << totalHP2     << "\n";   //new
    std::cerr<< "    L----total HP2-1 : " << totalHP2_1   << "\n";   //new
    std::cerr<< "    L----total HP2-2 : " << totalHP2_2   << "\n";   //new
    std::cerr<< "    L----total HP3   : " << totalHP3     << "\n";   //new
    std::cerr<< "    L----total HP4   : " << totalHP4     << "\n";   //new
    std::cerr<< "total untagged:        " << totalUnTagCount   << "\n";
    std::cerr<< "    L----total lower mapping quality:    " << totalLowerQuality   << "\n";   //new
    std::cerr<< "    L----total EmptyVariant:             " << totalEmptyVariant   << "\n";   //new
    std::cerr<< "    L----total start > last variant pos: " << totalOtherCase   << "\n";   //new
    std::cerr<< "    L----total untag:                    " << totalunTag_HP0   << "\n";   //new
    std::cerr<< "         L----total HighSimilarity:      " << totalHighSimilarity   << "\n";   //new
    std::cerr<< "         L----total CrossTwoBlock:       " << totalCrossTwoBlock   << "\n";   //new
    std::cerr<< "         L----total WithOut Variant:     " << totalWithOutVaraint   << "\n";   //new
    std::cerr<< "-------------------------------------------\n";
};

