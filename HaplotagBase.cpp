#include "HaplotagBase.h"

BamBaseCounter::BamBaseCounter(bool enableFilter){
    ChrVariantBase = new std::map<std::string, std::map<int, PosBase>>();
    applyFilter = enableFilter;
};

BamBaseCounter::~BamBaseCounter(){
    delete ChrVariantBase;
};

void BamBaseCounter::CountingBamBase(
    const std::string &BamFile, 
    const HaplotagParameters &params, 
    std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat, 
    std::vector<std::string> &chrVec, 
    std::map<std::string, int> &chrLength, 
    VCF_Info *vcfSet, 
    int genmoeType
){
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
    
    std::vector<int> last_pos;
    // get the last variant position of the reference
    germlineGetRefLastVarPos(last_pos, chrVec, vcfSet, genmoeType);
    // reference fasta parser
    FastaParser fastaParser(params.fastaFile, chrVec, last_pos, params.numThreads);

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
        std::map<int, MultiGenomeVar> currentVariants;

        #pragma omp critical
        {
            variantBase = &((*ChrVariantBase)[chr]);
            currentVariants = mergedChrVarinat[chr];
        }
        // fetch chromosome string
        std::string ref_string = fastaParser.chrString.at(chr);

        //inintial iterator
        std::map<int, MultiGenomeVar>::iterator firstVariantIter = currentVariants.begin();

        std::map<int, MultiGenomeVar>::reverse_iterator last = currentVariants.rbegin();

        std::string region = !params.region.empty() ? params.region : chr + ":1-" + std::to_string(chrLength[chr]);
        hts_itr_t* iter = sam_itr_querys(bam.idx, bam.bamHdr, region.c_str());

        while (sam_itr_multi_next(bam.in, iter, bam.aln) >= 0) {
            
            int flag = bam.aln->core.flag;

            // if ( bam.aln->core.qual < params.qualityThreshold ){
            //    // mapping quality is lower than threshold
            //    continue;
            // }

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
                StatisticBaseInfo(*bam.aln, *bam.bamHdr, chr, params, genmoeType, *variantBase, currentVariants, firstVariantIter, vcfSet, ref_string);
            }
        }
        // std::cerr<<"calculate Base Information ...\n";
        CalculateBaseInfo(chr, *variantBase, currentVariants);
        
        variantBase = nullptr;
        hts_itr_destroy(iter);
    }
    // if(tagResult != nullptr){
    //     (*tagResult).close();
    // }
    hts_tpool_destroy(threadPool.pool);
    std::cerr<< difftime(time(NULL), begin) << "s\n";
    return;
}

void BamBaseCounter::StatisticBaseInfo(
    const bam1_t &aln, 
    const bam_hdr_t &bamHdr,
    const std::string &chrName, 
    const HaplotagParameters &params, 
    int genmoeType, 
    std::map<int, PosBase> &variantBase, 
    std::map<int, MultiGenomeVar> &currentVariants,
    std::map<int, MultiGenomeVar>::iterator &firstVariantIter, 
    VCF_Info *vcfSet, 
    const std::string &ref_string
){

    int hp1Count = 0;
    int hp2Count = 0;
    //record variants on this read
    std::map<int, int> variantsHP;

    std::map<int,int> countPS;

    std::vector<int> tumVarPosVec;

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
    std::map<int, MultiGenomeVar>::iterator currentVariantIter = firstVariantIter;

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
                    // printf("into match cigar\n");

                    //statistically analyze SNP information exclusive to the tumor
                    if((*currentVariantIter).second.isExists(TUMOR)){
                        int curPos = (*currentVariantIter).first;
                        //std::cout << "curPos :" << curPos << "is tumor "<<(*currentVariantIter).second.isExistTumor  << " isNormal: "<< (*currentVariantIter).second.isExistNormal<< std::endl;

                        //detect ref bsae length(temp :tumor SNP)
                        int tumRefLength = (*currentVariantIter).second.Variant[TUMOR].allele.Ref.length();
                        int tumAltLength = (*currentVariantIter).second.Variant[TUMOR].allele.Alt.length();
                        
                        // the variant is SNP
                        if(tumRefLength == 1 && tumAltLength == 1){
                            //record tumor SNP position
                            tumVarPosVec.push_back(curPos);
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

                        }

                        // the indel(del) SNP position is the start position, and the deletion occurs at the next position
                        else if(tumRefLength > 1 && tumAltLength == 1){
                            // the indel SNP start position is at the end of the deletion, and the next cigar operator is deletion
                            if(curPos == (ref_pos + length - 1) && bam_cigar_op(cigar[i+1]) == 2 && i+1 < aln_core_n_cigar){

                            }
                        }

                    }       

                    if ( aln.core.qual >= params.somaticCallingMpqThreshold && (*currentVariantIter).second.isExists(NORMAL)){
                        // only judge the heterozygous SNP
                        if((*currentVariantIter).second.Variant[NORMAL].is_phased_hetero){
                            auto norVar = (*currentVariantIter).second.Variant[NORMAL].allele;
                            germlineJudgeSnpHap(chrName, vcfSet, norVar, base, ref_pos, length, i, aln_core_n_cigar, cigar, currentVariantIter, hp1Count, hp2Count, variantsHP, countPS);
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
            bool alreadyJudgeDel = false;
            while( currentVariantIter != currentVariants.end() && (*currentVariantIter).first < ref_pos + length){
                //statistically analyze SNP information exclusive to the tumor
                if((*currentVariantIter).second.isExists(TUMOR)){

                    int curPos = (*currentVariantIter).first;
                    
                    //record tumor SNP position
                    tumVarPosVec.push_back(curPos);

                    //detect ref bsae length(temp :tumor SNP)
                    int tumRefLength = (*currentVariantIter).second.Variant[TUMOR].allele.Ref.length();
                    int tumAltLength = (*currentVariantIter).second.Variant[TUMOR].allele.Alt.length();

                    // the variant is SNP
                    if(tumRefLength == 1 && tumAltLength == 1){
                        variantBase[curPos].delCount++;
                        variantBase[curPos].depth++;
                    }
                    // the variant is deletion
                    else if(tumRefLength > 1 && tumAltLength == 1){

                    }
                }

                // only execute at the first phased normal snp
                if ( aln.core.qual >= params.somaticCallingMpqThreshold && (*currentVariantIter).second.isExists(NORMAL) && !alreadyJudgeDel){
                    if((*currentVariantIter).second.Variant[NORMAL].is_phased_hetero){
                        // longphase v1.73 only execute once
                        alreadyJudgeDel = true;
                        auto norVar = (*currentVariantIter).second.Variant[NORMAL];
                        germlineJudgeDeletionHap(chrName, ref_string, ref_pos, length, query_pos, currentVariantIter, vcfSet, &aln, hp1Count, hp2Count, variantsHP, countPS);
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

    // printf("end cigar operation\n");

    // get the number of SVs occurring on different haplotypes in a read
    if( aln.core.qual >= params.somaticCallingMpqThreshold ){
        germlineJudgeSVHap(aln, vcfSet, hp1Count, hp2Count, genmoeType);
    }

    double min = 0.0;
    double max = 0.0;
    int hpResult = ReadHP::unTag;
    int pqValue = 0;
    int psValue = 0;
    double percentageThreshold = params.percentageThreshold;
    // determine the haplotype of the read
    hpResult = germlineDetermineReadHap(hp1Count, hp2Count, min, max, percentageThreshold, pqValue, psValue, countPS, nullptr, nullptr);

    //record read hp to tumor SNP position
    for(auto pos : tumVarPosVec){
        variantBase[pos].ReadHpCount[hpResult]++;
    }

}

void BamBaseCounter::CalculateBaseInfo(const std::string &chr, std::map<int, PosBase> &VariantBase, std::map<int, MultiGenomeVar> &currentVariants){
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
        std::string tumRefBase = currentVariants[(*currentPosIter).first].Variant[TUMOR].allele.Ref;
        std::string tumAltBase = currentVariants[(*currentPosIter).first].Variant[TUMOR].allele.Alt;

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
            baseInfo->nonDelAF = (float)AltCount / (float)(depth - baseInfo->delCount);
        }

        int filteredMpqDepth = baseInfo->filteredMpqDepth;
        
        if(filteredMpqAltCount != 0 && filteredMpqDepth != 0){
            baseInfo->filteredMpqVAF = (float)filteredMpqAltCount / (float)filteredMpqDepth;
        }

        if(depth != 0){
            baseInfo->lowMpqReadRatio = (float)(depth - filteredMpqDepth) / (float)depth;
        }
        
        float refAllele_threshold = 0.90; //(VAF <= 0.1)      
        if( !applyFilter ||
            //(baseInfo->max_ratio >= refAllele_threshold && depth > 1 && baseInfo->lowMpqReadRatio <= 0.1) ){
            // (baseInfo->VAF <= 0.1 && depth > 1 && baseInfo->lowMpqReadRatio <= 0.1) ){
            (baseInfo->VAF <= 0.1 && depth > 1) ){
            baseInfo->isHighRefAllelleFreq = true;
        }
        baseInfo->isHighRefAllelleFreq = true;
        
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

float BamBaseCounter::getNoDelAF(std::string chr, int pos){
    if((*ChrVariantBase)[chr].find(pos) == (*ChrVariantBase)[chr].end()){
        std::cerr << "ERROR (getNoDelAF) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return (*ChrVariantBase)[chr][pos].nonDelAF;
}

float BamBaseCounter::getFilterdMpqVAF(std::string chr, int pos){
    if((*ChrVariantBase)[chr].find(pos) == (*ChrVariantBase)[chr].end()){
        std::cerr << "ERROR (getFilterdMpqVAF) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return (*ChrVariantBase)[chr][pos].filteredMpqVAF;
}

int BamBaseCounter::getReadHpCountInNorBam(std::string chr, int pos, int Haplotype){
    if((*ChrVariantBase)[chr].find(pos) == (*ChrVariantBase)[chr].end()){
        std::cerr << "ERROR (getReadHpCountInNorBam) => can't find the position:" << " chr: " << chr << " pos: " << pos << std::endl;
        exit(1);
    }
    return (*ChrVariantBase)[chr][pos].ReadHpCount[Haplotype];
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

void GermlineJudgeBase::germlineGetRefLastVarPos(
    std::vector<int>& last_pos, 
    const std::vector<std::string>& chrVec, 
    VCF_Info* vcfSet, 
    int geneType
){
    for( auto chr : chrVec ){
        auto lastVariantIter = vcfSet[geneType].chrVariantPS[chr].rbegin();
        if( lastVariantIter != vcfSet[geneType].chrVariantPS[chr].rend() ){
            last_pos.push_back(lastVariantIter->second);
        }
        else{
            last_pos.push_back(0);
        }
    }
}


void GermlineJudgeBase::germlineJudgeSnpHap(
    const std::string& chrName,
    VCF_Info* vcfSet,
    RefAlt& norVar,
    const std::string& base,
    int& ref_pos,
    int& length,
    int& i,
    int& aln_core_n_cigar,
    uint32_t* cigar,
    std::map<int, MultiGenomeVar>::iterator &currentVariantIter,
    int& hp1Count,
    int& hp2Count,
    std::map<int, int>& variantsHP,
    std::map<int, int>& countPS
){
    int curPos = (*currentVariantIter).first;
    int refAlleleLen = norVar.Ref.length();
    int altAlleleLen = norVar.Alt.length();

    // currentVariant is SNP
    if( refAlleleLen == 1 && altAlleleLen == 1 ){
        // Detected that the base of the read is either REF or ALT. 
        if( (base == norVar.Ref) || (base == norVar.Alt) ){


            std::map<int, int>::iterator posPSiter = vcfSet[NORMAL].chrVariantPS[chrName].find((*currentVariantIter).first);

            if( posPSiter == vcfSet[NORMAL].chrVariantPS[chrName].end() ){
                std::cerr << "ERROR (germlineJudgeSnpHap) => can't find the position:" 
                          << " chr: " << chrName << "\t"
                          << " pos: " << curPos << "\t"
                          << " ref: " << norVar.Ref << "\t"
                          << " alt: " << norVar.Alt << "\n";
                exit(EXIT_SUCCESS);
            }
            else{
                if( base == vcfSet[NORMAL].chrVariantHP1[chrName][curPos]){
                    hp1Count++;
                    variantsHP[curPos]=0;
                }
                if( base == vcfSet[NORMAL].chrVariantHP2[chrName][curPos]){
                    hp2Count++;
                    variantsHP[curPos]=1;
                }
                countPS[vcfSet[NORMAL].chrVariantPS[chrName][curPos]]++;
            }
            
        }
    }
    // currentVariant is insertion
    else if( refAlleleLen == 1 && altAlleleLen != 1 && i+1 < aln_core_n_cigar){
        
        int hp1Length = vcfSet[NORMAL].chrVariantHP1[chrName][curPos].length();
        int hp2Length = vcfSet[NORMAL].chrVariantHP2[chrName][curPos].length();
        
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
        countPS[vcfSet[NORMAL].chrVariantPS[chrName][curPos]]++;
    } 
    // currentVariant is deletion
    else if( refAlleleLen != 1 && altAlleleLen == 1 && i+1 < aln_core_n_cigar) {

        int hp1Length = vcfSet[NORMAL].chrVariantHP1[chrName][curPos].length();
        int hp2Length = vcfSet[NORMAL].chrVariantHP2[chrName][curPos].length();
        
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
        countPS[vcfSet[NORMAL].chrVariantPS[chrName][curPos]]++;
    } 
}


void GermlineJudgeBase::germlineJudgeDeletionHap(
    const std::string& chrName,
    const std::string& ref_string,
    int& ref_pos,
    int& length,
    int& query_pos,
    std::map<int, MultiGenomeVar>::iterator &currentVariantIter,
    VCF_Info* vcfSet,
    const bam1_t* aln,
    int& hp1Count,
    int& hp2Count,
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
                int refAlleleLen = norVar.allele.Ref.length();
                int altAlleleLen = norVar.allele.Alt.length();
                
                // SNP
                if (refAlleleLen == 1 && altAlleleLen == 1) {
                    // get the next match
                    char base_chr = seq_nt16_str[bam_seqi(bam_get_seq(aln), query_pos)];
                    std::string base(1, base_chr);

                    if (base == vcfSet[NORMAL].chrVariantHP1[chrName][curPos]) {
                        hp1Count++;
                        variantsHP[curPos] = 0;
                    }
                    if (base == vcfSet[NORMAL].chrVariantHP2[chrName][curPos]) {
                        hp2Count++;
                        variantsHP[curPos] = 1;
                    }
                    countPS[vcfSet[NORMAL].chrVariantPS[chrName][curPos]]++;
                }
                
                // the read deletion contain VCF's deletion
                else if (refAlleleLen != 1 && altAlleleLen == 1) {

                    int hp1Length = vcfSet[NORMAL].chrVariantHP1[chrName][curPos].length();
                    int hp2Length = vcfSet[NORMAL].chrVariantHP2[chrName][curPos].length();
                    // hp1 occur deletion
                    if (hp1Length != 1 && hp2Length == 1) {
                        hp1Count++;
                        variantsHP[curPos] = 0;
                    }
                    // hp2 occur deletion
                    else if (hp1Length == 1 && hp2Length != 1) {
                        hp2Count++;
                        variantsHP[curPos] = 1;
                    }
                    countPS[vcfSet[NORMAL].chrVariantPS[chrName][curPos]]++;
                }
            }
        }
    }
}

void GermlineJudgeBase::germlineJudgeSVHap(const bam1_t &aln, VCF_Info* vcfSet, int& hp1Count, int& hp2Count, const int& tagGeneType){
    auto readIter = vcfSet[tagGeneType].readSVHapCount.find(bam_get_qname(&aln));
    if( readIter != vcfSet[tagGeneType].readSVHapCount.end() ){
        hp1Count += vcfSet[tagGeneType].readSVHapCount[bam_get_qname(&aln)][0];
        hp2Count += vcfSet[tagGeneType].readSVHapCount[bam_get_qname(&aln)][1];
    }
}

int GermlineJudgeBase::germlineDetermineReadHap(
    int& hp1Count, 
    int& hp2Count, 
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

    if(hp1Count > hp2Count){
        min = hp2Count;
        max = hp1Count;
    }
    else{
        min = hp1Count;
        max = hp2Count;
    }

    if( max/(max+min) < percentageThreshold){
        // no tag
        pqValue = 0;
        if(totalHighSimilarity != nullptr) (*totalHighSimilarity)++;
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

void GermlineJudgeBase::writeGermlineTagLog(std::ofstream& tagResult, const bam1_t& aln, const bam_hdr_t& bamHdr, int& hpResult, double& max, double& min, int& hp1Count, int& hp2Count, int& pqValue, const std::map<int, int>& variantsHP, const std::map<int, int>& countPS){
    // write tag log file
    std::string hpResultStr = ((hpResult == ReadHP::unTag) ? "." : std::to_string(hpResult));
    std::string psResultStr = ".";

    if (hpResultStr != ".") {
        auto psIter = countPS.begin();
        psResultStr = std::to_string(psIter->first);
    }

    (tagResult) << bam_get_qname(&aln) << "\t"
                << bamHdr.target_name[aln.core.tid] << "\t"
                << aln.core.pos << "\t"
                << max / (max + min) << "\t"
                << "H" << hpResultStr << "\t"
                << psResultStr << "\t"
                << hp1Count + hp2Count << "\t"
                << hp1Count << "\t"
                << hp2Count << "\t"
                << pqValue << "\t";

    // print position and HP
    for (const auto& v : variantsHP) {
        (tagResult) << " " << v.first << "," << v.second;
    }

    (tagResult) << "\t";

    // belong PS, number of variant
    for (const auto& v : countPS) {
        (tagResult) << " " << v.first << "," << v.second;
    }

    (tagResult) << "\n";
}

void SomaticJudgeBase::SomaticJudgeSnpHP(std::map<int, MultiGenomeVar>::iterator &currentVariantIter, VCF_Info *vcfSet, std::string chrName, std::string base, std::map<int, int> &hpCount, std::map<int, int> &norCountPS, std::map<int, int> &tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *tumorAllelePosVec, BamBaseCounter *NorBase, std::map<int, HP3_Info> *SomaticPos){
    int curPos = (*currentVariantIter).first;
    auto curVar = (*currentVariantIter).second;

    // normal & tumor SNP at the current position (base on normal phased SNPs)
    // both normal and tumor samples that do not exist in the high-confidence set
    if(curVar.isExists(NORMAL) && curVar.isExists(TUMOR)){

        // the tumor & normal SNP GT are phased heterozygous 
        if((curVar.Variant[NORMAL].is_phased_hetero) && (curVar.Variant[TUMOR].is_phased_hetero)){   
            if(curVar.Variant[NORMAL].allele.Ref == base || curVar.Variant[NORMAL].allele.Alt == base){
                //std::cerr<< "tag tumor normal SNP\n";
                std::map<int, int>::iterator NorPosPSiter = vcfSet[NORMAL].chrVariantPS[chrName].find(curPos);
                if( NorPosPSiter == vcfSet[NORMAL].chrVariantPS[chrName].end()){
                    std::cerr<< "Unable to locate the phase set of the current normal SNP\n"
                             << curPos << "\t"
                             << curVar.Variant[NORMAL].allele.Ref << "\t"
                             << curVar.Variant[NORMAL].allele.Alt  << "\n";
                    exit(EXIT_SUCCESS);
                }

                std::string norHP1 = vcfSet[NORMAL].chrVariantHP1[chrName][curPos];
                std::string norHP2 = vcfSet[NORMAL].chrVariantHP2[chrName][curPos];
                std::string tumHP1 = vcfSet[TUMOR].chrVariantHP1[chrName][curPos];
                std::string tumHP2 = vcfSet[TUMOR].chrVariantHP2[chrName][curPos];

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
                norCountPS[vcfSet[NORMAL].chrVariantPS[chrName][curPos]]++;
            }
        //the normal SNP GT is phased heterozgous & the tumor SNP GT is unphased heterozgous 
        }else if((curVar.Variant[NORMAL].is_phased_hetero) && (curVar.Variant[TUMOR].is_unphased_hetero)){   
            if(curVar.Variant[NORMAL].allele.Ref == base || curVar.Variant[NORMAL].allele.Alt == base){

                std::map<int, int>::iterator NorPosPSiter = vcfSet[NORMAL].chrVariantPS[chrName].find(curPos);
                if( NorPosPSiter == vcfSet[NORMAL].chrVariantPS[chrName].end()){
                    std::cerr<< "Unable to locate the phase set of the current normal SNP\n"
                             << curPos << "\t"
                             << curVar.Variant[NORMAL].allele.Ref << "\t"
                             << curVar.Variant[NORMAL].allele.Alt  << "\n";
                    exit(EXIT_SUCCESS);
                }

                std::string norHP1 = vcfSet[NORMAL].chrVariantHP1[chrName][curPos];
                std::string norHP2 = vcfSet[NORMAL].chrVariantHP2[chrName][curPos];

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
                norCountPS[vcfSet[NORMAL].chrVariantPS[chrName][curPos]]++;
            }

        //the normal SNP GT is phased heterozgous & the tumor SNP GT is homozygous 
        }else if((curVar.Variant[NORMAL].is_phased_hetero) && (curVar.Variant[TUMOR].is_homozygous)){   
            if(curVar.Variant[NORMAL].allele.Ref == base || curVar.Variant[NORMAL].allele.Alt == base){

                std::map<int, int>::iterator NorPosPSiter = vcfSet[NORMAL].chrVariantPS[chrName].find(curPos);
                if( NorPosPSiter == vcfSet[NORMAL].chrVariantPS[chrName].end()){
                    std::cerr<< "Unable to locate the phase set of the current normal SNP\n"
                             << curPos << "\t"
                             << curVar.Variant[NORMAL].allele.Ref << "\t"
                             << curVar.Variant[NORMAL].allele.Alt  << "\n";
                    exit(EXIT_SUCCESS);
                }

                std::string norHP1 = vcfSet[NORMAL].chrVariantHP1[chrName][curPos];
                std::string norHP2 = vcfSet[NORMAL].chrVariantHP2[chrName][curPos];

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
                norCountPS[vcfSet[NORMAL].chrVariantPS[chrName][curPos]]++;
            }
        }
    // only normal SNP at the current position
    }else if(curVar.isExists(NORMAL)){
        // the normal SNP GT is phased heterozgous SNP
        if((curVar.Variant[NORMAL].is_phased_hetero)){
            if(curVar.Variant[NORMAL].allele.Ref == base || curVar.Variant[NORMAL].allele.Alt == base){

                std::map<int, int>::iterator NorPosPSiter = vcfSet[NORMAL].chrVariantPS[chrName].find(curPos);
                if( NorPosPSiter == vcfSet[NORMAL].chrVariantPS[chrName].end()){
                    std::cerr<< "Unable to locate the phase set of the current normal SNP\n"
                             << curPos << "\t"
                             << curVar.Variant[NORMAL].allele.Ref << "\t"
                             << curVar.Variant[NORMAL].allele.Alt  << "\n";
                    exit(EXIT_SUCCESS);
                }

                std::string norHP1 = vcfSet[NORMAL].chrVariantHP1[chrName][curPos];
                std::string norHP2 = vcfSet[NORMAL].chrVariantHP2[chrName][curPos];

                if( base == norHP1){
                    hpCount[1]++;
                    if(variantsHP != nullptr) (*variantsHP)[curPos] = SnpHP::GERMLINE_H1;
                }
                if(base == norHP2){
                    hpCount[2]++;
                    if(variantsHP != nullptr) (*variantsHP)[curPos] = SnpHP::GERMLINE_H2;
                }
                norCountPS[vcfSet[NORMAL].chrVariantPS[chrName][curPos]]++;
            }
        }
    // only tumor SNP at the current position
    }else if(curVar.isExists(TUMOR)){
        //the tumor SNP GT is phased heterozygous
        if(curVar.Variant[TUMOR].is_phased_hetero == true){
            if(curVar.Variant[TUMOR].allele.Ref == base || curVar.Variant[TUMOR].allele.Alt == base){
                std::map<int, int>::iterator posPSiter = vcfSet[TUMOR].chrVariantPS[chrName].find(curPos);
                if( posPSiter == vcfSet[TUMOR].chrVariantPS[chrName].end() ){
                    std::cerr<< curPos << "\t"
                             << curVar.Variant[TUMOR].allele.Ref << "\t"
                             << curVar.Variant[TUMOR].allele.Alt << "\n";
                    exit(EXIT_SUCCESS);
                }else{
                    OnlyTumorSNPjudgeHP(chrName, curPos, curVar, base, vcfSet, hpCount, &tumCountPS, variantsHP, tumorAllelePosVec, NorBase, SomaticPos);
                }
            }
        //the tumor SNP GT is unphased heterozygous
        }else if(curVar.Variant[TUMOR].is_unphased_hetero == true){
            if(curVar.Variant[TUMOR].allele.Ref == base || curVar.Variant[TUMOR].allele.Alt == base){
                OnlyTumorSNPjudgeHP(chrName, curPos, curVar, base, vcfSet, hpCount, nullptr, variantsHP, tumorAllelePosVec, NorBase, SomaticPos);
            }           
        //the tumor SNP GT is homozygous
        }else if(curVar.Variant[TUMOR].is_homozygous == true){
            if(curVar.Variant[TUMOR].allele.Ref == base || curVar.Variant[TUMOR].allele.Alt == base){
                OnlyTumorSNPjudgeHP(chrName, curPos, curVar, base, vcfSet, hpCount, nullptr, variantsHP, tumorAllelePosVec, NorBase, SomaticPos);
            }
        }
    }
}

void SomaticJudgeBase::OnlyTumorSNPjudgeHP(const std::string &chrName, int &curPos, MultiGenomeVar &curVar, std::string base, VCF_Info *vcfSet, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *tumorAllelePosVec, BamBaseCounter *NorBase, std::map<int, HP3_Info> *SomaticPos){

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
    varReadHpResult[pos].readHpCounter[hpResult]++;
    
    if(hpResult != ReadHP::unTag){
        if(BaseHP == SnpHP::SOMATIC_H3){
            if(hpResult != ReadHP::H1_1 && hpResult != ReadHP::H2_1 && hpResult != ReadHP::H3){
                std::cerr << "Error(recordReadHp) => error read hp : BaseHP: " <<BaseHP << " readHP: " << hpResult << " pos: " << pos+1 << std::endl; 
                exit(1);
            }
            varReadHpResult[pos].somaticSnpH3count++;
            varReadHpResult[pos].somaticBaseReadHpCounter[hpResult]++;
        }
    }
}

void SomaticJudgeBase::recordDeriveHp(int &pos, int &deriveHP, float deriveHPsimilarity, std::map<int, ReadHpResult> &varReadHpResult){
    if(deriveHP != SnpHP::GERMLINE_H1 && deriveHP != SnpHP::GERMLINE_H2 && deriveHP != SnpHP::NONE_SNP){
        std::cerr << "Error(recordDeriveHp) => error derive hp : pos: " <<pos+1 << " deriveHP: " << deriveHP << std::endl; 
        exit(1);        
    }
    varReadHpResult[pos].deriveHP = deriveHP;
    if(deriveHPsimilarity != 0.0){
        varReadHpResult[pos].deriveHPsimilarVec.emplace_back(deriveHPsimilarity);
        if(deriveHPsimilarity != 1.0){
            // std::cout << "deriveHPsimilarity: " << deriveHPsimilarity << "\n";
            // std::cout << "deriveHPsimilarityVec: " << varReadHpResult[pos].deriveHPsimilarVec.back() << "\n";
        }
    }
} 

ReadHpDistriLog::ReadHpDistriLog(){

}

ReadHpDistriLog::~ReadHpDistriLog(){

}

void ReadHpDistriLog::mergeLocalReadHp(const std::string &chr, std::map<int, ReadHpResult> &localReadHpResult){
    std::map<int, ReadHpResult>::iterator localReadHpIter = localReadHpResult.begin();
    while(localReadHpIter != localReadHpResult.end()){
        int pos = (*localReadHpIter).first;
        chrVarReadHpResult[chr][pos] = (*localReadHpIter).second;
        localReadHpIter++;
    }
}

void ReadHpDistriLog::writeReadHpDistriLog(const HaplotagParameters &params, std::string logPosfix, const std::vector<std::string> &chrVec){
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
        std::map<int, ReadHpResult>::iterator curVarReadHpIter = chrVarReadHpResult[chr].begin();
        while(curVarReadHpIter != chrVarReadHpResult[chr].end()){
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

void ReadHpDistriLog::writeTagReadCoverRegionLog(const HaplotagParameters &params, std::string logPosfix, const std::vector<std::string> &chrVec, std::map<std::string, int> &chrLength){
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


BamFileRAII::BamFileRAII(const std::string& BamFile, const std::string& fastaFile, htsThreadPool& threadPool):
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
                RefAlt tmp;
                tmp.Ref = fields[3]; 
                tmp.Alt = fields[4];

                tmp.is_phased_hetero = true;
                tmp.is_unphased_hetero = false;
                tmp.is_homozygous= false;

                Info.chrVariant[chr][pos]=tmp;

                VarData varData;
                varData.allele.Ref = fields[3];
                varData.allele.Alt = fields[4];
                varData.is_phased_hetero = true;

                if(Info.gene_type == NORMAL){
                    mergedChrVarinat[chr][pos].Variant[NORMAL] = varData;
                }else if(Info.gene_type == TUMOR){
                    mergedChrVarinat[chr][pos].Variant[TUMOR] = varData;
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
                    mergedChrVarinat[chr][pos].Variant[Info.gene_type].HP1 = fields[3];
                    mergedChrVarinat[chr][pos].Variant[Info.gene_type].HP2 = fields[4];
                }
                else if( fields[9][modifu_start] == '1' && fields[9][modifu_start+2] == '0' ){
                    Info.chrVariantHP1[chr][pos]=fields[4];
                    Info.chrVariantHP2[chr][pos]=fields[3];
                    mergedChrVarinat[chr][pos].Variant[Info.gene_type].HP1 = fields[4];
                    mergedChrVarinat[chr][pos].Variant[Info.gene_type].HP2 = fields[3];
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

                    VarData varData;
                    varData.allele.Ref = fields[3];
                    varData.allele.Alt = fields[4];

                    varData.is_unphased_hetero = false;

                    if(Info.gene_type == NORMAL){
                        mergedChrVarinat[chr][pos].Variant[NORMAL] = varData;
                    }else if(Info.gene_type == TUMOR){
                        mergedChrVarinat[chr][pos].Variant[TUMOR] = varData;
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

                    VarData varData;
                    varData.allele.Ref = fields[3];
                    varData.allele.Alt = fields[4];
                    varData.is_phased_hetero = false;
                    varData.is_unphased_hetero = true;
                    varData.is_homozygous = false;

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