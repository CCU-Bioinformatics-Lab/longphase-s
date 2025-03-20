#include "HaplotagProcess.h"



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
    std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat, 
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

    // if(params.writeReadLog){
    //     tagResult=new std::ofstream(params.resultPrefix + "_NorBam_tag.log");
    //     (*tagResult) << "##snpFile:"                 << params.snpFile                    << "\n";
    //     (*tagResult) << "##svFile:"                  << params.svFile                     << "\n";
    //     (*tagResult) << "##bamFile:"                 << params.bamFile                    << "\n";
    //     (*tagResult) << "##resultPrefix:"            << params.resultPrefix               << "\n";
    //     (*tagResult) << "##numThreads:"              << params.numThreads                 << "\n";
    //     (*tagResult) << "##region:"                  << params.region                     << "\n";
    //     (*tagResult) << "##qualityThreshold:"        << params.qualityThreshold           << "\n";
    //     (*tagResult) << "##percentageThreshold:"     << params.percentageThreshold        << "\n";
    //     (*tagResult) << "##tagSupplementary:"        << params.tagSupplementary           << "\n";
    //     (*tagResult) << "#ReadID\t"
    //                  << "CHROM\t"
    //                  << "ReadStart\t"
    //                  << "Confidnet(%)\t";
    //     (*tagResult) << "Haplotype\t"
    //                  << "PhaseSet\t"
    //                  << "TotalAllele\t"
    //                  << "HP1Allele\t"
    //                  << "HP2Allele\t";
    //     (*tagResult) << "phasingQuality(PQ)\t"
    //                  << "(Variant,HP)\t"
    //                  << "(PhaseSet,Variantcount)\n";
    // }

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
        // fetch chromosome string
        std::string ref_string = fastaParser.chrString.at(chr);

        //inintial iterator
        std::map<int, RefAltSet>::iterator firstVariantIter = currentVariants.begin();

        std::map<int, RefAltSet>::reverse_iterator last = currentVariants.rbegin();

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
    std::map<int, RefAltSet> &currentVariants,
    std::map<int, RefAltSet>::iterator &firstVariantIter, 
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
                    // printf("into match cigar\n");

                    //statistically analyze SNP information exclusive to the tumor
                    if((*currentVariantIter).second.isExistTumor){
                        int curPos = (*currentVariantIter).first;
                        //std::cout << "curPos :" << curPos << "is tumor "<<(*currentVariantIter).second.isExistTumor  << " isNormal: "<< (*currentVariantIter).second.isExistNormal<< std::endl;

                        //detect ref bsae length(temp :tumor SNP)
                        int tumRefLength = (*currentVariantIter).second.Variant[Genome::TUMOR].Ref.length();
                        int tumAltLength = (*currentVariantIter).second.Variant[Genome::TUMOR].Alt.length();
                        
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

                    if ( aln.core.qual >= params.somaticCallingMpqThreshold && (*currentVariantIter).second.isExistNormal){
                        // only judge the heterozygous SNP
                        if((*currentVariantIter).second.Variant[Genome::NORMAL].is_phased_hetero){
                            auto norVar = (*currentVariantIter).second.Variant[Genome::NORMAL];
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
                if((*currentVariantIter).second.isExistTumor){

                    int curPos = (*currentVariantIter).first;
                    
                    //record tumor SNP position
                    tumVarPosVec.push_back(curPos);

                    //detect ref bsae length(temp :tumor SNP)
                    int tumRefLength = (*currentVariantIter).second.Variant[Genome::TUMOR].Ref.length();
                    int tumAltLength = (*currentVariantIter).second.Variant[Genome::TUMOR].Alt.length();

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
                if ( aln.core.qual >= params.somaticCallingMpqThreshold && (*currentVariantIter).second.isExistNormal && !alreadyJudgeDel){
                    if((*currentVariantIter).second.Variant[Genome::NORMAL].is_phased_hetero){
                        // longphase v1.73 only execute once
                        alreadyJudgeDel = true;
                        auto norVar = (*currentVariantIter).second.Variant[Genome::NORMAL];
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

    //write tag log file for testing the logic of the judge haplotype
    // if(params.writeReadLog && aln.core.qual >= params.somaticCallingMpqThreshold){
        // writeGermlineTagLog(*tagResult, aln, bamHdr, hpResult, max, min, hp1Count, hp2Count, pqValue, variantsHP, countPS);
    // }
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

void germlineJudgeBase::germlineGetRefLastVarPos(
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


void germlineJudgeBase::germlineJudgeSnpHap(
    const std::string& chrName,
    VCF_Info* vcfSet,
    RefAlt& norVar,
    const std::string& base,
    int& ref_pos,
    int& length,
    int& i,
    int& aln_core_n_cigar,
    uint32_t* cigar,
    std::map<int, RefAltSet>::iterator &currentVariantIter,
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


            std::map<int, int>::iterator posPSiter = vcfSet[Genome::NORMAL].chrVariantPS[chrName].find((*currentVariantIter).first);

            if( posPSiter == vcfSet[Genome::NORMAL].chrVariantPS[chrName].end() ){
                std::cerr << "ERROR (germlineJudgeSnpHap) => can't find the position:" 
                          << " chr: " << chrName << "\t"
                          << " pos: " << curPos << "\t"
                          << " ref: " << norVar.Ref << "\t"
                          << " alt: " << norVar.Alt << "\n";
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


void germlineJudgeBase::germlineJudgeDeletionHap(
    const std::string& chrName,
    const std::string& ref_string,
    int& ref_pos,
    int& length,
    int& query_pos,
    std::map<int, RefAltSet>::iterator &currentVariantIter,
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
                auto norVar = (*currentVariantIter).second.Variant[Genome::NORMAL];
                int refAlleleLen = norVar.Ref.length();
                int altAlleleLen = norVar.Alt.length();
                
                // SNP
                if (refAlleleLen == 1 && altAlleleLen == 1) {
                    // get the next match
                    char base_chr = seq_nt16_str[bam_seqi(bam_get_seq(aln), query_pos)];
                    std::string base(1, base_chr);

                    if (base == vcfSet[Genome::NORMAL].chrVariantHP1[chrName][curPos]) {
                        hp1Count++;
                        variantsHP[curPos] = 0;
                    }
                    if (base == vcfSet[Genome::NORMAL].chrVariantHP2[chrName][curPos]) {
                        hp2Count++;
                        variantsHP[curPos] = 1;
                    }
                    countPS[vcfSet[Genome::NORMAL].chrVariantPS[chrName][curPos]]++;
                }
                
                // the read deletion contain VCF's deletion
                else if (refAlleleLen != 1 && altAlleleLen == 1) {

                    int hp1Length = vcfSet[Genome::NORMAL].chrVariantHP1[chrName][curPos].length();
                    int hp2Length = vcfSet[Genome::NORMAL].chrVariantHP2[chrName][curPos].length();
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
                    countPS[vcfSet[Genome::NORMAL].chrVariantPS[chrName][curPos]]++;
                }
            }
        }
    }
}

void germlineJudgeBase::germlineJudgeSVHap(const bam1_t &aln, VCF_Info* vcfSet, int& hp1Count, int& hp2Count, const int& tagGeneType){
    auto readIter = vcfSet[tagGeneType].readSVHapCount.find(bam_get_qname(&aln));
    if( readIter != vcfSet[tagGeneType].readSVHapCount.end() ){
        hp1Count += vcfSet[tagGeneType].readSVHapCount[bam_get_qname(&aln)][0];
        hp2Count += vcfSet[tagGeneType].readSVHapCount[bam_get_qname(&aln)][1];
    }
}

int germlineJudgeBase::germlineDetermineReadHap(
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

void germlineJudgeBase::writeGermlineTagLog(std::ofstream& tagResult, const bam1_t& aln, const bam_hdr_t& bamHdr, int& hpResult, double& max, double& min, int& hp1Count, int& hp2Count, int& pqValue, const std::map<int, int>& variantsHP, const std::map<int, int>& countPS){
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

void SomaticJudgeBase::SomaticJudgeSnpHP(std::map<int, RefAltSet>::iterator &currentVariantIter, VCF_Info *vcfSet, std::string chrName, std::string base, std::map<int, int> &hpCount, std::map<int, int> &norCountPS, std::map<int, int> &tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *tumorAllelePosVec, BamBaseCounter *NorBase, std::map<int, HP3_Info> *SomaticPos){
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
                    OnlyTumorSNPjudgeHP(chrName, curPos, curVar, base, vcfSet, hpCount, &tumCountPS, variantsHP, tumorAllelePosVec, NorBase, SomaticPos);
                }
            }
        //the tumor SNP GT is unphased heterozygous
        }else if(curVar.Variant[Genome::TUMOR].is_unphased_hetero == true){
            if(curVar.Variant[Genome::TUMOR].Ref == base || curVar.Variant[Genome::TUMOR].Alt == base){
                OnlyTumorSNPjudgeHP(chrName, curPos, curVar, base, vcfSet, hpCount, nullptr, variantsHP, tumorAllelePosVec, NorBase, SomaticPos);
            }           
        //the tumor SNP GT is homozygous
        }else if(curVar.Variant[Genome::TUMOR].is_homozygous == true){
            if(curVar.Variant[Genome::TUMOR].Ref == base || curVar.Variant[Genome::TUMOR].Alt == base){
                OnlyTumorSNPjudgeHP(chrName, curPos, curVar, base, vcfSet, hpCount, nullptr, variantsHP, tumorAllelePosVec, NorBase, SomaticPos);
            }
        }
    }
}

void SomaticJudgeBase::OnlyTumorSNPjudgeHP(const std::string &chrName, int &curPos, RefAltSet &curVar, std::string base, VCF_Info *vcfSet, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *tumorAllelePosVec, BamBaseCounter *NorBase, std::map<int, HP3_Info> *SomaticPos){

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
    chrPosSomaticInfo = new std::map<std::string, std::map<int, HP3_Info>>();
    chrVarReadHpResult = new std::map<std::string, std::map<int, ReadHpResult>>();
    callerReadHpDistri = new readHpDistriLog();
    denseTumorSnpInterval = new std::map<std::string, std::map<int, std::pair<int, denseSnpInterval>>>();
    chrReadHpResultSet = new std::map<std::string, std::map<std::string, ReadVarHpCount>>();
    chrTumorPosReadCorrBaseHP = new std::map<std::string, std::map<int, std::map<std::string, int>>>();
}

SomaticVarCaller::~SomaticVarCaller(){
    releaseMemory();
}

void SomaticVarCaller::releaseMemory(){
    delete chrPosSomaticInfo;
    delete chrVarReadHpResult;
    delete callerReadHpDistri;
    delete denseTumorSnpInterval;
    delete chrReadHpResultSet;
    delete chrTumorPosReadCorrBaseHP;
}

void SomaticVarCaller::VariantCalling(const std::string BamFile, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat,const std::vector<std::string> &chrVec, std::map<std::string, int> &chrLength, HaplotagParameters &params, VCF_Info *vcfSet, BamBaseCounter &NorBase){
    std::cerr << "collecting data for the tumor sample... ";
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
        (*chrPosSomaticInfo)[chr] = std::map<int, HP3_Info>();
        (*chrVarReadHpResult)[chr] = std::map<int, ReadHpResult>();
        (*denseTumorSnpInterval)[chr] = std::map<int, std::pair<int, denseSnpInterval>>();
        (*chrReadHpResultSet)[chr] = std::map<std::string, ReadVarHpCount>();
        (*chrTumorPosReadCorrBaseHP)[chr] = std::map<int, std::map<std::string, int>>();
    }

    // somatic calling filter params
    SomaticFilterParaemter somaticParams;

    // setting somatic calling filter params
    InitialSomaticFilterParams(somaticParams, params.enableFilter);  

    // init data structure and get core n
    htsThreadPool threadPool = {NULL, 0};
    // creat thread pool
    if (!(threadPool.pool = hts_tpool_init(params.numThreads))) {
        fprintf(stderr, "Error creating thread pool\n");
        exit(1);
    }

    // record reference last variant pos
    std::vector<int> last_tumor_pos;
    for( auto chr : chrVec ){
        bool existLastTumorPos = false;

        for (auto lastVariantIter = mergedChrVarinat[chr].rbegin(); lastVariantIter != mergedChrVarinat[chr].rend(); ++lastVariantIter) {
            if ((*lastVariantIter).second.isExistTumor) {
                last_tumor_pos.push_back((*lastVariantIter).first);
                existLastTumorPos = true;
                break;
            }
        }
        
        if(!existLastTumorPos){
            last_tumor_pos.push_back(0);
        }
    }

    // reference fasta parser
    FastaParser fastaParser(params.fastaFile, chrVec, last_tumor_pos, params.numThreads);

    // loop all chromosome
    #pragma omp parallel for schedule(dynamic) num_threads(params.numThreads) 
    for(auto chr : chrVec ){
        
        // bam file resource allocation
        BamFileRAII bam(BamFile, params.fastaFile, threadPool);

        // record the position that tagged as HP3
        // chr, variant position
        std::map<int, HP3_Info> *somaticPosInfo = nullptr;

        // read ID, SNP HP count 
        std::map<std::string, ReadVarHpCount> *readHpResultSet = nullptr;
        // position, read ID, baseHP 
        std::map<int, std::map<std::string, int>> *tumorPosReadCorrBaseHP = nullptr;

        // records all variants within this chromosome.
        std::map<int, RefAltSet> currentChrVariants;


        #pragma omp critical
        {
            somaticPosInfo = &((*chrPosSomaticInfo)[chr]);
            readHpResultSet = &((*chrReadHpResultSet)[chr]);
            tumorPosReadCorrBaseHP = &((*chrTumorPosReadCorrBaseHP)[chr]);
            currentChrVariants = mergedChrVarinat[chr];
        }

        // since each read is sorted based on the start coordinates, to save time, 
        // firstVariantIter keeps track of the first variant that each read needs to check.
        std::map<int, RefAltSet>::iterator firstVariantIter = currentChrVariants.begin();
        // get the coordinates of the last variant
        // the tagging process will not be perform if the read's start coordinate are over than last variant.
        std::map<int, RefAltSet>::reverse_iterator lastVariant = currentChrVariants.rbegin();
        
        // fetch chromosome string
        std::string ref_string = fastaParser.chrString.at(chr);

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
                StatisticTumorVariantData(*bam.bamHdr, *bam.aln, chr, params, &NorBase, vcfSet, *somaticPosInfo, currentChrVariants, firstVariantIter, *readHpResultSet, *tumorPosReadCorrBaseHP, ref_string);
            }
        }
        hts_itr_destroy(Iter);
    }
    hts_tpool_destroy(threadPool.pool);
    std::cerr<< difftime(time(NULL), begin) << "s\n";

    double tumorPurity = predictTumorPurity(params, chrVec, NorBase);

    // set filter params with tumor purity
    // std::cerr << "[Debug] Setting tumorPurity = 1.0" << std::endl;
    // tumorPurity = 1.0;
    SetFilterParamsWithPurity(somaticParams, tumorPurity);
    
    std::cerr<< "calling somatic variants ... ";
    begin = time(NULL);
    // loop all chromosome
    #pragma omp parallel for schedule(dynamic) num_threads(params.numThreads) 
    for(auto chr : chrVec ){
        // record the position that tagged as HP3
        // chr, variant position
        std::map<int, HP3_Info> *somaticPosInfo = nullptr;
        // read ID, reads hpResult 
        std::map<int, ReadHpResult> *localCallerReadHpDistri = nullptr;
        
        // pos, (endPos, denseSnpInterval)
        std::map<int, std::pair<int, denseSnpInterval>> *localDenseTumorSnpInterval = nullptr;

        // read ID, SNP HP count 
        std::map<std::string, ReadVarHpCount> *readHpResultSet = nullptr;
        // position, read ID, baseHP 
        std::map<int, std::map<std::string, int>> *tumorPosReadCorrBaseHP = nullptr;

        // records all variants within this chromosome.
        std::map<int, RefAltSet> currentChrVariants;


        #pragma omp critical
        {
            somaticPosInfo = &((*chrPosSomaticInfo)[chr]);
            localCallerReadHpDistri = &((*chrVarReadHpResult)[chr]);
            localDenseTumorSnpInterval = &((*denseTumorSnpInterval)[chr]);
            readHpResultSet = &((*chrReadHpResultSet)[chr]);
            tumorPosReadCorrBaseHP = &((*chrTumorPosReadCorrBaseHP)[chr]);
            currentChrVariants = mergedChrVarinat[chr];
        }

        // since each read is sorted based on the start coordinates, to save time, 
        // firstVariantIter keeps track of the first variant that each read needs to check.
        std::map<int, RefAltSet>::iterator firstVariantIter = currentChrVariants.begin();
        // get the coordinates of the last variant
        // the tagging process will not be perform if the read's start coordinate are over than last variant.
        std::map<int, RefAltSet>::reverse_iterator lastVariant = currentChrVariants.rbegin();

        //get close somatic SNP interval
        getDenseTumorSnpInterval(*somaticPosInfo, *readHpResultSet, *tumorPosReadCorrBaseHP, *localDenseTumorSnpInterval);

        //calculate information and filter somatic SNPs
        SomaticFeatureFilter(somaticParams, currentChrVariants, chr, *somaticPosInfo, NorBase, tumorPurity);

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

        //merge local read hp result to global read hp result
        #pragma omp critical
        {
            callerReadHpDistri->mergeLocalReadHp(chr, *localCallerReadHpDistri);
        }
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
        WriteSomaticVarCallingLog(params ,somaticParams, chrVec, NorBase, mergedChrVarinat);

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

void SomaticVarCaller::InitialSomaticFilterParams(SomaticFilterParaemter &somaticParams, bool enableFilter){
    
    // Determine whether to apply the filter
    somaticParams.applyFilter = enableFilter;
    somaticParams.writeVarLog = true;

    somaticParams.tumorPurity = 1.0;

    // Below the mapping quality read ratio threshold
    somaticParams.LowMpqRatioThreshold = 0.1;

    somaticParams.MessyReadRatioThreshold = 1.0;
    somaticParams.MessyReadCountThreshold = 3;

    somaticParams.HapConsistency_ReadCount_Thr = 10;
    somaticParams.HapConsistency_VAF_Thr = 0.2;

    somaticParams.IntervalSnpCount_ReadCount_Thr = 10;
    somaticParams.IntervalSnpCount_VAF_Thr = 0.15;
}

void SomaticVarCaller::StatisticTumorVariantData(const  bam_hdr_t &bamHdr,const bam1_t &aln, const std::string &chr, HaplotagParameters &params, BamBaseCounter *NorBase, VCF_Info *vcfSet, std::map<int, HP3_Info> &somaticPosInfo, std::map<int, RefAltSet> &currentChrVariants, std::map<int, RefAltSet>::iterator &firstVariantIter, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, std::map<std::string, int>> &tumorPosReadCorrBaseHP, std::string &ref_string){
    
    std::map<int, int> hpCount;
    hpCount[1] = 0; 
    hpCount[2] = 0;
    hpCount[3] = 0;
    hpCount[4] = 0;

    //record variants on this read
    std::map<int, int> variantsHP;

    //record tumor-unique variants with low VAF in the normal.bam on this read
    std::vector<int> tumorAllelePosVec;

    //record tumor SNPs on this read
    std::vector<int> tumorSnpPosVec;

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
                        SomaticJudgeSnpHP(currentVariantIter, vcfSet , chr, base, hpCount, NorCountPS, TumCountPS, &variantsHP, &tumorAllelePosVec, NorBase, &somaticPosInfo);
                        if((*currentVariantIter).second.isExistTumor){
                            tumorSnpPosVec.push_back((*currentVariantIter).first);
                        }
                    }

                    //statistically analyze SNP information exclusive to the tumor
                    if((*currentVariantIter).second.isExistTumor){
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
                if((*currentVariantIter).second.isExistTumor){

                    int curPos = (*currentVariantIter).first;

                    //detect ref bsae length(temp :tumor SNP)
                    int tumRefLength = (*currentVariantIter).second.Variant[Genome::TUMOR].Ref.length();
                    int tumAltLength = (*currentVariantIter).second.Variant[Genome::TUMOR].Alt.length();

                    if(tumRefLength == 1 && tumAltLength == 1){
                        somaticPosInfo[curPos].base.delCount++;
                        somaticPosInfo[curPos].base.depth++;
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
    
    int pqValue = 0;   
    double normalHPsimilarity = 0.0;
    double tumorHPsimilarity = 0.0;
    //calculate germline read HP result for predict tumor purity
    //calculate somatic read HP result for somatic read consistency filter
    int hpResult =determineReadHP(hpCount, pqValue, NorCountPS, normalHPsimilarity, tumorHPsimilarity, params.percentageThreshold, nullptr, nullptr, nullptr);
   
    //classify read cases where tumor SNPs have low VAF in normal samples
    if(!tumorAllelePosVec.empty()){
        ClassifyReadsByCase(tumorAllelePosVec, NorCountPS, hpCount, params, somaticPosInfo);
        
        for(auto pos : tumorAllelePosVec){
            int baseHP = SnpHP::NONE_SNP;
            if(variantsHP.find(pos) != variantsHP.end()){
                baseHP = variantsHP[pos];
            }else{
                std::cerr << "Error (SomaticStatisticSomaticPosInfo) => can't find the position" << std::endl;
                std::cerr << "chr:" << chr << " pos: " << pos+1 << std::endl;
                std::cerr << "readID: " << bam_get_qname(&aln) << std::endl;
                exit(1);
            }
            if(baseHP != SnpHP::SOMATIC_H3){
                std::cerr << "Error (SomaticStatisticSomaticPosInfo) => baseHP is not HP3 : chr:" << chr << " pos: " << pos+1 << " baseHP: " << baseHP << std::endl;
                exit(1);
            }
            if(hpResult == ReadHP::H1_1 || hpResult == ReadHP::H2_1 || hpResult == ReadHP::H3 || hpResult == ReadHP::unTag){
                // record the somatic read HP for somatic read consistency filter
                somaticPosInfo[pos].somaticReadHpCount[hpResult]++;
            }else if(hpResult == ReadHP::H1 || hpResult == ReadHP::H2){
                std::cerr << "Error (SomaticStatisticSomaticPosInfo) => error somatic read HP : chr:" << chr << " pos: " << pos+1 << " hpResult: " << hpResult << std::endl;
                exit(1);      
            }
        }
    }

    //record variants HP count for each read in tumor-only position
    if(!tumorSnpPosVec.empty()){

        std::string readID = bam_get_qname(&aln);
        //read ID overide
        if(readHpResultSet.find(readID) != readHpResultSet.end()){
            readHpResultSet[readID].readIDcount++;
            readID = readID + "-" + std::to_string(readHpResultSet[readID].readIDcount);
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
        for(auto pos : tumorSnpPosVec){
            int tumorSnpBaseHP = SnpHP::NONE_SNP;

            if(variantsHP.find(pos) != variantsHP.end()){
                tumorSnpBaseHP = variantsHP[pos];
            }
            // record the base HP whatever it is germline or somatic for statistic read distribution at current position
            tumorPosReadCorrBaseHP[pos][readID] = tumorSnpBaseHP;
            // record the base HP for predict tumor purity
            somaticPosInfo[pos].allReadHpCount[hpResult]++;
        } 
    }
}

void SomaticVarCaller::OnlyTumorSNPjudgeHP(const std::string &chrName, int &curPos, RefAltSet &curVar, std::string base, VCF_Info *vcfSet, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *tumorAllelePosVec, BamBaseCounter *NorBase, std::map<int, HP3_Info> *SomaticPos){
    //the tumor SNP GT is phased heterozygous
    //all bases of the same type at the current position in normal.bam

    if(tumorAllelePosVec == nullptr){
        std::cerr << "ERROR (SomaticDetectJudgeHP) => tumorAllelePosVec pointer cannot be nullptr"<< std::endl;
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
        // if((*NorBase).getMaxFreqBase(chrName, curPos) == TumorRefBase){
        
            if(base == TumorAltBase){
                hpCount[3]++;
                if(variantsHP != nullptr) (*variantsHP)[curPos] = SnpHP::SOMATIC_H3;

                //record postions that tagged as HP3 for calculating the confidence of somatic positions
                (*tumorAllelePosVec).push_back(curPos);
                (*SomaticPos)[curPos].isNormalPosLowVAF = true;                            

            //base is not match to TumorRefBase & TumorAltBase (other HP)
            }else if(base != TumorRefBase && base != TumorAltBase){
                //hpCount[4]++;
                //if(variantsHP != nullptr) (*variantsHP)[curPos] = SnpHP::SOMATIC_H4;
                //(*tumorAllelePosVec).push_back(curPos);
                //(*SomaticPos)[curPos].isNormalPosLowVAF = true;  
            }

            if(tumCountPS != nullptr) (*tumCountPS)[vcfSet[Genome::TUMOR].chrVariantPS[chrName][curPos]]++;

        //max count base not match to tumorRefBase in normal.bam
        // }else if((*NorBase).getMaxFreqBase(chrName, curPos) != TumorRefBase){
        //     //temp 
        // }
    //exist more than one type of base at the current position in normal.bam 
    // }else{
    //     //temp
    }
}

double SomaticVarCaller::predictTumorPurity(const HaplotagParameters &params, const std::vector<std::string> &chrVec, BamBaseCounter &NorBase){
    std::cerr << "predicting tumor purity... ";
    std::time_t begin = time(NULL);

    int filter_consistencyRatioInNorBam_count = 0;
    int filter_consistencyRatioInNorBam_over_thr_count = 0;
    int filter_consistencyRatio_count = 0;
    int filter_readHpCountInNorBam_count = 0;
    int filter_percentageOfGermlineHpInNorBam_count = 0;

    // median, iqr
    std::vector<double> purityFeatureValueVec;
    std::vector<std::pair<std::string, int>> valueIndexVec;

    size_t initial_data_size = 0;
    for(auto chr : chrVec){
        std::map<int, HP3_Info>::iterator somaticPosIter = (*chrPosSomaticInfo)[chr].begin();
        while(somaticPosIter != (*chrPosSomaticInfo)[chr].end()){
            int H1readCount = (*somaticPosIter).second.allReadHpCount[ReadHP::H1];
            int H2readCount = (*somaticPosIter).second.allReadHpCount[ReadHP::H2];
            int tumDepth = (*somaticPosIter).second.base.depth;
            initial_data_size++;

            // the ratio of the read count of H1 and H2
            int germlineReadHpCount = H1readCount + H2readCount;
            double germlineReadHpConsistencyRatio = 0.0;
            if(H1readCount > 0 && H2readCount > 0){
                germlineReadHpConsistencyRatio = (H1readCount > H2readCount) ? ((double)H1readCount / (double)(H1readCount + H2readCount)) : ((double)H2readCount / (double)(H1readCount + H2readCount));
            }else if(H1readCount == 0 && H2readCount == 0){
                germlineReadHpConsistencyRatio = 0.0;
            }
            else{
                germlineReadHpConsistencyRatio = 1.0;
            }

            // the ratio of the germlineHP based on depth
            double percentageOfGermlineHp = 0.0;
            if(tumDepth > 0 && germlineReadHpCount > 0){
                percentageOfGermlineHp = (double)germlineReadHpCount / (double)tumDepth;            
            }


            //read hp count in the normal bam
            int norDepth = NorBase.getDepth(chr, (*somaticPosIter).first);
            int H1readCountInNorBam = NorBase.getReadHpCountInNorBam(chr, (*somaticPosIter).first, ReadHP::H1);
            int H2readCountInNorBam = NorBase.getReadHpCountInNorBam(chr, (*somaticPosIter).first, ReadHP::H2);
            int germlineReadHpCountInNorBam = H1readCountInNorBam + H2readCountInNorBam;

            double germlineReadHpConsistencyRatioInNorBam = 0.0;
            if(H1readCountInNorBam > 0 && H2readCountInNorBam > 0){
                germlineReadHpConsistencyRatioInNorBam = (H1readCountInNorBam > H2readCountInNorBam) ? ((double)H1readCountInNorBam / (double)germlineReadHpCountInNorBam) : ((double)H2readCountInNorBam / (double)germlineReadHpCountInNorBam);
            }else if(H1readCountInNorBam == 0 && H2readCountInNorBam == 0){
                germlineReadHpConsistencyRatioInNorBam = 0.0;
            }
            else{
                germlineReadHpConsistencyRatioInNorBam = 1.0;
            }

            double percentageOfGermlineHpInNorBam = 0.0;
            if(norDepth > 0 && germlineReadHpCountInNorBam > 0){
                percentageOfGermlineHpInNorBam = (double)germlineReadHpCountInNorBam / (double)norDepth;
            }

            bool includeInStatistics = true;
            // if(   germlineReadHpConsistencyRatioInNorBam == 0.0 
            //    || germlineReadHpConsistencyRatioInNorBam >= 0.7
            //    || germlineReadHpCountInNorBam <= 5
            //    || percentageOfGermlineHpInNorBam <= 0.6
            //    || germlineReadHpConsistencyRatio == 0.0
            // ){
            //     includeInStatistics = false;
            // }

            if(germlineReadHpConsistencyRatioInNorBam == 0.0){
                includeInStatistics = false;
                filter_consistencyRatioInNorBam_count++;
            }else if(germlineReadHpConsistencyRatio == 0.0){
                includeInStatistics = false;
                filter_consistencyRatio_count++;
            }else if(germlineReadHpConsistencyRatioInNorBam >= 0.7){
                includeInStatistics = false;
                filter_consistencyRatioInNorBam_over_thr_count++;
            }else if(germlineReadHpCountInNorBam <= 5){
                includeInStatistics = false;
                filter_readHpCountInNorBam_count++;
            }else if(percentageOfGermlineHpInNorBam <= 0.6){
                includeInStatistics = false;
                filter_percentageOfGermlineHpInNorBam_count++;
            }else if(includeInStatistics){
                purityFeatureValueVec.push_back(germlineReadHpConsistencyRatio);
                valueIndexVec.push_back(std::make_pair(chr, (*somaticPosIter).first));
                (*somaticPosIter).second.statisticPurity = true;
            }

            somaticPosIter++;
        }
    }
    std::cerr << "initial data size: " << initial_data_size << std::endl;


    std::cerr << "\n==========first filter==========" << std::endl;
    std::cerr << "[INFO] Data count: " << purityFeatureValueVec.size() << std::endl;
    std::cerr << "[INFO] consistencyRatioInNorBam == 0.0: " << filter_consistencyRatioInNorBam_count << std::endl;
    std::cerr << "[INFO] consistencyRatio == 0.0: " << filter_consistencyRatio_count << std::endl;
    std::cerr << "[INFO] consistencyRatioInNorBam <= 0.7: " << filter_consistencyRatioInNorBam_over_thr_count << std::endl;
    std::cerr << "[INFO] readHpCountInNorBam <= 5: " << filter_readHpCountInNorBam_count << std::endl;
    std::cerr << "[INFO] percentageOfGermlineHpInNorBam: " << filter_percentageOfGermlineHpInNorBam_count << std::endl;

    BoxPlotValue plotValue = statisticPurityData(purityFeatureValueVec);

    double median = plotValue.median;
    double iqr = plotValue.iqr;

    //purity prediction model
    double purity;

    //remove outliers for reduce noise 
    size_t test_outliers = 0;
    int iteration_times = 3;
    for(int i = 0; i < iteration_times; i++){
        std::printf("==========iteration %d==========\n", i+1);
        std::vector<double>::iterator vecIter = purityFeatureValueVec.begin();
        std::vector<std::pair<std::string, int>>::iterator valueIndexIter = valueIndexVec.begin();

        if(purityFeatureValueVec.size() != valueIndexVec.size()) {
            std::cerr << "Error: Vector sizes don't match!" << std::endl;
            exit(1);
        }

        while(vecIter != purityFeatureValueVec.end()){
            if(*vecIter < plotValue.lowerWhisker || *vecIter > plotValue.upperWhisker){
                //remove outliers for reduce noise
                vecIter = purityFeatureValueVec.erase(vecIter);

                //remove ouliers index
                std::string chr = valueIndexIter->first;
                int pos = valueIndexIter->second;
                (*chrPosSomaticInfo)[chr][pos].statisticPurity = false;
                valueIndexIter = valueIndexVec.erase(valueIndexIter);
                test_outliers++;
            }else{
                vecIter++;
                valueIndexIter++;
            }
        }

        plotValue = statisticPurityData(purityFeatureValueVec);
        std::printf("[INFO] wisker : lower = %f, upper = %f\n", plotValue.lowerWhisker, plotValue.upperWhisker);
        std::printf("[INFO] after remove outliers: %ld, data size: %ld\n", test_outliers, plotValue.data_size);
    }


    median = plotValue.median;
    iqr = plotValue.iqr;

    // purity prediction model
    // purity = -2.1031 * median + 13.1080 * iqr + 3.1823 * median * median + -11.9640 * median * iqr + -4.8450 * iqr * iqr + -0.0921;
    // purity = -2.5622 * median + 13.3930 * iqr + 3.4934 * median * median + -12.2045 * median * iqr + -4.9453 * iqr * iqr + 0.0548; //iteration = 2
    purity = -2.7106 * median + 13.6299 * iqr + 3.6107 * median * median + -12.5210 * median * iqr + -4.9128 * iqr * iqr + 0.0978; // iteration = 3

    if(purity > 1.15){
        std::cerr << "[ERROR] The value of purity exceeds the model's prediction : " << purity << std::endl;
        exit(1);
    }else if(purity > 1.0){
        purity = 1.0;
    }else if(purity < 0.0){
        std::cerr << "[ERROR] The value of purity exceeds the model's prediction : " << purity << std::endl;
        exit(1);
    }

    std::cerr<< difftime(time(NULL), begin) << "s\n";
    std::cerr << "[INFO] ===== tumor purity ===== : " << purity << std::endl;
    std::cerr << "[INFO] median: " << plotValue.median << std::endl;
    std::cerr << "[INFO] iqr: " << plotValue.iqr << std::endl;
    std::cerr << "[INFO] q1: " << plotValue.q1 << std::endl;
    std::cerr << "[INFO] q3: " << plotValue.q3 << std::endl;
    std::cerr << "[INFO] lowerWhisker: " << plotValue.lowerWhisker << std::endl;
    std::cerr << "[INFO] upperWhisker: " << plotValue.upperWhisker << std::endl;
    std::cerr << "[INFO] outliers: " << plotValue.outliers << std::endl;
    
    std::ofstream *purityLog = new std::ofstream(params.resultPrefix+"_purity.out");

    if(!purityLog->is_open()){
        std::cerr<< "Fail to open write file: " << params.resultPrefix+"_purity.out" << "\n";
        exit(1);
    }

    (*purityLog) << "#Initial data size: " << initial_data_size << std::endl;
    (*purityLog) << "#==========initial filter==========" << std::endl;
    (*purityLog) << "#consistencyRatioInNorBam: " << filter_consistencyRatioInNorBam_count << std::endl;
    (*purityLog) << "#consistencyRatio: " << filter_consistencyRatio_count << std::endl;
    (*purityLog) << "#consistencyRatioInNorBam_over_thr: " << filter_consistencyRatioInNorBam_over_thr_count << std::endl;
    (*purityLog) << "#readHpCountInNorBam: " << filter_readHpCountInNorBam_count << std::endl;
    (*purityLog) << "#percentageOfGermlineHpInNorBam: " << filter_percentageOfGermlineHpInNorBam_count << std::endl;
    (*purityLog) << "#==========wisker filter===========" << std::endl;
    (*purityLog) << "#iteration times: " << iteration_times << std::endl;
    (*purityLog) << "#remove outliers: " << test_outliers << std::endl;
    (*purityLog) << "#==========purity prediction===========" << std::endl;
    (*purityLog) << "Tumor purity: " << purity << std::endl;
    (*purityLog) << "Data size: " << plotValue.data_size << std::endl;
    (*purityLog) << "Group:" << std::endl;
    (*purityLog) << "Median: " << median << std::endl;
    (*purityLog) << "Q1: " << plotValue.q1 << std::endl;
    (*purityLog) << "Q3: " << plotValue.q3 << std::endl;
    (*purityLog) << "IQR: " << plotValue.iqr << std::endl;
    (*purityLog) << "Whiskers: " << plotValue.lowerWhisker << " to " << plotValue.upperWhisker << std::endl;
    (*purityLog) << "Outliers: " << plotValue.outliers << std::endl;

    (*purityLog).close();
    delete purityLog;
    purityLog = nullptr;

    return purity;
}

BoxPlotValue SomaticVarCaller::statisticPurityData(std::vector<double> &purityFeatureValueVec){
    BoxPlotValue plotValue;
    plotValue.data_size = purityFeatureValueVec.size();
    
    // Handle empty data case
    if (plotValue.data_size == 0) {
        plotValue.median = plotValue.q1 = plotValue.q3 = 0.0;
        plotValue.iqr = plotValue.lowerWhisker = plotValue.upperWhisker = 0.0;
        plotValue.outliers = 0;
        std::cerr << "[ERROR] (statisticPurityData) The data size is 0" << std::endl;
        return plotValue;
    }

    // Sort the original data
    std::sort(purityFeatureValueVec.begin(), purityFeatureValueVec.end());

    // Define percentile calculation function using linear interpolation
    auto percentile = [&](double p) -> double {
        if (p < 0.0 || p > 1.0) {
            throw std::invalid_argument("Percentile must be between 0 and 1");
        }
        
        double pos = p * (plotValue.data_size - 1);
        size_t idx = static_cast<size_t>(pos);
        double frac = pos - idx;
        
        if (idx + 1 >= plotValue.data_size) {
            return purityFeatureValueVec[plotValue.data_size - 1];
        }
        
        return purityFeatureValueVec[idx] * (1.0 - frac) + 
               purityFeatureValueVec[idx + 1] * frac;
    };

    try {
        // Calculate quartiles and median
        plotValue.q1 = percentile(0.25);      // First quartile
        plotValue.median = percentile(0.5);    // Median
        plotValue.q3 = percentile(0.75);      // Third quartile
        plotValue.iqr = plotValue.q3 - plotValue.q1;  // Interquartile range
        
        // Calculate whiskers (1.5 * IQR rule)
        plotValue.lowerWhisker = std::max(0.0, plotValue.q1 - 1.5 * plotValue.iqr);
        plotValue.upperWhisker = plotValue.q3 + 1.5 * plotValue.iqr;
        
        // Count outliers (points beyond whiskers)
        plotValue.outliers = 0;
        for (const auto& value : purityFeatureValueVec) {
            if (value < plotValue.lowerWhisker || value > plotValue.upperWhisker) {
                plotValue.outliers++;
            }
        }
            
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] Error in statistical calculations: " << e.what() << std::endl;
        exit(1);
        // Handle error case
    }

    return plotValue;
}

void SomaticVarCaller::ClassifyReadsByCase(std::vector<int> &tumorAllelePosVec, std::map<int, int> &NorCountPS, std::map<int, int> &hpCount, const HaplotagParameters &params, std::map<int, HP3_Info> &somaticPosInfo){
    
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

void SomaticVarCaller::SetFilterParamsWithPurity(SomaticFilterParaemter &somaticParams, double &tumorPurity){

    somaticParams.tumorPurity = tumorPurity;

    if (tumorPurity >= 0.9 && tumorPurity <= 1.0) {
        somaticParams.HapConsistency_ReadCount_Thr = 10;
        somaticParams.HapConsistency_VAF_Thr = 0.2;

        somaticParams.IntervalSnpCount_ReadCount_Thr = 10;
        somaticParams.IntervalSnpCount_VAF_Thr = 0.15;
    } else if (tumorPurity >= 0.7 && tumorPurity < 0.9) {
        somaticParams.HapConsistency_ReadCount_Thr = 10;
        somaticParams.HapConsistency_VAF_Thr = 0.15;

        somaticParams.IntervalSnpCount_ReadCount_Thr = 10;
        somaticParams.IntervalSnpCount_VAF_Thr = 0.12;
    } else if (tumorPurity >= 0.5 && tumorPurity < 0.7) {
        somaticParams.HapConsistency_ReadCount_Thr = 10;
        somaticParams.HapConsistency_VAF_Thr = 0.1;

        somaticParams.IntervalSnpCount_ReadCount_Thr = 8;
        somaticParams.IntervalSnpCount_VAF_Thr = 0.1;
    } else if (tumorPurity >= 0.3 && tumorPurity < 0.5) {
        somaticParams.HapConsistency_ReadCount_Thr = 6;
        somaticParams.HapConsistency_VAF_Thr = 0.08;

        somaticParams.IntervalSnpCount_ReadCount_Thr = 6;
        somaticParams.IntervalSnpCount_VAF_Thr = 0.08;
    } else if (tumorPurity >= 0.1 && tumorPurity < 0.3) {
        somaticParams.HapConsistency_ReadCount_Thr = 6;
        somaticParams.HapConsistency_VAF_Thr = 0.06;

        //not used interval snp count filter 
        somaticParams.IntervalSnpCount_ReadCount_Thr = 1;
        somaticParams.IntervalSnpCount_VAF_Thr = 0.00;
    } else if (tumorPurity > 0.0 && tumorPurity < 0.1) {
        somaticParams.HapConsistency_ReadCount_Thr = 5;
        somaticParams.HapConsistency_VAF_Thr = 0.04;

        //not used interval snp count filter 
        somaticParams.IntervalSnpCount_ReadCount_Thr = 1;
        somaticParams.IntervalSnpCount_VAF_Thr = 0.00;
    } else {
        std::cerr << "[Error] tumor purity is not in the range of 0.0 to 1.0: " << tumorPurity << std::endl;
        exit(1);
    }

    //messy read count threshold
    if ( 0.4 < tumorPurity && tumorPurity <= 1.0) {
        somaticParams.MessyReadCountThreshold = 3;
    }else{
        somaticParams.MessyReadCountThreshold = 2;
    }
}

void SomaticVarCaller::SomaticFeatureFilter(const SomaticFilterParaemter &somaticParams, std::map<int, RefAltSet> &currentChrVariants,const std::string &chr, std::map<int, HP3_Info> &somaticPosInfo, BamBaseCounter &NorBase, double& tumorPurity){
//calculate the information of the somatic positon 
    std::map<int, HP3_Info>::iterator somaticVarIter = somaticPosInfo.begin();
    while( somaticVarIter != somaticPosInfo.end()){

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
        (*somaticVarIter).second.base.nonDelAF = (float)tumAltCount / (float)(tumDepth - (*somaticVarIter).second.base.delCount);

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
            (*somaticVarIter).second.tumDelRatio = (float)tumDeletionCount / (float)tumDepth;
        }

        int somaticReadH1_1 = (*chrPosSomaticInfo)[chr][(*somaticVarIter).first].somaticReadHpCount[ReadHP::H1_1];
        int somaticReadH2_1 = (*chrPosSomaticInfo)[chr][(*somaticVarIter).first].somaticReadHpCount[ReadHP::H2_1];

        //filter threshold 
        // float OnlyHP3ReadRatioThreshold = 0.0;
        // float VAF_upper_threshold = 1.0;
        // float VAF_lower_threshold = 0.0;
        // float tumDeletionRatioThreshold = 1.0;
        // float tumLowMpqRatioThreshold = 1.0;

        //messy read filter parameter
        float messyReadRatioThreshold = 1.0;
        int readCountThreshold = 3;

        //haplotype consistency filter parameter
        float HapConsistency_VAF_Thr=0.2;
        int HapConsistency_ReadCount_Thr=10;

        //interval snp count filter parameter
        float IntervalSnpCount_VAF_Thr=0.15;
        int IntervalSnpCount_ReadCount_Thr=10;

        messyReadRatioThreshold = somaticParams.MessyReadRatioThreshold;
        readCountThreshold = somaticParams.MessyReadCountThreshold;

        HapConsistency_VAF_Thr = somaticParams.HapConsistency_VAF_Thr;
        HapConsistency_ReadCount_Thr = somaticParams.HapConsistency_ReadCount_Thr;

        IntervalSnpCount_VAF_Thr = somaticParams.IntervalSnpCount_VAF_Thr;
        IntervalSnpCount_ReadCount_Thr = somaticParams.IntervalSnpCount_ReadCount_Thr;

        //Stage 1 filter
        // if((*somaticVarIter).second.isNormalPosLowVAF == false && somaticParams.applyFilter){
        // // if((*somaticVarIter).second.isNormalPosLowVAF == false){
        //     somaticVarIter++;
        //     continue;
        // }

        // //stage 2 filter
        // if( (!somaticParams.applyFilter) ||
        //     (
        //     //((*somaticVarIter).second.OnlyHP3ReadRatio < OnlyHP3ReadRatioThreshold) && 
        //     //((*somaticVarIter).second.base.lowMpqReadRatio <= tumLowMpqRatioThreshold) &&
        //     ((*somaticVarIter).second.Mixed_HP_readRatio < messyReadRatioThreshold) && 
        //     ((*somaticVarIter).second.CaseReadCount > readCountThreshold) /*&&
        //     ((*somaticVarIter).second.base.VAF > VAF_lower_threshold) && 
        //     ((*somaticVarIter).second.base.VAF < VAF_upper_threshold)*/
        //     //((*somaticVarIter).second.tumDelRatio < tumDeletionRatioThreshold) &&
        //     ) ){
            
        //     // haplotype consistency filter
        //     if((somaticParams.applyFilter) && (*somaticVarIter).second.CaseReadCount <= 10 && (*somaticVarIter).second.base.VAF <= 0.2){
        //         if((somaticReadH1_1 > 0 && somaticReadH2_1 > 0)){
        //             somaticVarIter++;
        //             continue;
        //         }
        //     }     

        //     //z-score filter
        //     if((somaticParams.applyFilter) && (*somaticVarIter).second.CaseReadCount <= 10 && (*somaticVarIter).second.base.VAF <= 0.15){
        //         int intervalSnpCount = (*somaticVarIter).second.intervalSnpCount;
        //         float zScore = (*somaticVarIter).second.zScore;
        //         if(intervalSnpCount > 4 && zScore <= 2.0 && zScore >= 0.0){
        //             somaticVarIter++;
        //             continue;
        //         }
        //     }
        //     (*somaticVarIter).second.isHighConSomaticSNP = true;
        // }

        //filter out the somatic SNP
        (*somaticVarIter).second.isFilterOut = false;
        
        // 檢查所有過濾條件
        bool stage1_filtered = false;
        bool messy_read_filtered = false;
        bool read_count_filtered = false;

        float norVAF = NorBase.getVAF(chr, (*somaticVarIter).first);
        float norDepth = NorBase.getDepth(chr, (*somaticVarIter).first);
        
        //stage 1 filter
        if (!(*somaticVarIter).second.isNormalPosLowVAF || !(norVAF <= 0.1 && norDepth > 1)) {
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
        if ((*somaticVarIter).second.CaseReadCount <= HapConsistency_ReadCount_Thr && 
            (*somaticVarIter).second.base.VAF <= HapConsistency_VAF_Thr) {
            if (somaticReadH1_1 > 0 && somaticReadH2_1 > 0) {
                haplotype_filtered = true;
            }
        }
        
        // interval snp count filter check
        bool zscore_filtered = false;
        if ((*somaticVarIter).second.CaseReadCount <= IntervalSnpCount_ReadCount_Thr && 
            (*somaticVarIter).second.base.VAF <= IntervalSnpCount_VAF_Thr) {

            int intervalSnpCount = (*somaticVarIter).second.intervalSnpCount;
            float zScore = (*somaticVarIter).second.zScore;
            if (intervalSnpCount > 4 && zScore <= 2.0 && zScore >= 0.0) {
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

        somaticVarIter++;          
    }
}

void SomaticVarCaller::ShannonEntropyFilter(const std::string &chr, std::map<int, HP3_Info> &somaticPosInfo, std::map<int, RefAltSet> &currentChrVariants, std::string &ref_string){
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

void SomaticVarCaller::calculateIntervalZScore(bool &isStartPos, int &startPos, int &endPos, int &snpCount, denseSnpInterval &denseSnp, std::map<int, std::pair<int, denseSnpInterval>> &localDenseTumorSnpInterval){
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

void SomaticVarCaller::getDenseTumorSnpInterval(std::map<int, HP3_Info> &somaticPosInfo, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, std::map<std::string, int>> &tumorPosReadCorrBaseHP, std::map<int, std::pair<int, denseSnpInterval>> &localDenseTumorSnpInterval){
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
    denseSnpInterval denseSnp;

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
            std::map<int, int> NorCountPS = (*readTotalHPcountIter).second.NorCountPS;

            (*readTotalHPcountIter).second.hpResult = determineReadHP(hpCount, pqValue, NorCountPS, normalHPsimilarity, tumorHPsimilarity, params.percentageThreshold, nullptr, nullptr, nullptr);
            
            readTotalHPcountIter++;
        }
}

void SomaticVarCaller::StatisticSomaticPosReadHP(const std::string &chr, std::map<int, HP3_Info> &somaticPosInfo, std::map<int, std::map<std::string, int>> &tumorPosReadCorrBaseHP, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, ReadHpResult> &localReadHpDistri){
    std::map<int, HP3_Info>::iterator somaticVarIter = somaticPosInfo.begin();
    while(somaticVarIter != somaticPosInfo.end()){
        if((*somaticVarIter).second.isHighConSomaticSNP){
            int pos = (*somaticVarIter).first;

            if(tumorPosReadCorrBaseHP.find(pos) != tumorPosReadCorrBaseHP.end()){
                localReadHpDistri[pos] = ReadHpResult();

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
                    recordReadHp(pos, hpResult, baseHP, localReadHpDistri);

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

                localReadHpDistri[pos].existDeriveByH1andH2 = false;

                if(HP1_1ratio >= 1.0){
                    (*somaticVarIter).second.somaticReadDeriveByHP = SnpHP::GERMLINE_H1;
                }else if(HP2_1ratio >= 1.0){
                    (*somaticVarIter).second.somaticReadDeriveByHP = SnpHP::GERMLINE_H2;
                }else{
                    (*somaticVarIter).second.somaticReadDeriveByHP = SnpHP::NONE_SNP;

                    if((0 < HP1_1ratio && HP1_1ratio < 1.0)  || (0 < HP2_1ratio && HP2_1ratio < 1.0)){
                        localReadHpDistri[pos].existDeriveByH1andH2 = true;
                        // std::cerr << "ERROR (statistic all read HP) => somatic Variant had HP1-1 and HP2-1 reads :"<< std::endl;
                        // std::cerr << "chr: "<< chr << " pos: " << pos+1 <<std::endl;
                        // std::cerr << "HP1-1_ratio: "<< HP1_1ratio << " HP2-1_ratio: "<< HP2_1ratio << std::endl;
                        // std::cerr << "HP1-1: "<< deriveByHPfromBaseHp3["1-1"] << " HP2-1: "<< deriveByHPfromBaseHp3["2-1"] << " HP3: "<< readHpDistributed[pos].HP3read << " HP4: "<< readHpDistributed[pos].HP4read<< std::endl;
                    }
                }
                //record derive HP of current snp
                recordDeriveHp(pos, (*somaticVarIter).second.somaticReadDeriveByHP, 0.0, localReadHpDistri);

                //error: haven't exist somatic HP read in current position 
                if( localReadHpDistri[pos].readHpCounter[ReadHP::H3] == 0 && 
                    localReadHpDistri[pos].readHpCounter[ReadHP::H1_1] == 0 &&
                    localReadHpDistri[pos].readHpCounter[ReadHP::H2_1] == 0){
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


void SomaticVarCaller::WriteSomaticVarCallingLog(const HaplotagParameters &params, const SomaticFilterParaemter &somaticParams, const std::vector<std::string> &chrVec, BamBaseCounter &NorBase, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat){
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
                 << "##Messy read ratio threshold : " << somaticParams.MessyReadRatioThreshold << "\n"
                 << "##Messy read count threshold : " << somaticParams.MessyReadCountThreshold << "\n"
                 << "##Haplotag consistency filter threshold : " << somaticParams.HapConsistency_VAF_Thr << "\n"
                 << "##Haplotag consistency filter read count threshold : " << somaticParams.HapConsistency_ReadCount_Thr << "\n"
                 << "##Interval SNP count filter threshold : " << somaticParams.IntervalSnpCount_VAF_Thr << "\n"
                 << "##Interval SNP count filter read count threshold : " << somaticParams.IntervalSnpCount_ReadCount_Thr << "\n"
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
                tumMpqAltCount = (*chrPosSomaticInfo)[chr][(*somaticVarIter).first].base.MPQ_A_count;
                norAltCount = NorBase.getBaseAcount(chr, (*somaticVarIter).first);
            }else if(AltBase == "C"){
                tumMpqAltCount = (*chrPosSomaticInfo)[chr][(*somaticVarIter).first].base.MPQ_C_count;
                norAltCount = NorBase.getBaseCcount(chr, (*somaticVarIter).first);
            }else if(AltBase == "T"){
                tumMpqAltCount = (*chrPosSomaticInfo)[chr][(*somaticVarIter).first].base.MPQ_T_count;       
                norAltCount = NorBase.getBaseTcount(chr, (*somaticVarIter).first);
            }else if(AltBase == "G"){
                tumMpqAltCount = (*chrPosSomaticInfo)[chr][(*somaticVarIter).first].base.MPQ_G_count;
                norAltCount = NorBase.getBaseGcount(chr, (*somaticVarIter).first);
            }else{
                std::cerr << "Error(write tag HP3 log file) => can't match RefBase or AltBase : chr:" << chr << " pos: " << ((*somaticVarIter).first) + 1 << " RefBase:" << RefBase << " AltBase:" << AltBase;
                exit(1);
            }

            //calculating tumor VAF
            float tumVAF = (*somaticVarIter).second.base.VAF;
            float tumMpqVAF = (float)tumMpqAltCount / (float)tumMpqDepth;

            float norVAF =  0.0;
            if(norAltCount > 0) norVAF = (float)norAltCount / (float)norDepth;
            
            float norCountVAF = NorBase.getVAF(chr,(*somaticVarIter).first);
            float norMpqVAF = NorBase.getFilterdMpqVAF(chr,(*somaticVarIter).first);

            if(norVAF != norCountVAF){
                std::cerr << "Error(diff nor VAF) => can't find GTtype at chr: " << chr << " pos: " << ((*somaticVarIter).first) + 1;
                std::cerr << "norVAF: " << norVAF << " norCountVAF: " << norCountVAF << std::endl;
                exit(1); 
            }
            float norVAF_substract = (norMpqVAF - norCountVAF);

            float norMPQReadRatio = 0.0; 

            if(norDepth > 0) norMPQReadRatio = (float)(norDepth - norMpqDepth) / (float)norDepth;
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
            int tumDeletionCount = (*chrPosSomaticInfo)[chr][(*somaticVarIter).first].base.delCount;

            //the ratio of deletions occurring
            float norDeletionRatio;
            float tumDeletionRatio = (*somaticVarIter).second.tumDelRatio;

            if(norDeletionCount == 0){
                norDeletionRatio = 0.0;     
            }else{
                norDeletionRatio = (float)norDeletionCount / (float)norDepth;
            }

            // the distribution of reads HP at the current position
            int H1readCount = (*chrVarReadHpResult)[chr][(*somaticVarIter).first].readHpCounter[ReadHP::H1];
            int H2readCount = (*chrVarReadHpResult)[chr][(*somaticVarIter).first].readHpCounter[ReadHP::H2];
            int H1_1readCount = (*chrVarReadHpResult)[chr][(*somaticVarIter).first].readHpCounter[ReadHP::H1_1];
            int H2_1readCount = (*chrVarReadHpResult)[chr][(*somaticVarIter).first].readHpCounter[ReadHP::H2_1];
            int H3readCount = (*chrVarReadHpResult)[chr][(*somaticVarIter).first].readHpCounter[ReadHP::H3];

            // if (H1readCount != (*somaticVarIter).second.allReadHpCount[ReadHP::H1]){
            //     std::cerr << "Error(H1readCount) => chr: " << chr << " pos: " << ((*somaticVarIter).first) + 1 << " H1readCount: " << H1readCount << " allReadHpCount: " << (*somaticVarIter).second.allReadHpCount[ReadHP::H1];
            //     exit(1);
            // }else if (H2readCount != (*somaticVarIter).second.allReadHpCount[ReadHP::H2]){
            //     std::cerr << "Error(H2readCount) => chr: " << chr << " pos: " << ((*somaticVarIter).first) + 1 << " H2readCount: " << H2readCount << " allReadHpCount: " << (*somaticVarIter).second.allReadHpCount[ReadHP::H2];
            //     exit(1);
            // }else if (H1_1readCount != (*somaticVarIter).second.allReadHpCount[ReadHP::H1_1]){
            //     std::cerr << "Error(H1_1readCount) => chr: " << chr << " pos: " << ((*somaticVarIter).first) + 1 << " H1_1readCount: " << H1_1readCount << " allReadHpCount: " << (*somaticVarIter).second.allReadHpCount[ReadHP::H1_1];
            //     exit(1);
            // }else if (H2_1readCount != (*somaticVarIter).second.allReadHpCount[ReadHP::H2_1]){
            //     std::cerr << "Error(H2_1readCount) => chr: " << chr << " pos: " << ((*somaticVarIter).first) + 1 << " H2_1readCount: " << H2_1readCount << " allReadHpCount: " << (*somaticVarIter).second.allReadHpCount[ReadHP::H2_1];
            //     exit(1);
            // }else if (H3readCount != (*somaticVarIter).second.allReadHpCount[ReadHP::H3]){
            //     std::cerr << "Error(H3readCount) => chr: " << chr << " pos: " << ((*somaticVarIter).first) + 1 << " H3readCount: " << H3readCount << " allReadHpCount: " << (*somaticVarIter).second.allReadHpCount[ReadHP::H3];
            //     exit(1);
            // }
            
            // the ratio of the read count of H1 and H2
            int germlineReadHpCount = H1readCount + H2readCount;
            double germlineReadHpConsistencyRatio = 0.0;
            if(H1readCount > 0 && H2readCount > 0){         
                germlineReadHpConsistencyRatio = (H1readCount > H2readCount) ? ((double)H1readCount / (double)(H1readCount + H2readCount)) : ((double)H2readCount / (double)(H1readCount + H2readCount));
            }else if(H1readCount == 0 && H2readCount == 0){
                germlineReadHpConsistencyRatio = 0.0;
            }
            else{
                germlineReadHpConsistencyRatio = 1.0;
            }

            // the ratio of the germlineHP based on depth   
            double percentageOfGermlineHp = 0.0;
            if(tumDepth > 0 && germlineReadHpCount > 0){
                percentageOfGermlineHp = (double)germlineReadHpCount / (double)tumDepth;            
            }

            // the count of reads have somatic base at current position
            int somaticVarReadH1_1 = (*chrPosSomaticInfo)[chr][(*somaticVarIter).first].somaticReadHpCount[ReadHP::H1_1];
            int somaticVarReadH2_1 = (*chrPosSomaticInfo)[chr][(*somaticVarIter).first].somaticReadHpCount[ReadHP::H2_1];
            int somaticVarReadH3 = (*chrPosSomaticInfo)[chr][(*somaticVarIter).first].somaticReadHpCount[ReadHP::H3];
            int untaggedReadCount = (*chrPosSomaticInfo)[chr][(*somaticVarIter).first].somaticReadHpCount[ReadHP::unTag];

            // the count of reads based on germline HP
            int BaseOnH1ReadCount = H1readCount + H1_1readCount;
            int BaseOnH2ReadCount = H2readCount + H2_1readCount;


            // the ratio of the read count of H1 and H2 based on germline HP
            double baseOnGermlineReadHpConsistencyRatio = 0.0;
            if(BaseOnH1ReadCount > 0 && BaseOnH2ReadCount > 0){
                baseOnGermlineReadHpConsistencyRatio = (BaseOnH1ReadCount > BaseOnH2ReadCount) ? ((double)BaseOnH1ReadCount / (double)(BaseOnH1ReadCount + BaseOnH2ReadCount)) : ((double)BaseOnH2ReadCount / (double)(BaseOnH1ReadCount + BaseOnH2ReadCount));
            }else if(BaseOnH1ReadCount == 0 && BaseOnH2ReadCount == 0){
                baseOnGermlineReadHpConsistencyRatio = 0.0;
            }
            else{
                baseOnGermlineReadHpConsistencyRatio = 1.0;
            }

            double somaticReadHpConsistencyRatio = 0.0;
            if(H1_1readCount > 0 && H2_1readCount > 0){     
                somaticReadHpConsistencyRatio = (H1_1readCount > H2_1readCount) ? ((double)H1_1readCount / (double)(H1_1readCount + H2_1readCount)) : ((double)H2_1readCount / (double)(H1_1readCount + H2_1readCount));
            }else if(H1_1readCount == 0 && H2_1readCount == 0){
                somaticReadHpConsistencyRatio = 0.0;
            }
            else{
                somaticReadHpConsistencyRatio = 1.0;
            }

            //read hp count in the normal bam
            int H1readCountInNorBam = NorBase.getReadHpCountInNorBam(chr, (*somaticVarIter).first, ReadHP::H1);
            int H2readCountInNorBam = NorBase.getReadHpCountInNorBam(chr, (*somaticVarIter).first, ReadHP::H2);
            int germlineReadHpCountInNorBam = H1readCountInNorBam + H2readCountInNorBam;

            double germlineReadHpConsistencyRatioInNorBam = 0.0;
            if(H1readCountInNorBam > 0 && H2readCountInNorBam > 0){
                germlineReadHpConsistencyRatioInNorBam = (H1readCountInNorBam > H2readCountInNorBam) ? ((double)H1readCountInNorBam / (double)germlineReadHpCountInNorBam) : ((double)H2readCountInNorBam / (double)germlineReadHpCountInNorBam);
            }else if(H1readCountInNorBam == 0 && H2readCountInNorBam == 0){
                germlineReadHpConsistencyRatioInNorBam = 0.0;
            }
            else{
                germlineReadHpConsistencyRatioInNorBam = 1.0;
            }

            double percentageOfGermlineHpInNorBam = 0.0;
            if(norDepth > 0 && germlineReadHpCountInNorBam > 0){
                percentageOfGermlineHpInNorBam = (double)germlineReadHpCountInNorBam / (double)norDepth;
            }

            //difference of the germline consistency ratio in tumor and normal bam
            double germlineReadHpConsistencyRatioDifference = 0.0;
            germlineReadHpConsistencyRatioDifference = germlineReadHpConsistencyRatio - germlineReadHpConsistencyRatioInNorBam;
            
            //difference of the percentage of germline HP in tumor and normal bam
            double percentageOfGermlineHpDifference = 0.0;
            percentageOfGermlineHpDifference = percentageOfGermlineHp - percentageOfGermlineHpInNorBam;

            // If the filter is not applied, the altCount will not equal the sum of somaticRead 
            //(somaticRead H1_1 to H3 does not include SNPs present in both tumor and normal positions)
            // if(somaticReadH1_1 + somaticReadH2_1 + somaticReadH3 + somaticReadUnTag != tumMpqAltCount){
            //     std::cerr << "Error (SomaticFeatureFilter) => altCount is not equal to the sum of somaticRead : chr:" << chr << " pos: " << (*somaticVarIter).first + 1 << std::endl;
            //     std::cerr << "altCount: " << tumMpqAltCount << " somaticReadH1-1: " << somaticReadH1_1 << " somaticReadH2-1: " << somaticReadH2_1 << " somaticReadH3: " << somaticReadH3 << " somaticReadUnTag: " << somaticReadUnTag << std::endl;
            //     exit(1);
            // }

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

            float norNonDelAF = NorBase.getNoDelAF(chr, (*somaticVarIter).first);
            float tumNonDelAF = (*somaticVarIter).second.base.nonDelAF;
            
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
                        << germlineReadHpConsistencyRatio << "\t" //41
                        << somaticReadHpConsistencyRatio << "\t" //42
                        << baseOnGermlineReadHpConsistencyRatio << "\t" //43
                        << percentageOfGermlineHp << "\t" //44
                        << H1readCountInNorBam << "\t" //45
                        << H2readCountInNorBam << "\t" //46
                        << germlineReadHpCountInNorBam << "\t" //47
                        << germlineReadHpConsistencyRatioInNorBam << "\t" //48
                        << percentageOfGermlineHpInNorBam << "\t" //49
                        << germlineReadHpConsistencyRatioDifference << "\t" //50
                        << percentageOfGermlineHpDifference << "\t" //51
                        << somaticVarReadH1_1 << "\t" //52
                        << somaticVarReadH2_1 << "\t" //53
                        << somaticVarReadH3 << "\t" //54
                        << untaggedReadCount << "\t" //55
                        << (*somaticVarIter).second.MeanAltCountPerVarRead << "\t" //56
                        << zScore << "\t" //57
                        << (*somaticVarIter).second.intervalSnpCount << "\t" //58
                        << mergedChrVarinat[chr][(*somaticVarIter).first].isExistNormal << "\t" //59
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

void SomaticVarCaller::WriteOtherSomaticHpLog(const HaplotagParameters &params, const std::vector<std::string> &chrVec, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat){
    std::ofstream *OtherHpSomaticVarLog=NULL;
    std::string logPosfix = "_otherHpSomaticVar.log";
    OtherHpSomaticVarLog=new std::ofstream(params.resultPrefix + logPosfix);

    int totalOtherSomaticHpVar = 0;
    for(auto chr: chrVec){
        for(auto somaticVar: (*chrPosSomaticInfo)[chr]){
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
        for(auto somaticVar: (*chrPosSomaticInfo)[chr]){
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

std::map<std::string, std::map<int, HP3_Info>> SomaticVarCaller::getSomaticChrPosInfo(){
    return (*chrPosSomaticInfo);
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
        mergedChrVarinat[chr][pos].Variant[Genome::HIGH_CON_SOMATIC] = tmp;
        mergedChrVarinat[chr][pos].isExistHighConSomatic = true;
    }
}

void highConBenchmark::recordDelReadCount(const std::string &chr, std::map<int, RefAltSet>::iterator &currentVariantIter){
    if(!openTestingFunc) return;
    
    if(currentVariantIter->second.isExistHighConSomatic){
        int pos = (*currentVariantIter).first;
        posAltRefDelCount[chr][pos].delCount++;

        //record somatic position for record crossing high con snp read
        highConSomaticPos.push_back(std::make_pair(pos, SnpHP::NONE_SNP));
    }
}

void highConBenchmark::recordRefAltAlleleCount(const std::string &chr, std::string &base, std::map<int, RefAltSet>::iterator &currentVariantIter){
    if(!openTestingFunc) return;

    if(currentVariantIter->second.isExistHighConSomatic){
        int pos = currentVariantIter->first;
        std::string refAllele = currentVariantIter->second.Variant[Genome::HIGH_CON_SOMATIC].Ref;
        std::string altAllele = currentVariantIter->second.Variant[Genome::HIGH_CON_SOMATIC].Alt;

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

somaticReadLog highConBenchmark::createBasicSomaticReadLog(const std::string &chr, std::string &readID, std::string &hpResult, double &norHPsimilarity, float &deriveByHpSimilarity, std::map<int, int> &hpCount){
    somaticReadLog tmp;
    tmp.chr = chr;
    tmp.readID = readID;
    tmp.hpResult = hpResult;
    tmp.germlineVarSimilarity = norHPsimilarity;
    tmp.deriveByHpSimilarity = deriveByHpSimilarity;
    tmp.germlineSnpCount = hpCount[SnpHP::GERMLINE_H1] + hpCount[SnpHP::GERMLINE_H2];
    tmp.tumorSnpCount = hpCount[SnpHP::SOMATIC_H3];

    return tmp;
}

void highConBenchmark::recordCrossingHighConSnpRead(const std::string &chr, std::string &readID, std::string &hpResult, std::map<int, int> &variantsHP, std::map<int, int> &hpCount, double &norHPsimilarity, float &deriveByHpSimilarity, std::map<int, RefAltSet> &currentChrVariants){
    // if not open testing function, return
    if(!openTestingFunc) return;

    somaticReadLog tmp = createBasicSomaticReadLog(chr, readID, hpResult, norHPsimilarity, deriveByHpSimilarity, hpCount);

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

void highConBenchmark::recordTaggedSomaticRead(const std::string &chr, std::string &readID, std::string &hpResult, std::map<int, int> &variantsHP, std::map<int, int> &hpCount, double &norHPsimilarity, float &deriveByHpSimilarity, std::map<int, RefAltSet> &currentChrVariants){
    // if not open testing function, return
    if(!openTestingFunc) return;

    somaticReadLog tmp = createBasicSomaticReadLog(chr, readID, hpResult, norHPsimilarity, deriveByHpSimilarity, hpCount);

    bool readExistHighConSomatic = false;

    auto varIter = variantsHP.begin();
    while(varIter != variantsHP.end()){
        int pos = varIter->first;
        int snpHP = varIter->second;
        if(currentChrVariants.find(pos) != currentChrVariants.end()){
            if(currentChrVariants[pos].isExistHighConSomatic && (snpHP == SnpHP::SOMATIC_H3 || snpHP == SnpHP::SOMATIC_H4)){
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

void highConBenchmark::recordTaggedRead(const std::string &chr, std::string &readID, std::string &hpResult, std::map<int, int> &variantsHP, std::map<int, int> &hpCount, double &norHPsimilarity, float &deriveByHpSimilarity, std::map<int, RefAltSet> &currentChrVariants){
    // if not open testing function, return
    if(!openTestingFunc) return;
    
    if(hpResult != "."){
        somaticReadLog tmp = createBasicSomaticReadLog(chr, readID, hpResult, norHPsimilarity, deriveByHpSimilarity, hpCount);

        auto varIter = variantsHP.begin();
        while(varIter != variantsHP.end()){
            int pos = varIter->first;
            int snpHP = varIter->second;
            if(currentChrVariants.find(pos) != currentChrVariants.end()){
                if(currentChrVariants[pos].isExistHighConSomatic && (snpHP == SnpHP::SOMATIC_H3 || snpHP == SnpHP::SOMATIC_H4)){
                    tmp.somaticSnpHp[pos] = snpHP;
                }
            }
            varIter++;
        }

        totalReadVec.push_back(tmp);
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
        (*refAltCountLog) << "##High confidence VCF:"  << params.benchmarkVcf << "\n";
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
                              << mergedChrVarinat[chr][posIter.first].Variant[Genome::HIGH_CON_SOMATIC].Ref << "\t"
                              << mergedChrVarinat[chr][posIter.first].Variant[Genome::HIGH_CON_SOMATIC].Alt << "\t"
                              << posIter.second.refCount << "\t"
                              << posIter.second.altCount << "\t"
                              << posIter.second.delCount << "\n";
        }
    }

    refAltCountLog->close();
    delete refAltCountLog;
    refAltCountLog = nullptr;
}

void highConBenchmark::writeTaggedReadLog(HaplotagParameters &params, std::string logPosfix){
    // if not open testing function, return
    if(!openTestingFunc) return;
    writeReadLog(params, logPosfix, totalReadVec);
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

    int totalReads = totalReadVec.size();

    if(!somaticReadLog->is_open()){
        std::cerr<< "Fail to open write file: " << params.resultPrefix + logPosfix << "\n";
        exit(1);
    }else{
        (*somaticReadLog) << "#####################\n";
        (*somaticReadLog) << "# Somatic Reads Log #\n";
        (*somaticReadLog) << "#####################\n";
        (*somaticReadLog) << "##High confidence VCF: "  << params.benchmarkVcf << "\n";
        (*somaticReadLog) << "##MappingQualityThreshold: "  << params.qualityThreshold << "\n";
        (*somaticReadLog) << "##Tatal reads: "  << totalReads << "\n";
        (*somaticReadLog) << "##Tatal tagged somatic reads: "  << totalTaggedSomaticReads << "\n";
        (*somaticReadLog) << "##Tatal truth somatic reads: "  << totalTruthSomaticReads << "\n";
        (*somaticReadLog) << "##Truth somatic read ratio: "  << (float)totalTaggedSomaticReads / (float)totalTruthSomaticReads << "\n";
        (*somaticReadLog) << "#CHROM\t"
                          << "ReadID\t"
                          << "germlineVarSimilarity\t"
                          << "deriveByHpSimilarity\t"
                          << "germlineSnpCount\t"
                          << "tumorSnpCount\t"
                          << "Haplotype\t"
                          << "somaticVariant,HP\n";
    }

    for(auto somaticRead: somaticReadVec){

        (*somaticReadLog) << somaticRead.chr << "\t"
                          << somaticRead.readID << "\t"
                          << somaticRead.germlineVarSimilarity << "\t"
                          << somaticRead.deriveByHpSimilarity << "\t"
                          << somaticRead.germlineSnpCount << "\t"
                          << somaticRead.tumorSnpCount << "\t"
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
            if(chrVariantIter->second.isExistHighConSomatic){
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

    std::vector<int> last_pos;
    // record reference last variant pos
    germlineGetRefLastVarPos(last_pos, *chrVec, vcfSet, geneType);
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
            (*tagResult) << "#ReadID\t"
                         << "CHROM\t"
                         << "ReadStart\t"
                         << "Confidnet(%)\t";
            if(tagTumorMode)
            (*tagResult) << "deriveByHpSimilarity\t";
            (*tagResult) << "Haplotype\t"
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
        highConSomaticData.writeTaggedReadLog(params, "_totalRead.out");
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

                auto norVar = (*currentVariantIter).second.Variant[Genome::NORMAL];

                int offset = (*currentVariantIter).first - ref_pos;

                if( offset < 0){
                }
                else{
                    uint8_t *q = bam_get_seq(&aln);
                    char base_chr = seq_nt16_str[bam_seqi(q,query_pos + offset)];
                    std::string base(1, base_chr);

                    germlineJudgeSnpHap(chrName, vcfSet, norVar, base, ref_pos, length, i, aln_core_n_cigar, cigar, currentVariantIter, hp1Count, hp2Count, variantsHP, countPS);

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

            germlineJudgeDeletionHap(chrName, ref_string, ref_pos, length, query_pos, currentVariantIter, vcfSet, &aln, hp1Count, hp2Count, variantsHP, countPS);
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

    // get the number of SVs occurring on different haplotypes in a read
    germlineJudgeSVHap(aln, vcfSet, hp1Count, hp2Count, tagGeneType);

    double min,max;
    int hpResult = ReadHP::unTag;

    // determine the haplotype of the read
    hpResult = germlineDetermineReadHap(hp1Count, hp2Count, min, max, percentageThreshold, pqValue, psValue, countPS, &totalHighSimilarity, &totalWithOutVaraint);
     
    //write tag log file
    if(tagResult != nullptr){
        writeGermlineTagLog(*tagResult, aln, bamHdr, hpResult, max, min, hp1Count, hp2Count, pqValue, variantsHP, countPS);
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

    // poition, <BaseHP,deriveByHp>
    std::map<int, std::pair<int , int>> somaticVarDeriveHP;

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
                            int BaseHp = SnpHP::NONE_SNP;
                            // if(base == (*currentVariantIter).second.Variant[Genome::TUMOR].Alt){
                            //     BaseHp = SnpHP::SOMATIC_H3;
                            // }
                            if(variantsHP.find((*currentVariantIter).first) != variantsHP.end()){
                                if(variantsHP[(*currentVariantIter).first] == SnpHP::SOMATIC_H3){
                                    BaseHp = SnpHP::SOMATIC_H3;
                                }
                            }
                            somaticVarDeriveHP[(*currentVariantIter).first] = std::make_pair(BaseHp, deriveByHp);
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
    std::string readID_test = bam_get_qname(&aln);

    //Record read HP result before correction hp result
    if(!somaticVarDeriveHP.empty()){
        for(auto somaticVarIter : somaticVarDeriveHP){
            int pos = somaticVarIter.first;
            int baseHP = somaticVarIter.second.first;
            int deriveHP = somaticVarIter.second.second;
            //if current read have somatic SNP then record the somatic SNP
            // if(variantsHP.find(pos) != variantsHP.end()){
            //     baseHP = variantsHP[pos];
            // }
            recordReadHp(pos, hpResult, baseHP, (*beforeCorrReadHpResult)[chrName]);
            recordDeriveHp(pos, deriveHP, 0.0, (*beforeCorrReadHpResult)[chrName]);
        }
    }

    //correction read HP result by deriveByHp
    float deriveByHpSimilarity = 0.0;
 
    if(hpResult == ReadHP::H3){
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
    }

    if(hpCount[1] == 0 && hpCount[2] == 0){
        if(hpResult == ReadHP::H3 && hpCount[3] != 0){
            totalreadOnlyH3Snp++;
        }
    }

    //Record read HP result after correction hp result
    if(!somaticVarDeriveHP.empty()){
        for(auto somaticVarIter : somaticVarDeriveHP){
            int pos = somaticVarIter.first;
            int baseHP = somaticVarIter.second.first;
            int deriveHP = somaticVarIter.second.second;
            //if current read have somatic SNP then record the somatic SNP
            // if(variantsHP.find(pos) != variantsHP.end()){
            // baseHP = variantsHP[pos];
            // }
            recordReadHp(pos, hpResult, baseHP, (*afterCorrReadHpResult)[chrName]);
            recordDeriveHp(pos, deriveHP, deriveByHpSimilarity , (*afterCorrReadHpResult)[chrName]);

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
            // if(tumCountPS.size() != 0){
            //     psIter = tumCountPS.begin();
            //     psResultStr = std::to_string((*psIter).first);
            //     psValue = (*psIter).first;    
            // }
            if(norCountPS.size() != 0){
                psIter = norCountPS.begin();
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
    highConSomaticData.recordTaggedRead(chrName, readID , hpResultStr, variantsHP, hpCount, norHPsimilarity, deriveByHpSimilarity, currentChrVariants);
    highConSomaticData.recordTaggedSomaticRead(chrName, readID , hpResultStr, variantsHP, hpCount, norHPsimilarity, deriveByHpSimilarity, currentChrVariants);
    highConSomaticData.recordCrossingHighConSnpRead(chrName, readID, hpResultStr, variantsHP, hpCount, norHPsimilarity, deriveByHpSimilarity, currentChrVariants);

    //write tag log file
    if(tagResult!=NULL){

        (*tagResult)<< readID                                       << "\t"
                    << bamHdr.target_name[aln.core.tid]             << "\t"
                    << aln.core.pos                                 << "\t"
                    << norHPsimilarity                              << "\t"
                    << deriveByHpSimilarity                         << "\t"
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
            (*tagResult)<< " " << (v.first + 1) << "," << v.second ;
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

void HaplotagProcess::OnlyTumorSNPjudgeHP(const std::string &chrName, int &curPos, RefAltSet &curVar, std::string base, VCF_Info *vcfSet, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *tumorAllelePosVec, BamBaseCounter *NorBase, std::map<int, HP3_Info> *SomaticPos){

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
    std::cerr<< "tumor SNP file:        " << params.tumorSnpFile        << "\n";  
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
    std::cerr<< "enable somatic variant filter:   " << (params.enableFilter ? "true" : "false") << "\n"; 
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
        if(params.benchmarkVcf != ""){
            std::time_t begin = time(NULL);
            std::cerr<< "loading high confidence SNP ... ";
            highConSomaticData.setTestingFunc(true);
            highConSomaticData.loadHighConSomatic(params.benchmarkVcf, vcfSet[Genome::HIGH_CON_SOMATIC], *mergedChrVarinat);
            std::cerr<< difftime(time(NULL), begin) << "s\n";
            highConSomaticData.displaySomaticVarCount(vcfSet[Genome::HIGH_CON_SOMATIC].chrVec, *mergedChrVarinat);
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

    int tumor_snp_count = 0;
    int normal_snp_count = 0;
    int overlap_snp_count = 0;
    for(auto& chrIter : (*chrVec)){
        auto chrVarIter = (*mergedChrVarinat)[chrIter].begin();
        while(chrVarIter != (*mergedChrVarinat)[chrIter].end()){
            if((*chrVarIter).second.isExistTumor){
                tumor_snp_count++;
            }
            if((*chrVarIter).second.isExistNormal){
                normal_snp_count++;
            }
            if((*chrVarIter).second.isExistTumor && (*chrVarIter).second.isExistNormal){
                overlap_snp_count++;
            }
            chrVarIter++;
        }
    }
    std::cerr << "Normal SNP count: " << normal_snp_count << std::endl;
    std::cerr << "Tumor SNP count: " << tumor_snp_count << std::endl;
    std::cerr << "Overlap SNP count: " << overlap_snp_count << std::endl;

    //somatic SNPs calling
    if(tagTumorMode){

        //Count each base numbers at tumor SNP position in the Normal.bam
        BamBaseCounter *NorBase = new BamBaseCounter(params.enableFilter);
        NorBase->CountingBamBase(params.bamFile, params, (*mergedChrVarinat), *chrVec, *chrLength, vcfSet, Genome::NORMAL);

        //record the HP3 confidence of each read
        SomaticVarCaller *SomaticVar = new SomaticVarCaller();
        SomaticVar->VariantCalling(params.tumorBamFile, (*mergedChrVarinat), *chrVec, *chrLength, params, vcfSet, *NorBase);
        (*chrPosReadCase) = SomaticVar->getSomaticChrPosInfo();

        delete NorBase;
        delete SomaticVar;
        NorBase = nullptr;
        SomaticVar = nullptr;
        return;
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



