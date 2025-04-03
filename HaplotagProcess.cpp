#include "HaplotagProcess.h"


void HaplotagProcess::tagRead(HaplotagParameters &params, const Genome& geneType){

    // input file management
    std::string openBamFile = params.bamFile;

    if(geneType == TUMOR){
        openBamFile = params.tumorBamFile;
    }else if(geneType == NORMAL){
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
    getLastVarPos(last_pos, *chrVec, *mergedChrVarinat, geneType);
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
        std::map<int, MultiGenomeVar>::reverse_iterator last = currentChrVariants.rbegin();
        
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
                    haplotype = somaticJudgeHaplotype(*bamHdr, *aln, chr, params.percentageThreshold, tagResult, pqValue, psValue, geneType, chr_reference);
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

    std::map<int, int> hpCount;
    hpCount[SnpHP::GERMLINE_H1] = 0;
    hpCount[SnpHP::GERMLINE_H2] = 0;

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
    std::map<int, MultiGenomeVar>::iterator currentVariantIter = firstVariantIter;

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

                auto norVar = (*currentVariantIter).second.Variant[NORMAL];

                int offset = (*currentVariantIter).first - ref_pos;

                if( offset < 0){
                }
                else{
                    uint8_t *q = bam_get_seq(&aln);
                    char base_chr = seq_nt16_str[bam_seqi(q,query_pos + offset)];
                    std::string base(1, base_chr);

                    germlineJudgeSnpHap(chrName, norVar, base, ref_pos, length, i, aln_core_n_cigar, cigar, currentVariantIter, hpCount, variantsHP, countPS);

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

            germlineJudgeDeletionHap(chrName, ref_string, ref_pos, length, query_pos, currentVariantIter, &aln, hpCount, variantsHP, countPS);
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
    germlineJudgeSVHap(aln, vcfSet, hpCount, tagGeneType);

    double min,max;
    int hpResult = ReadHP::unTag;

    // determine the haplotype of the read
    hpResult = germlineDetermineReadHap(hpCount, min, max, percentageThreshold, pqValue, psValue, countPS, &totalHighSimilarity, &totalWithOutVaraint);
     
    //write tag log file
    if(tagResult != nullptr){
        writeGermlineTagLog(*tagResult, aln, bamHdr, hpResult, max, min, hpCount, pqValue, variantsHP, countPS);
    }

    return hpResult;
}

int HaplotagProcess::somaticJudgeHaplotype(const bam_hdr_t &bamHdr,const bam1_t &aln, const std::string &chrName, double percentageThreshold, std::ofstream *tagResult, int &pqValue, int &psValue, const int tagGeneType, std::string &ref_string){

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
    std::map<int, MultiGenomeVar>::iterator currentVariantIter = firstVariantIter;

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
                    SomaticJudgeSnpHP(currentVariantIter, chrName, base, hpCount, norCountPS, tumCountPS, &variantsHP, nullptr, &((*chrPosReadCase)[chrName]));
                    if((*chrPosReadCase)[chrName].find((*currentVariantIter).first) != (*chrPosReadCase)[chrName].end()){

                        //record the somatic snp derive by which germline hp in this read
                        if((*chrPosReadCase)[chrName][(*currentVariantIter).first].isHighConSomaticSNP){
                            int deriveByHp = (*chrPosReadCase)[chrName][(*currentVariantIter).first].somaticReadDeriveByHP;
                            int BaseHp = SnpHP::NONE_SNP;
                            // if(base == (*currentVariantIter).second.Variant[TUMOR].Alt){
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
    auto readIter = vcfSet[NORMAL].readSVHapCount.find(bam_get_qname(&aln));
    if( readIter != vcfSet[NORMAL].readSVHapCount.end() ){
        hpCount[1] += vcfSet[NORMAL].readSVHapCount[bam_get_qname(&aln)][0];
        hpCount[2] += vcfSet[NORMAL].readSVHapCount[bam_get_qname(&aln)][1];
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

void HaplotagProcess::OnlyTumorSNPjudgeHP(const std::string &chrName, int &curPos, MultiGenomeVar &curVar, std::string base, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *tumorAllelePosVec, std::map<int, HP3_Info> *SomaticPos){

    if(SomaticPos == nullptr){
        std::cerr << "ERROR (SomaticTaggingJudgeHP) => SomaticPos pointer cannot be nullptr"<< std::endl;
        exit(1);
    }

    if((*SomaticPos).find(curPos) != (*SomaticPos).end()){
        if((*SomaticPos)[curPos].isHighConSomaticSNP){
            std::string TumorRefBase = curVar.Variant[TUMOR].allele.Ref;
            std::string TumorAltBase = curVar.Variant[TUMOR].allele.Alt; 

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
}


HaplotagProcess::HaplotagProcess():
totalAlignment(0),totalSupplementary(0),totalSecondary(0),totalUnmapped(0),totalTagCount(0),totalUnTagCount(0),processBegin(time(NULL))
{
    //initialize variable
    chrVec = nullptr;
    chrLength = nullptr;

    vcfSet[NORMAL].gene_type = NORMAL;
    vcfSet[TUMOR].gene_type = TUMOR;

    mergedChrVarinat = new std::map<std::string, std::map<int, MultiGenomeVar>>();;
    chrPosReadCase = new std::map<std::string, std::map<int, HP3_Info>>();

    beforeCorrReadHpResult = new std::map<std::string, std::map<int, ReadHpResult>>();
    afterCorrReadHpResult = new std::map<std::string, std::map<int, ReadHpResult>>();

    hpBeforeInheritance = new ReadHpDistriLog();
    hpAfterInheritance = new ReadHpDistriLog();

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


void HaplotagProcess::taggingProcess(HaplotagParameters &params)
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
    Genome tagGeneType;

    if(tagTumorMode){
        tagGeneType = TUMOR;
    }else{
        tagGeneType = NORMAL;
    }

    VcfParser vcfParser(tagTumorMode);

    if(tagTumorMode){
        //load seqc high con file for benchmarking
        if(params.benchmarkVcf != ""){
            std::time_t begin = time(NULL);
            std::cerr<< "loading high confidence SNP ... ";
            highConSomaticData.setTestingFunc(true);
            highConSomaticData.loadHighConSomatic(params.benchmarkVcf, vcfSet[HIGH_CON_SOMATIC], *mergedChrVarinat);
            std::cerr<< difftime(time(NULL), begin) << "s\n";
            highConSomaticData.displaySomaticVarCount(vcfSet[HIGH_CON_SOMATIC].chrVec, *mergedChrVarinat);
        }
        //load tumor snp vcf
        if(params.tumorSnpFile != ""){
            std::time_t begin = time(NULL);
            std::cerr<< "parsing tumor SNP VCF ... ";
            vcfParser.setParseSnpFile(true);
            vcfParser.variantParser(params.tumorSnpFile, vcfSet[TUMOR], *mergedChrVarinat);
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
    vcfParser.variantParser(params.snpFile, vcfSet[NORMAL], *mergedChrVarinat);
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
    if(tagGeneType == NORMAL){
        chrVec = &(vcfSet[NORMAL].chrVec);
        chrLength = &(vcfSet[NORMAL].chrLength); 
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

        chrVec = &(vcfSet[TUMOR].chrVec);
        chrLength = &(vcfSet[TUMOR].chrLength); 
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

    //somatic SNPs calling
    if(tagTumorMode){

        //Count each base numbers at tumor SNP position in the Normal.bam
        BamBaseCounter *NorBase = new BamBaseCounter(*chrVec);
        // NorBase->CountingBamBase(params.bamFile, params, (*mergedChrVarinat), *chrVec, *chrLength, vcfSet, NORMAL);
        NorBase->extractNormalData(params.bamFile, params, (*mergedChrVarinat), *chrVec, *chrLength, vcfSet, NORMAL);

        //record the HP3 confidence of each read
        SomaticVarCaller *SomaticVar = new SomaticVarCaller(*chrVec, params);
        SomaticVar->VariantCalling(params.tumorBamFile, (*mergedChrVarinat), *chrVec, *chrLength, params, *NorBase);
        (*chrPosReadCase) = SomaticVar->getSomaticChrPosInfo();

        delete NorBase;
        delete SomaticVar;
        NorBase = nullptr;
        SomaticVar = nullptr;
        // return;
    }

    // tag read
    begin = time(NULL);

    if(tagGeneType == TUMOR){
        std::cerr<< "somatic tagging start ...\n";
    }else{
        std::cerr<< "tag read start ...\n";
    }
    tagRead(params, tagGeneType);

    std::cerr<< "tag read " << difftime(time(NULL), begin) << "s\n";

    return;
};



