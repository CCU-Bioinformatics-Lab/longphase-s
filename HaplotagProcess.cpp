#include "HaplotagProcess.h"


void HaplotagProcess::tagRead(HaplotagParameters &params, const Genome& geneType){

    // input file management
    std::string openBamFile = params.bamFile;

    if(geneType == TUMOR){
        openBamFile = params.tumorBamFile;
    }else if(geneType == NORMAL){
        openBamFile = params.bamFile;
    }

    std::vector<int> last_pos;
    // record reference last variant pos
    getLastVarPos(last_pos, *chrVec, *mergedChrVarinat, geneType);
    // reference fasta parser
    FastaParser fastaParser(params.fastaFile, *chrVec, last_pos, params.numThreads);

    // init data structure and get core n
    htsThreadPool threadPool = {NULL, 0};
    // creat thread pool
    if (!(threadPool.pool = hts_tpool_init(params.numThreads))) {
        fprintf(stderr, "Error creating thread pool\n");
    }

    bool writeOutputBam = true;
    BamFileRAII bamRAII(openBamFile, params.fastaFile, threadPool, params, writeOutputBam);

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

        int result = 0;
        
        // tagging will be attempted for reads within the specified coordinate region.
        std::string region = !params.region.empty() ? params.region : chr + ":1-" + std::to_string((*chrLength)[chr]);
        hts_itr_t *iter = sam_itr_querys(bamRAII.idx, bamRAII.bamHdr, region.c_str());
        // iter all reads
        while ((sam_itr_multi_next(bamRAII.in, iter, bamRAII.aln)) >= 0) {
            totalAlignment++;
            int flag = bamRAII.aln->core.flag;
            
            if ( bamRAII.aln->core.qual < params.qualityThreshold ){
                // mapping quality is lower than threshold
                totalUnTagCount++;
                totalLowerQuality++;
                continue;
            }
            
            if( (flag & 0x4) != 0 ){
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
            else if(int(bamRAII.aln->core.pos) <= (*last).first){
                if( (flag & 0x800) != 0 ){
                    totalSupplementary++;
                }

                int pqValue = 0;
                int psValue = 0; 
                int haplotype = ReadHP::unTag;

                if(tagTumorMode){
                    haplotype = somaticJudgeHaplotype(*bamRAII.bamHdr, *bamRAII.aln, chr, params.percentageThreshold, tagResult, pqValue, psValue, geneType, chr_reference, params);
                }else{
                    haplotype = judgeHaplotype(*bamRAII.bamHdr, *bamRAII.aln, chr, params.percentageThreshold, tagResult, pqValue, psValue, geneType, chr_reference, params);
                }

                initFlag(bamRAII.aln, "HP");
                initFlag(bamRAII.aln, "PS");
                initFlag(bamRAII.aln, "PQ");

                if (haplotype != ReadHP::unTag){

                    totalHpCount[haplotype]++;
                    totalTagCount++;

                    if(tagTumorMode){
                        std::string haplotype_str = "";
                        haplotype_str = convertHpResultToString(haplotype);

                        bam_aux_append(bamRAII.aln, "HP", 'Z', (haplotype_str.size() + 1), (const uint8_t*) haplotype_str.c_str());
                        if(psValue != -1) bam_aux_append(bamRAII.aln, "PS", 'i', sizeof(psValue), (uint8_t*) &psValue);
                        bam_aux_append(bamRAII.aln, "PQ", 'i', sizeof(pqValue), (uint8_t*) &pqValue);
                    }else{
                        bam_aux_append(bamRAII.aln, "HP", 'i', sizeof(haplotype), (uint8_t*) &haplotype);
                        bam_aux_append(bamRAII.aln, "PS", 'i', sizeof(psValue), (uint8_t*) &psValue);
                        bam_aux_append(bamRAII.aln, "PQ", 'i', sizeof(pqValue), (uint8_t*) &pqValue);
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
            result = sam_write1(bamRAII.out, bamRAII.bamHdr, bamRAII.aln);
        }

        if (tagTumorMode && params.writeReadLog){
            hpBeforeInheritance->mergeLocalReadHp(chr, (*beforeCorrReadHpResult)[chr]);
            hpAfterInheritance->mergeLocalReadHp(chr, (*afterCorrReadHpResult)[chr]);
        }
        hts_itr_destroy(iter);  
        std::cerr<< difftime(time(NULL), begin) << "s\n";
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

    //must destroy bamRAII before destroying threadPool
    bamRAII.destroy();
    hts_tpool_destroy(threadPool.pool);

    return;
}

void HaplotagProcess::initFlag(bam1_t *aln, std::string flag){

    uint8_t *hpTag = bam_aux_get(aln, flag.c_str() );

    if( hpTag != NULL )
        bam_aux_del(aln, hpTag);

    return;
}

int HaplotagProcess::judgeHaplotype(const bam_hdr_t &bamHdr,const bam1_t &aln, std::string chrName, double percentageThreshold, std::ofstream *tagResult, int &pqValue, int &psValue, const int tagGeneType, std::string &ref_string, const HaplotagParameters &params){

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
    CigarParser* cigarParser = new GermlineHaplotagCigarParser(ref_pos, query_pos);
    cigarParser->parsingCigar(aln, bamHdr, chrName, params, firstVariantIter, currentChrVariants, ref_string, hpCount, variantsHP, countPS);
    delete cigarParser;

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

int HaplotagProcess::somaticJudgeHaplotype(
    const bam_hdr_t &bamHdr,
    const bam1_t &aln,
    const std::string &chrName,
    double percentageThreshold,
    std::ofstream *tagResult,
    int &pqValue,
    int &psValue,
    const int tagGeneType,
    std::string &ref_string,
    const HaplotagParameters &params
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

    /// Create a CIGAR parser using polymorphism design
    // Use GermlineHaplotagCigarParser to process CIGAR strings for normal samples    
    CigarParser* cigarParser = new SomaticJudgeHpCigarParser(ref_pos, query_pos, tumCountPS, somaticVarDeriveHP, highConSomaticData);
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
    hpResult = determineReadHP(hpCount, pqValue, norCountPS, norHPsimilarity, tumHPsimilarity, percentageThreshold, &totalHighSimilarity, &totalCrossTwoBlock, &totalWithOutVaraint);

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

    if(hpCount[1] == 0 && hpCount[2] == 0 && hpCount[3] != 0){
        if(hpResult == ReadHP::H3){
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
                    << hpCount[1]+hpCount[2]+hpCount[3]+hpCount[4]  << "\t"
                    << hpCount[1]                                   << "\t"
                    << hpCount[2]                                   << "\t"
                    << hpCount[3]                                   << "\t"
                    << hpCount[4]                                   << "\t"
                    << pqValue                                      << "\t\t";

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

void HaplotagProcess::OnlyTumorSNPjudgeHP(const std::string &chrName, int &curPos, MultiGenomeVar &curVar, std::string base, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *tumorAllelePosVec){

    if(curVar.isSomaticVariant){
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


HaplotagProcess::HaplotagProcess(HaplotagParameters &params):
params(params),tagResult(nullptr)
,totalAlignment(0),totalSupplementary(0),totalSecondary(0),totalUnmapped(0),totalTagCount(0),totalUnTagCount(0)
,processBegin(time(NULL))
{
    //initialize variable
    chrVec = nullptr;
    chrLength = nullptr;

    vcfSet[NORMAL].gene_type = NORMAL;
    vcfSet[TUMOR].gene_type = TUMOR;

    mergedChrVarinat = new std::map<std::string, std::map<int, MultiGenomeVar>>();;

    beforeCorrReadHpResult = new std::map<std::string, std::map<int, ReadHpResult>>();
    afterCorrReadHpResult = new std::map<std::string, std::map<int, ReadHpResult>>();

    hpBeforeInheritance = new ReadHpDistriLog();
    hpAfterInheritance = new ReadHpDistriLog();

    if(params.writeReadLog){
        tagResult = new std::ofstream(params.resultPrefix + ".tag.log");
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
    delete mergedChrVarinat;

    delete beforeCorrReadHpResult;
    delete afterCorrReadHpResult;

    delete hpBeforeInheritance;
    delete hpAfterInheritance;

    if(tagResult){
        tagResult->close();
        delete tagResult;
    }
};


void HaplotagProcess::taggingProcess()
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

        //somatic variant calling
        SomaticVarCaller *somaticVarCaller = new SomaticVarCaller(*chrVec, params);
        somaticVarCaller->VariantCalling(params, *chrVec, *chrLength, (*mergedChrVarinat), vcfSet, Genome::TUMOR);
        somaticVarCaller->getSomaticFlag(*chrVec, *mergedChrVarinat);

        delete somaticVarCaller;
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
    printTaggingResult();

    return;
};

void HaplotagProcess::printTaggingResult(){
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
}

GermlineHaplotagCigarParser::GermlineHaplotagCigarParser(int& ref_pos, int& query_pos)
:CigarParser(ref_pos, query_pos)
{

}

GermlineHaplotagCigarParser::~GermlineHaplotagCigarParser(){

}

void GermlineHaplotagCigarParser::processMatchOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, std::string& base){
    auto norVar = (*currentVariantIter).second.Variant[NORMAL];
    germlineJudgeSnpHap(*chrName, norVar, base, ref_pos, length, i, aln_core_n_cigar, cigar, currentVariantIter, *hpCount, *variantsHP, *norCountPS);
}

void GermlineHaplotagCigarParser::processDeletionOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, bool& alreadyJudgeDel){
    // only execute at the first phased normal snp
    if ((*currentVariantIter).second.isExists(NORMAL) && !alreadyJudgeDel){
        if((*currentVariantIter).second.Variant[NORMAL].is_phased_hetero){
            // longphase v1.73 only execute once
            alreadyJudgeDel = true;
            germlineJudgeDeletionHap(*chrName, *ref_string, ref_pos, length, query_pos, currentVariantIter, aln, *hpCount, *variantsHP, *norCountPS);
        }
    }
}

SomaticJudgeHpCigarParser::SomaticJudgeHpCigarParser(
    int& ref_pos,
    int& query_pos,
    std::map<int, int>& tumCountPS,
    std::map<int, std::pair<int, int>>& somaticVarDeriveHP,
    SomaticReadVerifier& highConSomaticData
):CigarParser(ref_pos, query_pos),
    tumCountPS(tumCountPS),
    somaticVarDeriveHP(somaticVarDeriveHP),
    highConSomaticData(highConSomaticData)
{

}

SomaticJudgeHpCigarParser::~SomaticJudgeHpCigarParser(){

}

void SomaticJudgeHpCigarParser::OnlyTumorSNPjudgeHP(const std::string &chrName, int &curPos, MultiGenomeVar &curVar, std::string base, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *tumorAllelePosVec){

    if(curVar.isSomaticVariant){
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

void SomaticJudgeHpCigarParser::processMatchOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, std::string& base){
    SomaticJudgeSnpHP(currentVariantIter, *chrName, base, *hpCount, *norCountPS, tumCountPS, variantsHP, nullptr);

    //record the somatic snp derive by which germline hp in this read
    if(currentVariantIter->second.isSomaticVariant){
        int deriveByHp = currentVariantIter->second.somaticReadDeriveByHP;
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
    
    highConSomaticData.recordRefAltAlleleCount(*chrName, base, currentVariantIter);
}

void SomaticJudgeHpCigarParser::processDeletionOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, bool& alreadyJudgeDel){
    highConSomaticData.recordDelReadCount(*chrName, currentVariantIter);
}

