#include "HaplotagStragtegy.h"

/**
 * @brief Judges SNP haplotype for germline samples
 * @param chrName Chromosome name
 * @param norVar Normal variant data
 * @param base Base nucleotide
 * @param ref_pos Reference position
 * @param length Length of the operation
 * @param i CIGAR operation index
 * @param aln_core_n_cigar Number of CIGAR operations
 * @param cigar CIGAR array
 * @param currentVariantIter Iterator to current variant
 * @param hpCount Haplotype count map
 * @param variantsHP Variant haplotype map
 * @param countPS Phase set count map
 * 
 * Determines haplotype assignment for SNPs in germline samples
 */
void GermlineHaplotagStrategy::judgeSnpHap(
    const std::string& chrName,
    VarData& norVar,
    const std::string& base,
    int& ref_pos,
    int& length,
    int& i,
    int& aln_core_n_cigar,
    uint32_t* cigar,
    std::map<int, MultiGenomeVar>::iterator &currentVariantIter,
    std::map<int, int>& hpCount,
    std::map<int, int>& variantsHP,
    std::map<int, int>& countPS
){
    int curPos = (*currentVariantIter).first;

    // currentVariant is SNP
    if( norVar.variantType == VariantType::SNP ){
        // Detected that the base of the read is either REF or ALT. 
        if( (base == norVar.allele.Ref) || (base == norVar.allele.Alt) ){


            if(!norVar.isExistPhasedSet()){
                std::cerr << "[ERROR] (judgeSnpHap) => can't find the position:" 
                          << " chr: " << chrName << "\t"
                          << " pos: " << curPos << "\t"
                          << " ref: " << norVar.allele.Ref << "\t"
                          << " alt: " << norVar.allele.Alt << "\n";
                exit(EXIT_SUCCESS);
            }
            else{
                if( base == norVar.HP1){
                    hpCount[SnpHP::GERMLINE_H1]++;
                    variantsHP[curPos]=0;
                }
                if( base == norVar.HP2){
                    hpCount[SnpHP::GERMLINE_H2]++;
                    variantsHP[curPos]=1;
                }
                countPS[norVar.PhasedSet]++;
            }
            
        }
    }
    // currentVariant is insertion
    else if( norVar.variantType == VariantType::INSERTION && i+1 < aln_core_n_cigar){
        
        int hp1Length = norVar.HP1.length();
        int hp2Length = norVar.HP2.length();
        
        if ( ref_pos + length - 1 == (*currentVariantIter).first && bam_cigar_op(cigar[i+1]) == 1 ) {
            // hp1 occur insertion
            if( hp1Length != 1 && hp2Length == 1 ){
                hpCount[SnpHP::GERMLINE_H1]++;
                variantsHP[curPos]=0;
            }
            // hp2 occur insertion
            else if( hp1Length == 1 && hp2Length != 1 ){
                hpCount[SnpHP::GERMLINE_H2]++;
                variantsHP[curPos]=1;
            }
        }
        else {
            // hp1 occur insertion
            if( hp1Length != 1 && hp2Length == 1 ){
                hpCount[SnpHP::GERMLINE_H2]++;
                variantsHP[curPos]=1;
            }
            // hp2 occur insertion
            else if( hp1Length == 1 && hp2Length != 1 ){
                hpCount[SnpHP::GERMLINE_H1]++;
                variantsHP[curPos]=0;
            }
        }
        countPS[norVar.PhasedSet]++;
    } 
    // currentVariant is deletion
    else if( norVar.variantType == VariantType::DELETION && i+1 < aln_core_n_cigar) {

        int hp1Length = norVar.HP1.length();
        int hp2Length = norVar.HP2.length();
        
        if ( ref_pos + length - 1 == (*currentVariantIter).first && bam_cigar_op(cigar[i+1]) == 2 ) {
            // hp1 occur deletion
            if( hp1Length != 1 && hp2Length == 1 ){
                hpCount[SnpHP::GERMLINE_H1]++;
                variantsHP[curPos]=0;
            }
            // hp2 occur deletion
            else if( hp1Length == 1 && hp2Length != 1 ){
                hpCount[SnpHP::GERMLINE_H2]++;
                variantsHP[curPos]=1;
            }
        }
        else {
            // hp2 occur deletion
            if( hp1Length != 1 && hp2Length == 1 ){
                hpCount[SnpHP::GERMLINE_H2]++;
                variantsHP[curPos]=1;
            }
            // hp1 occur deletion
            else if( hp1Length == 1 && hp2Length != 1 ){
                hpCount[SnpHP::GERMLINE_H1]++;
                variantsHP[curPos]=0;
            }
        }
        countPS[norVar.PhasedSet]++;
    } 
}

/**
 * @brief Judges deletion haplotype for germline samples
 * @param chrName Chromosome name
 * @param ref_string Reference sequence string
 * @param ref_pos Reference position
 * @param length Length of deletion
 * @param query_pos Query position
 * @param currentVariantIter Iterator to current variant
 * @param aln BAM alignment record
 * @param hpCount Haplotype count map
 * @param variantsHP Variant haplotype map
 * @param countPS Phase set count map
 * 
 * Determines haplotype assignment for deletions in germline samples
 */
void GermlineHaplotagStrategy::judgeDeletionHap(
    const std::string& chrName,
    const std::string& ref_string,
    int& ref_pos,
    int& length,
    int& query_pos,
    std::map<int, MultiGenomeVar>::iterator &currentVariantIter,
    const bam1_t* aln,
    std::map<int, int>& hpCount,
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
                
                // SNP
                if (norVar.variantType == VariantType::SNP) {
                    // get the next match
                    char base_chr = seq_nt16_str[bam_seqi(bam_get_seq(aln), query_pos)];
                    std::string base(1, base_chr);

                    if (base == norVar.HP1) {
                        hpCount[SnpHP::GERMLINE_H1]++;
                        variantsHP[curPos] = 0;
                    }
                    if (base == norVar.HP2) {
                        hpCount[SnpHP::GERMLINE_H2]++;
                        variantsHP[curPos] = 1;
                    }
                    countPS[norVar.PhasedSet]++;
                }
                
                // the read deletion contain VCF's deletion
                else if (norVar.variantType == VariantType::DELETION) {

                    int hp1Length = norVar.HP1.length();
                    int hp2Length = norVar.HP2.length();
                    // hp1 occur deletion
                    if (hp1Length != 1 && hp2Length == 1) {
                        hpCount[SnpHP::GERMLINE_H1]++;
                        variantsHP[curPos] = 0;
                    }
                    // hp2 occur deletion
                    else if (hp1Length == 1 && hp2Length != 1) {
                        hpCount[SnpHP::GERMLINE_H2]++;
                        variantsHP[curPos] = 1;
                    }
                    countPS[norVar.PhasedSet]++;
                }
            }
        }
    }
}

/**
 * @brief Judges SV haplotype
 * @param aln BAM alignment record
 * @param vcfSet VCF information map
 * @param hpCount Haplotype count map
 * @param genomeSample Genome sample type
 * 
 * Determines haplotype assignment for structural variants
 */
void GermlineHaplotagStrategy::judgeSVHap(const bam1_t &aln, std::map<Genome, VCF_Info> &vcfSet, std::map<int, int>& hpCount, const int& genomeSample){
    auto readIter = vcfSet[Genome::NORMAL].readSVHapCount.find(bam_get_qname(&aln));
    if( readIter != vcfSet[Genome::NORMAL].readSVHapCount.end()){
        hpCount[SnpHP::GERMLINE_H1] += vcfSet[Genome::NORMAL].readSVHapCount[bam_get_qname(&aln)][0];
        hpCount[SnpHP::GERMLINE_H2] += vcfSet[Genome::NORMAL].readSVHapCount[bam_get_qname(&aln)][1];
    }
}

/**
 * @brief Judges read haplotype for germline samples
 * @param hpCount Haplotype count map
 * @param min Minimum haplotype count
 * @param max Maximum haplotype count
 * @param percentageThreshold Percentage threshold for confidence
 * @param pqValue Phasing quality value
 * @param psValue Phase set value
 * @param countPS Phase set count map
 * @param totalHighSimilarity Counter for high similarity reads
 * @param totalWithOutVaraint Counter for reads without variants
 * @return Haplotype result (ReadHP enum)
 * 
 * Determines the final haplotype assignment for a read based on variant counts
 */
int GermlineHaplotagStrategy::judgeReadHap(
    std::map<int, int>& hpCount, 
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

    if(hpCount[SnpHP::GERMLINE_H1] > hpCount[SnpHP::GERMLINE_H2]){
        min = hpCount[SnpHP::GERMLINE_H2];
        max = hpCount[SnpHP::GERMLINE_H1];
    }
    else{
        min = hpCount[SnpHP::GERMLINE_H1];
        max = hpCount[SnpHP::GERMLINE_H2];
    }

    if( max/(max+min) < percentageThreshold){
        // no tag
        pqValue = 0;
        if(totalHighSimilarity != nullptr) (*totalHighSimilarity)++;
    }
    else{
        if(hpCount[SnpHP::GERMLINE_H1] > hpCount[SnpHP::GERMLINE_H2]){
            hpResult = ReadHP::H1;
        }
        if(hpCount[SnpHP::GERMLINE_H1] < hpCount[SnpHP::GERMLINE_H2]){
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

/**
 * @brief Judges somatic SNP haplotype
 * @param currentVariantIter Iterator to current variant
 * @param chrName Chromosome name
 * @param base Base nucleotide
 * @param hpCount Haplotype count map
 * @param norCountPS Normal phase set count map
 * @param tumCountPS Tumor phase set count map
 * @param variantsHP Variant haplotype map
 * @param tumorAllelePosVec Vector to store tumor allele positions
 * 
 * Determines haplotype assignment for SNPs in somatic samples
 */
void SomaticJudgeHapStrategy::judgeSomaticSnpHap(std::map<int, MultiGenomeVar>::iterator &currentVariantIter, std::string chrName, std::string base, std::map<int, int> &hpCount, std::map<int, int> &norCountPS, std::map<int, int> &tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *tumorAllelePosVec){
    int curPos = (*currentVariantIter).first;
    auto& curVar = (*currentVariantIter).second;

    // normal & tumor SNP at the current position (base on normal phased SNPs)
    // both normal and tumor samples that do not exist in the high-confidence set
    if(curVar.isExists(NORMAL) && curVar.isExists(TUMOR)){
        // the tumor & normal SNP GT are phased heterozygous 
        if((curVar.Variant[NORMAL].GT == GenomeType::PHASED_HETERO) && 
            (curVar.Variant[TUMOR].GT == GenomeType::PHASED_HETERO)){ 
            judgeNormalSnpHap(chrName, curPos, curVar, base, hpCount, norCountPS, variantsHP);

        }
        //the normal SNP GT is phased heterozgous & the tumor SNP GT is unphased heterozgous 
        else if((curVar.Variant[NORMAL].GT == GenomeType::PHASED_HETERO) && 
                (curVar.Variant[TUMOR].GT == GenomeType::UNPHASED_HETERO)){   
            judgeNormalSnpHap(chrName, curPos, curVar, base, hpCount, norCountPS, variantsHP);

        }
        //the normal SNP GT is phased heterozgous & the tumor SNP GT is homozygous 
        else if((curVar.Variant[NORMAL].GT == GenomeType::PHASED_HETERO) && 
                (curVar.Variant[TUMOR].GT == GenomeType::UNPHASED_HOMO)){ 
            judgeNormalSnpHap(chrName, curPos, curVar, base, hpCount, norCountPS, variantsHP);
        }
    // only normal SNP at the current position
    }else if(curVar.isExists(NORMAL)){
        // the normal SNP GT is phased heterozgous SNP
        if((curVar.Variant[NORMAL].GT == GenomeType::PHASED_HETERO)){
            judgeNormalSnpHap(chrName, curPos, curVar, base, hpCount, norCountPS, variantsHP);
        }
    // only tumor SNP at the current position
    }else if(curVar.isExists(TUMOR)){
        //the tumor SNP GT is phased heterozygous
        if(curVar.Variant[TUMOR].GT == GenomeType::PHASED_HETERO){
            if(curVar.Variant[TUMOR].allele.Ref == base || curVar.Variant[TUMOR].allele.Alt == base){
                if(!curVar.Variant[TUMOR].isExistPhasedSet()){
                    std::cerr<< curPos << "\t"
                             << curVar.Variant[TUMOR].allele.Ref << "\t"
                             << curVar.Variant[TUMOR].allele.Alt << "\n";
                    exit(EXIT_SUCCESS);
                }else{
                    judgeTumorOnlySnpHap(chrName, curPos, curVar, base, hpCount, &tumCountPS, variantsHP, tumorAllelePosVec);
                }
            }
        //the tumor SNP GT is unphased heterozygous
        }else if(curVar.Variant[TUMOR].GT == GenomeType::UNPHASED_HETERO){
            if(curVar.Variant[TUMOR].allele.Ref == base || curVar.Variant[TUMOR].allele.Alt == base){
                judgeTumorOnlySnpHap(chrName, curPos, curVar, base, hpCount, nullptr, variantsHP, tumorAllelePosVec);
            }           
        //the tumor SNP GT is homozygous
        }else if(curVar.Variant[TUMOR].GT == GenomeType::UNPHASED_HOMO){
            if(curVar.Variant[TUMOR].allele.Ref == base || curVar.Variant[TUMOR].allele.Alt == base){
                judgeTumorOnlySnpHap(chrName, curPos, curVar, base, hpCount, nullptr, variantsHP, tumorAllelePosVec);
            }
        }
    }
}

/**
 * @brief Judges normal SNP haplotype in somatic context
 * @param chrName Chromosome name
 * @param curPos Current position
 * @param curVar Current variant data
 * @param base Base nucleotide
 * @param hpCount Haplotype count map
 * @param norCountPS Normal phase set count map
 * @param variantsHP Variant haplotype map
 * 
 * Determines haplotype assignment for normal SNPs in somatic samples
 */
void SomaticJudgeHapStrategy::judgeNormalSnpHap(
    const std::string& chrName, 
    int& curPos,
    MultiGenomeVar& curVar,
    std::string& base,
    std::map<int, int>& hpCount, 
    std::map<int, int>& norCountPS,
    std::map<int, int> *variantsHP
){
    if(curVar.Variant[NORMAL].allele.Ref == base || curVar.Variant[NORMAL].allele.Alt == base){

        if(!curVar.Variant[NORMAL].isExistPhasedSet()){
            std::cerr<< "Unable to locate the phase set of the current normal SNP\n"
                        << curPos << "\t"
                        << curVar.Variant[NORMAL].allele.Ref << "\t"
                        << curVar.Variant[NORMAL].allele.Alt  << "\n";
            exit(EXIT_SUCCESS);
        }

        std::string& norHP1 = curVar.Variant[NORMAL].HP1;
        std::string& norHP2 = curVar.Variant[NORMAL].HP2;

        if( base == norHP1){
            hpCount[1]++;
            if(variantsHP != nullptr) (*variantsHP)[curPos] = SnpHP::GERMLINE_H1;
        }
        if(base == norHP2){
            hpCount[2]++;
            if(variantsHP != nullptr) (*variantsHP)[curPos] = SnpHP::GERMLINE_H2;
        }
        norCountPS[curVar.Variant[NORMAL].PhasedSet]++;
    }
}

/**
 * @brief Judges somatic read haplotype
 * @param hpCount Haplotype count map
 * @param pqValue Phasing quality value
 * @param norCountPS Normal phase set count map
 * @param norHPsimilarity Normal haplotype similarity
 * @param tumHPsimilarity Tumor haplotype similarity
 * @param percentageThreshold Percentage threshold for confidence
 * @param totalHighSimilarity Counter for high similarity reads
 * @param totalCrossTwoBlock Counter for reads crossing two blocks
 * @param totalWithOutVaraint Counter for reads without variants
 * @return Haplotype result (ReadHP enum)
 * 
 * Determines the final haplotype assignment for a read in somatic samples
 */
int SomaticJudgeHapStrategy::judgeSomaticReadHap(
    std::map<int, int> &hpCount,
    int &pqValue,
    std::map<int, int> &norCountPS,
    double &norHPsimilarity,
    double &tumHPsimilarity,
    double percentageThreshold,
    int *totalHighSimilarity,
    int *totalCrossTwoBlock,
    int *totalWithOutVaraint
){
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
                        std::cerr << "[ERROR]: Unexpected haplotype : tumor Max HP= "<< maxTumorHP << std::endl;
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
                        std::cerr << "[ERROR]: Unexpected haplotype : tumor Max HP= "<< maxTumorHP << std::endl;
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

/**
 * @brief Judges tumor-only SNP haplotype for somatic data extraction
 * @param chrName Chromosome name
 * @param curPos Current position
 * @param curVar Current variant data
 * @param base Base nucleotide
 * @param hpCount Haplotype count map
 * @param tumCountPS Tumor phase set count map
 * @param variantsHP Variant haplotype map
 * @param tumorAllelePosVec Vector to store tumor allele positions
 * 
 * Determines haplotype assignment for tumor-only SNPs in somatic data extraction
 */
void ExtractSomaticDataStragtegy::judgeTumorOnlySnpHap(const std::string &chrName, int &curPos, MultiGenomeVar &curVar, std::string base, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *tumorAllelePosVec){

    if(tumorAllelePosVec == nullptr){
        std::cerr << "[ERROR] (SomaticDetectJudgeHP) => tumorAllelePosVec pointer cannot be nullptr"<< std::endl;
        exit(1);
    }

    std::string& TumorRefBase = curVar.Variant[TUMOR].allele.Ref;
    std::string& TumorAltBase = curVar.Variant[TUMOR].allele.Alt; 
    
    if(base == TumorAltBase){
        hpCount[3]++;
        if(variantsHP != nullptr) (*variantsHP)[curPos] = SnpHP::SOMATIC_H3;

        //record postions that tagged as HP3 for calculating the confidence of somatic positions
        (*tumorAllelePosVec).push_back(curPos);
    }else if(base != TumorRefBase && base != TumorAltBase){
        //other Haploype
    }

    if(tumCountPS != nullptr) (*tumCountPS)[curVar.Variant[TUMOR].PhasedSet]++;
}

/**
 * @brief Judges tumor-only SNP haplotype for somatic haplotagging
 * @param chrName Chromosome name
 * @param curPos Current position
 * @param curVar Current variant data
 * @param base Base nucleotide
 * @param hpCount Haplotype count map
 * @param tumCountPS Tumor phase set count map
 * @param variantsHP Variant haplotype map
 * @param tumorAllelePosVec Vector to store tumor allele positions
 * 
 * Determines haplotype assignment for tumor-only SNPs in somatic haplotagging
 */
void SomaticHaplotagStrategy::judgeTumorOnlySnpHap(const std::string &chrName, int &curPos, MultiGenomeVar &curVar, std::string base, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *tumorAllelePosVec){

    if(curVar.isSomaticVariant){
        std::string& TumorAltBase = curVar.Variant[TUMOR].allele.Alt; 

        if(base == TumorAltBase){
            hpCount[3]++;
            if(variantsHP != nullptr) (*variantsHP)[curPos] = SnpHP::SOMATIC_H3;
        }

        if(curVar.Variant[TUMOR].GT == GenomeType::PHASED_HETERO){
            if(tumCountPS != nullptr) (*tumCountPS)[curVar[TUMOR].PhasedSet]++;
        }
    }
}