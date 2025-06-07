#ifndef HAPLOTAG_TYPE_H
#define HAPLOTAG_TYPE_H

#include "Util.h"
#include "ParsingBam.h"
#include <cmath>
#include <iomanip>
#include <omp.h>
#include <climits>

// Option identifiers enum
enum HaplotagOption {
    OPT_HELP = 1,
    TAG_SUP,
    SV_FILE,
    REGION,
    LOG,
    MOD_FILE,
    CRAM
};

enum SomaticHaplotagOption{
    TUM_SNP = 50,
    TUM_BAM,
    BENCHMARK_VCF,
    BENCHMARK_BED,
    DISABLE_FILTER,
    TUMOR_PURITY
};

struct ParsingBamConfig{
    int numThreads;
    int qualityThreshold;
    double percentageThreshold;
    std::string resultPrefix;
    std::string region;
    std::string command;
    std::string version;
    std::string outputFormat;

    bool tagSupplementary;
    bool writeReadLog;
};

enum Genome
{
    NORMAL = 0,
    TUMOR = 1,
    TRUTH_SOMATIC = 2
};

enum GenomeType{
    NONE_GT = 0,
    PHASED_HETERO = 1,
    UNPHASED_HETERO = 2,
    UNPHASED_HOMO = 3
};

enum Nitrogenous
{
    UNKOWN = 0,
    A = 1,
    C = 2,
    G = 3,
    T = 4,
};

enum VariantType
{
    NONE_VAR = 0,
    SNP = 1,
    INSERTION = 2,
    DELETION = 3,
    MNP = 4
};


enum SnpHP
{
    NONE_SNP = 0,
    GERMLINE_H1 = 1,
    GERMLINE_H2 = 2,
    SOMATIC_H3 = 3,
    SOMATIC_H4 = 4,
    SOMATIC_H5 = 5
};

enum ReadHP
{
    unTag = 0,
    H1 = 1,
    H2 = 2,
    H3 = 3,
    H4 = 4,
    H1_1 = 5,
    H1_2 = 6,
    H2_1 = 7,
    H2_2 = 8,
};

struct VarData{
    const static int NONE_PHASED_SET = -1;

    RefAlt allele;
    //phased set (-1 means no phase set)
    int PhasedSet;

    //haplotype
    std::string HP1;
    std::string HP2;

    VariantType variantType;

    GenomeType GT;

    bool isExistPhasedSet(){
        return PhasedSet != NONE_PHASED_SET;
    }

    void setVariantType(){
        if(allele.Ref.length() == 1 && allele.Alt.length() == 1){
            variantType = VariantType::SNP;
        }else if(allele.Ref.length() == 1 && allele.Alt.length() > 1){
            variantType = VariantType::INSERTION;
        }else if(allele.Ref.length() > 1 && allele.Alt.length() == 1){
            variantType = VariantType::DELETION;
        }else if((allele.Ref.length() > 1) && (allele.Ref.length() == allele.Alt.length())){
            variantType = VariantType::MNP;
        }else{
            throw std::runtime_error("(loadVariantType)Invalid allele: " + allele.Ref + " " + allele.Alt);
        }
    }

    VarData(): PhasedSet(NONE_PHASED_SET), HP1(""), HP2(""), variantType(VariantType::NONE_VAR), GT(GenomeType::NONE_GT){}
};

struct MultiGenomeVar{
    // record the variants from the normal and tumor VCF files (normal, tumor, truth_somatic)
    std::map<Genome, VarData> Variant;
    bool isSomaticVariant;
    int somaticReadDeriveByHP;
    bool isInBedRegion;

    bool isExists(Genome sample){
        return Variant.find(sample) != Variant.end();
    }

    VarData& operator[](Genome sample){
        return Variant[sample];
    }
    
    MultiGenomeVar(): isSomaticVariant(false), somaticReadDeriveByHP(0), isInBedRegion(true){}
};

//record each type of base in specific position
struct PosBase{
    int A_count;
    int C_count;
    int G_count;
    int T_count;
    int unknow;
    int depth;
    int delCount;

    int MPQ_A_count;
    int MPQ_C_count;
    int MPQ_G_count;
    int MPQ_T_count;
    int MPQ_unknow;
    int filteredMpqDepth;

    float VAF;
    //Non-deletion Adjusted AF
    float nonDelVAF;
    float filteredMpqVAF;
    float lowMpqReadRatio;
    float delRatio;

    //germline haplotype imbalance ratio
    double germlineHaplotypeImbalanceRatio;
    double percentageOfGermlineHp;
    //snp position, read hp count
    std::map<int, int> ReadHpCount;
    
    PosBase(): A_count(0), C_count(0), G_count(0), T_count(0), unknow(0), depth(0), delCount(0)
             , MPQ_A_count(0), MPQ_C_count(0), MPQ_G_count(0), MPQ_T_count(0), MPQ_unknow(0), filteredMpqDepth(0) 
             , VAF(0.0), nonDelVAF(0.0), filteredMpqVAF(0.0), lowMpqReadRatio(0.0), delRatio(0.0)
             , germlineHaplotypeImbalanceRatio(0.0), percentageOfGermlineHp(0.0)
             , ReadHpCount(std::map<int, int>()){}

    int getBaseCount(const std::string& base) const {
        if (base == "A") return A_count;
        if (base == "T") return T_count;
        if (base == "C") return C_count;
        if (base == "G") return G_count;
        throw std::runtime_error("(getBaseCount)Invalid base: " + base);
    }

    int getMpqBaseCount(const std::string& base) const {
        if (base == "A") return MPQ_A_count;
        if (base == "T") return MPQ_T_count;
        if (base == "C") return MPQ_C_count;
        if (base == "G") return MPQ_G_count;
        throw std::runtime_error("(getMpqBaseCount)Invalid base: " + base);
    }
};

struct SomaticData{
    //Case ratio 
    int totalCleanHP3Read;
    int pure_H1_1_read;
    int pure_H2_1_read;
    int pure_H3_read;
    int Mixed_HP_read;
    int unTag;

    int CaseReadCount;
    float pure_H1_1_readRatio;
    float pure_H2_1_readRatio;
    float pure_H3_readRatio;
    float Mixed_HP_readRatio;

    PosBase base;
    std::string GTtype;
    int somaticHp4Base;
    int somaticHp5Base;
    int somaticHp4BaseCount;
    int somaticHp5BaseCount;

    bool isHighConSomaticSNP;

    int somaticReadDeriveByHP;
    double shannonEntropy;
    int homopolymerLength;

    bool statisticPurity;

    // imbalance ratio for predict purity
    double allelicImbalanceRatio;
    double somaticHaplotypeImbalanceRatio;

    //interval snp filter information
    float MeanAltCountPerVarRead;
    float zScore;
    int intervalSnpCount;
    int minDistance;
    bool inDenseTumorInterval;
    
    // filter out by somatic feature filter
    bool isFilterOut;

    //readHp, count
    std::map<int, int> somaticReadHpCount;

    SomaticData(): totalCleanHP3Read(0), pure_H1_1_read(0), pure_H2_1_read(0), pure_H3_read(0), Mixed_HP_read(0), unTag(0)
             , CaseReadCount(0), pure_H1_1_readRatio(0.0), pure_H2_1_readRatio(0.0), pure_H3_readRatio(0.0), Mixed_HP_readRatio(0.0)
             , base(), GTtype(""), somaticHp4Base(Nitrogenous::UNKOWN), somaticHp5Base(Nitrogenous::UNKOWN), somaticHp4BaseCount(0), somaticHp5BaseCount(0)
             , isHighConSomaticSNP(false), somaticReadDeriveByHP(0), shannonEntropy(0.0), homopolymerLength(0)
             , statisticPurity(false), allelicImbalanceRatio(0.0), somaticHaplotypeImbalanceRatio(0.0)
             , MeanAltCountPerVarRead(0.0), zScore(0.0), intervalSnpCount(0), minDistance(0), inDenseTumorInterval(false)
             , isFilterOut(false), somaticReadHpCount(std::map<int, int>()){}
};

//record vcf information
struct VCF_Info 
{
    std::vector<std::string> chrVec;
    std::map<std::string, int> chrLength;
        
    // The number of SVs occurring on different haplotypes in a read
    // read id, sv haplotype, count
    std::map<std::string, std::map<int, int>> readSVHapCount;

    Genome sample;
};



struct ReadHpResult{
    std::map<int, int> readHpCounter;
    std::map<int, int> somaticBaseReadHpCounter;
    std::vector<float> deriveHPsimilarVec;
    int somaticSnpH3count;
    bool existDeriveByH1andH2;
    int deriveHP;
    int coverRegionStartPos;
    int coverRegionEndPos;
    ReadHpResult(): somaticSnpH3count(0), existDeriveByH1andH2(false), deriveHP(0),
                    coverRegionStartPos(INT_MAX), coverRegionEndPos(INT_MIN){}
};
        

// ReadHP related utility functions
namespace ReadHapUtil {
    inline std::string readHapIntToString(int hpResult) {
        switch(hpResult) {
            case ReadHP::unTag: return ".";
            case ReadHP::H1:    return "1";
            case ReadHP::H2:    return "2";
            case ReadHP::H3:    return "3";
            case ReadHP::H4:    return "4";
            case ReadHP::H1_1:  return "1-1";
            case ReadHP::H2_1:  return "2-1";
            case ReadHP::H1_2:  return "1-2";
            case ReadHP::H2_2:  return "2-2";
            default:
                std::cerr << "[ERROR] (ReadHapUtil::toString) => Unsupported HP result: " << static_cast<int>(hpResult) << std::endl;
                exit(1);
        }
    }
}

#endif
