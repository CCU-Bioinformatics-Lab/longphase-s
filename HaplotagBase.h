#ifndef HAPLOTAG_BASE_H
#define HAPLOTAG_BASE_H

#include "Util.h"
#include "ParsingBam.h"
#include <cmath>
#include <iomanip>
#include <omp.h>
#include <climits>

struct HaplotagParameters
{
    int numThreads;
    int qualityThreshold;
    
    double percentageThreshold;
    // Filter the read mapping quality below than threshold
    int somaticCallingMpqThreshold;  
    
    std::string snpFile;
    std::string tumorSnpFile;   
    std::string svFile;
    std::string modFile;
    std::string bamFile;
    std::string tumorBamFile;  
    std::string fastaFile;
    std::string resultPrefix;
    std::string region;
    std::string command;
    std::string version;
    std::string outputFormat;

    std::string benchmarkVcf;
    bool enableFilter;

    bool tagSupplementary;
    bool writeReadLog;
    bool tagTumorSnp; 
};

enum Genome
{
    NORMAL = 0,
    TUMOR = 1,
    HIGH_CON_SOMATIC = 2
};

enum Nucleotide
{
    UNKOWN = 0,
    A = 1,
    C = 2,
    G = 3,
    T = 4,
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

class BamFileRAII {
    public:
        samFile* in;
        bam_hdr_t* bamHdr;
        hts_idx_t* idx;
        bam1_t* aln;

        BamFileRAII(const std::string& BamFile, const std::string& fastaFile, htsThreadPool &threadPool);
        ~BamFileRAII();
};


struct VarData{
    RefAlt allele;
    //phased set (-1 means no phase set)
    int PhasedSet;

    //haplotype
    std::string HP1;
    std::string HP2;

    bool is_phased_hetero;
    bool is_homozygous;
    bool is_unphased_hetero;

    bool isExistPhasedSet(){
        return PhasedSet != -1;
    }

    VarData(): PhasedSet(-1), HP1(""), HP2(""), is_phased_hetero(false), is_homozygous(false), is_unphased_hetero(false){}
};

struct MultiGenomeVar{
    // record the variants from the normal and tumor VCF files (normal:0, tumor:1, HighCon:2)
    std::map<Genome, VarData> Variant;

    bool isExists(Genome type){
        return Variant.find(type) != Variant.end();
    }

    VarData& operator[](Genome type){
        return Variant[type];
    }
};

// struct CoreStruct{
//     std::map<std::string, std::map<int, MultiGenomeVar>> mergedChrVariant;

//     // The number of SVs occurring on different haplotypes in a read
//     //chr, genome, read, haplotype
//     std::map<Genome, std::map<std::string, std::map<int, int>>> readSVHapCount;
//     std::map<Genome, std::map<std::string, int>> psIndex;

//     std::map<Genome, std::vector<std::string>> chrVec;
//     std::map<Genome, std::map<std::string, int>> chrLength;
// };


//record each type of base in specific position
struct PosBase{
    int A_count;
    int C_count;
    int G_count;
    int T_count;
    int unknow;
    int depth;
    int delCount;

    int max_count;
    int second_max_count;

    std::string max_base;
    std::string second_max_base;
 
    float max_ratio;
    float second_max_ratio;

    //if only one type of Base in this position : 1 
    bool isHighRefAllelleFreq;

    int MPQ_A_count;
    int MPQ_C_count;
    int MPQ_G_count;
    int MPQ_T_count;
    int MPQ_unknow;
    int filteredMpqDepth;

    float VAF;
    //Non-deletion Adjusted AF
    float nonDelAF;
    float filteredMpqVAF;
    float lowMpqReadRatio;

    //snp position, read hp count
    std::map<int, int> ReadHpCount;
    
    PosBase(): A_count(0), C_count(0), G_count(0), T_count(0), unknow(0), depth(0), delCount(0)
             , max_count(0), second_max_count(INT_MAX), max_base(" "), second_max_base(" ")
             , max_ratio(0), second_max_ratio(0), isHighRefAllelleFreq(false)
             , MPQ_A_count(0), MPQ_C_count(0), MPQ_G_count(0), MPQ_T_count(0), MPQ_unknow(0), filteredMpqDepth(0) 
             ,VAF(0.0), nonDelAF(0.0), filteredMpqVAF(0.0), lowMpqReadRatio(0.0), ReadHpCount(std::map<int, int>()){}
};

struct HP3_Info{
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

    float tumDelRatio;

    PosBase base;
    std::string GTtype;
    int somaticHp4Base;
    int somaticHp5Base;
    int somaticHp4BaseCount;
    int somaticHp5BaseCount;

    bool isHighConSomaticSNP;

    //Corresponding normal position has a very low VAF
    bool isNormalPosLowVAF;

    int somaticReadDeriveByHP;
    double shannonEntropy;
    int homopolymerLength;
    //readHp, count
    std::map<int, int> somaticReadHpCount;
    //for predict tumor purity
    std::map<int, int> allReadHpCount;
    bool statisticPurity;
    float MeanAltCountPerVarRead;
    float zScore;
    int intervalSnpCount;
    bool inDenseTumorInterval;
    bool isFilterOut;

    HP3_Info(): totalCleanHP3Read(0), pure_H1_1_read(0), pure_H2_1_read(0), pure_H3_read(0), Mixed_HP_read(0), unTag(0)
             , CaseReadCount(0), pure_H1_1_readRatio(0.0), pure_H2_1_readRatio(0.0), pure_H3_readRatio(0.0), Mixed_HP_readRatio(0.0)
             , tumDelRatio(0.0), base(), GTtype(""), somaticHp4Base(Nucleotide::UNKOWN), somaticHp5Base(Nucleotide::UNKOWN), somaticHp4BaseCount(0), somaticHp5BaseCount(0)
             , isHighConSomaticSNP(false), isNormalPosLowVAF(false), somaticReadDeriveByHP(0), shannonEntropy(0.0), homopolymerLength(0), statisticPurity(false), MeanAltCountPerVarRead(0.0), zScore(0.0), intervalSnpCount(0), inDenseTumorInterval(false)
             , isFilterOut(false){}
};

//record vcf information
struct VCF_Info 
{
    std::vector<std::string> chrVec;
    std::map<std::string, int> chrLength;

    // chr, variant position (0-base), allele haplotype
    std::map<std::string, std::map<int, RefAlt>> chrVariant;

    std::map<std::string, int > psIndex;
    // chr, variant position (0-base), phased set
    std::map<std::string, std::map<int, int>> chrVariantPS;
    
    // chr, variant position (0-base), haplotype
    
    // // The number of SVs occurring on different haplotypes in a read
    std::map<std::string, std::map<int, int>> readSVHapCount;

    Genome gene_type;
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
    ReadHpResult(): somaticSnpH3count(0), existDeriveByH1andH2(false), deriveHP(0), coverRegionStartPos(INT_MAX), coverRegionEndPos(INT_MIN){}
};

//-------

class ReadHpDistriLog{
    private :
        struct coverRegionInfo{
            int startPos;
            int endPos;
            int length;

            coverRegionInfo(): startPos(0), endPos(0), length(0){}
        };
        // chr, variant position (0-base), reads HP 
        std::map<std::string, std::map<int, ReadHpResult>> chrVarReadHpResult;
    protected:

    public:
        ReadHpDistriLog();
        ~ReadHpDistriLog();
        //use for multi-thread SomaticVarCaller
        void mergeLocalReadHp(const std::string &chr, std::map<int, ReadHpResult> &localReadHpResult);
        void writeReadHpDistriLog(const HaplotagParameters &params, std::string logPosfix, const std::vector<std::string> &chrVec);
        void writePosCoverRegionLog(const HaplotagParameters &params, std::string logPosfix, const std::vector<std::string> &chrVec);
        void writeTagReadCoverRegionLog(const HaplotagParameters &params, std::string logPosfix, const std::vector<std::string> &chrVec, std::map<std::string, int> &chrLength);
        void removeNotDeriveByH1andH2pos(const std::vector<std::string> &chrVec);
};

class GermlineJudgeBase{
    private:

    protected:
        void germlineJudgeSnpHap(const std::string& chrName, VCF_Info* vcfSet, VarData& norVar, const std::string& base, int& ref_pos, int& length, int& i, int& aln_core_n_cigar
        ,uint32_t* cigar, std::map<int, MultiGenomeVar>::iterator& currentVariantIter, int& hp1Count, int& hp2Count, std::map<int, int>& variantsHP, std::map<int, int>& countPS);

        void germlineJudgeDeletionHap(const std::string& chrName, const std::string& ref_string, int& ref_pos, int& length, int& query_pos, std::map<int, MultiGenomeVar>::iterator &currentVariantIter, VCF_Info* vcfSet, const bam1_t* aln, int& hp1Count, int& hp2Count, std::map<int, int>& variantsHP, std::map<int, int>& countPS);
        void germlineJudgeSVHap(const bam1_t &aln, VCF_Info* vcfSet, int& hp1Count, int& hp2Count, const int& tagGeneType);
        int germlineDetermineReadHap(int& hp1Count, int& hp2Count, double& min, double& max, double& percentageThreshold, int& pqValue, int& psValue, std::map<int, int>& countPS, int* totalHighSimilarity, int* totalWithOutVaraint);
        void germlineGetRefLastVarPos(std::vector<int>& last_pos, const std::vector<std::string>& chrVec, VCF_Info* vcfSet, int geneType);
        void writeGermlineTagLog(std::ofstream& tagResult, const bam1_t& aln, const bam_hdr_t& bamHdr, int& hpResult, double& max, double& min, int& hp1Count, int& hp2Count, int& pqValue, const std::map<int, int>& variantsHP, const std::map<int, int>& countPS);
    public:
};


class BamBaseCounter : public GermlineJudgeBase{
    private:
        // chr, variant position (0-base), base count & depth
        std::map<std::string, std::map<int, PosBase>> *ChrVariantBase;
        bool applyFilter;

        // Test if the logic of judgeHaplotype is consistent (only single chromosome)
        std::ofstream *tagResult;

        void StatisticBaseInfo(const bam1_t &aln, const bam_hdr_t &bamHdr, const std::string &chrName, const HaplotagParameters &params, int genmoeType
                         , std::map<int, PosBase> &VariantBase, std::map<int, MultiGenomeVar> &currentVariants ,std::map<int, MultiGenomeVar>::iterator &firstVariantIter, VCF_Info *vcfSet, const std::string &ref_string);
        void CalculateBaseInfo(const std::string &chr, std::map<int, PosBase> &VariantBase, std::map<int, MultiGenomeVar> &currentVariants);

    public:
        BamBaseCounter(bool enableFilter);
        ~BamBaseCounter();

        void CountingBamBase(const std::string &BamFile, const HaplotagParameters &params, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat, std::vector<std::string> &chrVec, std::map<std::string, int> &chrLength, VCF_Info *vcfSet, int genmoeType);

        bool isHighRefAllelleFreq(std::string chr, int pos);
        std::string getMaxFreqBase(std::string chr, int pos);
        float getMaxBaseRatio(std::string chr, int pos);
        float getSecondMaxBaseRatio(std::string chr, int pos);
        float getVAF(std::string chr, int pos);
        float getNoDelAF(std::string chr, int pos);
        float getFilterdMpqVAF(std::string chr, int pos);
        float getLowMpqReadRatio(std::string chr, int pos);
        int getReadHpCountInNorBam(std::string chr, int pos, int Haplotype);

        int getBaseAcount(std::string chr, int pos);
        int getBaseCcount(std::string chr, int pos);
        int getBaseTcount(std::string chr, int pos);
        int getBaseGcount(std::string chr, int pos);
        int getDepth(std::string chr, int pos);
        int getMpqDepth(std::string chr, int pos);

        int getVarDeletionCount(std::string chr, int pos);
        void displayPosInfo(std::string chr, int pos);
};

class SomaticJudgeBase{
    private :

    protected:
        void SomaticJudgeSnpHP(std::map<int, MultiGenomeVar>::iterator &currentVariantIter, VCF_Info *vcfSet, std::string chrName, std::string base, std::map<int, int> &hpCount
        , std::map<int, int> &NorCountPS, std::map<int, int> &tumCountPS, std::map<int, int> *variantsHP
        , std::vector<int> *readPosHP3, BamBaseCounter *NorBase, std::map<int, HP3_Info> *SomaticPos);

        virtual void OnlyTumorSNPjudgeHP(const std::string &chrName, int &curPos, MultiGenomeVar &curVar, std::string base, VCF_Info *vcfSet, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *readPosHP3, BamBaseCounter *NorBase, std::map<int, HP3_Info> *SomaticPos)=0;
        int determineReadHP(std::map<int, int> &hpCount, int &pqValue,std::map<int, int> &norCountPS, double &norHPsimilarity, double &tumHPsimilarity,  double percentageThreshold, int *totalHighSimilarity, int *totalCrossTwoBlock, int *totalWithOutVaraint);

        int convertStrNucToInt(std::string &base);
        std::string convertIntNucToStr(int base);
        void recordReadHp(int &pos, int &hpResult, int &BaseHP, std::map<int, ReadHpResult> &varReadHpResult);
        void recordDeriveHp(int &pos, int &deriveHP, float deriveHPsimilarity, std::map<int, ReadHpResult> &varReadHpResult);
    public:
};
        
class VcfParser{
    private:
        bool parseSnpFile;
        bool parseSVFile;
        bool parseMODFile;
        bool tagTumorMode;
        bool integerPS;
        void compressParser(std::string &variantFile, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
        void unCompressParser(std::string &variantFile, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
        virtual void parserProcess(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
    
    protected:

    public:
        VcfParser();
        VcfParser(bool tagTumorMode);
        virtual ~VcfParser();
        void setParseSnpFile(bool parseSnpFile);
        void setParseSVFile(bool parseSVFile);
        void setParseMODFile(bool parseMODFile);
        bool getParseSnpFile();
        void reset();
        void variantParser(std::string &variantFile, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
};

#endif

