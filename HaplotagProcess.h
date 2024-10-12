#ifndef HAPLOTAGPROCESS_H
#define HAPLOTAGPROCESS_H

#include "Util.h"
#include "ParsingBam.h"
#include <cmath>
#include <iomanip>
#include <omp.h>

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

    std::string seqcHighCon;

    bool tagSupplementary;
    bool writeReadLog;
    bool tagTumorSnp; 
};

struct SomaticFilterParaemter
{
    // Determine whether to apply the filter
    bool applyFilter;
    // Determine whether to write log file
    bool writeVarLog;  

    // Below the mapping quality read ratio threshold
    float LowMpqRatioThreshold; 

    // Phased heterozygous SNPs filter threshold
    float Hetero_OnlyHP3ReadRatioThreshold;
    float Hetero_MessyReadRatioThreshold;
    int Hetero_readCountThreshold;
    float Hetero_VAF_upper_threshold;
    float Hetero_VAF_lower_threshold;
    float Hetero_tumDeletionRatio;

    //Unphased Heterozygous SNPs filter threshold
    float Unphased_Hetero_OnlyHP3ReadRatioThreshold;
    float Unphased_Hetero_MessyReadRatioThreshold;
    int Unphased_Hetero_readCountThreshold;
    float Unphased_Hetero_VAF_upper_threshold;
    float Unphased_Hetero_VAF_lower_threshold;
    float Unphased_Hetero_tumDeletionRatio;

    //Homozygous SNPs filter threshold
    float Homo_OnlyHP3ReadRatioThreshold; //not used
    float Homo_MessyReadRatioThreshold;
    int Homo_readCountThreshold;
    float Homo_VAF_upper_threshold ; //not used
    float Homo_VAF_lower_threshold;
    float Homo_tumDeletionRatio;
};

enum Genome
{
    NORMAL = 0,
    TUMOR = 1,
    SEQC_HIGH_CON = 2
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

// record the variants from the normal and tumor VCF files (normal:0, tumor:1, seqcHighCon:2)
struct RefAltSet{
    RefAlt Variant[3];
    bool isExistNormal;
    bool isExistTumor;
    bool isExistSeqcHighCon;
    RefAltSet(): isExistNormal(false), isExistTumor(false), isExistSeqcHighCon(false){}
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
    
    std::map<std::string, std::map<int, std::string>> chrVariantHP1;
    std::map<std::string, std::map<int, std::string>> chrVariantHP2;
    
    // // The number of SVs occurring on different haplotypes in a read
    std::map<std::string, std::map<int, int>> readSVHapCount;


    int gene_type;
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
    float filteredMpqVAF;
    float lowMpqReadRatio;
    
    PosBase(): A_count(0), C_count(0), G_count(0), T_count(0), unknow(0), depth(0), delCount(0)
             , max_count(0), second_max_count(INT_MAX), max_base(" "), second_max_base(" ")
             , max_ratio(0), second_max_ratio(0), isHighRefAllelleFreq(false)
             , MPQ_A_count(0), MPQ_C_count(0), MPQ_G_count(0), MPQ_T_count(0), MPQ_unknow(0), filteredMpqDepth(0) 
             ,VAF(0.0), filteredMpqVAF(0.0), lowMpqReadRatio(0.0){}
};

struct ReadHpResult{
    std::map<int, int> hpResultCounter;
    int somaticSnpH3count;
    bool existDeriveByH1andH2;
    int coverRegionStartPos;
    int coverRegionEndPos;
    ReadHpResult(): somaticSnpH3count(0), existDeriveByH1andH2(false), coverRegionStartPos(INT_MAX), coverRegionEndPos(INT_MIN){}
};

struct ReadVarHpCount{
    int HP1;
    int HP2;
    int HP3;
    int HP4;

    // to modify the read ID to avoid supplementary read IDs overriding the primary read HP count.
    int readIDcount;
    int hpResult;
    int startPos;
    int endPos;
    int readLength;
    std::map<int, int> NorCountPS;
    
    ReadVarHpCount(): HP1(0), HP2(0), HP3(0), HP4(0), readIDcount(0), hpResult(0), startPos(0), endPos(0), readLength(0){}
};

struct coverRegionInfo{
    int startPos;
    int endPos;
    int length;

    coverRegionInfo(): startPos(0), endPos(0), length(0){}
};

struct somaticReadLog{
    std::string chr;
    std::string readID;
    std::string hpResult;
    //pos, hp
    std::map<int, int> somaticSnpHp;
    somaticReadLog(): chr(""), readID(""), hpResult(""){}
};

struct HP3_Info{
    //Case ratio 
    int totalCleanHP3Read;
    int HP1withHP3Read;
    int HP2withHP3Read;
    int OnlyHP3Read;
    int MessyHPRead;
    int unTag;

    int CaseReadCount;
    float HP1withHP3ReadRatio;
    float HP2WithHP3ReadRatio;
    float OnlyHP3ReadRatio;
    float MessyHPReadRatio;

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

    HP3_Info(): totalCleanHP3Read(0), HP1withHP3Read(0), HP2withHP3Read(0), OnlyHP3Read(0), MessyHPRead(0), unTag(0)
             , CaseReadCount(0), HP1withHP3ReadRatio(0.0), HP2WithHP3ReadRatio(0.0), OnlyHP3ReadRatio(0.0), MessyHPReadRatio(0.0)
             , tumDelRatio(0.0), base(), GTtype(""), somaticHp4Base(Nucleotide::UNKOWN), somaticHp5Base(Nucleotide::UNKOWN), somaticHp4BaseCount(0), somaticHp5BaseCount(0)
             , isHighConSomaticSNP(false), isNormalPosLowVAF(false), somaticReadDeriveByHP(0){}
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

class readHpDistriLog{
    private :
        // chr, variant position (0-base), reads HP 
        std::map<std::string, std::map<int, ReadHpResult>> chrVarReadHpResult;
    protected:

    public:
        readHpDistriLog();
        ~readHpDistriLog();
        //use for multi-thread SomaticVarCaller
        void mergeLocalReadHp(const std::string &chr, std::map<int, ReadHpResult> &localReadHpResult);
        void writeReadHpDistriLog(HaplotagParameters &params, std::string logPosfix, const std::vector<std::string> &chrVec);
        void writePosCoverRegionLog(HaplotagParameters &params, std::string logPosfix, const std::vector<std::string> &chrVec);
        void writeTagReadCoverRegionLog(HaplotagParameters &params, std::string logPosfix, const std::vector<std::string> &chrVec, std::map<std::string, int> &chrLength);
        void removeNotDeriveByH1andH2pos(const std::vector<std::string> &chrVec);
};

class BamBaseCounter{
    private:
        // chr, variant position (0-base), base count & depth
        std::map<std::string, std::map<int, PosBase>> *ChrVariantBase;

        void StatisticBaseInfo(const bam1_t &aln, const std::string &chrName, const HaplotagParameters &params, int genmoeType
                         , std::map<int, PosBase> &VariantBase, std::map<int, RefAltSet> &currentVariants ,std::map<int, RefAltSet>::iterator &firstVariantIter);
        void CalculateBaseInfo(const std::string &chr, std::map<int, PosBase> &VariantBase, std::map<int, RefAltSet> &currentVariants);

    public:
        BamBaseCounter();
        ~BamBaseCounter();

        void CountingBamBase(const std::string &BamFile, const HaplotagParameters &params, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat, std::vector<std::string> &chrVec, std::map<std::string, int> &chrLength, int genmoeType);

        bool isHighRefAllelleFreq(std::string chr, int pos);
        std::string getMaxFreqBase(std::string chr, int pos);
        float getMaxBaseRatio(std::string chr, int pos);
        float getSecondMaxBaseRatio(std::string chr, int pos);
        float getVAF(std::string chr, int pos);
        float getFilterdMpqVAF(std::string chr, int pos);
        float getLowMpqReadRatio(std::string chr, int pos);

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
        void SomaticJudgeSnpHP(std::map<int, RefAltSet>::iterator &currentVariantIter, VCF_Info *vcfSet, std::string chrName, std::string base, std::map<int, int> &hpCount
        , std::map<int, int> &NorCountPS, std::map<int, int> &tumCountPS, std::map<int, int> *variantsHP
        , std::vector<int> *readPosHP3, BamBaseCounter *NorBase, std::map<int, HP3_Info> *SomaticPos);

        virtual void OnlyTumorSNPjudgeHP(const std::string &chrName, int &curPos, RefAltSet &curVar, std::string base, VCF_Info *vcfSet, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *readPosHP3, BamBaseCounter *NorBase, std::map<int, HP3_Info> *SomaticPos)=0;
        int determineReadHP(std::map<int, int> &hpCount, int &pqValue,std::map<int, int> &norCountPS, double &norHPsimilarity, double &tumHPsimilarity,  double percentageThreshold, int *totalHighSimilarity, int *totalCrossTwoBlock, int *totalWithOutVaraint);

        int convertStrNucToInt(std::string &base);
        std::string convertIntNucToStr(int base);
        void recordReadHp(int &pos, int &hpResult, int &BaseHP, std::map<int, ReadHpResult> &varReadHpResult);
    public:
};

class SomaticVarCaller: public SomaticJudgeBase{
    private:
        // record the position that tagged as HP3
        // chr, variant position
        std::map<std::string, std::map<int, HP3_Info>> *SomaticChrPosInfo;
        
        // chr, variant position (0-base), reads HP 
        std::map<std::string, std::map<int, ReadHpResult>> *chrVarReadHpResult;

        readHpDistriLog *callerReadHpDistri;

        void InitialSomaticFilterParams(SomaticFilterParaemter &somaticParams);

        void SetSomaticFilterParams(const SomaticFilterParaemter &somaticParams, std::string GTtype, float &OnlyHP3ReadRatioThreshold
                                  , float &messyReadRatioThreshold,int &readCountThreshold, float &VAF_upper_threshold, float &VAF_lower_threshold
                                  , float &tumDeletionRatioThreshold, float &tumLowMpqRatioThreshold);
    
        void StatisticSomaticPosInfo(const bam_hdr_t &bamHdr,const bam1_t &aln, const std::string &chr, HaplotagParameters &params, BamBaseCounter *NorBase, VCF_Info *vcfSet
                                   , std::map<int, HP3_Info> &posReadCase, std::map<int, RefAltSet> &currentChrVariants, std::map<int, RefAltSet>::iterator &firstVariantIter
                                   , std::map<std::string, ReadVarHpCount> &readTotalHPcount, std::map<int, std::map<std::string, int>> &somaticPosReadID);
        void ClassifyReadsByCase(std::vector<int> &readPosHP3, std::map<int, int> &NorCountPS, std::map<int, int> &hpCount, const HaplotagParameters &params, std::map<int, HP3_Info> &somaticPosInfo);

        void SomaticFeatureFilter(const SomaticFilterParaemter &somaticParams, std::map<int, RefAltSet> &currentChrVariants,const std::string &chr, std::map<int, HP3_Info> &somaticPosInfo);
        void FindOtherSomaticSnpHP(const std::string &chr, std::map<int, HP3_Info> &somaticPosInfo, std::map<int, RefAltSet> &currentChrVariants);
        void CalibrateReadHP(const std::string &chr, const SomaticFilterParaemter &somaticParams, std::map<int, HP3_Info> &somaticPosInfo, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, std::map<std::string, int>> &somaticPosReadHPCount);
        void CalculateChrReadHP(const HaplotagParameters &params, const std::string &chr, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, std::map<std::string, int>> &somaticPosReadHPCount);
        void StatisticSomaticPosReadHP(const std::string &chr, std::map<int, HP3_Info> &somaticPosInfo, std::map<int, std::map<std::string, int>> &somaticPosReadHPCount, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, ReadHpResult> &readHpDistributed);
        void WriteSomaticVarCallingLog(const HaplotagParameters &params, const SomaticFilterParaemter &somaticParams, const std::vector<std::string> &chrVec, BamBaseCounter &NorBase
                                     , std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat);
        void WriteOtherSomaticHpLog(const HaplotagParameters &params, const std::vector<std::string> &chrVec, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat);

        void releaseMemory();

    protected:
        void OnlyTumorSNPjudgeHP(const std::string &chrName, int &curPos, RefAltSet &curVar, std::string base, VCF_Info *vcfSet, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *readPosHP3, BamBaseCounter *NorBase, std::map<int, HP3_Info> *SomaticPos);
    public:
        SomaticVarCaller();
        virtual ~SomaticVarCaller();

        void VariantCalling(const std::string BamFile, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat,const std::vector<std::string> &chrVec, std::map<std::string, int> &chrLength, HaplotagParameters &params, VCF_Info *vcfSet, BamBaseCounter &NorBase);
        std::map<std::string, std::map<int, HP3_Info>> getSomaticChrPosInfo();


};

class VcfParser{
    private:
        bool parseSnpFile;
        bool parseSVFile;
        bool parseMODFile;
        bool tagTumorMode;
        bool integerPS;
        void compressParser(std::string &variantFile, VCF_Info &Info, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat);
        void unCompressParser(std::string &variantFile, VCF_Info &Info, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat);
        virtual void parserProcess(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat);
    
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
        void variantParser(std::string &variantFile, VCF_Info &Info, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat);
};

struct RefAltDelCount{
    int refCount;
    int altCount;
    int delCount;
};

class highConBenchmark: public VcfParser{
    private:
        bool openTestingFunc;

        // store data
        std::map<std::string, std::map<int, RefAltDelCount>> posAltRefDelCount;
        std::vector<std::pair<int, int>> highConSomaticPos;
        std::vector<somaticReadLog> readsCrossingHighConSnpVec;
        std::vector<somaticReadLog> taggedSomaticReadVec;
        void parserProcess(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat);
    public:

        highConBenchmark();
        ~highConBenchmark();
        void setTestingFunc(bool openTestingFunc);
        void loadHighConSomatic(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat);
        
        void recordDelReadCount(const std::string &chr, std::map<int, RefAltSet>::iterator &currentVariantIter);
        void recordRefAltAlleleCount(const std::string &chr, std::string &base, std::map<int, RefAltSet>::iterator &currentVariantIter);
        void recordCrossingHighConSnpRead(const std::string &chr, std::string &readID, std::string &hpResult, std::map<int, int> &variantsHP, std::map<int, RefAltSet> &currentChrVariants);
        void recordTaggedSomaticRead(const std::string &chr, std::string &readID, std::string &hpResult, std::map<int, int> &variantsHP, std::map<int, RefAltSet> &currentChrVariants);
        
        void writePosAlleleCountLog(std::vector<std::string> &chrVec, HaplotagParameters &params, std::string logPosfix, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat);
        void writeTaggedSomaticReadLog(HaplotagParameters &params, std::string logPosfix);
        void writeCrossHighConSnpReadLog(HaplotagParameters &params, std::string logPosfix);
        void writeReadLog(HaplotagParameters &params, std::string logPosfix, std::vector<somaticReadLog> &somaticReadVec);
        
        void displaySomaticVarCount(std::vector<std::string> &chrVec, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat);
};

class HaplotagProcess: public SomaticJudgeBase
{
    private:
        std::vector<std::string> *chrVec;
        std::map<std::string, int> *chrLength;

        // chr, variant position (0-base), allele haplotype set
        std::map<std::string, std::map<int, RefAltSet>> *mergedChrVarinat;

        // variant position (0-base), allele haplotype set
        std::map<int, RefAltSet> currentChrVariants;
        std::map<int, RefAltSet>::iterator firstVariantIter;

        // chr, variant position (0-base), somatic SNP information
        std::map<std::string, std::map<int, HP3_Info>> *chrPosReadCase;  

        // record the VCF files of the normal and tumor datasets (normal:0, tumor:1, seqcHighCon:2)
        VCF_Info vcfSet[3];

        void tagRead(HaplotagParameters &params, const int geneType);

        void initFlag(bam1_t *aln, std::string flag);
        
        int judgeHaplotype(const bam_hdr_t &bamHdr,const bam1_t &aln, std::string chrName, double percentageThreshold, std::ofstream *tagResult, int &pqValue, int &psValue, const int tagGeneType, std::string &ref_string);
        int SomaticJudgeHaplotype(const bam_hdr_t &bamHdr,const bam1_t &aln,const std::string &chrName, double percentageThreshold, std::ofstream *tagResult, int &pqValue, int &psValue, const int tagGeneType, std::string &ref_string);
        std::string convertHpResultToString(int hpResult);
        
        int totalAlignment;
        int totalSupplementary;
        int totalSecondary;
        int totalUnmapped;
        int totalTagCount;
        int totalUnTagCount;

        //--------------------verification parameter---------------------
        bool tagTumorMode;

        // reads HP count
        std::map<int, int> totalHpCount;
        // reads untag count
        int totalLowerQuality;
        int totalOtherCase;
        int totalunTag_HP0;
        int totalreadOnlyH3Snp;
        int totalHighSimilarity;
        int totalCrossTwoBlock;
        int totalEmptyVariant;
        int totalWithOutVaraint;

        // chr, variant position (0-base), reads HP 
        std::map<std::string, std::map<int, ReadHpResult>> *beforeCorrReadHpResult;
        std::map<std::string, std::map<int, ReadHpResult>> *afterCorrReadHpResult;

        readHpDistriLog *hpBeforeInheritance;
        readHpDistriLog *hpAfterInheritance;
        
        highConBenchmark highConSomaticData;
        //---------------------------------------------------------------
        
        std::time_t processBegin;
    protected:
        void OnlyTumorSNPjudgeHP(const std::string &chrName, int &curPos, RefAltSet &curVar, std::string base, VCF_Info *vcfSet, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *readPosHP3, BamBaseCounter *NorBase, std::map<int, HP3_Info> *SomaticPos);
    public:
        HaplotagProcess();
        void TaggingProcess(HaplotagParameters &params);
        virtual ~HaplotagProcess();

};


#endif