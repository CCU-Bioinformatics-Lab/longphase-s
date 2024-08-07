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
    std::string tumorSnpFile;   //new
    std::string svFile;
    std::string modFile;
    std::string bamFile;
    std::string tumorBamFile;   //new
    std::string fastaFile;
    std::string resultPrefix;
    std::string region;
    std::string command;
    std::string version;
    std::string outputFormat;
    
    bool tagSupplementary;
    bool writeReadLog;
    bool tagTumorSnp;   //new
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
    float Homo_OnlyHP3ReadRatioThreshold;
    float Homo_MessyReadRatioThreshold;
    int Homo_readCountThreshold;
    float Homo_VAF_upper_threshold ; //not used
    float Homo_VAF_lower_threshold;
    float Homo_tumDeletionRatio;
};

enum Genome
{
    NORMAL = 0,
    TUMOR = 1
};

// record the variants from the normal and tumor VCF files (normal:0, tumor:1)
struct RefAltSet{
    RefAlt Variant[2];
    bool isExistNormal;
    bool isExistTumor;
    RefAltSet(): isExistNormal(false), isExistTumor(false){}
};

//record vcf information
struct VCF_Info 
{
    std::vector<std::string> chrVec;
    std::map<std::string, int> chrLength;

    // chr, variant position (0-base), allele haplotype
    std::map<std::string, std::map<int, RefAlt>> chrVariant;

    // chr, variant position (0-base), phased set
    std::map<std::string, int > psIndex;
    std::map<std::string, std::map<int, int>> chrVariantPS;
    
    std::map<std::string, std::map<int, std::string>> chrVariantHP1;
    std::map<std::string, std::map<int, std::string>> chrVariantHP2;
    
    // The number of SVs occurring on different haplotypes in a read
    std::map<std::string, std::map<int, int> > readSVHapCount;

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

    bool isHighConSomaticSNP;

    //Corresponding normal position has a very low VAF
    bool isNormalPosLowVAF;


    HP3_Info(): totalCleanHP3Read(0), HP1withHP3Read(0), HP2withHP3Read(0), OnlyHP3Read(0), MessyHPRead(0), unTag(0)
             , CaseReadCount(0), HP1withHP3ReadRatio(0.0), HP2WithHP3ReadRatio(0.0), OnlyHP3ReadRatio(0.0), MessyHPReadRatio(0.0)
             , tumDelRatio(0.0), base(), GTtype("")
             , isHighConSomaticSNP(false), isNormalPosLowVAF(false){}
};

class BamBaseCounter{
    private:
        // chr, variant position (0-base), base count & depth
        std::map<std::string, std::map<int, PosBase>> ChrVariantBase;

        void StatisticBaseInfo(const bam1_t &aln, std::string chrName, const HaplotagParameters &params, int genmoeType
                         , std::map<int, PosBase> &VariantBase, std::map<int, RefAltSet> &currentVariants ,std::map<int, RefAltSet>::iterator &firstVariantIter);
        void CalculateBaseInfo(std::string chr, std::map<int, PosBase> &VariantBase, std::map<int, RefAltSet> &currentVariants);

    public:
        BamBaseCounter(std::string BamFile, const HaplotagParameters &params, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat, std::vector<std::string> &chrVec, std::map<std::string, int> &chrLength, int genmoeType);
        ~BamBaseCounter();


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
        virtual void OnlyTumorSNPjudgeHP(RefAltSet &curVar, int &curPos, BamBaseCounter *NorBase, VCF_Info *vcfSet, std::string chrName, std::string base, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int,std::string> *variantsHP, std::vector<int> *readPosHP3, std::map<int, HP3_Info> *SomaticPos)=0;
    public:
        void SomaticJudgeSnpHP(std::map<int, RefAltSet>::iterator &currentVariantIter, VCF_Info *vcfSet, std::string chrName, std::string base, std::map<int, int> &hpCount
        , std::map<int, int> &NorCountPS, std::map<int, int> &tumCountPS, std::map<int,std::string> *variantsHP
        , std::vector<int> *readPosHP3, BamBaseCounter *NorBase, std::map<int, HP3_Info> *SomaticPos);
};

class SomaticDetectJudgeHP : public SomaticJudgeBase{
    protected:
        void OnlyTumorSNPjudgeHP(RefAltSet &curVar, int &curPos, BamBaseCounter *NorBase, VCF_Info *vcfSet, std::string chrName, std::string base, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int,std::string> *variantsHP, std::vector<int> *readPosHP3, std::map<int, HP3_Info> *SomaticPos);
    public:
};


class SomaticTaggingJudgeHP : public SomaticJudgeBase{
    protected:
        void OnlyTumorSNPjudgeHP(RefAltSet &curVar, int &curPos, BamBaseCounter *NorBase, VCF_Info *vcfSet, std::string chrName, std::string base, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int,std::string> *variantsHP, std::vector<int> *readPosHP3, std::map<int, HP3_Info> *SomaticPos);
    public:
};

class SomaticVarCaller{
    private:
        // record the position that tagged as HP3
        // chr, variant position
        std::map<std::string, std::map<int, HP3_Info>> SomaticChrPosInfo;

        void InitialSomaticFilterParams(SomaticFilterParaemter &somaticParams);

        void SetSomaticFilterParams(const SomaticFilterParaemter &somaticParams, std::string GTtype, float &OnlyHP3ReadRatioThreshold
                                  , float &messyReadRatioThreshold,int &readCountThreshold, float &VAF_upper_threshold, float &VAF_lower_threshold
                                  , float &tumDeletionRatioThreshold, float &tumLowMpqRatioThreshold);
    
        void StatisticSomaticPosInfo(const bam_hdr_t &bamHdr,const bam1_t &aln, std::string chr, HaplotagParameters &params, BamBaseCounter *NorBase, VCF_Info *vcfSet
                                   , std::map<int, HP3_Info> &posReadCase, std::map<int, RefAltSet> &currentChrVariants, std::map<int, RefAltSet>::iterator &firstVariantIter);
        void ClassifyReadsByCase(std::vector<int> &readPosHP3, std::string chr, std::map<int, int> &countPS, std::map<int, int> &hpCount, const HaplotagParameters &params, std::map<int, HP3_Info> &somaticPosInfo);

        void SomaticFeatureFilter(const SomaticFilterParaemter &somaticParams, std::map<int, RefAltSet> &currentChrVariants, std::string chr, std::map<int, HP3_Info> &somaticPosInfo);
        void WriteSomaticVarCallingLog(const HaplotagParameters &params, const SomaticFilterParaemter &somaticParams, std::ofstream *tagHP3Log, const std::vector<std::string> &chrVec, BamBaseCounter &NorBase
                                     , std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat);

    protected:

    public:
        SomaticVarCaller(std::string BamFile, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat, std::vector<std::string> &chrVec, std::map<std::string, int> &chrLength, HaplotagParameters &params, VCF_Info *vcfSet, BamBaseCounter &NorBase);
        ~SomaticVarCaller();

        bool isHighConfidneceSomaticSNP(std::string chr, int pos);
        std::map<int, HP3_Info>* getSomaticPosInfo(std::string chr);
        std::map<std::string, std::map<int, HP3_Info>> getSomaticChrPosInfo();


};

class HaplotagProcess
{
    void variantParser(std::string variantFile, VCF_Info &Info);
    void compressParser(std::string &variantFile, VCF_Info &Info);
    void unCompressParser(std::string &variantFile, VCF_Info &Info);
    void parserProcess(std::string &input, VCF_Info &Info);
    
    void tagRead(HaplotagParameters &params, const int geneType);
    
    std::vector<std::string> *chrVec;
    std::map<std::string, int> *chrLength;

    // chr, variant position (0-base), allele haplotype set
    std::map<std::string, std::map<int, RefAltSet>> mergedChrVarinat;

    // variant position (0-base), allele haplotype set
    std::map<int, RefAltSet> currentChrVariants;
    std::map<int, RefAltSet>::iterator firstVariantIter;

    // chr, variant position
    std::map<std::string, std::map<int, HP3_Info>>  chrPosReadCase;  

    // record the VCF files of the normal and tumor datasets (normal:0, tumor:1)
    VCF_Info vcfSet[2];

    void initFlag(bam1_t *aln, std::string flag);
    
    int judgeHaplotype(const bam_hdr_t &bamHdr,const bam1_t &aln, std::string chrName, double percentageThreshold, std::ofstream *tagResult, int &pqValue, int &psValue, const int tagGeneType, std::string &ref_string);
    int SomaticJudgeHaplotype(const bam_hdr_t &bamHdr,const bam1_t &aln, std::string chrName, double percentageThreshold, std::ofstream *tagResult, int &pqValue, int &psValue, const int tagGeneType, std::string &ref_string);
    
    int totalAlignment;
    int totalSupplementary;
    int totalSecondary;
    int totalUnmapped;
    int totalTagCount;
    int totalUnTagCount;

    //--------------------verification parameter---------------------
    bool tagTumorMode;

    int totalLowerQuality;
    int totalOtherCase;
    int totalunTag_HP0;
    int totalHP1;
    int totalHP2;
    int totalHP3;
    int totalHP4;
    int totalHighSimilarity;
    int totalCrossTwoBlock;
    int totalEmptyVariant;
    int totalWithOutVaraint;

    //---------------------------------------------------------------
    
    std::time_t processBegin;
    bool integerPS;
    bool parseSnpFile;
    bool parseSVFile;
    bool parseMODFile;
    
    public:
        HaplotagProcess(HaplotagParameters params);
        ~HaplotagProcess();

};


#endif