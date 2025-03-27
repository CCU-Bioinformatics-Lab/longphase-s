#ifndef SOMATIC_VAR_CALLER_H
#define SOMATIC_VAR_CALLER_H

#include "Util.h" 
#include "HaplotagBase.h"
#include "TumorPurityPredictor.h"
#include "HaplotagProcess.h"

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


struct DenseSnpInterval{
    int snpCount;
    std::map<int, double> snpAltMean;
    std::map<int, double> snpZscore;
    double totalAltMean;
    double StdDev;
    DenseSnpInterval(): snpCount(0), totalAltMean(0.0), StdDev(0.0){}
};

struct SomaticFilterParaemter
{
    // Determine whether to apply the filter
    bool applyFilter;
    // Determine whether to write log file
    bool writeVarLog;  
    
    // Tumor purity
    float tumorPurity;

    // Maximum normal VAF threshold
    float norVAF_maxThr;
    int norDepth_minThr;

    // Messy read ratio threshold
    float MessyReadRatioThreshold;
    int ReadCount_minThr;

    // Haplotag consistency filter threshold
    float HapConsistency_VAF_maxThr;
    int HapConsistency_ReadCount_maxThr;
    int HapConsistency_somaticRead_minThr;

    // Interval SNP count filter threshold
    float IntervalSnpCount_VAF_maxThr;
    int IntervalSnpCount_ReadCount_maxThr;
    int IntervalSnpCount_minThr;
    float zScore_maxThr;

    // Below the mapping quality read ratio threshold
    float LowMpqRatioThreshold; 
};

class SomaticVarCaller: public SomaticJudgeBase{
    private:
        // record the position that tagged as HP3
        // chr, tumor SNP pos, somatic info
        std::map<std::string, std::map<int, HP3_Info>> *chrPosSomaticInfo;
        
        // chr, variant position (0-base), reads HP 
        std::map<std::string, std::map<int, ReadHpResult>> *chrVarReadHpResult;
        
        ReadHpDistriLog *callerReadHpDistri;

        //  chr, startPos, endPos
        std::map<std::string, std::map<int, std::pair<int, DenseSnpInterval>>> *denseTumorSnpInterval;

        // chr, read ID, SNP,read HP
        std::map<std::string, std::map<std::string, ReadVarHpCount>> *chrReadHpResultSet;
        // chr, position, read ID, baseHP 
        std::map<std::string, std::map<int, std::map<std::string, int>>> *chrTumorPosReadCorrBaseHP;

        void InitialSomaticFilterParams(SomaticFilterParaemter &somaticParams, bool enableFilter);

        void SetFilterParamsWithPurity(SomaticFilterParaemter &somaticParams, double &tumorPurity);
    
        void extractTumorVariantData(const bam_hdr_t &bamHdr,const bam1_t &aln, const std::string &chr, const HaplotagParameters &params, BamBaseCounter *NorBase, VCF_Info *vcfSet
                                   , std::map<int, HP3_Info> &posReadCase, std::map<int, MultiGenomeVar> &currentChrVariants, std::map<int, MultiGenomeVar>::iterator &firstVariantIter
                                   , std::map<std::string, ReadVarHpCount> &readTotalHPcount, std::map<int, std::map<std::string, int>> &somaticPosReadID, std::string &ref_string);
        
        void ClassifyReadsByCase(std::vector<int> &readPosHP3, std::map<int, int> &NorCountPS, std::map<int, int> &hpCount, const HaplotagParameters &params, std::map<int, HP3_Info> &somaticPosInfo);

        double calculateStandardDeviation(const std::map<int, double>& data, double mean);
        void calculateZScores(const std::map<int, double>& data, double mean, double stdDev, std::map<int, double> &zScores);
        void calculateIntervalZScore(bool &isStartPos, int &startPos, int &pos, int &snpCount, DenseSnpInterval &denseSnp, std::map<int, std::pair<int, DenseSnpInterval>> &localDenseTumorSnpInterval);
        void getDenseTumorSnpInterval(std::map<int, HP3_Info> &somaticPosInfo, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, std::map<std::string, int>> &somaticPosReadHPCount, std::map<int, std::pair<int, DenseSnpInterval>> &closeSomaticSnpInterval);

        void SomaticFeatureFilter(const SomaticFilterParaemter &somaticParams, std::map<int, MultiGenomeVar> &currentChrVariants,const std::string &chr, std::map<int, HP3_Info> &somaticPosInfo, BamBaseCounter &NorBase, double& tumorPurity);
        
        void CalibrateReadHP(const std::string &chr, const SomaticFilterParaemter &somaticParams, std::map<int, HP3_Info> &somaticPosInfo, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, std::map<std::string, int>> &somaticPosReadHPCount);
        void CalculateReadSetHP(const HaplotagParameters &params, const std::string &chr, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, std::map<std::string, int>> &somaticPosReadHPCount);
        
        void StatisticSomaticPosReadHP(const std::string &chr, std::map<int, HP3_Info> &somaticPosInfo, std::map<int, std::map<std::string, int>> &somaticPosReadHPCount, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, ReadHpResult> &readHpDistributed);
        
        void WriteSomaticVarCallingLog(const HaplotagParameters &params, const SomaticFilterParaemter &somaticParams, const std::vector<std::string> &chrVec, BamBaseCounter &NorBase
                                     , std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
        void WriteOtherSomaticHpLog(const HaplotagParameters &params, const std::vector<std::string> &chrVec, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
        void WriteDenseTumorSnpIntervalLog(const HaplotagParameters &params, const std::vector<std::string> &chrVec);
        
        // temporary function
        void ShannonEntropyFilter(const std::string &chr, std::map<int, HP3_Info> &somaticPosInfo, std::map<int, MultiGenomeVar> &currentChrVariants, std::string &ref_string);
        double entropyComponent(int count, int total);
        double calculateShannonEntropy(int nA, int nC, int nT, int nG);
        double calculateMean(const std::map<int, double>& data);
        void FindOtherSomaticSnpHP(const std::string &chr, std::map<int, HP3_Info> &somaticPosInfo, std::map<int, MultiGenomeVar> &currentChrVariants);

        void releaseMemory();

    protected:
        void OnlyTumorSNPjudgeHP(const std::string &chrName, int &curPos, MultiGenomeVar &curVar, std::string base, VCF_Info *vcfSet, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *readPosHP3, BamBaseCounter *NorBase, std::map<int, HP3_Info> *SomaticPos);
    public:
        SomaticVarCaller();
        virtual ~SomaticVarCaller();

        void VariantCalling(const std::string BamFile, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat,const std::vector<std::string> &chrVec, std::map<std::string, int> &chrLength,const HaplotagParameters &params, VCF_Info *vcfSet, BamBaseCounter &NorBase);
        std::map<std::string, std::map<int, HP3_Info>> getSomaticChrPosInfo();


};

#endif

