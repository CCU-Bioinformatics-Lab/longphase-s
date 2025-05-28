#ifndef SOMATIC_VAR_CALLER_H
#define SOMATIC_VAR_CALLER_H

#include "Util.h" 
#include "HaplotagType.h"
#include "HaplotagParsingBam.h"
#include "TumorPurityPredictor.h"
#include "HaplotagLogging.h"


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
    std::map<int, int> norCountPS;
    
    ReadVarHpCount(): HP1(0), HP2(0), HP3(0), HP4(0), readIDcount(0), hpResult(0), startPos(0), endPos(0), readLength(0){}
};

struct DenseSnpData{
    double snpAltMean;
    double snpZscore;
    int minDistance;
    DenseSnpData(): snpAltMean(0.0), snpZscore(0.0), minDistance(0){}
};


struct DenseSnpInterval{
    std::map<int, double> snpAltMean;
    std::map<int, double> snpZscore;
    std::map<int, int> minDistance;
    int snpCount;
    double totalAltMean;
    double StdDev;
    DenseSnpInterval(): snpCount(0), totalAltMean(0.0), StdDev(0.0){}

    void clear(){
        snpAltMean.clear();
        snpZscore.clear();
        minDistance.clear();
        snpCount = 0;
        totalAltMean = 0.0;
        StdDev = 0.0;
    }
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


class ExtractNorDataChrProcessor : public ChromosomeProcessor{
    private:
        // store base information
        std::map<int, PosBase> *variantBase;
    protected:
        //override processRead
        void processRead(
            bam1_t &aln, 
            const bam_hdr_t &bamHdr,
            const std::string &chrName, 
            const HaplotagParameters &params, 
            const Genome& genmoeType, 
            std::map<int, MultiGenomeVar> &currentVariants,
            std::map<int, MultiGenomeVar>::iterator &firstVariantIter, 
            std::map<Genome, VCF_Info> &vcfSet, 
            const std::string &ref_string
        ) override;

        //override postProcess
        void postProcess(
            const std::string &chr,
            std::map<int, MultiGenomeVar> &currentVariants
        ) override;

    public:
        ExtractNorDataChrProcessor(std::map<std::string, std::map<int, PosBase>> &chrPosNorBase, const std::string &chr);
        virtual ~ExtractNorDataChrProcessor() override;
};

class ExtractNorDataCigarParser : public CigarParser{
    private:
        //specific data members
        std::map<int, PosBase>& variantBase;
        std::vector<int>& tumVarPosVec;
        const int& mappingQualityThr;
    protected:

        void processMatchOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, std::string& base) override;
        void processDeletionOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, bool& alreadyJudgeDel) override;

    public:
        ExtractNorDataCigarParser(
            std::map<int, PosBase>& variantBase,
            std::vector<int>& tumVarPosVec,
            int& ref_pos, 
            int& query_pos,
            const int& mappingQualityThr
        );
        ~ExtractNorDataCigarParser() override;
};

class ExtractNorDataBamParser : public HaplotagBamParser{
    private:
        // chr, variant position (0-base), base count & depth
        std::map<std::string, std::map<int, PosBase>>& chrPosNorBase;

        // Test if the logic of judgeHaplotype is consistent (only single chromosome)
        std::ofstream *tagResult;

    protected:
        // Factory method to create a new chromosome processor
        std::unique_ptr<ChromosomeProcessor> createProcessor(const std::string &chr) override{
            return std::unique_ptr<ChromosomeProcessor>(new ExtractNorDataChrProcessor(chrPosNorBase, chr));
        };

    public:
        ExtractNorDataBamParser(std::map<std::string, std::map<int, PosBase>>& chrPosNorBase);
        ~ExtractNorDataBamParser();
        void displayPosInfo(std::string chr, int pos);
};




class ExtractTumDataChrProcessor : public ChromosomeProcessor, public SomaticJudgeBase{
    private:
        // record the position that tagged as HP3
        // chr, variant position
        std::map<int, SomaticData> *somaticPosInfo;
        // read ID, SNP HP count 
        std::map<std::string, ReadVarHpCount> *readHpResultSet;
        // position, read ID, baseHP 
        std::map<int, std::map<std::string, int>> *tumorPosReadCorrBaseHP;

    protected:

        //override processRead
        void processRead(
            bam1_t &aln, 
            const bam_hdr_t &bamHdr,
            const std::string &chrName, 
            const HaplotagParameters &params, 
            const Genome& genmoeType, 
            std::map<int, MultiGenomeVar> &currentVariants,
            std::map<int, MultiGenomeVar>::iterator &firstVariantIter, 
            std::map<Genome, VCF_Info> &vcfSet, 
            const std::string &ref_string
        ) override;


        void classifyReadsByCase(std::vector<int> &readPosHP3, std::map<int, int> &norCountPS, std::map<int, int> &hpCount, std::map<int, SomaticData> &somaticPosInfo);

        void onlyTumorSNPjudgeHP(const std::string &chrName, int &curPos, MultiGenomeVar &curVar, std::string base, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *tumorAllelePosVec){};

        void postProcess(const std::string &chr, std::map<int, MultiGenomeVar> &currentVariants) override;
    public:
        ExtractTumDataChrProcessor(
            std::map<std::string, std::map<int, SomaticData>>& chrPosSomaticInfo,
            std::map<std::string, std::map<std::string, ReadVarHpCount>>& chrReadHpResultSet,
            std::map<std::string, std::map<int, std::map<std::string, int>>>& chrTumorPosReadCorrBaseHP,
            const std::string &chr
        );
        virtual ~ExtractTumDataChrProcessor() override;
};
        

class ExtractTumDataCigarParser : public CigarParser, public SomaticJudgeBase{
    private:
        std::map<int, SomaticData>& somaticPosInfo;

        //record tumor-unique variants on current read
        std::vector<int>& tumorAllelePosVec;

        //record tumor SNPs on current read
        std::vector<int>& tumorSnpPosVec;

        //record PS count( PS value, count)
        std::map<int, int>& tumCountPS;

        const int& mappingQualityThr;
    protected:

        void processMatchOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, std::string& base) override;
        void processDeletionOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, bool& alreadyJudgeDel) override;

        void onlyTumorSNPjudgeHP(const std::string &chrName, int &curPos, MultiGenomeVar &curVar, std::string base, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *tumorAllelePosVec) override;

    public:
        ExtractTumDataCigarParser(
            std::map<int, SomaticData> &somaticPosInfo, 
            std::vector<int>& tumorAllelePosVec, 
            std::vector<int>& tumorSnpPosVec, 
            std::map<int, int>& tumCountPS,
            int& ref_pos,
            int& query_pos,
            const int& mappingQualityThr
        );
        ~ExtractTumDataCigarParser() override;

};

class ExtractTumDataBamParser : public HaplotagBamParser{
    private:
        std::map<std::string, std::map<int, SomaticData>>& chrPosSomaticInfo;
        std::map<std::string, std::map<std::string, ReadVarHpCount>>& chrReadHpResultSet;
        std::map<std::string, std::map<int, std::map<std::string, int>>>& chrTumorPosReadCorrBaseHP;
    protected:

        std::unique_ptr<ChromosomeProcessor> createProcessor(const std::string &chr) override{
            return std::unique_ptr<ChromosomeProcessor>(new ExtractTumDataChrProcessor(chrPosSomaticInfo, chrReadHpResultSet, chrTumorPosReadCorrBaseHP, chr));
        };
    public:
        ExtractTumDataBamParser(
            std::map<std::string, std::map<int, SomaticData>>& chrPosSomaticInfo,
            std::map<std::string, std::map<std::string, ReadVarHpCount>>& chrReadHpResultSet,
            std::map<std::string, std::map<int, std::map<std::string, int>>>& chrTumorPosReadCorrBaseHP
        );
        ~ExtractTumDataBamParser();
};


class SomaticVarCaller: public SomaticJudgeBase, public GermlineJudgeBase{
    private:
        // somatic calling filter params
        SomaticFilterParaemter somaticParams;

        // record the position that tagged as HP3
        // chr, tumor SNP pos, somatic info
        std::map<std::string, std::map<int, SomaticData>> *chrPosSomaticInfo;

        // chr, variant position (0-base), normal variant base
        std::map<std::string, std::map<int, PosBase>> *chrPosNorBase;
        
        ReadHpDistriLog *callerReadHpDistri;

        //  chr, startPos, endPos
        std::map<std::string, std::map<int, std::pair<int, DenseSnpInterval>>> *denseTumorSnpInterval;

        // chr, read ID, SNP,read HP
        std::map<std::string, std::map<std::string, ReadVarHpCount>> *chrReadHpResultSet;
        // chr, position, read ID, baseHP 
        std::map<std::string, std::map<int, std::map<std::string, int>>> *chrTumorPosReadCorrBaseHP;

        void InitialSomaticFilterParams(bool enableFilter);

        void SetFilterParamsWithPurity(SomaticFilterParaemter &somaticParams, double &tumorPurity);

        double calculateStandardDeviation(const std::map<int, double>& data, double mean);
        void calculateZScores(const std::map<int, double>& data, double mean, double stdDev, std::map<int, double> &zScores);
        void calculateIntervalData(bool &isStartPos, int &startPos, int &pos, DenseSnpInterval &denseSnp, std::map<int, std::pair<int, DenseSnpInterval>> &localDenseTumorSnpInterval);
        void getDenseTumorSnpInterval(std::map<int, SomaticData> &somaticPosInfo, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, std::map<std::string, int>> &somaticPosReadHPCount, std::map<int, std::pair<int, DenseSnpInterval>> &closeSomaticSnpInterval);

        void somaticFeatureFilter(const SomaticFilterParaemter &somaticParams, std::map<int, MultiGenomeVar> &currentChrVariants,const std::string &chr, std::map<int, SomaticData> &somaticPosInfo, double& tumorPurity);
        
        void calibrateReadHP(const std::string &chr, std::map<int, SomaticData> &somaticPosInfo, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, std::map<std::string, int>> &somaticPosReadHPCount);
        void calculateReadSetHP(const std::string &chr, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, std::map<std::string, int>> &somaticPosReadHPCount, const double& percentageThreshold);
        
        void statisticSomaticPosReadHP(
            const std::string &chr,
            std::map<int, SomaticData> &somaticPosInfo,
            std::map<int, std::map<std::string, int>> &somaticPosReadHPCount,
            std::map<std::string, ReadVarHpCount> &readHpResultSet,
            chrReadHpResult &localReadHpDistri
        );
        
        void writeSomaticVarCallingLog(const HaplotagParameters &params, const SomaticFilterParaemter &somaticParams, const std::vector<std::string> &chrVec
                                     , std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
        void writeOtherSomaticHpLog(const std::string logFileName, const std::vector<std::string> &chrVec, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
        void writeDenseTumorSnpIntervalLog(const std::string logFileName, const std::vector<std::string> &chrVec);
        
        // temporary function
        void shannonEntropyFilter(const std::string &chr, std::map<int, SomaticData> &somaticPosInfo, std::map<int, MultiGenomeVar> &currentChrVariants, std::string &ref_string);
        double entropyComponent(int count, int total);
        double calculateShannonEntropy(int nA, int nC, int nT, int nG);
        double calculateMean(const std::map<int, double>& data);

        void findOtherSomaticSnpHP(const std::string &chr, std::map<int, SomaticData> &somaticPosInfo, std::map<int, MultiGenomeVar> &currentChrVariants);
        int convertStrNucToInt(std::string &base);
        std::string convertIntNucToStr(int base);
        
        void releaseMemory();

    protected:
        void onlyTumorSNPjudgeHP(
            const std::string &chrName,
            int &curPos,
            MultiGenomeVar &curVar,
            std::string base,
            std::map<int, int> &hpCount,
            std::map<int, int> *tumCountPS,
            std::map<int, int> *variantsHP,
            std::vector<int> *readPosHP3
        ){};
    public:
        SomaticVarCaller(const std::vector<std::string> &chrVec, const HaplotagParameters &params);
        virtual ~SomaticVarCaller();

        void variantCalling(
            const HaplotagParameters &params,
            const std::vector<std::string> &chrVec,
            const std::map<std::string, int> &chrLength,
            std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat,
            std::map<Genome, VCF_Info> &vcfSet,
            const Genome& genmoeType
        );
        void getSomaticFlag(const std::vector<std::string> &chrVec, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);

};

#endif

