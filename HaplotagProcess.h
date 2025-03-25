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

    std::string benchmarkVcf;
    bool enableFilter;

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

enum PeakTrend{
    NONE = 0,
    UP = 1,
    DOWN = 2,
    FLAG = 3
};

// record the variants from the normal and tumor VCF files (normal:0, tumor:1, seqcHighCon:2)
struct RefAltSet{
    RefAlt Variant[3];
    bool isExistNormal;
    bool isExistTumor;
    bool isExistHighConSomatic;
    RefAltSet(): isExistNormal(false), isExistTumor(false), isExistHighConSomatic(false){}
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


struct BoxPlotValue {
    size_t data_size;
    double median;
    double q1;
    double q3;
    double iqr;
    double lowerWhisker;
    double upperWhisker;
    int outliers;
    BoxPlotValue(): data_size(0), median(0.0), q1(0.0), q3(0.0), iqr(0.0), lowerWhisker(0.0), upperWhisker(0.0), outliers(0){}
};

struct PurityData{
    std:: string chr;
    int pos;
    double germlineReadHpConsistencyRatio;
    int germlineReadHpCountInNorBam;

    static bool compareByGermlineReadHpConsisRatio(const PurityData& a, const PurityData& b){
        return a.germlineReadHpConsistencyRatio < b.germlineReadHpConsistencyRatio;
    }
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


struct DenseSnpInterval{
    int snpCount;
    std::map<int, double> snpAltMean;
    std::map<int, double> snpZscore;
    double totalAltMean;
    double StdDev;
    DenseSnpInterval(): snpCount(0), totalAltMean(0.0), StdDev(0.0){}
};

struct SomaticReadLog{
    std::string chr;
    std::string readID;
    std::string hpResult;
    //pos, hp
    std::map<int, int> somaticSnpHp;
    float germlineVarSimilarity;
    float deriveByHpSimilarity;
    int germlineSnpCount;
    int tumorSnpCount;
    SomaticReadLog(): chr(""), readID(""), hpResult(""), germlineVarSimilarity(0.0), deriveByHpSimilarity(0.0), germlineSnpCount(0), tumorSnpCount(0){}
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

struct HistogramData {
    double count;     
    double percentage; 
    HistogramData(double c = 0, double p = 0.0) : count(c), percentage(p) {}
};

struct Peak{
    size_t histo_index;
    double height;
    PeakTrend left_trend;
    PeakTrend right_trend;

    bool is_main_peak;

    //sort by index in ascending order
    static bool compareByIndex(const Peak& a, const Peak& b){
        return a.histo_index < b.histo_index;
    }
    //sort by height in descending order
    static bool compareByHight(const Peak& a, const Peak& b){
        return a.height > b.height;
    }

    Peak(size_t index = 0, double height = 0) 
        : histo_index(index), height(height), 
          left_trend(PeakTrend::NONE), 
          right_trend(PeakTrend::NONE),
          is_main_peak(false){} 
};

class Histogram {
    private:
        const size_t MAX_HISTOGRAM_SIZE = 1000000;

        std::vector<HistogramData> histogram;

        size_t total_snp_count;
        double max_height;
        std::pair<size_t, size_t> data_range;

        // Gaussian filter parameters
        static constexpr double DEFAULT_SIGMA = 1.0;
        static constexpr int KERNEL_SIZE_MULTIPLIER = 6;
    
    public:
        Histogram();
        ~Histogram();
        void buildHistogram(const std::vector<PurityData>& purity_data);
        void calculateStatistics();    
        const std::vector<HistogramData>& getHistogram() const { return histogram; }
        size_t getTotalSnpCount() const { return total_snp_count; }
        double getMaxHeight() const { return max_height; }
        const std::pair<size_t, size_t>& getDataRange() const { return data_range; }

        // New Gaussian filter functions
        void applyGaussianFilter(double sigma = DEFAULT_SIGMA);
        std::vector<double> createGaussianKernel(double sigma);
        Histogram getSmoothedHistogram(double sigma);
};


class PeakProcessor {
    private:
        struct Valley {
            size_t index;
            double height;
            double percentage;
        };

        struct MainPeakInfo {
            bool found;
            size_t index;
            Peak peak;
            MainPeakInfo() : found(false), index(-1) ,peak() {}
        };

        struct SaddlePointInfo {
            bool found;
            size_t index;
            Peak peak;
            Peak next_peak;
            Peak pre_peak;
            SaddlePointInfo() : found(false), index(-1) ,peak() ,next_peak() ,pre_peak() {}
        };

        static constexpr double THRESHOLD_PERCENTAGE_LIMIT = 0.3;

        std::vector<Peak> peaksVec;
        int mainPeakCount;

        MainPeakInfo mainPeak;
        SaddlePointInfo saddlePoint;

        Valley lowestValley;
        int threshold;
        double thresholdPercentage;

    public:
        std::vector<std::string> exec_log;

        void findPeaks(const std::vector<HistogramData>& histogram, const double min_peak_count);
        void removeClosePeaks(const size_t minDistance);
        void determineTrends();  
        void findMainPeakCandidates();
        bool findFirstPriorityMainPeak();
        bool findSaddlePoint();
        bool findLowestValley(const std::vector<HistogramData>& histogram, size_t start_index, size_t end_index, Valley& result);
        void SetThresholdByValley(const std::vector<HistogramData>& histogram);
        Peak getPeak(size_t histo_index, int offset);
        std::string transformTrend(const PeakTrend &trend);
        int getThreshold();

        void writePeakValleyLog(const HaplotagParameters &params,
                                const std::vector<HistogramData>& histogram,
                                const std::vector<HistogramData>& smoothed_histogram,
                                size_t &total_snp_count,
                                const std::pair<size_t, size_t>& data_range,
                                double &max_height,
                                const double &MIN_PEAK_RATIO,
                                double &peak_threshold,
                                double &sigma);

        PeakProcessor();
        ~PeakProcessor();
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
        void germlineJudgeSnpHap(const std::string& chrName, VCF_Info* vcfSet, RefAlt& norVar, const std::string& base, int& ref_pos, int& length, int& i, int& aln_core_n_cigar
        ,uint32_t* cigar, std::map<int, RefAltSet>::iterator& currentVariantIter, int& hp1Count, int& hp2Count, std::map<int, int>& variantsHP, std::map<int, int>& countPS);

        void germlineJudgeDeletionHap(const std::string& chrName, const std::string& ref_string, int& ref_pos, int& length, int& query_pos, std::map<int, RefAltSet>::iterator &currentVariantIter, VCF_Info* vcfSet, const bam1_t* aln, int& hp1Count, int& hp2Count, std::map<int, int>& variantsHP, std::map<int, int>& countPS);
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
                         , std::map<int, PosBase> &VariantBase, std::map<int, RefAltSet> &currentVariants ,std::map<int, RefAltSet>::iterator &firstVariantIter, VCF_Info *vcfSet, const std::string &ref_string);
        void CalculateBaseInfo(const std::string &chr, std::map<int, PosBase> &VariantBase, std::map<int, RefAltSet> &currentVariants);

    public:
        BamBaseCounter(bool enableFilter);
        ~BamBaseCounter();

        void CountingBamBase(const std::string &BamFile, const HaplotagParameters &params, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat, std::vector<std::string> &chrVec, std::map<std::string, int> &chrLength, VCF_Info *vcfSet, int genmoeType);

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

class TumorPurityPredictor{
    private:
        struct FilterCounts {
            int consistencyRatioInNorBam = 0;
            int consistencyRatioInNorBamMaxThr = 0;
            int consistencyRatio = 0;
            int readHpCountInNorBam = 0;
            int percentageOfGermlineHp = 0;
            int peakValley = 0;
            size_t outliers = 0;
        };

        static constexpr float GERMLINE_HP_CONSISTENCY_RATIO_MIN_THR = 0.0;
        static constexpr float GERMLINE_HP_CONSISTENCY_RATIO_IN_NOR_BAM_MIN_THR = 0.0;
        static constexpr float GERMLINE_HP_CONSISTENCY_RATIO_IN_NOR_BAM_MAX_THR = 0.7;
        static constexpr float GERMLINE_HP_PERCENTAGE_IN_NOR_BAM_MAX_THR = 0.7;
        static const int GERMLINE_HP_READ_COUNT_IN_NOR_BAM_MIN_THR = 5;

        const HaplotagParameters& params;    
        const std::vector<std::string>& chrVec;    
        BamBaseCounter& norBase;
        std::map<std::string, std::map<int, HP3_Info>>& chrPosSomaticInfo;

        size_t initial_data_size;
        FilterCounts filterCounts;

    public:
        TumorPurityPredictor(
            const HaplotagParameters& params,
            const std::vector<std::string>& chrVec,
            BamBaseCounter& norBase,
            std::map<std::string, std::map<int, HP3_Info>>& chrPosSomaticInfo); 

        ~TumorPurityPredictor();
        double predictTumorPurity();
        void buildPurityFeatureValueVec(std::vector<PurityData> &purityFeatureValueVec);

        int findPeakValleythreshold(const HaplotagParameters& params, const std::vector<PurityData> &purityFeatureValueVec);
        void peakValleyFilter(std::vector<PurityData> &purityFeatureValueVec, int &germlineReadHpCountThreshold);

        BoxPlotValue statisticPurityData(std::vector<PurityData> &purityFeatureValueVec);
        void removeOutliers(std::vector<PurityData> &purityFeatureValueVec, BoxPlotValue &plotValue);

        void writePurityLog(const HaplotagParameters &params, double &purity, BoxPlotValue &plotValue, size_t &iteration_times, int &germlineReadHpCountThreshold);
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
        void recordDeriveHp(int &pos, int &deriveHP, float deriveHPsimilarity, std::map<int, ReadHpResult> &varReadHpResult);
    public:
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
                                   , std::map<int, HP3_Info> &posReadCase, std::map<int, RefAltSet> &currentChrVariants, std::map<int, RefAltSet>::iterator &firstVariantIter
                                   , std::map<std::string, ReadVarHpCount> &readTotalHPcount, std::map<int, std::map<std::string, int>> &somaticPosReadID, std::string &ref_string);
        
        void ClassifyReadsByCase(std::vector<int> &readPosHP3, std::map<int, int> &NorCountPS, std::map<int, int> &hpCount, const HaplotagParameters &params, std::map<int, HP3_Info> &somaticPosInfo);

        double calculateStandardDeviation(const std::map<int, double>& data, double mean);
        void calculateZScores(const std::map<int, double>& data, double mean, double stdDev, std::map<int, double> &zScores);
        void calculateIntervalZScore(bool &isStartPos, int &startPos, int &pos, int &snpCount, DenseSnpInterval &denseSnp, std::map<int, std::pair<int, DenseSnpInterval>> &localDenseTumorSnpInterval);
        void getDenseTumorSnpInterval(std::map<int, HP3_Info> &somaticPosInfo, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, std::map<std::string, int>> &somaticPosReadHPCount, std::map<int, std::pair<int, DenseSnpInterval>> &closeSomaticSnpInterval);

        void SomaticFeatureFilter(const SomaticFilterParaemter &somaticParams, std::map<int, RefAltSet> &currentChrVariants,const std::string &chr, std::map<int, HP3_Info> &somaticPosInfo, BamBaseCounter &NorBase, double& tumorPurity);
        
        void CalibrateReadHP(const std::string &chr, const SomaticFilterParaemter &somaticParams, std::map<int, HP3_Info> &somaticPosInfo, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, std::map<std::string, int>> &somaticPosReadHPCount);
        void CalculateReadSetHP(const HaplotagParameters &params, const std::string &chr, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, std::map<std::string, int>> &somaticPosReadHPCount);
        
        void StatisticSomaticPosReadHP(const std::string &chr, std::map<int, HP3_Info> &somaticPosInfo, std::map<int, std::map<std::string, int>> &somaticPosReadHPCount, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, ReadHpResult> &readHpDistributed);
        
        void WriteSomaticVarCallingLog(const HaplotagParameters &params, const SomaticFilterParaemter &somaticParams, const std::vector<std::string> &chrVec, BamBaseCounter &NorBase
                                     , std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat);
        void WriteOtherSomaticHpLog(const HaplotagParameters &params, const std::vector<std::string> &chrVec, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat);
        void WriteDenseTumorSnpIntervalLog(const HaplotagParameters &params, const std::vector<std::string> &chrVec);
        
        // temporary function
        void ShannonEntropyFilter(const std::string &chr, std::map<int, HP3_Info> &somaticPosInfo, std::map<int, RefAltSet> &currentChrVariants, std::string &ref_string);
        double entropyComponent(int count, int total);
        double calculateShannonEntropy(int nA, int nC, int nT, int nG);
        double calculateMean(const std::map<int, double>& data);
        void FindOtherSomaticSnpHP(const std::string &chr, std::map<int, HP3_Info> &somaticPosInfo, std::map<int, RefAltSet> &currentChrVariants);

        void releaseMemory();

    protected:
        void OnlyTumorSNPjudgeHP(const std::string &chrName, int &curPos, RefAltSet &curVar, std::string base, VCF_Info *vcfSet, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *readPosHP3, BamBaseCounter *NorBase, std::map<int, HP3_Info> *SomaticPos);
    public:
        SomaticVarCaller();
        virtual ~SomaticVarCaller();

        void VariantCalling(const std::string BamFile, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat,const std::vector<std::string> &chrVec, std::map<std::string, int> &chrLength,const HaplotagParameters &params, VCF_Info *vcfSet, BamBaseCounter &NorBase);
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

class HighConBenchmark: public VcfParser{
    private:
        bool openTestingFunc;

        // store data
        std::map<std::string, std::map<int, RefAltDelCount>> posAltRefDelCount;
        std::vector<std::pair<int, int>> highConSomaticPos;
        std::vector<SomaticReadLog> totalReadVec;
        std::vector<SomaticReadLog> readsCrossingHighConSnpVec;
        std::vector<SomaticReadLog> taggedSomaticReadVec;
        void parserProcess(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat);
    public:

        HighConBenchmark();
        ~HighConBenchmark();
        void setTestingFunc(bool openTestingFunc);
        void loadHighConSomatic(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat);
        
        void recordDelReadCount(const std::string &chr, std::map<int, RefAltSet>::iterator &currentVariantIter);
        void recordRefAltAlleleCount(const std::string &chr, std::string &base, std::map<int, RefAltSet>::iterator &currentVariantIter);
        void recordCrossingHighConSnpRead(const std::string &chr, std::string &readID, std::string &hpResult, std::map<int, int> &variantsHP, std::map<int, int> &hpCount, double &norHPsimilarity, float &deriveByHpSimilarity, std::map<int, RefAltSet> &currentChrVariants);
        void recordTaggedSomaticRead(const std::string &chr, std::string &readID, std::string &hpResult, std::map<int, int> &variantsHP, std::map<int, int> &hpCount, double &norHPsimilarity, float &deriveByHpSimilarity, std::map<int, RefAltSet> &currentChrVariants);
        void recordTaggedRead(const std::string &chr, std::string &readID, std::string &hpResult, std::map<int, int> &variantsHP, std::map<int, int> &hpCount, double &norHPsimilarity, float &deriveByHpSimilarity, std::map<int, RefAltSet> &currentChrVariants);

        SomaticReadLog createBasicSomaticReadLog(const std::string &chr, std::string &readID, std::string &hpResult, double &norHPsimilarity, float &deriveByHpSimilarity, std::map<int, int> &hpCount);
        
        void writePosAlleleCountLog(std::vector<std::string> &chrVec, HaplotagParameters &params, std::string logPosfix, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat);
        void writeTaggedSomaticReadLog(HaplotagParameters &params, std::string logPosfix);
        void writeCrossHighConSnpReadLog(HaplotagParameters &params, std::string logPosfix);
        void writeTaggedReadLog(HaplotagParameters &params, std::string logPosfix);
        void writeReadLog(HaplotagParameters &params, std::string logPosfix, std::vector<SomaticReadLog> &somaticReadVec);
        
        void displaySomaticVarCount(std::vector<std::string> &chrVec, std::map<std::string, std::map<int, RefAltSet>> &mergedChrVarinat);
};

class HaplotagProcess: public SomaticJudgeBase, public GermlineJudgeBase
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
        int somaticJudgeHaplotype(const bam_hdr_t &bamHdr,const bam1_t &aln,const std::string &chrName, double percentageThreshold, std::ofstream *tagResult, int &pqValue, int &psValue, const int tagGeneType, std::string &ref_string);
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

        ReadHpDistriLog *hpBeforeInheritance;
        ReadHpDistriLog *hpAfterInheritance;
        
        HighConBenchmark highConSomaticData;
        //---------------------------------------------------------------
        
        std::time_t processBegin;
    protected:
        void OnlyTumorSNPjudgeHP(const std::string &chrName, int &curPos, RefAltSet &curVar, std::string base, VCF_Info *vcfSet, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *readPosHP3, BamBaseCounter *NorBase, std::map<int, HP3_Info> *SomaticPos);
    public:
        HaplotagProcess();
        void taggingProcess(HaplotagParameters &params);
        virtual ~HaplotagProcess();

};


#endif