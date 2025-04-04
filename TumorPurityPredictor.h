#ifndef TUMOR_PURITY_PREDICTOR_H
#define TUMOR_PURITY_PREDICTOR_H

#include "Util.h" 
#include "HaplotagBase.h"


enum PeakTrend{
    NONE = 0,
    UP = 1,
    DOWN = 2,
    FLAG = 3
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
        void SetThresholdByValley(const std::vector<HistogramData>& histogram, const double &max_height);
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
        std::map<std::string, std::map<int, PosBase>>& chrPosNorBase;
        std::map<std::string, std::map<int, HP3_Info>>& chrPosSomaticInfo;
        std::map<std::string, std::map<int, HP3_Info>> chrPosSomaticFlag;

        size_t initial_data_size;
        FilterCounts filterCounts;

    public:
        TumorPurityPredictor(
            const HaplotagParameters& params,
            const std::vector<std::string>& chrVec,
            std::map<std::string, std::map<int, PosBase>>& chrPosNorBase,
            std::map<std::string, std::map<int, HP3_Info>>& chrPosSomaticInfo); 

        ~TumorPurityPredictor();
        double predictTumorPurity();
        void buildPurityFeatureValueVec(std::vector<PurityData> &purityFeatureValueVec);

        int findPeakValleythreshold(const HaplotagParameters& params, const std::vector<PurityData> &purityFeatureValueVec);
        void peakValleyFilter(std::vector<PurityData> &purityFeatureValueVec, int &germlineReadHpCountThreshold);

        BoxPlotValue statisticPurityData(std::vector<PurityData> &purityFeatureValueVec);
        void removeOutliers(std::vector<PurityData> &purityFeatureValueVec, BoxPlotValue &plotValue);
        void markStatisticFlag(std::map<std::string, std::map<int, HP3_Info>>& chrPosSomaticInfo);

        void writePurityLog(const HaplotagParameters &params, double &purity, BoxPlotValue &plotValue, size_t &iteration_times, int &germlineReadHpCountThreshold);
};

#endif
