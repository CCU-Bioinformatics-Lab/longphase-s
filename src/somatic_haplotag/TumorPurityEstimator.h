#ifndef TUMOR_PURITY_ESTIMATOR_H
#define TUMOR_PURITY_ESTIMATOR_H

#include "../shared/Util.h" 
#include "../haplotag/HaplotagType.h"

/**
 * @brief Peak trend enumeration for histogram analysis
 */
enum PeakTrend{
    NONE_PT = 0,
    UP = 1,
    DOWN = 2,
    FLAG = 3
};

/**
 * @brief Box plot statistical values
 */
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

/**
 * @brief Tumor purity estimation data structure
 */
struct PurityData{
    std:: string chr;
    int pos;
    double germlineReadHpImbalanceRatio;
    int germlineReadHpCountInNorBam;

    /**
     * @brief Comparison function for sorting by germline read haplotype consistency ratio
     */
    static bool compareByGermlineReadHpConsisRatio(const PurityData& a, const PurityData& b){
        return a.germlineReadHpImbalanceRatio < b.germlineReadHpImbalanceRatio;
    }
};

/**
 * @brief Histogram data structure
 */
struct HistogramData {
    double count;     
    double percentage; 
    HistogramData(double c = 0, double p = 0.0) : count(c), percentage(p) {}
};

/**
 * @brief Peak structure for histogram analysis
 */
struct Peak{
    size_t histo_index;     // Histogram index
    double height;          // Peak height
    PeakTrend left_trend;   // Left trend
    PeakTrend right_trend;  // Right trend
    bool is_main_peak;      // Whether this is a main peak

    /**
     * @brief Sort peaks by index in ascending order
     */
    static bool compareByIndex(const Peak& a, const Peak& b){
        return a.histo_index < b.histo_index;
    }
    
    /**
     * @brief Sort peaks by height in descending order
     */
    static bool compareByHight(const Peak& a, const Peak& b){
        return a.height > b.height;
    }

    Peak(size_t index = 0, double height = 0) 
        : histo_index(index), height(height), 
          left_trend(PeakTrend::NONE_PT), 
          right_trend(PeakTrend::NONE_PT),
          is_main_peak(false){} 
};

/**
 * @brief Histogram class for data visualization and analysis
 */
class Histogram {
    private:
        const size_t MAX_HISTOGRAM_SIZE = 1000000;  ///< Maximum histogram size

        std::vector<HistogramData> histogram;       ///< Histogram data

        size_t total_snp_count;
        double max_height;
        std::pair<size_t, size_t> data_range;

        // Gaussian filter parameters
        static constexpr double DEFAULT_SIGMA = 1.0;           // Default sigma
        static constexpr int KERNEL_SIZE_MULTIPLIER = 6;       // Kernel size multiplier
    
    public:
        Histogram();
        ~Histogram();
        
        /**
         * @brief Build histogram from purity data
         */
        void buildHistogram(const std::vector<PurityData>& purity_data);
    
        void calculateStatistics();    
        
        const std::vector<HistogramData>& getHistogram() const { return histogram; }
        size_t getTotalSnpCount() const { return total_snp_count; }
        double getMaxHeight() const { return max_height; }
        const std::pair<size_t, size_t>& getDataRange() const { return data_range; }

        // Gaussian filter functions
        /**
         * @brief Apply Gaussian filter to histogram
         */
        void applyGaussianFilter(double sigma = DEFAULT_SIGMA);
        
        /**
         * @brief Create Gaussian kernel
         */
        std::vector<double> createGaussianKernel(double sigma);
        
        /**
         * @brief Get smoothed histogram
         */
        Histogram getSmoothedHistogram(double sigma);
};

/**
 * @brief Peak processor for histogram peak analysis
 */
class PeakProcessor {
    private:
        /**
         * @brief Valley structure for peak analysis
         */
        struct Valley {
            size_t index;
            double height;
            double percentage;
        };

        /**
         * @brief Main peak information
         */
        struct MainPeakInfo {
            bool found;        ///< Whether main peak is found
            size_t index;      ///< Main peak index
            Peak peak;         ///< Main peak data
            MainPeakInfo() : found(false), index(-1) ,peak() {}
        };

        /**
         * @brief Saddle point information
         */
        struct SaddlePointInfo {
            bool found;        ///< Whether saddle point is found
            size_t index;      ///< Saddle point index
            Peak peak;         ///< Saddle point peak
            Peak next_peak;    ///< Next peak
            Peak pre_peak;     ///< Previous peak
            SaddlePointInfo() : found(false), index(-1) ,peak() ,next_peak() ,pre_peak() {}
        };

        static constexpr double THRESHOLD_PERCENTAGE_LIMIT = 0.3;

        std::vector<Peak> peaksVec;    ///< Vector of peaks
        int mainPeakCount;             ///< Number of main peaks

        MainPeakInfo mainPeak;         ///< Main peak information
        SaddlePointInfo saddlePoint;   ///< Saddle point information

        Valley lowestValley;           ///< Lowest valley
        int threshold;                 ///< Threshold value
        double thresholdPercentage;    ///< Threshold percentage

    public:
        std::vector<std::string> exec_log;  ///< Execution log

        /**
         * @brief Find peaks in histogram
         */
        void findPeaks(const std::vector<HistogramData>& histogram, const double min_peak_count);
        
        /**
         * @brief Remove close peaks
         */
        void removeClosePeaks(const size_t minDistance);
        
        /**
         * @brief Determine peak trends
         */
        void determineTrends();  
        
        /**
         * @brief Find main peak candidates
         */
        void findMainPeakCandidates();
        
        /**
         * @brief Find first priority main peak
         */
        bool findFirstPriorityMainPeak();
        
        /**
         * @brief Find saddle point
         */
        bool findSaddlePoint();
        
        /**
         * @brief Find lowest valley between two peaks
         */
        bool findLowestValley(const std::vector<HistogramData>& histogram, size_t start_index, size_t end_index, Valley& result);
        
        /**
         * @brief Set threshold by valley analysis
         */
        void SetThresholdByValley(const std::vector<HistogramData>& histogram, const double &max_height);
        
        /**
         * @brief Get peak with offset
         */
        Peak getPeak(size_t histo_index, int offset);
        
        /**
         * @brief Transform trend to string
         */
        std::string transformTrend(const PeakTrend &trend);
        
        /**
         * @brief Get threshold value
         */
        int getThreshold();

        /**
         * @brief Write peak valley analysis log
         */
        void writePeakValleyLog(const std::string& resultPrefix,
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

/**
 * @brief Tumor purity estimator class
 */
class TumorPurityEstimator{
    private:
        /**
         * @brief Filter counts for different filtering stages
         */
        struct FilterCounts {
            int imbalanceRatioInNorBam = 0;
            int imbalanceRatioInNorBamMaxThr = 0;
            int imbalanceRatio = 0;
            int readHpCountInNorBam = 0;
            int percentageOfGermlineHp = 0;
            int peakValley = 0;
            size_t outliers = 0;
        };

        static constexpr float GERMLINE_HP_IMBALANCE_RATIO_MIN_THR = 0.0;
        static constexpr float GERMLINE_HP_IMBALANCE_RATIO_IN_NOR_BAM_MIN_THR = 0.0;
        static constexpr float GERMLINE_HP_IMBALANCE_RATIO_IN_NOR_BAM_MAX_THR = 0.7;
        static constexpr float GERMLINE_HP_PERCENTAGE_IN_NOR_BAM_MAX_THR = 0.7;
        static const int GERMLINE_HP_READ_COUNT_IN_NOR_BAM_MIN_THR = 5;
        
        const std::vector<std::string>& chrVec;    
        std::map<std::string, std::map<int, PosBase>>& chrPosNorBase;
        std::map<std::string, std::map<int, SomaticData>>& chrPosSomaticInfo;
        std::map<std::string, std::map<int, SomaticData>> chrPosSomaticFlag;

        const bool writeLog;
        const std::string& resultPrefix;   

        size_t initial_data_size;
        FilterCounts filterCounts;


        /**
         * @brief Build purity feature value vector
         */
        void buildPurityFeatureValueVec(std::vector<PurityData> &purityFeatureValueVec);

        /**
         * @brief Find peak valley threshold
         */
        int findPeakValleythreshold(const std::vector<PurityData> &purityFeatureValueVec);
        
        /**
         * @brief Apply peak valley filter
         */
        void peakValleyFilter(std::vector<PurityData> &purityFeatureValueVec, int &germlineReadHpCountThreshold);

        /**
         * @brief Calculate box plot statistics for purity data
         */
        BoxPlotValue statisticPurityData(std::vector<PurityData> &purityFeatureValueVec);
        
        /**
         * @brief Remove outliers from purity data
         */
        void removeOutliers(std::vector<PurityData> &purityFeatureValueVec, BoxPlotValue &plotValue);

        /**
         * @brief Write purity estimation results
         */
        void writePurityResult(double &purity, BoxPlotValue &plotValue, size_t &iteration_times, int &germlineReadHpCountThreshold);

    public:
        TumorPurityEstimator(
            const std::vector<std::string>& chrVec,
            std::map<std::string, std::map<int, PosBase>>& chrPosNorBase,
            std::map<std::string, std::map<int, SomaticData>>& chrPosSomaticInfo,
            const bool writeLog,
            const std::string& resultPrefix
        ); 
        ~TumorPurityEstimator();
        
        /**
         * @brief Estimate tumor purity
         */
        double estimateTumorPurity();

        /**
         * @brief Mark statistical flag at somatic positions
         */
        void markStatisticFlag(std::map<std::string, std::map<int, SomaticData>>& chrPosSomaticInfo);
};

#endif
