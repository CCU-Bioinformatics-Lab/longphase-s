#ifndef SOMATIC_BENCHMARK_H
#define SOMATIC_BENCHMARK_H

#include <iomanip>
#include "Util.h"
#include "HaplotagType.h"
#include "HaplotagVcfParser.h"

/**
 * @brief Structure to store somatic read logging information
 * 
 * Contains comprehensive information about somatic reads including
 * haplotype assignments, similarity scores, and variant positions
 */
struct SomaticReadLog{
    std::string chr;
    std::string readID;
    int hpResult;
    //pos, hp
    std::map<int, int> somaticSnpHp;   /** Map of variant positions to haplotype assignments */
    float germlineVarSimilarity;       /** Similarity score with germline variants */
    float deriveByHpSimilarity;        /** Similarity score derived by haplotype */
    int germlineSnpCount;              /** Count of germline SNPs */
    int tumorSnpCount;                 /** Count of tumor SNPs */
    
    /**
     * @brief Default constructor with initialization
     */
    SomaticReadLog(): chr(""), readID(""), hpResult(ReadHP::unTag), germlineVarSimilarity(0.0), deriveByHpSimilarity(0.0), germlineSnpCount(0), tumorSnpCount(0){}
};

/**
 * @brief Structure to store somatic read metrics and statistics
 * 
 * Contains various metrics for evaluating somatic variant detection performance
 * including allele counts, truth somatic positions, and read vectors
 */
struct SomaticReadMetrics{

    /**
     * @brief Structure to store reference, alternate, and deletion allele counts
     */
    struct RefAltDelCount{
        int refCount;   /** Reference allele count */
        int altCount;   /** Alternate allele count */
        int delCount;   /** Deletion count */
    };

    std::map<int, RefAltDelCount> posAltRefDelCount;           /** Map of positions to allele counts */
    std::vector<std::pair<int, int>> truthSomaticPosVec;       /** Vector of truth somatic positions */
    std::vector<SomaticReadLog> totalReadVec;                  /** Vector of all reads */
    std::vector<SomaticReadLog> coverTruthSomaticPosReadVec;   /** Vector of reads covering truth somatic positions */
    std::vector<SomaticReadLog> taggedSomaticReadVec;          /** Vector of tagged somatic reads */
};

/**
 * @brief Class for verifying and recording somatic read information
 * 
 * This class handles the verification and recording of somatic read data
 * including allele counts, haplotype assignments, and read statistics.
 * It provides methods to track reads that cross truth somatic somatic variants
 * and record tagged reads for performance evaluation.
 * 
 * Used by SomaticReadBenchmark for comprehensive somatic variant analysis
 */
class SomaticReadVerifier{
    private:
        bool openTestingFunc;           /** Flag to enable testing functionality */
        SomaticReadMetrics *metrics;    /** Pointer to metrics for data collection */

        /**
         * @brief Create basic somatic read log entry
         * @param chr Chromosome name
         * @param readID Read identifier
         * @param hpResult Haplotype result
         * @param norHPsimilarity Normal haplotype similarity
         * @param deriveByHpSimilarity Derived haplotype similarity
         * @param hpCount Haplotype count map
         * @return SomaticReadLog object
         */
        SomaticReadLog createBasicSomaticReadLog(const std::string &chr, std::string &readID, int &hpResult, double &norHPsimilarity, float &deriveByHpSimilarity, std::map<int, int> &hpCount);
    
    public:
        /**
         * @brief Constructor
         * @param openTestingFunc Flag to enable testing functionality
         * @param metrics Pointer to metrics for data collection
         */
        SomaticReadVerifier(bool openTestingFunc, SomaticReadMetrics *metrics);
        
        /**
         * @brief Destructor
         */
        ~SomaticReadVerifier();

        /**
         * @brief Record deletion read count for a variant position
         * @param chr Chromosome name
         * @param currentVariantIter Iterator to current variant
         */
        void recordDelReadCount(const std::string &chr, std::map<int, MultiGenomeVar>::iterator &currentVariantIter);
        
        /**
         * @brief Record reference/alternate allele count for a variant position
         * @param chr Chromosome name
         * @param base Base at the position
         * @param currentVariantIter Iterator to current variant
         */
        void recordRefAltAlleleCount(const std::string &chr, std::string &base, std::map<int, MultiGenomeVar>::iterator &currentVariantIter);
        
        /**
         * @brief Record reads crossing truth somatic somatic SNPs
         * @param chr Chromosome name
         * @param readID Read identifier
         * @param hpResult Haplotype result
         * @param variantsHP Map of variant positions to haplotypes
         * @param hpCount Haplotype count map
         * @param norHPsimilarity Normal haplotype similarity
         * @param deriveByHpSimilarity Derived haplotype similarity
         * @param currentChrVariants Current chromosome variants
         */
        void recordCrossingTruthSomaticSnpRead(const std::string &chr, std::string &readID, int &hpResult, std::map<int, int> &variantsHP, std::map<int, int> &hpCount, double &norHPsimilarity, float &deriveByHpSimilarity, std::map<int, MultiGenomeVar> &currentChrVariants);
        
        /**
         * @brief Record tagged reads for performance evaluation
         * @param chr Chromosome name
         * @param readID Read identifier
         * @param hpResult Haplotype result
         * @param variantsHP Map of variant positions to haplotypes
         * @param hpCount Haplotype count map
         * @param norHPsimilarity Normal haplotype similarity
         * @param deriveByHpSimilarity Derived haplotype similarity
         * @param currentChrVariants Current chromosome variants
         */
        void recordTaggedRead(const std::string &chr, std::string &readID, int &hpResult, std::map<int, int> &variantsHP, std::map<int, int> &hpCount, double &norHPsimilarity, float &deriveByHpSimilarity, std::map<int, MultiGenomeVar> &currentChrVariants);
};

/**
 * @brief Main benchmark class for somatic variant analysis
 * 
 * This class extends VcfParser to provide comprehensive benchmarking
 * capabilities for somatic variant detection. It handles truth VCF parsing,
 * BED region processing, and generates detailed performance metrics.
 * 
 * Key functionalities:
 * - Parse truth somatic VCF files for benchmarking
 * - Process BED regions for targeted analysis
 * - Mark variants in/out of BED regions
 * - Generate comprehensive performance reports
 * - Calculate precision, recall, and F1 scores
 * 
 * Used for evaluating the accuracy of somatic variant detection algorithms
 */
class SomaticReadBenchmark: public VcfParser{
    private:

        /**
         * @brief Structure to store BED region information
         */
        struct BedRegion {
            int start;
            int end;
        };
        
        bool openTestingFunc;           /** Flag to enable testing functionality */
        bool loadedBedFile;             /** Flag indicating if BED file is loaded */
        std::map<std::string, SomaticReadMetrics> chrMetrics;  /** Map of chromosome names to metrics */
        std::map<std::string, std::vector<BedRegion>> bedRegions;  /** Map of chromosome names to BED regions */
        std::string benchmarkVcf;       /** Path to benchmark VCF file */
        std::string benchmarkBed;       /** Path to benchmark BED file */
        int mappingQualityThreshold;    /** Mapping quality threshold for analysis */
        std::map<Genome, int> variantInBedRegionCount;   /** Count of variants in BED regions */
        std::map<Genome, int> variantOutBedRegionCount;  /** Count of variants outside BED regions */

        /**
         * @brief Override parserProcess to handle truth somatic VCF parsing
         * @param input VCF line content
         * @param Info VCF metadata and sample information
         * @param mergedChrVarinat Output container for parsed variants
         */
        void parserProcess(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat) override;
        
        /**
         * @brief Set chromosome somatic read vector pointer for multi-threaded access
         * @param chr Chromosome name
         * @param somaticReadVecMap Map of chromosome names to read vector pointers
         * @param somaticReadVec Somatic read vector
         */
        void setChrSomaticReadVecPtr(
            const std::string &chr,
            std::map<std::string, std::vector<SomaticReadLog>*> &somaticReadVecMap,
            std::vector<SomaticReadLog> &somaticReadVec
        );
        
        /**
         * @brief Write read log to output file
         * @param chrVec Vector of chromosome names
         * @param outputFileName Output file name
         * @param somaticReadVecMap Map of chromosome names to read vector pointers
         */
        void writeReadLog(
            const std::vector<std::string>& chrVec,
            std::string outputFileName,
            std::map<std::string, std::vector<SomaticReadLog>*> &somaticReadVecMap
        );

        /**
         * @brief Process a single BED file line
         * @param line BED file line content
         * @return True if line is valid, false otherwise
         */
        bool processBedLine(const std::string& line);

        /**
         * @brief Calculate recall metric
         * @param TP True positives
         * @param TP_FN True positives plus false negatives
         * @return Recall value
         */
        static float calculateRecall(int TP, int TP_FN){
            if(TP_FN == 0 || TP == 0) return 0.0;
            return (float)TP / (float)TP_FN;
        }

        /**
         * @brief Calculate precision metric
         * @param TP True positives
         * @param TP_FP True positives plus false positives
         * @return Precision value
         */
        static float calculatePrecision(int TP, int TP_FP){
            if(TP_FP == 0 || TP == 0) return 0.0;
            return (float)TP / (float)TP_FP;
        }

        /**
         * @brief Calculate F1 score metric
         * @param recall Recall value
         * @param precision Precision value
         * @return F1 score value
         */
        static float calculateF1Score(float recall, float precision){
            if(recall == 0.0 || precision == 0.0) return 0.0;
            return 2 * recall * precision / (recall + precision);
        }

    public:

        /**
         * @brief Constructor
         * @param benchmarkVcf Path to benchmark VCF file
         * @param benchmarkBed Path to benchmark BED file
         * @param mappingQualityThreshold Mapping quality threshold
         */
        SomaticReadBenchmark(std::string benchmarkVcf, std::string benchmarkBed, int mappingQualityThreshold);
        
        /**
         * @brief Destructor
         */
        ~SomaticReadBenchmark();
        
        /**
         * @brief Enable or disable testing functionality
         * @param openTestingFunc True to enable testing
         */
        void setEnabled(bool openTestingFunc);
        
        /**
         * @brief Check if testing functionality is enabled
         * @return True if testing is enabled
         */
        bool isEnabled();
        
        /**
         * @brief Check if BED file is loaded
         * @return True if BED file is loaded
         */
        bool isLoadBedFile();
        
        /**
         * @brief Initialize chromosome key for multi-threaded access
         * @param chr Chromosome name
         */
        void loadChrKey(const std::string &chr);
        
        /**
         * @brief Load truth somatic VCF file
         * @param input Input VCF file path
         * @param Info VCF metadata and sample information
         * @param mergedChrVarinat Output container for parsed variants
         */
        void loadTruthSomaticVCF(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
        
        /**
         * @brief Parse benchmark BED file
         * @param bedFile Path to BED file
         */
        void parseBedFile(const std::string& bedFile);

        /**
         * @brief Mark variants in BED regions
         * @param chrVec Vector of chromosome names
         * @param mergedChrVarinat Variant data container
         */
        void markVariantsInBedRegions(
            std::vector<std::string> &chrVec,
            std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
        
        /**
         * @brief Remove variants outside BED regions
         * @param mergedChrVarinat Variant data container
         */
        void removeVariantsOutBedRegion(
            std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat
        );

        /**
         * @brief Write BED region log
         * @param chrVec Vector of chromosome names
         * @param mergedChrVarinat Variant data container
         * @param outPrefix Output file prefix
         */
        void writeBedRegionLog(const std::vector<std::string>& chrVec, 
                        const std::map<std::string, std::map<int, MultiGenomeVar>>& mergedChrVarinat,
                        const std::string& outPrefix);

        /**
         * @brief Get metrics pointer for multi-threaded parallel processing
         * @param chr Chromosome name
         * @return Pointer to chromosome metrics
         */
        SomaticReadMetrics* getMetricsPtr(const std::string &chr);

        /**
         * @brief Write position allele count log
         * @param chrVec Vector of chromosome names
         * @param outputFileName Output file name
         * @param mergedChrVarinat Variant data container
         */
        void writePosAlleleCountLog(
            std::vector<std::string>& chrVec,
            std::string outputFileName,
            std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat
        );
        
        /**
         * @brief Write tagged somatic read report
         * @param chrVec Vector of chromosome names
         * @param outputFileName Output file name
         */
        void writeTaggedSomaticReadReport(
            const std::vector<std::string>& chrVec,
            std::string outputFileName
        );
        
        /**
         * @brief Write total truth somatic read report
         * @param chrVec Vector of chromosome names
         * @param outputFileName Output file name
         */
        void writeTotalTruthSomaticReadReport(
            const std::vector<std::string>& chrVec,
            std::string outputFileName
        );
        
        /**
         * @brief Write tagged read report
         * @param chrVec Vector of chromosome names
         * @param outputFileName Output file name
         */
        void writeTaggedReadReport(
            const std::vector<std::string>& chrVec,
            std::string outputFileName
        );
        
        /**
         * @brief Display somatic variant count
         * @param chrVec Vector of chromosome names
         * @param mergedChrVarinat Variant data container
         */
        void displaySomaticVarCount(std::vector<std::string> &chrVec, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
        
        /**
         * @brief Display BED region count
         * @param chrVec Vector of chromosome names
         */
        void displayBedRegionCount(std::vector<std::string> &chrVec);
};

#endif
