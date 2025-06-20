#ifndef SOMATIC_VAR_CALLER_H
#define SOMATIC_VAR_CALLER_H

#include "Util.h" 
#include "HaplotagType.h"
#include "HaplotagParsingBam.h"
#include "HaplotagStrategy.h"
#include "TumorPurityEstimator.h"
#include "HaplotagLogging.h"

/**
 * @brief Context information for somatic variant calling
 * 
 * Contains all input file paths and configuration needed for somatic variant calling
 * including normal/tumor BAM files, SNP files, and reference FASTA
 */
struct CallerContext
{
    std::string normalBamFile;
    std::string tumorBamFile;
    std::string normalSnpFile;
    std::string tumorSnpFile;
    std::string fastaFile;
    CallerContext(std::string normalBamFile, std::string tumorBamFile, std::string normalSnpFile, std::string tumorSnpFile, std::string fastaFile)
    :normalBamFile(normalBamFile), tumorBamFile(tumorBamFile), normalSnpFile(normalSnpFile), tumorSnpFile(tumorSnpFile), fastaFile(fastaFile){}
};

/**
 * @brief Configuration parameters for somatic variant calling
 * 
 * Controls filtering behavior, tumor purity estimation, and logging options
 */
struct CallerConfig
{
    bool enableFilter = true;           /** Whether to enable somatic variant filtering */
    bool estimateTumorPurity = true;    /** Whether to estimate tumor purity automatically */
    double tumorPurity = 0.0;          /** Manual tumor purity value (0.0-1.0) */
    bool writeCallingLog = false;       /** Whether to write detailed calling logs */
    
    CallerConfig() = default;
    
    /**
     * @brief Constructor for CallerConfig
     * @param filter Enable filtering
     * @param isEstimate Enable tumor purity estimation
     * @param purity Manual tumor purity value
     * @param writeLog Enable detailed logging
     */
    CallerConfig(bool filter, bool isEstimate, double purity, bool writeLog)
        : enableFilter(filter), estimateTumorPurity(isEstimate), tumorPurity(purity), writeCallingLog(writeLog) {}
};

/**
 * @brief Filter parameters for somatic variant calling
 * 
 * Contains all threshold values used for filtering somatic variants
 * based on tumor purity, read quality, and haplotype consistency
 */
struct SomaticVarFilterParams
{  
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

/**
 * @brief Haplotype count information for a single read
 * 
 * Tracks the number of variants of each haplotype type found on a read
 * and stores read metadata for analysis
 */
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

/**
 * @brief Data structure for dense SNP interval analysis
 * 
 * Stores statistical information about SNPs in dense intervals
 * including mean values, z-scores, and distances
 */
struct DenseSnpData{
    // Mean alternative allele count
    double snpAltMean;
    // Z-score for statistical analysis
    double snpZscore;
    // Minimum distance to neighboring SNP
    int minDistance;
    
    DenseSnpData(): snpAltMean(0.0), snpZscore(0.0), minDistance(0){}
};

/**
 * @brief Container for dense SNP interval analysis
 * 
 * Manages collections of SNPs within dense intervals and their statistical properties
 * including mean calculations, z-scores, and interval metadata
 */
struct DenseSnpInterval{
    std::map<int, double> snpAltMean;
    std::map<int, double> snpZscore;
    std::map<int, int> minDistance;
    int snpCount;
    double totalAltMean;
    double StdDev;
    DenseSnpInterval(): snpCount(0), totalAltMean(0.0), StdDev(0.0){}

    /**
     * @brief Clears all data in the interval
     * 
     * Resets all maps and counters to initial state
     */
    void clear(){
        snpAltMean.clear();
        snpZscore.clear();
        minDistance.clear();
        snpCount = 0;
        totalAltMean = 0.0;
        StdDev = 0.0;
    }
};

/**
 * @brief Namespace for tumor-normal analysis utilities
 * 
 * Contains common functions used for analyzing both tumor and normal samples
 */
namespace tumor_normal_analysis{
    /**
     * @brief Calculate common base information for tumor SNPs
     * 
     * Computes VAF, depth ratios, and haplotype imbalance metrics
     * @param baseInfo Base information structure to update
     * @param tumorAltBase Alternative base in tumor sample
     */
    void calculateBaseCommonInfo(PosBase& baseInfo, std::string& tumorAltBase);
};

namespace statisticsUtils{
    double calculateMean(const std::map<int, double>& data);
    double calculateStandardDeviation(const std::map<int, double>& data, double mean);
    void calculateZScores(const std::map<int, double>& data, double mean, double stdDev, std::map<int, double> &zScores);
};

/**
 * @brief Chromosome processor for extracting normal sample data
 * 
 * Processes reads from normal BAM file to extract base counts, depths,
 * and haplotype information for somatic variant calling
 */
class ExtractNorDataChrProcessor : public ChromosomeProcessor{
    private:
        // store base information
        std::map<int, PosBase> *variantBase;
        GermlineHaplotagStrategy judger;
    protected:
        //override processRead
        void processRead(
            bam1_t &aln, 
            const bam_hdr_t &bamHdr,
            const std::string &ref_string,
            std::map<int, MultiGenomeVar> &currentVariants,
            std::map<int, MultiGenomeVar>::iterator &firstVariantIter,
            ChrProcContext& ctx
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

/**
 * @brief CIGAR parser for normal sample data extraction
 * 
 * Specialized CIGAR parser that handles normal sample variant processing
 * and haplotype determination for germline variants
 */
class ExtractNorDataCigarParser : public CigarParser{
    private:
        //specific data members
        std::map<int, PosBase>& variantBase;
        std::vector<int>& tumVarPosVec;
        const int& mappingQualityThr;
        GermlineHaplotagStrategy judger;
    protected:

        void processMatchOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, std::string& base) override;
        void processDeletionOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, bool& alreadyJudgeDel) override;

    public:
        ExtractNorDataCigarParser(
            CigarParserContext& ctx,
            std::map<int, PosBase>& variantBase,
            std::vector<int>& tumVarPosVec,
            int& ref_pos, 
            int& query_pos,
            const int& mappingQualityThr
        );
        ~ExtractNorDataCigarParser() override;
};

/**
 * @brief BAM parser for normal sample data extraction
 * 
 * Manages parsing of normal BAM file to extract base information
 * and haplotype data for somatic variant calling
 */
class ExtractNorDataBamParser : public HaplotagBamParser{
    private:
        // chr, variant position (0-base), base count & depth
        std::map<std::string, std::map<int, PosBase>>& chrPosNorBase;

        // Test if the logic of judgeHaplotype is consistent (only single chromosome)
        std::ofstream *tagResult;

    protected:
        /**
         * @brief Factory method to create chromosome processor
         * @param chr Chromosome name
         * @return Unique pointer to chromosome processor
         */
        std::unique_ptr<ChromosomeProcessor> createProcessor(const std::string &chr) override{
            return std::unique_ptr<ChromosomeProcessor>(new ExtractNorDataChrProcessor(chrPosNorBase, chr));
        };

    public:
        ExtractNorDataBamParser(
            const ParsingBamConfig& config,
            const ParsingBamControl& control,
            std::map<std::string, std::map<int, PosBase>>& chrPosNorBase
        );
        ~ExtractNorDataBamParser();
        void displayPosInfo(std::string chr, int pos);
};

/**
 * @brief Chromosome processor for extracting tumor sample data
 * 
 * Processes reads from tumor BAM file to extract somatic variant information,
 * haplotype data, and read classification for somatic calling
 */
class ExtractTumDataChrProcessor : public ChromosomeProcessor{
    private:
        ExtractSomaticDataStragtegy somaticJudger;

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
            const std::string &ref_string,
            std::map<int, MultiGenomeVar> &currentVariants,
            std::map<int, MultiGenomeVar>::iterator &firstVariantIter,
            ChrProcContext& ctx
        ) override;


        void classifyReadsByCase(std::vector<int> &readPosHP3, std::map<int, int> &norCountPS, std::map<int, int> &hpCount, std::map<int, SomaticData> &somaticPosInfo);

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
        

class ExtractTumDataCigarParser : public CigarParser{
    private:
        ExtractSomaticDataStragtegy somaticJudger;
        GermlineHaplotagStrategy judger;

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

    public:
        ExtractTumDataCigarParser(
            CigarParserContext& ctx,
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

/**
 * @brief BAM parser for tumor sample data extraction
 * 
 * Manages parsing of tumor BAM file to extract somatic variant information
 * and haplotype data for somatic calling
 */
class ExtractTumDataBamParser : public HaplotagBamParser{
    private:
        std::map<std::string, std::map<int, SomaticData>>& chrPosSomaticInfo;
        std::map<std::string, std::map<std::string, ReadVarHpCount>>& chrReadHpResultSet;
        std::map<std::string, std::map<int, std::map<std::string, int>>>& chrTumorPosReadCorrBaseHP;
    protected:
        /**
         * @brief Factory method to create chromosome processor
         * @param chr Chromosome name
         * @return Unique pointer to chromosome processor
         */
        std::unique_ptr<ChromosomeProcessor> createProcessor(const std::string &chr) override{
            return std::unique_ptr<ChromosomeProcessor>(new ExtractTumDataChrProcessor(chrPosSomaticInfo, chrReadHpResultSet, chrTumorPosReadCorrBaseHP, chr));
        };
        
    public:
        ExtractTumDataBamParser(
            const ParsingBamConfig& config,
            const ParsingBamControl& control,
            std::map<std::string, std::map<int, SomaticData>>& chrPosSomaticInfo,
            std::map<std::string, std::map<std::string, ReadVarHpCount>>& chrReadHpResultSet,
            std::map<std::string, std::map<int, std::map<std::string, int>>>& chrTumorPosReadCorrBaseHP
        );
        ~ExtractTumDataBamParser();
};

/**
 * @brief Main class for somatic variant calling
 * 
 * Orchestrates the entire somatic variant calling pipeline including:
 * - Data extraction from normal and tumor BAM files
 * - Tumor purity estimation
 * - Somatic variant filtering and calling
 * - Statistical analysis and reporting
 */
class SomaticVarCaller{
    private:
        // Maximum distance for dense SNP intervals
        static constexpr int INTERVAL_SNP_MAX_DISTANCE = 5000;

        const CallerConfig &callerCfg;
        const ParsingBamConfig &bamCfg;
        const std::vector<std::string>& chrVec;

        // somatic calling filter params
        SomaticVarFilterParams somaticParams;

        ExtractSomaticDataStragtegy somaticJudger;
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

        void initialSomaticFilterParams();

        /**
         * @brief Set filter parameters based on tumor purity
         * 
         * Adjusts filter thresholds based on estimated tumor purity
         * @param somaticParams Filter parameters to update
         * @param tumorPurity Estimated tumor purity value
         */
        void SetFilterParamsWithPurity(SomaticVarFilterParams &somaticParams, double &tumorPurity);

        
        /**
         * @brief Calculate interval data statistics
         * @param isStartPos Flag indicating start position
         * @param startPos Start position of interval
         * @param pos Current position
         * @param denseSnp Dense SNP interval data
         * @param localDenseTumorSnpInterval Local dense tumor SNP intervals
         */
        void calculateIntervalData(bool &isStartPos, int &startPos, int &pos, DenseSnpInterval &denseSnp, std::map<int, std::pair<int, DenseSnpInterval>> &localDenseTumorSnpInterval);
        
        /**
         * @brief Get dense tumor SNP intervals
         * @param somaticPosInfo Somatic position information
         * @param readHpResultSet Read haplotype results
         * @param somaticPosReadHPCount Somatic position read haplotype counts
         * @param closeSomaticSnpInterval Close somatic SNP intervals
         */
        void getDenseTumorSnpInterval(std::map<int, SomaticData> &somaticPosInfo, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, std::map<std::string, int>> &somaticPosReadHPCount, std::map<int, std::pair<int, DenseSnpInterval>> &closeSomaticSnpInterval);

        /**
         * @brief Apply somatic feature filtering
         * @param somaticParams Filter parameters
         * @param currentChrVariants Current chromosome variants
         * @param chr Chromosome name
         * @param somaticPosInfo Somatic position information
         * @param tumorPurity Tumor purity value
         */
        void somaticFeatureFilter(const SomaticVarFilterParams &somaticParams, std::map<int, MultiGenomeVar> &currentChrVariants,const std::string &chr, std::map<int, SomaticData> &somaticPosInfo, double& tumorPurity);
        
        /**
         * @brief Calibrate read haplotype assignments
         * @param chr Chromosome name
         * @param somaticPosInfo Somatic position information
         * @param readHpResultSet Read haplotype results
         * @param somaticPosReadHPCount Somatic position read haplotype counts
         */
        void calibrateReadHP(const std::string &chr, std::map<int, SomaticData> &somaticPosInfo, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, std::map<std::string, int>> &somaticPosReadHPCount);
        
        /**
         * @brief Calculate read set haplotype results
         * @param chr Chromosome name
         * @param readHpResultSet Read haplotype results
         * @param somaticPosReadHPCount Somatic position read haplotype counts
         * @param percentageThreshold Percentage threshold for haplotype determination
         */
        void calculateReadSetHP(const std::string &chr, std::map<std::string, ReadVarHpCount> &readHpResultSet, std::map<int, std::map<std::string, int>> &somaticPosReadHPCount, const double& percentageThreshold);
        
        /**
         * @brief Statistic somatic position read haplotype
         * @param chr Chromosome name
         * @param somaticPosInfo Somatic position information
         * @param somaticPosReadHPCount Somatic position read haplotype counts
         * @param readHpResultSet Read haplotype results
         * @param localReadHpDistri Local read haplotype distribution
         */
        void statisticSomaticPosReadHP(
            const std::string &chr,
            std::map<int, SomaticData> &somaticPosInfo,
            std::map<int, std::map<std::string, int>> &somaticPosReadHPCount,
            std::map<std::string, ReadVarHpCount> &readHpResultSet,
            chrReadHpResult &localReadHpDistri
        );
        
        /**
         * @brief Write somatic variant calling log
         * @param ctx Caller context
         * @param somaticParams Somatic filter parameters
         * @param chrVec Vector of chromosome names
         * @param chrMultiVariants Multi-genome chromosome variants
         */
        void writeSomaticVarCallingLog(const CallerContext &ctx, const SomaticVarFilterParams &somaticParams, const std::vector<std::string> &chrVec
                                     , std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants);
        
        
        /**
         * @brief Write dense tumor SNP interval log
         * @param logFileName Log file name
         * @param chrVec Vector of chromosome names
         */
        void writeDenseTumorSnpIntervalLog(const std::string logFileName, const std::vector<std::string> &chrVec);

        /**
         * @brief Release allocated memory
         * 
         * Deletes all dynamically allocated data structures
         */
        void releaseMemory();
        
        /**
         * @brief Find other somatic SNP haplotypes
          * @note Temporary function for analysis
         * Analyzes and identifies additional somatic SNP haplotypes (HP4, HP5)
         * @param chr Chromosome name
         * @param somaticPosInfo Somatic position information
         * @param currentChrVariants Current chromosome variants
         */
        void findOtherSomaticSnpHP(const std::string &chr, std::map<int, SomaticData> &somaticPosInfo, std::map<int, MultiGenomeVar> &currentChrVariants);

        /**
         * @brief Write other somatic haplotype log
         * @note Temporary function for analysis
         * @param logFileName Log file name
         * @param chrVec Vector of chromosome names
         * @param chrMultiVariants Multi-genome chromosome variants
         */
        void writeOtherSomaticHpLog(const std::string logFileName, const std::vector<std::string> &chrVec, std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants);

    protected:

    public:
        /**
         * @brief Constructor for SomaticVarCaller
         * @param callerCfg Caller configuration
         * @param bamCfg BAM parsing configuration
         * @param chrVec Vector of chromosome names
         */
        SomaticVarCaller(const CallerConfig &callerCfg, const ParsingBamConfig &bamCfg, const std::vector<std::string> &chrVec);
        virtual ~SomaticVarCaller();

        /**
         * @brief Main variant calling function
         * 
         * Orchestrates the entire somatic variant calling pipeline
         * @param ctx Caller context with input files
         * @param chrVec Vector of chromosome names
         * @param chrLength Map of chromosome names to lengths
         * @param chrMultiVariants Multi-genome chromosome variants
         * @param vcfSet VCF information by genome type
         */
        void variantCalling(
            const CallerContext &ctx,
            const std::vector<std::string> &chrVec,
            const std::map<std::string, int> &chrLength,
            std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants,
            std::map<Genome, VCF_Info> &vcfSet
        );

        /**
         * @brief Extract somatic data from BAM files
         * 
         * Processes normal and tumor BAM files to extract variant information
         * @param normalBamFile Path to normal BAM file
         * @param tumorBamFile Path to tumor BAM file
         * @param fastaFile Path to reference FASTA file
         * @param config Parsing configuration
         * @param chrVec Vector of chromosome names
         * @param chrLength Map of chromosome names to lengths
         * @param chrMultiVariants Multi-genome chromosome variants
         * @param vcfSet VCF information by genome type
         */
        void extractSomaticData(
            const std::string &normalBamFile,
            const std::string &tumorBamFile,
            const std::string &fastaFile,
            const ParsingBamConfig &config,
            const std::vector<std::string> &chrVec,
            const std::map<std::string, int> &chrLength,
            std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants,
            std::map<Genome, VCF_Info> &vcfSet           
        );
        
        /**
         * @brief Run tumor purity estimation
         * @param writeReadLog Whether to write read logs
         * @param resultPrefix Result file prefix
         * @return Estimated tumor purity value
         */
        double runTumorPurityEstimator(bool writeReadLog, const std::string resultPrefix);
        
        /**
         * @brief Get somatic variant flags
         * 
         * Marks variants as somatic based on calling results
         * @param chrVec Vector of chromosome names
         * @param chrMultiVariants Multi-genome chromosome variants
         */
        void getSomaticFlag(const std::vector<std::string> &chrVec, std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants);

        /**
         * @brief Display calling SNP count
         * 
         * Prints the number of somatic SNPs called for debugging
         */
        void displayCallingSnpCount();
};

#endif

