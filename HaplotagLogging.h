#ifndef HAPLOTAG_LOGGING_H
#define HAPLOTAG_LOGGING_H

#include "HaplotagType.h"

/**
 * @brief Structure to store read haplotype results for a chromosome
 * 
 * Contains haplotype distribution data for each variant position,
 * including read counts, derived haplotype information, and coverage regions
 */
struct chrReadHpResult{
    /** Map of variant positions to read haplotype results */
    std::map<int, ReadHpResult> posReadHpResult;

    /**
     * @brief Record read haplotype assignment for a specific position
     * @param pos Variant position (0-based)
     * @param hpResult Assigned haplotype for the read
     * @param BaseHP Base haplotype type (germline or somatic)
     */
    void recordReadHp(int &pos, int &hpResult, int &BaseHP);
    
    /**
     * @brief Record derived haplotype information for a position
     * @param pos Variant position (0-based)
     * @param deriveHP Derived haplotype (H1, H2, or NONE)
     * @param deriveHPsimilarity Similarity score for derived haplotype
     */
    void recordDeriveHp(int &pos, int &deriveHP, float deriveHPsimilarity);
    
    /**
     * @brief Record alignment coverage region for a variant position
     * @param curVarPos Current variant position
     * @param startPos Start position of coverage region
     * @param endPos End position of coverage region
     */
    void recordAlignCoverRegion(int& curVarPos, int &startPos, int &endPos);
};

/**
 * @brief Main logging class for read haplotype distribution analysis
 * 
 * This class manages the collection and output of haplotype distribution data
 * across multiple chromosomes. It provides both single-thread and multi-thread
 * safe interfaces for recording read haplotype assignments and generating
 * comprehensive distribution reports.
 * 
 * Used by:
 * - SomaticHaplotagProcess: For tracking haplotype distributions before/after inheritance
 * - SomaticVarCaller: For somatic variant calling analysis
 * 
 * Key functionalities:
 * - Record read haplotype assignments for each variant position
 * - Track derived haplotype information and similarity scores
 * - Generate comprehensive distribution reports
 * - Calculate coverage regions and ratios
 * - Support multi-threaded access with thread-safe pointers
 */
class ReadHpDistriLog{
    private :
        /**
         * @brief Structure to store coverage region information
         */
        struct coverRegionInfo{
            int startPos;  /** Start position of coverage region */
            int endPos;    /** End position of coverage region */
            int length;    /** Length of coverage region */

            coverRegionInfo(): startPos(0), endPos(0), length(0){}
        };
        /** Map of chromosome names to read haplotype results */
        std::map<std::string, chrReadHpResult> chrVarReadHpResult;
    protected:

    public:
        ReadHpDistriLog();
        ~ReadHpDistriLog();

        /**
         * @brief Initialize chromosome key for multi-threaded access
         * @param chr Chromosome name
         */
        void loadChrKey(const std::string &chr);
        
        /**
         * @brief Get thread-safe pointer to chromosome haplotype results
         * @param chr Chromosome name
         * @return Pointer to chromosome haplotype results
         */
        chrReadHpResult* getChrHpResultsPtr (const std::string &chr);

        /**
         * @brief Record read haplotype for a chromosome (single-thread use only)
         * @param chr Chromosome name
         * @param pos Variant position (0-based)
         * @param hpResult Assigned haplotype for the read
         * @param BaseHP Base haplotype type
         */
        void recordChrReadHp(const std::string &chr, int &pos, int &hpResult, int &BaseHP);
        
        /**
         * @brief Record derived haplotype for a chromosome (single-thread use only)
         * @param chr Chromosome name
         * @param pos Variant position (0-based)
         * @param deriveHP Derived haplotype
         * @param deriveHPsimilarity Similarity score
         */
        void recordChrDeriveHp(const std::string &chr, int &pos, int &deriveHP, float deriveHPsimilarity);
        
        /**
         * @brief Record alignment coverage region for a chromosome (single-thread use only)
         * @param chr Chromosome name
         * @param pos Variant position (0-based)
         * @param startPos Start position of coverage region
         * @param endPos End position of coverage region
         * TODO: [PENDING] This function needs review and validation
         */
        void recordChrAlignCoverRegion(const std::string &chr, int &pos, int &startPos, int &endPos);

        /**
         * @brief Write comprehensive read haplotype distribution report
         * @param logFileName Output file name
         * @param chrVec Vector of chromosome names to process
         */
        void writeReadHpDistriLog(const std::string logFileName, const std::vector<std::string> &chrVec);
        
        /**
         * @brief Write position coverage region report
         * @param logFileName Output file name
         * @param chrVec Vector of chromosome names to process
         */
        void writePosCoverRegionLog(const std::string logFileName, const std::vector<std::string> &chrVec);
        
        /**
         * @brief Write tagged read coverage region report with merged regions
         * @param logFileName Output file name
         * @param chrVec Vector of chromosome names to process
         * @param chrLength Map of chromosome names to their lengths
         */
        void writeTagReadCoverRegionLog(const std::string logFileName, const std::vector<std::string> &chrVec, std::map<std::string, int> &chrLength);
        
        /**
         * @brief Remove positions not derived by both H1 and H2 haplotypes
         * @param chrVec Vector of chromosome names to process
         */
        void removeNotDeriveByH1andH2pos(const std::vector<std::string> &chrVec);
};

/**
 * @brief Template class for haplotype read logging with parameterized types
 * 
 * This template class provides a flexible framework for logging haplotype
 * read data with different parameter types and log data structures.
 * It handles file I/O operations and provides virtual methods for
 * customizing header and data writing behavior.
 * 
 * Template parameters:
 * - ParamType: Type of parameters used for configuration
 * - LogType: Type of log data structure
 * 
 * Used by:
 * - GermlineTagLog: For germline haplotype tagging logs
 */
template<typename ParamType, typename LogType>
class HaplotagReadLog{
    protected:
        /** Reference to parameters for configuration */
        const ParamType& params;

        /**
         * @brief Check if output stream is valid
         * @throws std::runtime_error if stream is not valid
         */
        void checkStreamStatus() {
            if (!tagReadLog || !tagReadLog->is_open()) {
                throw std::runtime_error("Output stream is not valid");
            }
        }

        /**
         * @brief Pure virtual method to add parameter-specific messages to header
         */
        virtual void addParamsMessage() = 0;

        /**
         * @brief Pure virtual method to write basic column headers
         */
        virtual void writeBasicColumns() = 0;

    public:
        /** Output file stream for logging */
        std::ofstream* tagReadLog;

        /**
         * @brief Constructor
         * @param params Reference to parameters
         * @param fileName Output file name
         */
        HaplotagReadLog(const ParamType& params, std::string fileName) : params(params) {
            tagReadLog = new std::ofstream(fileName);
            if(!tagReadLog->is_open()){
                std::cerr<< "Fail to open write file: " << fileName << "\n";
                exit(1);
            }
        };

        /**
         * @brief Virtual destructor - cleanup file stream
         */
        virtual ~HaplotagReadLog(){
            tagReadLog->close();
            if(tagReadLog) delete tagReadLog;
        };

        /**
         * @brief Write header information to log file
         * 
         * Calls addParamsMessage() and writeBasicColumns() to generate
         * complete header with parameter information and column names
         */
        virtual void writeHeader() {
            addParamsMessage();
            writeBasicColumns();
        }

        /**
         * @brief Pure virtual method to write log data
         * @param data Log data to write
         */
        virtual void writeTagReadLog(LogType& data) = 0;
};

#endif
