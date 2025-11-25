#ifndef HAPLOTAG_PARSING_BAM_H
#define HAPLOTAG_PARSING_BAM_H

#include "HaplotagType.h"
#include "HaplotagStrategy.h"
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/khash.h>
#include <htslib/kbitset.h>
#include <htslib/thread_pool.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <zlib.h>

class HaplotagBamParser;
class ChromosomeProcessor;
class CigarParser;

/**
 * @brief Enumeration for BAM parsing modes
 * 
 * Defines whether to use single-threaded or multi-threaded processing
 * for BAM file parsing and haplotagging operations
 */
enum ParsingBamMode{
    SINGLE_THREAD = 0,
    MULTI_THREAD = 1
};

/**
 * @brief Configuration structure for BAM parsing parameters
 * 
 * Contains all the parameters needed for BAM file processing and haplotagging
 * including thread count, quality thresholds, output settings, and processing options
 */
struct ParsingBamConfig{
    int numThreads;              /** Number of threads for parallel processing */
    int qualityThreshold;        /** Mapping quality threshold for read filtering */
    double percentageThreshold;  /** Percentage threshold for variant calling */
    std::string resultPrefix;    /** Output file prefix */
    std::string region;          /** Genomic region to process (optional) */
    std::string command;         /** Command line string for BAM header */
    std::string version;         /** Program version string */
    std::string outputFormat;    /** Output format (bam/cram) */

    bool tagSupplementary;       /** Whether to tag supplementary alignments */
    bool writeReadLog;           /** Whether to write detailed read logs */
};

/**
 * @brief Control structure for BAM parsing behavior
 * 
 * Controls the processing mode and output behavior for BAM file operations
 * including thread mode selection and output file generation
 */
struct ParsingBamControl{
    ParsingBamMode mode = ParsingBamMode::MULTI_THREAD;  /** Processing mode (single/multi-thread) */
    bool writeOutputBam = false;                         /** Whether to write output BAM file */
    bool mappingQualityFilter = false;                   /** Whether to apply mapping quality filtering */
};

/**
 * @brief Context structure for BAM parser operations
 * 
 * Contains all the data structures and files needed for BAM processing
 * including input files, chromosome information, variant data, and genome type
 */
struct BamParserContext{
    const std::string &BamFile;
    const std::string &fastaFile;
    const std::vector<std::string> &chrVec; 
    const std::map<std::string, int> &chrLength; 
    std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants; 
    std::map<Genome, VCF_Info> &vcfSet;
    const Genome genomeSample;
    
    BamParserContext(
        const std::string &BamFile, 
        const std::string &fastaFile,
        const std::vector<std::string> &chrVec, 
        const std::map<std::string, int> &chrLength, 
        std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants, 
        std::map<Genome, VCF_Info> &vcfSet,
        const Genome genomeSample
    ): 
    BamFile(BamFile),
    fastaFile(fastaFile),
    chrVec(chrVec),
    chrLength(chrLength),
    chrMultiVariants(chrMultiVariants),
    vcfSet(vcfSet),
    genomeSample(genomeSample)
    {}
};

/**
 * @brief Context structure for chromosome processing
 * 
 * Contains data specific to processing a single chromosome
 */
struct ChrProcContext{
    const std::string &chrName; 
    const int& chrLength;
    const ParsingBamConfig &params; 
    const Genome genomeSample;
    std::map<Genome, VCF_Info> &vcfSet; 

    ChrProcContext(
        const std::string &chrName,
        const int& chrLength,
        const ParsingBamConfig &params,
        const Genome genomeSample,
        std::map<Genome, VCF_Info> &vcfSet
    ):
    chrName(chrName),
    chrLength(chrLength),
    params(params),
    genomeSample(genomeSample),
    vcfSet(vcfSet)
    {}
};

/**
 * @brief Context structure for CIGAR parsing operations
 * 
 * Contains all the data needed for parsing CIGAR strings and processing alignments
 */
struct CigarParserContext{
    const bam1_t& aln;
    const bam_hdr_t& bamHdr;
    const std::string& chrName;
    const ParsingBamConfig& params;
    std::map<int, MultiGenomeVar>::iterator& firstVariantIter;
    std::map<int, MultiGenomeVar>& currentVariants;
    const std::string& ref_string;
    
    CigarParserContext(
        const bam1_t& aln,
        const bam_hdr_t& bamHdr,
        const std::string& chrName,
        const ParsingBamConfig& params,
        std::map<int, MultiGenomeVar>::iterator& firstVariantIter,
        std::map<int, MultiGenomeVar>& currentVariants,
        const std::string& ref_string
    ):
    aln(aln),
    bamHdr(bamHdr),
    chrName(chrName),
    params(params),
    firstVariantIter(firstVariantIter),
    currentVariants(currentVariants),
    ref_string(ref_string)
    {}
};

/**
 * @brief RAII wrapper for BAM file operations
 * 
 * Manages BAM file resources and provides safe file handling with automatic cleanup.
 * Handles input/output BAM files, headers, indexes, and alignment records.
 * Ensures proper resource management and error handling.
 */
class BamFileRAII {
    private:
        bool writeOutputBam;     /** Whether to write output BAM file */
        bool isReleased;         /** Whether resources have been released */

        /**
         * @brief Template function to check for null pointers
         * @param ptr Pointer to check
         * @param errorMessage Error message to display if pointer is null
         * 
         * Throws runtime_error if the pointer is null
         */
        template<typename T>
        void checkNullPointer(const T* ptr, const std::string& errorMessage) const;
        
    public:
        samFile* in;             /** Input BAM file handle */
        samFile* out;            /** Output BAM file handle */
        bam_hdr_t* bamHdr;       /** BAM header information */
        hts_idx_t* idx;          /** BAM index for random access */
        bam1_t* aln;             /** BAM alignment record */

        BamFileRAII(
            const std::string& BamFile,
            const std::string& fastaFile,
            htsThreadPool &threadPool,
            const std::string& version,
            const std::string& command,
            const std::string& resultPrefix,
            const std::string& outputFormat,
            const bool writeOutputBam = false
        );
        
        /**
         * @brief Destructor for BamFileRAII
         * 
         * Ensures proper cleanup of BAM file resources
         */
        ~BamFileRAII();

        /**
         * @brief Validates the state of BAM file resources
         * @return true if all resources are valid, exits program if not
         * 
         * Checks that all BAM file pointers are valid and not null
         */
        bool validateState();
        
        /**
         * @brief Writes the current alignment to the output BAM file
         * 
         * Writes the alignment record to the output file and handles errors
         */
        void samWriteBam();

        /**
         * @brief Destroys all BAM file resources
         * 
         * Frees memory and closes file handles
         */
        void destroy();
};

/**
 * @brief Base class for BAM file parsing and haplotagging
 * 
 * Provides the framework for parsing BAM files and performing haplotagging operations.
 * Supports both single-threaded and multi-threaded processing modes.
 * Uses factory pattern to create chromosome-specific processors.
 * 
 * Key functionalities:
 * - Parse BAM files with configurable parameters
 * - Support parallel processing across chromosomes
 * - Handle different processing modes (single/multi-thread)
 * - Manage BAM file resources and error handling
 */
class HaplotagBamParser{
    private:
        /**
         * @brief Processes BAM file using parallel processing
         * @param ctx BAM parser context
         * @param fastaParser Reference sequence parser
         * @param threadPool Thread pool for parallel execution
         * 
         * Uses OpenMP to process chromosomes in parallel
         */
        void processBamParallel(
            BamParserContext& ctx,
            const FastaParser &fastaParser,
            htsThreadPool &threadPool
        );

        /**
         * @brief Processes BAM file using single-threaded processing with output
         * @param ctx BAM parser context
         * @param fastaParser Reference sequence parser
         * @param threadPool Thread pool for parallel execution
         * 
         * Processes chromosomes sequentially and writes output BAM file
         */
        void processBamWithOutput(
            BamParserContext& ctx,
            const FastaParser &fastaParser,
            htsThreadPool &threadPool
        );
    protected: 
        const ParsingBamConfig& config;
        const ParsingBamControl& control;

        /**
         * @brief Factory method to create a chromosome processor
         * @param chr Chromosome name
         * @return Unique pointer to chromosome processor
         * 
         * Pure virtual function that derived classes must implement
         */
        virtual std::unique_ptr<ChromosomeProcessor> createProcessor(const std::string &chr) = 0;

        /**
         * @brief Gets the last variant position for each chromosome
         * @param last_pos Vector to store last variant positions
         * @param chrVec Vector of chromosome names
         * @param chrMultiVariants Map of variants by chromosome
         * @param geneSample Genome type (NORMAL/TUMOR)
         * 
         * Determines the last variant position for efficient processing
         */
        void getLastVarPos(
            std::vector<int>& last_pos,
            const std::vector<std::string>& chrVec,
            std::map<std::string,std::map<int, MultiGenomeVar>> &chrMultiVariants,
            const Genome& geneSample
        );

    public:
        /**
         * @brief Constructor for HaplotagBamParser
         * @param config Configuration parameters for BAM parsing
         * @param control Control parameters for processing behavior
         */
        HaplotagBamParser(
            const ParsingBamConfig &config, 
            const ParsingBamControl &control
        );
        
        /**
         * @brief Virtual destructor for HaplotagBamParser
         */
        virtual ~HaplotagBamParser();

        /**
         * @brief Main function for parsing BAM files and performing haplotagging
         * @param ctx Context containing all necessary data structures and files
         * 
         * Orchestrates the entire BAM parsing and haplotagging process
         */
        void parsingBam(BamParserContext& ctx);
};

/**
 * @brief Base class for processing individual chromosomes
 * 
 * Handles the processing of reads within a single chromosome.
 * Provides virtual methods for different types of read processing.
 * Supports both output BAM generation and mapping quality filtering.
 * 
 * Key functionalities:
 * - Process reads from a single chromosome
 * - Handle different read types (mapped, unmapped, secondary, supplementary)
 * - Apply mapping quality filtering
 * - Generate output BAM files
 * - Support post-processing operations
 */
class ChromosomeProcessor{
    private:
        bool writeOutputBam;         /** Whether to write output BAM file */
        bool mappingQualityFilter;   /** Whether to filter by mapping quality */
    protected:

        virtual void processLowMappingQuality(){};
        virtual void processUnmappedRead(){};
        virtual void processSecondaryAlignment(){};
        virtual void processSupplementaryAlignment(){};
        virtual void processEmptyVariants(){};
        virtual void processOtherCase(){};

        /**
         * @brief Process a single read
         * @param aln BAM alignment record
         * @param bamHdr BAM header information
         * @param ref_string Reference sequence string
         * @param currentVariants Current chromosome variants
         * @param firstVariantIter Iterator to first variant
         * @param ctx Chromosome processing context
         * 
         * Pure virtual method that derived classes must implement
         */
        virtual void processRead(
            bam1_t &aln, 
            const bam_hdr_t &bamHdr,
            const std::string &ref_string,
            std::map<int, MultiGenomeVar> &currentVariants,
            std::map<int, MultiGenomeVar>::iterator &firstVariantIter,
            ChrProcContext& ctx
        ) = 0;

        /**
         * @brief Post-processing operations
         * @param chr Chromosome name
         * @param currentVariants Current chromosome variants
         * 
         * Virtual method for post-processing operations after chromosome completion
         */
        virtual void postProcess(
            const std::string &chr,
            std::map<int, MultiGenomeVar> &currentVariants
        ){};
        
    public:
        /**
         * @brief Constructor for ChromosomeProcessor
         * @param writeOutputBam Whether to write output BAM file
         * @param mappingQualityFilter Whether to filter by mapping quality
         */
        ChromosomeProcessor( bool writeOutputBam=false, bool mappingQualityFilter=false);
        
        /**
         * @brief Virtual destructor for ChromosomeProcessor
         */
        virtual ~ChromosomeProcessor();

        /**
         * @brief Processes a single chromosome
         * @param ctx Chromosome processing context
         * @param bamRAII BAM file RAII wrapper
         * @param fastaParser Reference sequence parser
         * @param chrMultiVariants Map of variants by chromosome
         * 
         * Main function for processing all reads in a chromosome
         */
        void processSingleChrom(
            ChrProcContext& ctx,
            BamFileRAII& bamRAII,
            const FastaParser& fastaParser,
            std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants        
        );
};

/**
 * @brief Base class for parsing CIGAR strings
 * 
 * Provides the framework for parsing CIGAR operations and processing alignments.
 * Handles different CIGAR operations (match, insertion, deletion, etc.) and
 * determines haplotype assignments based on variant positions.
 * 
 * Key functionalities:
 * - Parse CIGAR strings from BAM alignments
 * - Handle different CIGAR operations (M, I, D, S, H, N, P, =, X)
 * - Track reference and query positions
 * - Count haplotype assignments and phase sets
 * - Support base nucleotide counting with quality filtering
 */
class CigarParser{
    private:
    protected:
        // Common data members that derived classes might need
        CigarParserContext ctx;

        std::map<int, int>* hpCount;
        std::map<int, int>* norCountPS;
        std::map<int, int>* variantsHP;

        // position relative to reference
        int& ref_pos;
        // position relative to read
        int& query_pos;

        std::map<int, MultiGenomeVar>::iterator currentVariantIter;

        /**
         * @brief Process match operation (M, =, X)
         * @param length Length of the match operation
         * @param cigar CIGAR array
         * @param i Current CIGAR index
         * @param aln_core_n_cigar Total number of CIGAR operations
         * @param base Base nucleotide at current position
         * 
         * Virtual method for handling match operations
         */
        virtual void processMatchOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, std::string& base, bool& isAlt, int& offset){};
        
        /**
         * @brief Process insertion operation (I)
         * @param length Length of the insertion
         * 
         * Virtual method for handling insertion operations
         */
        virtual void processInsertionOperation(int& length){};
        
        /**
         * @brief Process deletion operation (D)
         * @param length Length of the deletion
         * @param cigar CIGAR array
         * @param i Current CIGAR index
         * @param aln_core_n_cigar Total number of CIGAR operations
         * @param alreadyJudgeDel Whether deletion has already been judged
         * 
         * Virtual method for handling deletion operations
         */
        virtual void processDeletionOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, bool& alreadyJudgeDel){};
        
        /**
         * @brief Process skipped operation (N)
         * @param length Length of the skipped region
         * 
         * Virtual method for handling skipped operations
         */
        virtual void processSkippedOperation(int& length){};
        
        /**
         * @brief Process soft clipping operation (S)
         * @param length Length of the soft clipping
         * 
         * Virtual method for handling soft clipping operations
         */
        virtual void processSoftClippingOperation(int& length){};
        
        /**
         * @brief Process hard clipping operation (H)
         * 
         * Virtual method for handling hard clipping operations
         */
        virtual void processHardClippingOperation(){};
        
        /**
         * @brief Count base nucleotides and apply mapping quality filtering
         * @param posBase Position base information structure
         * @param base Base nucleotide
         * @param aln BAM alignment record
         * @param mpqThreshold Mapping quality threshold
         * 
         * Updates base counts and filtered depth based on mapping quality
         */
        void countBaseNucleotide(PosBase& posBase, std::string& base, const bam1_t& aln, const float& mpqThreshold, bool isAlt, HaplotagVariantType::VariantType variantType);
        
        /**
         * @brief Count deletion bases
         * @param posBase Position base information structure
         * 
         * Updates deletion count and total depth
         */
        void countDeletionBase(PosBase& posBase);

        bool IsAltIndel(
            int& ref_pos, 
            int& length, 
            int& i, 
            int& aln_core_n_cigar, 
            uint32_t* cigar,
            std::map<int, MultiGenomeVar>::iterator& currentVariantIter, 
            std::string& base, 
            HaplotagVariantType::VariantType variantType,
            Genome sample
        );

    public:
        /**
         * @brief Constructor for CigarParser
         * @param ctx CIGAR parsing context
         * @param ref_pos Reference position (will be modified)
         * @param query_pos Query position (will be modified)
         */
        CigarParser(CigarParserContext& ctx, int& ref_pos, int& query_pos);

        /**
         * @brief Virtual destructor for CigarParser
         */
        virtual ~CigarParser();

        /**
         * @brief Parses CIGAR string and processes alignment operations
         * @param hpCount Map to count haplotype assignments
         * @param variantsHP Map to record variant haplotype assignments
         * @param norCountPS Map to count phase set assignments
         * 
         * Main function for parsing CIGAR operations and determining haplotype assignments
         */
        void parsingCigar(
            std::map<int, int>& hpCount,
            std::map<int, int>& variantsHP,
            std::map<int, int>& norCountPS
        );
};

#endif
