#ifndef HAPLOTAG_VCF_PARSER_H
#define HAPLOTAG_VCF_PARSER_H

#include "HaplotagType.h"

/**
 * @brief VCF parser operation modes
 * 
 * VCF_PARSER_LOAD_NODE: Parse and load variant data from VCF files
 * VCF_PARSER_WRITE_NODE: Write processed variant data to output VCF files
 */
enum VcfParserMode{
    VCF_PARSER_LOAD_NODE = 0,
    VCF_PARSER_WRITE_NODE
};

/**
 * @brief VCF file parser for haplotype tagging and somatic variant analysis
 * 
 * This class handles parsing of VCF files containing phased variants, structural variants (SV),
 * and modification variants (MOD). It supports both compressed (.vcf.gz) and uncompressed (.vcf) files.
 * 
 * Key functionalities:
 * - Parse phased SNP variants with haplotype information (PS field)
 * - Process structural variants and track read-haplotype associations
 * - Handle modification variants for methylation analysis
 * - Write processed results to output VCF files
 * - Support both normal and tumor sample processing
 * 
 * Used by:
 * - HaplotagProcess: For normal sample haplotype analysis
 * - SomaticHaplotagProcess: For somatic variant analysis in tumor samples
 * - PurityEstimationProcess: For tumor purity estimation
 * - SomaticReadBenchmark: For benchmarking somatic variant detection (inherits from VcfParser)
 */
class VcfParser{
    private:
        // Current operation mode (load or write)
        VcfParserMode mode;
        // Sample type being processed (NORMAL or TUMOR)
        Genome tagSample;
        // Flag to enable SNP variant parsing
        bool parseSnpFile;
        // Flag to enable structural variant parsing
        bool parseSVFile;
        // Flag to enable modification variant parsing
        bool parseMODFile;
        // Flag indicating if PS field is integer type (vs string)
        bool integerPS;
        // Mapping of string PS values to integer indices for string-type PS fields
        std::map<std::string, int> psIndex;

        // Output VCF file stream for writing results
        std::ofstream* resultVcf;
        // Flag to control command line header writing
        bool writeCommandline;

        // Command line string to write in VCF header
        std::string commandline;
        // Version string to write in VCF header
        std::string version;

        /**
         * @brief Parse compressed VCF file (.vcf.gz)
         */
        void compressParser(std::string &variantFile, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants);
        
        /**
         * @brief Parse uncompressed VCF file (.vcf)
         */
        void unCompressParser(std::string &variantFile, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants);
        
        /**
         * @brief Process a single VCF line based on current mode
         */
        void processLine(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants);
        
        /**
         * @brief Parse and load variant data from VCF line
         */
        virtual void parserProcess(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants);
        
        /**
         * @brief Write processed variant data to output VCF
         */
        virtual void writeProcess(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants);
    protected:

    public:
        VcfParser(Genome tagSample=Genome::NORMAL);
        virtual ~VcfParser();
        
        /**
         * @brief Set parser operation mode
         * @param mode VCF_PARSER_LOAD_NODE or VCF_PARSER_WRITE_NODE
         */
        void setMode(VcfParserMode mode);
        
        /**
         * @brief Set command line string for VCF header
         * @param commandline Command line string
         */
        void setCommandLine(std::string &commandline);
        
        /**
         * @brief Set version string for VCF header
         * @param version Version string
         */
        void setVersion(std::string &version);
        
        /**
         * @brief Enable/disable SNP variant parsing
         * @param parseSnpFile True to enable SNP parsing
         */
        void setParseSnpFile(bool parseSnpFile);
        
        /**
         * @brief Enable/disable structural variant parsing
         * @param parseSVFile True to enable SV parsing
         */
        void setParseSVFile(bool parseSVFile);
        
        /**
         * @brief Enable/disable modification variant parsing
         * @param parseMODFile True to enable MOD parsing
         */
        void setParseMODFile(bool parseMODFile);
        
        /**
         * @brief Get current SNP parsing status
         * @return True if SNP parsing is enabled
         */
        bool getParseSnpFile();
        
        /**
         * @brief Reset parser state and clear internal data
         */
        void reset();
        
        /**
         * @brief Parse VCF file and load variant data
         */
        void parsingVCF(std::string &variantFile, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants);
        
        /**
         * @brief Write processed variant data to output VCF file
         */
        void writingResultVCF(
            std::string &variantFile,
            VCF_Info &Info,
            std::map<std::string, std::map<int, MultiGenomeVar>> &chrMultiVariants,
            const std::string &outputPrefix
        );
};

#endif
