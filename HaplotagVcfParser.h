#ifndef HAPLOTAG_VCF_PARSER_H
#define HAPLOTAG_VCF_PARSER_H

#include "HaplotagType.h"

enum VcfParserMode{
    VCF_PARSER_LOAD_NODE = 0,
    VCF_PARSER_WRITE_NODE
};

class VcfParser{
    private:
        
        VcfParserMode mode;
        Genome tagSample;
        bool parseSnpFile;
        bool parseSVFile;
        bool parseMODFile;
        bool integerPS;
        std::map<std::string, int> psIndex;

        std::ofstream* resultVcf;
        bool writeCommandline;

        std::string commandline;
        std::string version;

        void compressParser(std::string &variantFile, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
        void unCompressParser(std::string &variantFile, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
        void processLine(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
        virtual void parserProcess(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
        virtual void writeProcess(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
    protected:

    public:
        VcfParser(Genome tagSample=Genome::NORMAL);
        virtual ~VcfParser();
        void setMode(VcfParserMode mode);
        void setCommandLine(std::string &commandline);
        void setVersion(std::string &version);
        void setParseSnpFile(bool parseSnpFile);
        void setParseSVFile(bool parseSVFile);
        void setParseMODFile(bool parseMODFile);
        bool getParseSnpFile();
        void reset();
        void parsingVCF(std::string &variantFile, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
        void writingResultVCF(
            std::string &variantFile,
            VCF_Info &Info,
            std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat,
            const std::string &outputPrefix
        );
};

#endif
