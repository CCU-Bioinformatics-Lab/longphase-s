#ifndef HAPLOTAG_VCF_PARSER_H
#define HAPLOTAG_VCF_PARSER_H

#include "HaplotagType.h"


class VcfParser{
    private:
        Genome tagGeneType;
        bool parseSnpFile;
        bool parseSVFile;
        bool parseMODFile;
        bool integerPS;
        std::map<std::string, int> psIndex;

        void compressParser(std::string &variantFile, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
        void unCompressParser(std::string &variantFile, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
        virtual void parserProcess(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
    
    protected:

    public:
        VcfParser();
        VcfParser(Genome tagGeneType);
        virtual ~VcfParser();
        void setParseSnpFile(bool parseSnpFile);
        void setParseSVFile(bool parseSVFile);
        void setParseMODFile(bool parseMODFile);
        bool getParseSnpFile();
        void reset();
        void variantParser(std::string &variantFile, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
};

#endif
