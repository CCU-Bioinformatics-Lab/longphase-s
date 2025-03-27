#ifndef SOMATIC_BENCHMARK_H
#define SOMATIC_BENCHMARK_H

#include "Util.h"
#include "HaplotagBase.h"

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

class SomaticReadVerifier: public VcfParser{
    private:
        struct RefAltDelCount{
            int refCount;
            int altCount;
            int delCount;
        };
        
        bool openTestingFunc;

        // store data
        std::map<std::string, std::map<int, RefAltDelCount>> posAltRefDelCount;
        std::vector<std::pair<int, int>> highConSomaticPos;
        std::vector<SomaticReadLog> totalReadVec;
        std::vector<SomaticReadLog> readsCrossingHighConSnpVec;
        std::vector<SomaticReadLog> taggedSomaticReadVec;
        void parserProcess(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
    public:

        SomaticReadVerifier();
        ~SomaticReadVerifier();
        void setTestingFunc(bool openTestingFunc);
        void loadHighConSomatic(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
        
        void recordDelReadCount(const std::string &chr, std::map<int, MultiGenomeVar>::iterator &currentVariantIter);
        void recordRefAltAlleleCount(const std::string &chr, std::string &base, std::map<int, MultiGenomeVar>::iterator &currentVariantIter);
        void recordCrossingHighConSnpRead(const std::string &chr, std::string &readID, std::string &hpResult, std::map<int, int> &variantsHP, std::map<int, int> &hpCount, double &norHPsimilarity, float &deriveByHpSimilarity, std::map<int, MultiGenomeVar> &currentChrVariants);
        void recordTaggedSomaticRead(const std::string &chr, std::string &readID, std::string &hpResult, std::map<int, int> &variantsHP, std::map<int, int> &hpCount, double &norHPsimilarity, float &deriveByHpSimilarity, std::map<int, MultiGenomeVar> &currentChrVariants);
        void recordTaggedRead(const std::string &chr, std::string &readID, std::string &hpResult, std::map<int, int> &variantsHP, std::map<int, int> &hpCount, double &norHPsimilarity, float &deriveByHpSimilarity, std::map<int, MultiGenomeVar> &currentChrVariants);

        SomaticReadLog createBasicSomaticReadLog(const std::string &chr, std::string &readID, std::string &hpResult, double &norHPsimilarity, float &deriveByHpSimilarity, std::map<int, int> &hpCount);
        
        void writePosAlleleCountLog(std::vector<std::string> &chrVec, HaplotagParameters &params, std::string logPosfix, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
        void writeTaggedSomaticReadLog(HaplotagParameters &params, std::string logPosfix);
        void writeCrossHighConSnpReadLog(HaplotagParameters &params, std::string logPosfix);
        void writeTaggedReadLog(HaplotagParameters &params, std::string logPosfix);
        void writeReadLog(HaplotagParameters &params, std::string logPosfix, std::vector<SomaticReadLog> &somaticReadVec);
        
        void displaySomaticVarCount(std::vector<std::string> &chrVec, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
};

#endif
