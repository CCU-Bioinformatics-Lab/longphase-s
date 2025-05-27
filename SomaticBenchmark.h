#ifndef SOMATIC_BENCHMARK_H
#define SOMATIC_BENCHMARK_H

#include "Util.h"
#include "HaplotagType.h"
#include "HaplotagVcfParser.h"

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


struct SomaticReadMetrics{

    struct RefAltDelCount{
        int refCount;
        int altCount;
        int delCount;
    };

    std::map<int, RefAltDelCount> posAltRefDelCount;
    std::vector<std::pair<int, int>> highConSomaticPos;
    std::vector<SomaticReadLog> totalReadVec;
    std::vector<SomaticReadLog> readsCrossingHighConSnpVec;
    std::vector<SomaticReadLog> taggedSomaticReadVec;
};


class SomaticReadVerifier{
    private:
        bool openTestingFunc;

        SomaticReadMetrics *metrics;

        SomaticReadLog createBasicSomaticReadLog(const std::string &chr, std::string &readID, std::string &hpResult, double &norHPsimilarity, float &deriveByHpSimilarity, std::map<int, int> &hpCount);
    
    public:
        SomaticReadVerifier(bool openTestingFunc, SomaticReadMetrics *metrics);
        ~SomaticReadVerifier();

        void recordDelReadCount(const std::string &chr, std::map<int, MultiGenomeVar>::iterator &currentVariantIter);
        void recordRefAltAlleleCount(const std::string &chr, std::string &base, std::map<int, MultiGenomeVar>::iterator &currentVariantIter);
        void recordCrossingHighConSnpRead(const std::string &chr, std::string &readID, std::string &hpResult, std::map<int, int> &variantsHP, std::map<int, int> &hpCount, double &norHPsimilarity, float &deriveByHpSimilarity, std::map<int, MultiGenomeVar> &currentChrVariants);
        void recordTaggedRead(const std::string &chr, std::string &readID, std::string &hpResult, std::map<int, int> &variantsHP, std::map<int, int> &hpCount, double &norHPsimilarity, float &deriveByHpSimilarity, std::map<int, MultiGenomeVar> &currentChrVariants);
};

class SomaticReadBenchmark: public VcfParser{
    private:

        struct BedRegion {
            int start;
            int end;
        };
        
        bool openTestingFunc;
        // chr, metrics
        std::map<std::string, SomaticReadMetrics> chrMetrics;

        std::map<std::string, std::vector<BedRegion>> bedRegions;

        void parserProcess(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
        void setChrSomaticReadVecPtr(
            const std::string &chr,
            std::map<std::string, std::vector<SomaticReadLog>*> &somaticReadVecMap,
            std::vector<SomaticReadLog> &somaticReadVec
        );
        void writeReadLog(
            const std::vector<std::string>& chrVec,
            HaplotagParameters &params,
            std::string logPosfix,
            std::map<std::string, std::vector<SomaticReadLog>*> &somaticReadVecMap
        );

        void processBedLine(const std::string& line);

    public:

        SomaticReadBenchmark();
        ~SomaticReadBenchmark();
        void setTestingFunc(bool openTestingFunc);
        bool getTestingFunc();
        void loadChrKey(const std::string &chr);
        void loadHighConSomatic(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
        
        // parse benchmark bed file
        void parseBedFile(const std::string& bedFile);

        // mark variants in bed regions
        void markVariantsInBedRegions(
            std::vector<std::string> &chrVec,
            std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
        
        // remove variants out bed regions
        void removeVariantsOutBedRegion(
            std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat
        );

        void writeBedRegionLog(const std::vector<std::string>& chrVec, 
                        const std::map<std::string, std::map<int, MultiGenomeVar>>& mergedChrVarinat,
                        const std::string& outPrefix);

        // get metrics pointer for multiple threads parallel processing
        SomaticReadMetrics* getMetricsPtr(const std::string &chr);

        void writePosAlleleCountLog(
            std::vector<std::string>& chrVec,
            HaplotagParameters &params,
            std::string logPosfix,
            std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat
        );
        void writeTaggedSomaticReadLog(
            const std::vector<std::string>& chrVec,
            HaplotagParameters &params,
            std::string logPosfix
        );
        void writeTaggedTruthSomaticReadLog(
            const std::vector<std::string>& chrVec,
            HaplotagParameters &params,
            std::string logPosfix
        );
        void writeTaggedReadLog(
            const std::vector<std::string>& chrVec,
            HaplotagParameters &params,
            std::string logPosfix
        );
        
        void displaySomaticVarCount(std::vector<std::string> &chrVec, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
        void displayBedRegionCount(std::vector<std::string> &chrVec);
};

#endif
