#ifndef SOMATICHAPLOTAGPROCESS_H
#define SOMATICHAPLOTAGPROCESS_H

#include "../haplotag/HaplotagType.h"
#include "../haplotag/HaplotagProcess.h"
#include "../haplotag/HaplotagParsingBam.h"
#include "../haplotag/HaplotagVcfParser.h"
#include "../haplotag/HaplotagLogging.h"
#include "SomaticVarCaller.h"
#include "SomaticBenchmark.h"

struct SomaticHaplotagParameters
{
    HaplotagParameters basic;
    // Somatic haplotag parameters
    std::string tumorSnvFile;   
    std::string tumorBamFile;  

    // somatic benchmark parameters
    std::string benchmarkVcf;
    std::string benchmarkBedFile;
    std::string metricsSuffix;
    bool writeBenchmarkLog;

    // somatic calling config
    CallerConfig callerCfg;
    bool writeSomaticCallingVcf;
    bool enableTwoSiteSplit;
};

struct TwoSiteSomaticPS
{
    int ps;
    int siteA;
    int siteB;
};

enum TwoSiteReadAllele
{
    READ_ALLELE_UNK = 0,
    READ_ALLELE_REF = 1,
    READ_ALLELE_ALT = 2
};

enum PairAlleleState
{
    PAIR_UNKNOWN = 0,
    ALT_REF = 1,
    REF_ALT = 2,
    ALT_ALT = 3,
    REF_REF = 4,
    PARTIAL = 5
};

// Tumor specific header
class SomaticTagLog : public GermlineTagLog {
    private:
        const SomaticHaplotagParameters& sParams;
    public:
        SomaticTagLog(const SomaticHaplotagParameters& sParams) : GermlineTagLog(sParams.basic), sParams(sParams){}
        ~SomaticTagLog(){};

        void addParamsMessage() override;

        void writeBasicColumns() override;

        void writeTagReadLog(TagReadLog& data) override;
};

class SomaticHaplotagCigarParser: public CigarParser{
    private:
        std::map<int, int>& tumCountPS;
        std::map<int, std::pair<int, int>>& somaticVarDeriveHP;
        std::map<int, int>& somaticReadAlleleState;
        SomaticReadVerifier& somaticReadCounter;
        SomaticHaplotagStrategy somaticJudger;
        int resolveSomaticReadAllele(MultiGenomeVar& currentVar, const std::string& base) const;
    protected:

        void processMatchOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, std::string& base, bool& isAlt, int& offset) override;
        void processDeletionOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, bool& alreadyJudgeDel) override;
    public:
        SomaticHaplotagCigarParser(
            CigarParserContext& ctx,
            int& ref_pos,
            int& query_pos,
            std::map<int, int>& tumCountPS,
            std::map<int, std::pair<int, int>>& somaticVarDeriveHP,
            std::map<int, int>& somaticReadAlleleState,
            SomaticReadVerifier& somaticReadCounter
        );
        ~SomaticHaplotagCigarParser() override;
};

class SomaticHaplotagChrProcessor: public GermlineHaplotagChrProcessor{
    private:

        SomaticReadVerifier* somaticReadCounter;

        chrReadHpResult* localHpBeforeInheritance;
        chrReadHpResult* localHpAfterInheritance;

        std::string chr;
        std::time_t begin;

        SomaticHaplotagStrategy somaticJudger;
        bool enableTwoSiteSplit;
        const std::map<int, TwoSiteSomaticPS>* twoSiteSomaticPS;

        PairAlleleState resolvePairAlleleState(
            const TwoSiteSomaticPS& twoSite,
            const std::map<int, int>& somaticReadAlleleState,
            int& alleleA,
            int& alleleB
        ) const;

        int resolveTwoSiteSplitLabel(
            int rawHpResult,
            const std::map<int, int>& norCountPS,
            const std::map<int, int>& somaticReadAlleleState,
            const bam1_t& aln
        );
    protected:
        int judgeHaplotype(
            const bam_hdr_t &bamHdr,
            const bam1_t &aln,
            std::string chrName,
            double percentageThreshold,
            GermlineTagLog *tagResult,
            int &pqValue,
            int &psValue,
            const int tagSample,
            const std::string &ref_string,
            const ParsingBamConfig &params,
            std::map<int, MultiGenomeVar>::iterator &firstVariantIter,
            std::map<int, MultiGenomeVar> &currentChrVariants,
            std::map<Genome, VCF_Info> &vcfSet
        ) override;

        virtual void addAuxiliaryTags(bam1_t *aln, int& haplotype, int& pqValue, int& psValue) override;

        int inheritHaplotype(
            float &deriveByHpSimilarity,
            double percentageThreshold,
            std::map<int, std::pair<int , int>>& somaticVarDeriveHP,
            std::map<int, int>& hpCount,
            int &hpResult,
            const bam1_t &aln
        );

    public:
        SomaticHaplotagChrProcessor(
            bool writeOutputBam,
            bool mappingQualityFilter,
            ReadStatistics& readStats,
            GermlineTagLog *tagResult,
            SomaticReadBenchmark& highConSomaticData,
            ReadHpDistriLog& hpBeforeInheritance,
            ReadHpDistriLog& hpAfterInheritance,
            const std::string &chr,
            bool enableTwoSiteSplit,
            const std::map<int, TwoSiteSomaticPS>* twoSiteSomaticPS
        );
        ~SomaticHaplotagChrProcessor() override;
};

class SomaticHaplotagBamParser: public GermlineHaplotagBamParser{
    private:
        SomaticReadBenchmark& somaticBenchmark;
        ReadHpDistriLog& hpBeforeInheritance;
        ReadHpDistriLog& hpAfterInheritance;
        const std::map<std::string, std::map<int, TwoSiteSomaticPS>>& chrTwoSiteSomaticPS;
        std::map<int, TwoSiteSomaticPS> emptyTwoSiteSomaticPS;
        const SomaticHaplotagParameters& sParams;
    protected:
        std::unique_ptr<ChromosomeProcessor> createProcessor(const std::string &chr) override {
            const std::map<int, TwoSiteSomaticPS>* chrTwoSitePS = &emptyTwoSiteSomaticPS;
            auto twoSitePSIter = chrTwoSiteSomaticPS.find(chr);
            if(twoSitePSIter != chrTwoSiteSomaticPS.end()){
                chrTwoSitePS = &twoSitePSIter->second;
            }

            return std::unique_ptr<ChromosomeProcessor>(new SomaticHaplotagChrProcessor(
                control.writeOutputBam,
                control.mappingQualityFilter,
                readStats,
                tagResult,
                somaticBenchmark,
                hpBeforeInheritance,
                hpAfterInheritance,
                chr,
                sParams.enableTwoSiteSplit,
                chrTwoSitePS
            ));
        };

        // create somatic tag read log
        virtual GermlineTagLog* createTagReadLog() override {
            return new SomaticTagLog(sParams);
        }

    public:
        SomaticHaplotagBamParser(
            const ParsingBamConfig &config,
            const ParsingBamControl &control,
            ReadStatistics& readStats,
            SomaticReadBenchmark& highConSomaticData,
            ReadHpDistriLog& hpBeforeInheritance,
            ReadHpDistriLog& hpAfterInheritance,
            const std::map<std::string, std::map<int, TwoSiteSomaticPS>>& chrTwoSiteSomaticPS,
            const SomaticHaplotagParameters& sParams
        );
        ~SomaticHaplotagBamParser() override;
};

class SomaticHaplotagProcess: public HaplotagProcess{
    private:
        SomaticHaplotagParameters& sParams;
        virtual void printParamsMessage() override;
    protected:

        ReadHpDistriLog *hpBeforeInheritance;
        ReadHpDistriLog *hpAfterInheritance;
        
        SomaticReadBenchmark somaticBenchmark;
        std::map<std::string, std::map<int, TwoSiteSomaticPS>> chrTwoSiteSomaticPS;
        int totalTwoSiteSomaticPS;

        void parseVariantFiles(VcfParser& vcfParser) override;
        void setChrVecAndChrLength() override;
        void buildTwoSiteSomaticPS();
        
        // calculate SNP counts
        void displaySnpCounts();

        void postprocessForHaplotag() override;

        GermlineHaplotagBamParser* createHaplotagBamParser(
            const ParsingBamConfig &config,
            const ParsingBamControl &control,
            ReadStatistics& readStats
        ) override{
            return new SomaticHaplotagBamParser(config, control, readStats, somaticBenchmark, *hpBeforeInheritance, *hpAfterInheritance, chrTwoSiteSomaticPS, sParams);
        }

        virtual std::string getNormalSnpParsingMessage() override{
            return "parsing normal SNP VCF ... ";
        };

        virtual std::string getTagReadStartMessage() override{
            return "somatic tagging start ...\n";
        };

    public:
        SomaticHaplotagProcess(SomaticHaplotagParameters &params);
        ~SomaticHaplotagProcess() override;

        virtual void pipelineProcess() override;
};

#endif
