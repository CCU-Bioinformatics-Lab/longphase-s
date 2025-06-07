#ifndef SOMATICHAPLOTAGPROCESS_H
#define SOMATICHAPLOTAGPROCESS_H

#include "HaplotagType.h"
#include "HaplotagProcess.h"
#include "HaplotagParsingBam.h"
#include "HaplotagVcfParser.h"
#include "SomaticVarCaller.h"
#include "SomaticBenchmark.h"
#include "HaplotagLogging.h"

struct SomaticHaplotagParameters
{
    HaplotagParameters basic;
    // Somatic haplotag parameters
    std::string tumorSnpFile;   
    std::string tumorBamFile;  
    std::string benchmarkVcf;
    std::string benchmarkBedFile;

    std::string metricsSuffix;
    
    double tumorPurity;
    bool predictTumorPurity;

    bool enableFilter;
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
        SomaticReadVerifier& somaticReadCounter;
        SomaticHaplotagStrategy somaticJudger;
    protected:

        void processMatchOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, std::string& base) override;
        void processDeletionOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, bool& alreadyJudgeDel) override;
    public:
        SomaticHaplotagCigarParser(
            CigarParserContext& ctx,
            int& ref_pos,
            int& query_pos,
            std::map<int, int>& tumCountPS,
            std::map<int, std::pair<int, int>>& somaticVarDeriveHP,
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
            int &hpResult
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
            const std::string &chr
        );
        ~SomaticHaplotagChrProcessor() override;
};

class SomaticHaplotagBamParser: public GermlineHaplotagBamParser{
    private:
        SomaticReadBenchmark& somaticBenchmark;
        ReadHpDistriLog& hpBeforeInheritance;
        ReadHpDistriLog& hpAfterInheritance;
        const SomaticHaplotagParameters& sParams;
    protected:
        std::unique_ptr<ChromosomeProcessor> createProcessor(const std::string &chr) override {
            return std::unique_ptr<ChromosomeProcessor>(new SomaticHaplotagChrProcessor(
                control.writeOutputBam,
                control.mappingQualityFilter,
                readStats,
                tagResult,
                somaticBenchmark,
                hpBeforeInheritance,
                hpAfterInheritance,
                chr
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

        void parseVariantFiles(VcfParser& vcfParser) override;
        void setChrVecAndChrLength() override;
        
        // calculate SNP counts
        void calculateSnpCounts();

        void postprocessForHaplotag() override;

        GermlineHaplotagBamParser* createHaplotagBamParser(
            const ParsingBamConfig &config,
            const ParsingBamControl &control,
            ReadStatistics& readStats
        ) override{
            return new SomaticHaplotagBamParser(config, control, readStats, somaticBenchmark, *hpBeforeInheritance, *hpAfterInheritance, sParams);
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

        virtual void taggingProcess() override;
};

#endif

