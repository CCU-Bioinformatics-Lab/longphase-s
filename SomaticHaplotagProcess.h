#ifndef SOMATICHAPLOTAGPROCESS_H
#define SOMATICHAPLOTAGPROCESS_H

#include "HaplotagType.h"
#include "HaplotagProcess.h"
#include "HaplotagParsingBam.h"
#include "HaplotagVcfParser.h"
#include "SomaticVarCaller.h"
#include "SomaticBenchmark.h"
#include "HaplotagLogging.h"

class SomaticHaplotagParamsMessage : public HaplotagParamsMessage{
    public:
        SomaticHaplotagParamsMessage(const HaplotagParameters& params):HaplotagParamsMessage(params){}

        virtual void addParamsMessage() override {
            addCommonParamsMessage();
            insertAfterKey("tumor SNP file", params.tumorSnpFile, "phased SNP file");
            insertAfterKey("input tumor bam file", params.tumorBamFile, "input bam file");
            insertAfterKey("somatic calling mapping quality", params.qualityThreshold, "percentage threshold");
            insertAfterKey("enable somatic variant filter", params.enableFilter, "somatic calling mapping quality");

            // sort the entries by order
            sortEntries();
        } 
};

// Tumor specific header
class SomaticTagLog : public GermlineTagLog {
    public:
        SomaticTagLog(const HaplotagParameters& params) : GermlineTagLog(params){}
        ~SomaticTagLog(){};

        void addBasicEntries() override {
            // add basic entries
            addCommonBasicEntries();
            // add tumor specific entries
            insertAfterKey("tumorSnpFile", params.tumorSnpFile, "snpFile");
            insertAfterKey("tumorBamFile", params.tumorBamFile, "bamFile");
            insertAfterKey("somaticCallingThreshold", params.qualityThreshold, "region");
        }

        void writeBasicColumns() override {
            *tagReadLog << "#ReadID\t"
                        << "CHROM\t"
                        << "ReadStart\t"
                        << "Confidnet(%)\t"
                        << "deriveByHpSimilarity\t"
                        << "Haplotype\t"
                        << "PhaseSet\t"
                        << "TotalAllele\t"
                        << "HP1Allele\t"
                        << "HP2Allele\t"
                        << "HP3Allele\t"
                        << "HP4Allele\t"
                        << "phasingQuality(PQ)\t"
                        << "(Variant,HP)\t"
                        << "(PhaseSet,Variantcount)\n";
        }

        void writeTagReadLog(TagReadLog& data) override;
};

class SomaticHaplotagCigarParser: public CigarParser, public SomaticJudgeBase{
    private:
        std::map<int, int>& tumCountPS;
        std::map<int, std::pair<int, int>>& somaticVarDeriveHP;
        SomaticReadVerifier& somaticReadCounter;
    protected:

        void OnlyTumorSNPjudgeHP(
            const std::string &chrName,
            int &curPos,
            MultiGenomeVar &curVar,
            std::string base,
            std::map<int, int> &hpCount,
            std::map<int, int> *tumCountPS,
            std::map<int, int> *variantsHP,
            std::vector<int> *tumorAllelePosVec
        ) override;

        void processMatchOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, std::string& base) override;
        void processDeletionOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, bool& alreadyJudgeDel) override;
    public:
        SomaticHaplotagCigarParser(
            int& ref_pos
          , int& query_pos
          , std::map<int, int>& tumCountPS
          , std::map<int, std::pair<int, int>>& somaticVarDeriveHP
          , SomaticReadVerifier& somaticReadCounter
        );
        ~SomaticHaplotagCigarParser() override;
};

class SomaticHaplotagChrProcessor: public GermlineHaplotagChrProcessor, public SomaticJudgeBase{
    private:

        SomaticReadVerifier* somaticReadCounter;

        chrReadHpResult* localHpBeforeInheritance;
        chrReadHpResult* localHpAfterInheritance;

        std::string chr;
        std::time_t begin;
    protected:
        int judgeHaplotype(
            const bam_hdr_t &bamHdr,
            const bam1_t &aln,
            std::string chrName,
            double percentageThreshold,
            GermlineTagLog *tagResult,
            int &pqValue,
            int &psValue,
            const int tagGeneType,
            const std::string &ref_string,
            const HaplotagParameters &params,
            std::map<int, MultiGenomeVar>::iterator &firstVariantIter,
            std::map<int, MultiGenomeVar> &currentChrVariants,
            std::map<Genome, VCF_Info> &vcfSet
        ) override;

        virtual void addAuxiliaryTags(bam1_t *aln, int& haplotype, int& pqValue, int& psValue) override;
        std::string convertHpResultToString(int hpResult);

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
        SomaticReadBenchmark& highConSomaticData;
        ReadHpDistriLog& hpBeforeInheritance;
        ReadHpDistriLog& hpAfterInheritance;
    protected:
        std::unique_ptr<ChromosomeProcessor> createProcessor(const std::string &chr) override {
            return std::unique_ptr<ChromosomeProcessor>(new SomaticHaplotagChrProcessor(
                writeOutputBam,
                mappingQualityFilter,
                readStats,
                tagResult,
                highConSomaticData,
                hpBeforeInheritance,
                hpAfterInheritance,
                chr
            ));
        };

        // create somatic tag read log
        virtual GermlineTagLog* createTagReadLog(const HaplotagParameters& params) override {
            return new SomaticTagLog(params);
        }

    public:
        SomaticHaplotagBamParser(
            ParsingBamMode mode,
            bool writeOutputBam,
            bool mappingQualityFilter,
            ReadStatistics& readStats,
            SomaticReadBenchmark& highConSomaticData,
            ReadHpDistriLog& hpBeforeInheritance,
            ReadHpDistriLog& hpAfterInheritance
        );
        ~SomaticHaplotagBamParser() override;
};

class SomaticHaplotagProcess: public HaplotagProcess, public SomaticJudgeBase{
    protected:

        SomaticHaplotagParamsMessage somaticParamsMessage;

        ReadHpDistriLog *hpBeforeInheritance;
        ReadHpDistriLog *hpAfterInheritance;
        
        SomaticReadBenchmark highConSomaticData;

        void parseVariantFiles(VcfParser& vcfParser) override;
        void setChrVecAndChrLength() override;
        
        // calculate SNP counts
        void calculateSnpCounts();

        void postprocessForHaplotag() override;

        GermlineHaplotagBamParser* createHaplotagBamParser(
            ParsingBamMode mode,
            bool writeOutputBam,
            bool mappingQualityFilter,
            ReadStatistics& readStats
        ) override{
            return new SomaticHaplotagBamParser(mode, writeOutputBam, mappingQualityFilter, readStats, highConSomaticData, *hpBeforeInheritance, *hpAfterInheritance);
        }

        virtual std::string getNormalSnpParsingMessage() override{
            return "parsing normal SNP VCF ... ";
        };

        virtual std::string getTagReadStartMessage() override{
            return "somatic tagging start ...\n";
        };

    public:
        SomaticHaplotagProcess(HaplotagParameters &params);
        ~SomaticHaplotagProcess() override;

        virtual void taggingProcess() override;
};

#endif

