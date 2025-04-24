#ifndef HAPLOTAGPROCESS_H
#define HAPLOTAGPROCESS_H

#include "Util.h"
#include "ParsingBam.h"
#include "HaplotagBase.h"
#include "SomaticVarCaller.h"
#include "SomaticBenchmark.h"

struct ReadStatistics {
    // reads HP count
    std::map<int, int> totalHpCount;

    int totalAlignment;
    int totalSupplementary;
    int totalSecondary;
    int totalUnmapped;
    int totalTagCount;
    int totalUnTagCount;
    // reads untag count
    int totalLowerQuality;
    int totalOtherCase;
    int totalunTag_HP0;
    int totalreadOnlyH3Snp;
    int totalHighSimilarity;
    int totalCrossTwoBlock;
    int totalEmptyVariant;
    int totalWithOutVaraint;
    ReadStatistics ():totalAlignment(0), totalSupplementary(0), totalSecondary(0),
                        totalUnmapped(0), totalTagCount(0), totalUnTagCount(0),
                        totalLowerQuality(0), totalOtherCase(0), totalunTag_HP0(0),
                        totalreadOnlyH3Snp(0), totalHighSimilarity(0), totalCrossTwoBlock(0),
                        totalEmptyVariant(0), totalWithOutVaraint(0){}
};

class GermlineHaplotagCigarParser: public CigarParser{
    protected:
        void processMatchOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, std::string& base) override;
        void processDeletionOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, bool& alreadyJudgeDel) override;
    public:
        GermlineHaplotagCigarParser(int& ref_pos, int& query_pos);
        ~GermlineHaplotagCigarParser() override;
};

class GermlineHaplotagChrProcessor: public ChromosomeProcessor{
    private:

    protected:
        ReadStatistics& readStats;
        std::ofstream *tagResult;
        virtual void processLowMappingQuality() override;
        virtual void processUnmappedRead() override;
        virtual void processSecondaryAlignment() override;
        virtual void processSupplementaryAlignment() override;
        virtual void processEmptyVariants() override;
        virtual void processOtherCase() override;

        virtual void processRead(
            bam1_t &aln, 
            const bam_hdr_t &bamHdr,
            const std::string &chrName, 
            const HaplotagParameters &params, 
            const Genome& genmoeType, 
            std::map<int, MultiGenomeVar> &currentVariants,
            std::map<int, MultiGenomeVar>::iterator &firstVariantIter, 
            std::map<Genome, VCF_Info> &vcfSet, 
            const std::string &ref_string
        ) override;

        virtual int judgeHaplotype(
            const bam_hdr_t &bamHdr,
            const bam1_t &aln,
            std::string chrName,
            double percentageThreshold,
            std::ofstream *tagResult,
            int &pqValue,
            int &psValue,
            const int tagGeneType,
            const std::string &ref_string,
            const HaplotagParameters &params,
            std::map<int, MultiGenomeVar>::iterator &firstVariantIter,
            std::map<int, MultiGenomeVar> &currentChrVariants,
            std::map<Genome, VCF_Info> &vcfSet
        );

        void initFlag(bam1_t *aln, std::string flag);

        virtual void addAuxiliaryTags(bam1_t *aln, int& haplotype, int& pqValue, int& psValue);

    public:
        GermlineHaplotagChrProcessor(
            bool writeOutputBam,
            bool mappingQualityFilter,
            ReadStatistics& readStats,
            std::ofstream *tagResult
        );
        ~GermlineHaplotagChrProcessor() override;
};

class GermlineHaplotagBamParser: public HaplotagBamParser{
    private:
    protected:
        ReadStatistics& readStats;
        std::ofstream *tagResult;
        std::unique_ptr<ChromosomeProcessor> createProcessor(const std::string &chr) override{
            return std::unique_ptr<ChromosomeProcessor>(new GermlineHaplotagChrProcessor(writeOutputBam, mappingQualityFilter, readStats, tagResult));
        };
    public:
        GermlineHaplotagBamParser(
            bool& writeOutputBam,
            bool& mappingQualityFilter,
            ReadStatistics& readStats,
            std::ofstream *tagResult
        );
        ~GermlineHaplotagBamParser() override;
};

class SomaticHaplotagCigarParser: public CigarParser, public SomaticJudgeBase{
    private:
        std::map<int, int>& tumCountPS;
        std::map<int, std::pair<int, int>>& somaticVarDeriveHP;
        SomaticReadVerifier& highConSomaticData;
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
          , SomaticReadVerifier& highConSomaticData
        );
        ~SomaticHaplotagCigarParser() override;
};

class SomaticHaplotagChrProcessor: public GermlineHaplotagChrProcessor, public SomaticJudgeBase{
    private:
        SomaticReadVerifier& highConSomaticData;
        ReadHpDistriLog& hpBeforeInheritance;
        ReadHpDistriLog& hpAfterInheritance;
    protected:
        int judgeHaplotype(
            const bam_hdr_t &bamHdr,
            const bam1_t &aln,
            std::string chrName,
            double percentageThreshold,
            std::ofstream *tagResult,
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

        void postProcess(
            const std::string &chr,
            std::map<int, MultiGenomeVar> &currentVariants
        ) override;

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
            std::ofstream *tagResult,
            SomaticReadVerifier& highConSomaticData,
            ReadHpDistriLog& hpBeforeInheritance,
            ReadHpDistriLog& hpAfterInheritance
        );
        ~SomaticHaplotagChrProcessor() override;
};

class SomaticHaplotagBamParser: public GermlineHaplotagBamParser{
    private:
        SomaticReadVerifier& highConSomaticData;
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
                hpAfterInheritance
            ));
        };
    public:
        SomaticHaplotagBamParser(
            bool& writeOutputBam,
            bool& mappingQualityFilter,
            ReadStatistics& readStats,
            std::ofstream *tagResult,
            SomaticReadVerifier& highConSomaticData,
            ReadHpDistriLog& hpBeforeInheritance,
            ReadHpDistriLog& hpAfterInheritance
        );
        ~SomaticHaplotagBamParser() override;
};

class HaplotagProcess: public SomaticJudgeBase
{
    private:
        HaplotagParameters &params;

        std::vector<std::string> *chrVec;
        std::map<std::string, int> *chrLength;

        // chr, variant position (0-base), allele haplotype set
        std::map<std::string, std::map<int, MultiGenomeVar>> *mergedChrVarinat;

        // variant position (0-base), allele haplotype set
        std::map<int, MultiGenomeVar> currentChrVariants;
        std::map<int, MultiGenomeVar>::iterator firstVariantIter;

        // record the VCF files
        std::map<Genome, VCF_Info> vcfSet;

        //--------------------verification parameter---------------------
        bool tagTumorMode;

        // write tag read detail information
        std::ofstream *tagResult;

        ReadStatistics readStats;

        ReadHpDistriLog *hpBeforeInheritance;
        ReadHpDistriLog *hpAfterInheritance;
        
        SomaticReadVerifier highConSomaticData;
        
        std::time_t processBegin;

        void tagRead(HaplotagParameters &params, const Genome& geneType);
        void printAlignmentStaristics();
    public:
        HaplotagProcess(HaplotagParameters &params);
        void taggingProcess();
        virtual ~HaplotagProcess();

};


#endif