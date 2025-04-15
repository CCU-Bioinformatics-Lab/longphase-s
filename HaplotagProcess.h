#ifndef HAPLOTAGPROCESS_H
#define HAPLOTAGPROCESS_H

#include "Util.h"
#include "ParsingBam.h"
#include "HaplotagBase.h"
#include "SomaticVarCaller.h"
#include "SomaticBenchmark.h"

class GermlineHaplotagCigarParser: public CigarParser{
    protected:
        void processMatchOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, std::string& base) override;
        void processDeletionOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, bool& alreadyJudgeDel) override;
    public:
        GermlineHaplotagCigarParser(int& ref_pos, int& query_pos);
        ~GermlineHaplotagCigarParser() override;
};


class SomaticJudgeHpCigarParser: public CigarParser, public SomaticJudgeBase{
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
        SomaticJudgeHpCigarParser(
            int& ref_pos
          , int& query_pos
          , std::map<int, int>& tumCountPS
          , std::map<int, std::pair<int, int>>& somaticVarDeriveHP
          , SomaticReadVerifier& highConSomaticData
        );
        ~SomaticJudgeHpCigarParser() override;
};
class HaplotagProcess: public SomaticJudgeBase, public GermlineJudgeBase
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

        // record the VCF files of the normal and tumor datasets (normal:0, tumor:1, seqcHighCon:2)
        VCF_Info vcfSet[3];

        void tagRead(HaplotagParameters &params, const Genome& geneType);

        void initFlag(bam1_t *aln, std::string flag);
        
        int judgeHaplotype(const bam_hdr_t &bamHdr,const bam1_t &aln, std::string chrName, double percentageThreshold, std::ofstream *tagResult, int &pqValue, int &psValue, const int tagGeneType, std::string &ref_string, const HaplotagParameters &params);
        int somaticJudgeHaplotype(
            const bam_hdr_t &bamHdr,
            const bam1_t &aln,
            const std::string &chrName,
            double percentageThreshold,
            std::ofstream *tagResult,
            int &pqValue,
            int &psValue,
            const int tagGeneType,
            std::string &ref_string,
            const HaplotagParameters &params
        );
        std::string convertHpResultToString(int hpResult);

        void printTaggingResult();
        

        //--------------------verification parameter---------------------
        bool tagTumorMode;

        // write tag read detail information
        std::ofstream *tagResult;

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

        // chr, variant position (0-base), reads HP 
        std::map<std::string, std::map<int, ReadHpResult>> *beforeCorrReadHpResult;
        std::map<std::string, std::map<int, ReadHpResult>> *afterCorrReadHpResult;

        ReadHpDistriLog *hpBeforeInheritance;
        ReadHpDistriLog *hpAfterInheritance;
        
        SomaticReadVerifier highConSomaticData;
        
        std::time_t processBegin;
    protected:
        void OnlyTumorSNPjudgeHP(const std::string &chrName, int &curPos, MultiGenomeVar &curVar, std::string base, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *readPosHP3) override;
    public:
        HaplotagProcess(HaplotagParameters &params);
        void taggingProcess();
        virtual ~HaplotagProcess();

};


#endif