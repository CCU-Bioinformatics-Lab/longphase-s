#ifndef HAPLOTAGPROCESS_H
#define HAPLOTAGPROCESS_H

#include "Util.h"
#include "ParsingBam.h"
#include "HaplotagBase.h"
#include "SomaticVarCaller.h"
#include "SomaticBenchmark.h"

class HaplotagProcess: public SomaticJudgeBase, public GermlineJudgeBase
{
    private:
        std::vector<std::string> *chrVec;
        std::map<std::string, int> *chrLength;

        // chr, variant position (0-base), allele haplotype set
        std::map<std::string, std::map<int, MultiGenomeVar>> *mergedChrVarinat;

        // variant position (0-base), allele haplotype set
        std::map<int, MultiGenomeVar> currentChrVariants;
        std::map<int, MultiGenomeVar>::iterator firstVariantIter;

        // chr, variant position (0-base), somatic SNP information
        std::map<std::string, std::map<int, HP3_Info>> *chrPosReadCase;  

        // record the VCF files of the normal and tumor datasets (normal:0, tumor:1, seqcHighCon:2)
        VCF_Info vcfSet[3];

        void tagRead(HaplotagParameters &params, const Genome& geneType);

        void initFlag(bam1_t *aln, std::string flag);
        
        int judgeHaplotype(const bam_hdr_t &bamHdr,const bam1_t &aln, std::string chrName, double percentageThreshold, std::ofstream *tagResult, int &pqValue, int &psValue, const int tagGeneType, std::string &ref_string);
        int somaticJudgeHaplotype(const bam_hdr_t &bamHdr,const bam1_t &aln,const std::string &chrName, double percentageThreshold, std::ofstream *tagResult, int &pqValue, int &psValue, const int tagGeneType, std::string &ref_string);
        std::string convertHpResultToString(int hpResult);
        
        int totalAlignment;
        int totalSupplementary;
        int totalSecondary;
        int totalUnmapped;
        int totalTagCount;
        int totalUnTagCount;

        //--------------------verification parameter---------------------
        bool tagTumorMode;

        // reads HP count
        std::map<int, int> totalHpCount;
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
        //---------------------------------------------------------------
        
        std::time_t processBegin;
    protected:
        void OnlyTumorSNPjudgeHP(const std::string &chrName, int &curPos, MultiGenomeVar &curVar, std::string base, VCF_Info *vcfSet, std::map<int, int> &hpCount, std::map<int, int> *tumCountPS, std::map<int, int> *variantsHP, std::vector<int> *readPosHP3, BamBaseCounter *NorBase, std::map<int, HP3_Info> *SomaticPos);
    public:
        HaplotagProcess();
        void taggingProcess(HaplotagParameters &params);
        virtual ~HaplotagProcess();

};


#endif