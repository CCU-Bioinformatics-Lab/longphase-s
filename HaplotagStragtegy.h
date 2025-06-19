#ifndef HAPLOTAG_STRATEGY_H
#define HAPLOTAG_STRATEGY_H

#include "HaplotagType.h"

/**
 * @brief Strategy class for germline haplotagging decisions
 * 
 * Implements the logic for determining haplotype assignments in germline samples
 */
class GermlineHaplotagStrategy{
    private:
    protected:
    public:
        void judgeSnpHap(
            const std::string& chrName,
            VarData& norVar, const std::string& base,
            int& ref_pos,
            int& length,
            int& i,
            int& aln_core_n_cigar,
            uint32_t* cigar,
            std::map<int, MultiGenomeVar>::iterator& currentVariantIter,
            std::map<int, int>& hpCount,
            std::map<int, int>& variantsHP,
            std::map<int, int>& countPS
        );

        void judgeDeletionHap(
            const std::string& chrName,
            const std::string& ref_string,
            int& ref_pos,
            int& length,
            int& query_pos,
            std::map<int, MultiGenomeVar>::iterator &currentVariantIter,
            const bam1_t* aln, std::map<int, int>& hpCount,
            std::map<int, int>& variantsHP,
            std::map<int, int>& countPS
        );

        void judgeSVHap(
            const bam1_t &aln,
            std::map<Genome, VCF_Info> &vcfSet,
            std::map<int, int>& hpCount,
            const int& tagSample
        );

        int judgeReadHap(
            std::map<int, int>& hpCount,
            double& min,
            double& max,
            double& percentageThreshold,
            int& pqValue,
            int& psValue,
            std::map<int, int>& countPS,
            int* totalHighSimilarity,
            int* totalWithOutVaraint
        );
};

/**
 * @brief Strategy class for somatic haplotagging decisions
 * 
 * Implements the logic for determining haplotype assignments in tumor samples
 */
class SomaticJudgeHapStrategy{
    private :

    protected:
        void judgeNormalSnpHap(
            const std::string& chrName, 
            int& curPos,
            MultiGenomeVar& curVar,
            std::string& base,
            std::map<int, int>& hpCount, 
            std::map<int, int>& norCountPS,
            std::map<int, int> *variantsHP
        );

        virtual void judgeTumorOnlySnpHap(
            const std::string &chrName,
            int &curPos, MultiGenomeVar &curVar,
            std::string base,
            std::map<int, int> &hpCount,
            std::map<int, int> *tumCountPS,
            std::map<int, int> *variantsHP,
            std::vector<int> *tumorAllelePosVec
        ) = 0;

    public:
        void judgeSomaticSnpHap(
            std::map<int, MultiGenomeVar>::iterator &currentVariantIter,
            std::string chrName,
            std::string base,
            std::map<int, int> &hpCount,
            std::map<int, int> &norCountPS,
            std::map<int, int> &tumCountPS,
            std::map<int, int> *variantsHP,
            std::vector<int> *tumorAllelePosVec
        );


        int judgeSomaticReadHap(
            std::map<int, int> &hpCount,
            int &pqValue,
            std::map<int, int> &norCountPS,
            double &norHPsimilarity,
            double &tumHPsimilarity,
            double percentageThreshold,
            int *totalHighSimilarity,
            int *totalCrossTwoBlock,
            int *totalWithOutVaraint
        );

};

/**
 * @brief Strategy class for extracting somatic data
 */
class ExtractSomaticDataStragtegy: public SomaticJudgeHapStrategy{
    private:
    protected:
        virtual void judgeTumorOnlySnpHap(
            const std::string &chrName,
            int &curPos, MultiGenomeVar &curVar,
            std::string base,
            std::map<int, int> &hpCount,
            std::map<int, int> *tumCountPS,
            std::map<int, int> *variantsHP,
            std::vector<int> *tumorAllelePosVec
        ) override;
    public:
};

/**
 * @brief Strategy class for somatic haplotagging
 * 
 * Implements the complete somatic haplotagging logic
 */
class SomaticHaplotagStrategy: public SomaticJudgeHapStrategy{
    private:
    protected:
        virtual void judgeTumorOnlySnpHap(
            const std::string &chrName,
            int &curPos, MultiGenomeVar &curVar,
            std::string base,
            std::map<int, int> &hpCount,
            std::map<int, int> *tumCountPS,
            std::map<int, int> *variantsHP,
            std::vector<int> *tumorAllelePosVec
        ) override;
    public:
};

/**
 * @brief Namespace containing base analysis utility functions
 * 
 * Provides common calculation functions for variant analysis
 */
namespace base_analysis{

    inline float calculateVAF(int altCount, int depth){
        return (depth == 0 || altCount == 0) ? 0.0 : (float)altCount / (float)depth;
    }

    inline float calculateLowMpqReadRatio(int depth, int filteredMpqDepth){
        return depth == 0 ? 0.0 : (float)(depth - filteredMpqDepth) / (float)depth;
    }

    inline float calculateDelRatio(int delCount, int depth){
        return (depth == 0 || delCount == 0) ? 0.0 : (float)delCount / (float)depth;
    }

    inline double calculateHaplotypeImbalanceRatio(int& H1readCount, int& H2readCount, int& totalReadCount) {
        if(H1readCount > 0 && H2readCount > 0) {
            return (H1readCount > H2readCount) ? 
                ((double)H1readCount / (double)totalReadCount) : 
                ((double)H2readCount / (double)totalReadCount);
        } else if(H1readCount == 0 && H2readCount == 0) {
            return 0.0;
        } else {
            return 1.0;
        }
    }

    inline double calculatePercentageOfGermlineHp(int& totalGermlineReadCount, int& depth) {
        return (depth == 0 || totalGermlineReadCount == 0) ? 0.0 : (double)totalGermlineReadCount / (double)depth;
    }   
}

#endif