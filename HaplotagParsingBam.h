#ifndef HAPLOTAG_PARSING_BAM_H
#define HAPLOTAG_PARSING_BAM_H

#include "HaplotagType.h"

class HaplotagBamParser;
class ChromosomeProcessor;
class CigarParser;


class GermlineJudgeBase{
    private:

    protected:
        void germlineJudgeSnpHap(
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

        void germlineJudgeDeletionHap(
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
        void germlineJudgeSVHap(
            const bam1_t &aln,
            std::map<Genome, VCF_Info> &vcfSet,
            std::map<int, int>& hpCount,
            const int& tagGeneType
        );

        int germlineDetermineReadHap(
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
    public:
};

class SomaticJudgeBase{
    private :

    protected:
        void SomaticJudgeSnpHP(
            std::map<int, MultiGenomeVar>::iterator &currentVariantIter,
            std::string chrName,
            std::string base,
            std::map<int, int> &hpCount,
            std::map<int, int> &norCountPS,
            std::map<int, int> &tumCountPS,
            std::map<int, int> *variantsHP,
            std::vector<int> *tumorAllelePosVec
        );

        virtual void OnlyTumorSNPjudgeHP(
            const std::string &chrName,
            int &curPos, MultiGenomeVar &curVar,
            std::string base,
            std::map<int, int> &hpCount,
            std::map<int, int> *tumCountPS,
            std::map<int, int> *variantsHP,
            std::vector<int> *tumorAllelePosVec
        );

        int determineReadHP(
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

    public:
};

enum ParsingBamMode{
    SINGLE_THREAD = 0,
    MULTI_THREAD = 1
};

class BamFileRAII {
    private:
        bool writeOutputBam;
        bool isReleased;

        template<typename T>
        void checkNullPointer(const T* ptr, const std::string& errorMessage) const;
        
    public:
        samFile* in;
        samFile* out;
        bam_hdr_t* bamHdr;
        hts_idx_t* idx;
        bam1_t* aln;

        BamFileRAII(const std::string& BamFile
                  , const std::string& fastaFile
                  , htsThreadPool &threadPool
                  , const HaplotagParameters& params
                  , const bool writeOutputBam = false);
        ~BamFileRAII();

        bool validateState();
        void samWriteBam();

        void destroy();
};

class HaplotagBamParser{
    private:
        ParsingBamMode mode;

        void processBamParallel(
            const std::string &BamFile, 
            const HaplotagParameters &params, 
            const std::vector<std::string> &chrVec, 
            const std::map<std::string, int> &chrLength, 
            const FastaParser &fastaParser,
            htsThreadPool &threadPool,
            std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat, 
            std::map<Genome, VCF_Info> &vcfSet, 
            const Genome& genmoeType
        );

        void processBamWithOutput(
            const std::string &BamFile, 
            const HaplotagParameters &params, 
            const std::vector<std::string> &chrVec, 
            const std::map<std::string, int> &chrLength, 
            const FastaParser &fastaParser,
            htsThreadPool &threadPool,
            std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat, 
            std::map<Genome, VCF_Info> &vcfSet, 
            const Genome& genmoeType
        );
    protected: 

        bool writeOutputBam;
        bool mappingQualityFilter;
        // Factory method to create a chromosome processor
        virtual std::unique_ptr<ChromosomeProcessor> createProcessor(const std::string &chr) = 0;

        void getLastVarPos(
            std::vector<int>& last_pos,
            const std::vector<std::string>& chrVec,
            std::map<std::string,std::map<int, MultiGenomeVar>> &mergedChrVarinat,
            const Genome& geneType
        );

    public:
        HaplotagBamParser(
            ParsingBamMode mode = ParsingBamMode::MULTI_THREAD,
            bool writeOutputBam = false, 
            bool mappingQualityFilter = false
        );
        virtual ~HaplotagBamParser();
        void parsingBam(
            const std::string &BamFile, 
            const HaplotagParameters &params, 
            const std::vector<std::string> &chrVec, 
            const std::map<std::string, int> &chrLength, 
            std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat, 
            std::map<Genome, VCF_Info> &vcfSet,
            const Genome& genmoeType
        );

};

class ChromosomeProcessor : public GermlineJudgeBase{
    private:
        bool writeOutputBam;
        bool mappingQualityFilter;
    protected:

        virtual void processLowMappingQuality(){};
        virtual void processUnmappedRead(){};
        virtual void processSecondaryAlignment(){};
        virtual void processSupplementaryAlignment(){};
        virtual void processEmptyVariants(){};
        virtual void processOtherCase(){};

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
        ) = 0;

        virtual void postProcess(
            const std::string &chr,
            std::map<int, MultiGenomeVar> &currentVariants
        ){};
        
        void calculateBaseCommonInfo(PosBase& baseInfo, std::string& tumorAltBase);

        float calculateVAF(int altCount, int depth);
        float calculateLowMpqReadRatio(int depth, int filteredMpqDepth);
        float calculateDelRatio(int delCount, int depth);

        double calculateHaplotypeImbalanceRatio(int& H1readCount, int& H2readCount, int& totalReadCount);
        double calculatePercentageOfGermlineHp(int& totalReadCount, int& depth);

    public:
        ChromosomeProcessor( bool writeOutputBam=false, bool mappingQualityFilter=false);
        virtual ~ChromosomeProcessor();
        
        void processSingleChromosome(
            const std::string& chr,
            const std::map<std::string, int>& chrLength,
            const HaplotagParameters& params, 
            const FastaParser& fastaParser,
            std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat, 
            BamFileRAII& bam,
            const Genome& genmoeType,
            std::map<Genome, VCF_Info> &vcfSet
        );

};

class CigarParser : public GermlineJudgeBase{
    private:
    protected:
        // Common data members that derived classes might need
        const bam1_t* aln;
        const bam_hdr_t* bamHdr;
        const std::string* chrName;
        const HaplotagParameters* params;
        const std::string* ref_string;
        std::map<int, int>* hpCount;
        std::map<int, int>* norCountPS;
        std::map<int, int>* variantsHP;

        // position relative to reference
        int& ref_pos;
        // position relative to read
        int& query_pos;

        std::map<int, MultiGenomeVar>::iterator currentVariantIter;

        // Virtual methods that derived classes must implement
        virtual void processMatchOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, std::string& base){};
        virtual void processInsertionOperation(int& length){};
        virtual void processDeletionOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, bool& alreadyJudgeDel){};
        virtual void processSkippedOperation(int& length){};
        virtual void processSoftClippingOperation(int& length){};
        virtual void processHardClippingOperation(){};
        
        //count base nucleotide
        void countBaseNucleotide(PosBase& posBase, std::string& base, const bam1_t& aln, const float& mpqThreshold);
        void countDeletionBase(PosBase& posBase);

    public:

        CigarParser(int& ref_pos, int& query_pos);

        virtual ~CigarParser();
        void parsingCigar(
            const bam1_t& aln,
            const bam_hdr_t& bamHdr,
            const std::string& chrName,
            const HaplotagParameters& params,
            std::map<int, MultiGenomeVar>::iterator& firstVariantIter,
            std::map<int, MultiGenomeVar>& currentVariants,
            const std::string& ref_string,
            std::map<int, int>& hpCount,
            std::map<int, int>& variantsHP,
            std::map<int, int>& norCountPS
        );
};

#endif
