#ifndef HAPLOTAGPROCESS_H
#define HAPLOTAGPROCESS_H

#include "HaplotagType.h"
#include "HaplotagParsingBam.h"
#include "HaplotagStragtegy.h"
#include "HaplotagVcfParser.h"
#include "HaplotagLogging.h"

struct HaplotagParameters
{
    std::string snpFile;
    std::string svFile;
    std::string modFile;
    std::string bamFile;
    std::string fastaFile;
    
    ParsingBamConfig bamCfg;
};

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

struct TagReadLog{
    const bam1_t& aln;
    const bam_hdr_t& bamHdr;
    double& norHPsimilarity;
    const std::string& hpResultStr;
    const std::string& psResultStr;
    std::map<int, int>& hpCount;
    int& pqValue;
    const std::map<int, int>& variantsHP;
    const std::map<int, int>& norCountPS;
    
    //Somatic Read Log
    const std::map<int, int>* tumCountPS;
    double deriveByHpSimilarity;
};


class GermlineTagLog : public HaplotagReadLog<HaplotagParameters, TagReadLog>{
    protected:
        virtual void addParamsMessage() override;

        virtual void writeBasicColumns() override;

    public:
        GermlineTagLog(const HaplotagParameters& params);
        ~GermlineTagLog();

        virtual void writeHeader() {
            addParamsMessage();

            writeBasicColumns();
        }

        virtual void writeTagReadLog(TagReadLog& data);
};

class GermlineHaplotagCigarParser: public CigarParser{
    private:
    protected:
        GermlineHaplotagStrategy judger;
        void processMatchOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, std::string& base) override;
        void processDeletionOperation(int& length, uint32_t* cigar, int& i, int& aln_core_n_cigar, bool& alreadyJudgeDel) override;
    public:
        GermlineHaplotagCigarParser(CigarParserContext& ctx, int& ref_pos, int& query_pos);
        ~GermlineHaplotagCigarParser() override;
};

class GermlineHaplotagChrProcessor: public ChromosomeProcessor{
    private:
    protected:
        ReadStatistics& readStats;
        // std::ofstream *tagResult;
        GermlineTagLog* tagResult;

        ReadStatistics localReadStats;

        GermlineHaplotagStrategy judger;
        
        virtual void processLowMappingQuality() override;
        virtual void processUnmappedRead() override;
        virtual void processSecondaryAlignment() override;
        virtual void processSupplementaryAlignment() override;
        virtual void processEmptyVariants() override;
        virtual void processOtherCase() override;

        virtual void processRead(
            bam1_t &aln, 
            const bam_hdr_t &bamHdr,
            const std::string &ref_string,
            std::map<int, MultiGenomeVar> &currentVariants,
            std::map<int, MultiGenomeVar>::iterator &firstVariantIter,
            ChrProcContext& ctx
        ) override;

        virtual int judgeHaplotype(
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
        );

        void initFlag(bam1_t *aln, std::string flag);

        virtual void addAuxiliaryTags(bam1_t *aln, int& haplotype, int& pqValue, int& psValue);

       void postProcess(
            const std::string &chr,
            std::map<int, MultiGenomeVar> &currentVariants
        ) override;
        
    public:
        GermlineHaplotagChrProcessor(
            bool writeOutputBam,
            bool mappingQualityFilter,
            ReadStatistics& readStats,
            GermlineTagLog *tagResult
        );
        ~GermlineHaplotagChrProcessor() override;
};

class GermlineHaplotagBamParser: public HaplotagBamParser{
    private:
        const HaplotagParameters& params;
    protected:
        ReadStatistics& readStats;
        GermlineTagLog *tagResult;

        std::unique_ptr<ChromosomeProcessor> createProcessor(const std::string &chr) override{
            return std::unique_ptr<ChromosomeProcessor>(
                new GermlineHaplotagChrProcessor(
                    control.writeOutputBam,
                    control.mappingQualityFilter,
                    readStats,
                    tagResult
                )
            );
        };

        // create tag read log
        virtual GermlineTagLog* createTagReadLog(){
            return new GermlineTagLog(params);
        };

    public:
        GermlineHaplotagBamParser(
            const ParsingBamConfig &config,
            const ParsingBamControl &control,
            ReadStatistics& readStats,
            const HaplotagParameters& params
        );
        ~GermlineHaplotagBamParser() override;

        void createTagLog();
};

class HaplotagProcess
{
    private:
        HaplotagParameters &params;

        virtual void printParamsMessage();
    protected:

        std::vector<std::string> *chrVec;
        std::map<std::string, int> *chrLength;

        // chr, variant position (0-base), allele haplotype set
        std::map<std::string, std::map<int, MultiGenomeVar>> *mergedChrVarinat;

        // variant position (0-base), allele haplotype set
        std::map<int, MultiGenomeVar> currentChrVariants;
        std::map<int, MultiGenomeVar>::iterator firstVariantIter;

        // record the VCF files
        std::map<Genome, VCF_Info> vcfSet;

        ReadStatistics readStats;
        
        std::time_t processBegin;

        // load SNP, SV, MOD vcf file
        virtual void parseVariantFiles(VcfParser& vcfParser);
        // decide which genome sample chrVec and chrLength belong to
        virtual void setChrVecAndChrLength();
        // update chromosome processing based on region
        void setProcessingChromRegion();

        virtual void tagRead(HaplotagParameters &params, std::string& tagBamFile, const Genome& geneSample);
        virtual void postprocessForHaplotag(){};
        virtual void printExecutionReport();

        virtual GermlineHaplotagBamParser* createHaplotagBamParser(
            const ParsingBamConfig &config,
            const ParsingBamControl &control,
            ReadStatistics& readStats
        ){
            return new GermlineHaplotagBamParser(config, control, readStats, params);
        }

        virtual std::string getNormalSnpParsingMessage(){
            return "parsing SNP VCF ... ";
        };

        virtual std::string getTagReadStartMessage(){
            return "tag read start ...\n";
        };

    public:
        HaplotagProcess(HaplotagParameters &params);
        virtual void pipelineProcess();
        virtual ~HaplotagProcess();

};

#endif