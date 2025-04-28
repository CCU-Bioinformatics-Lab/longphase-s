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

// Normal header
class GermlineTagLog : public MessageManager {
    protected:
        const HaplotagParameters& params;

        void checkStreamStatus() {
            if (!tagReadLog || !tagReadLog->is_open()) {
                throw std::runtime_error("Output stream is not valid");
            }
        }

        void addCommonBasicEntries(){
            addEntry("snpFile", params.snpFile);       
            addEntry("svFile", params.svFile);       
            addEntry("bamFile", params.bamFile);    
            addEntry("resultPrefix", params.resultPrefix);
            addEntry("numThreads", params.numThreads);  
            addEntry("region", params.region);  
            addEntry("qualityThreshold", params.qualityThreshold);         
            addEntry("percentageThreshold", params.percentageThreshold); 
            addEntry("tagSupplementary", params.tagSupplementary); 
        }
        virtual void addBasicEntries(){
            addCommonBasicEntries();
        }

        virtual void writeBasicColumns(){
            *tagReadLog << "#ReadID\t"
                    << "CHROM\t"
                    << "ReadStart\t"
                    << "Confidnet(%)\t"
                    << "Haplotype\t"
                    << "PhaseSet\t"
                    << "TotalAllele\t"
                    << "HP1Allele\t"
                    << "HP2Allele\t"
                    << "phasingQuality(PQ)\t"
                    << "(Variant,HP)\t"
                    << "(PhaseSet,Variantcount)\n";
        }

    public:
        std::ofstream* tagReadLog;

        GermlineTagLog(const HaplotagParameters& params);
        ~GermlineTagLog();

        virtual void writeHeader() {
            addBasicEntries();
            // Sort entries by order
            std::sort(entries.begin(), entries.end(), 
                    [](const LogEntry& a, const LogEntry& b) {
                        return a.order < b.order;
                    });

            // Write all entries
            for (const auto& entry : entries) {
                *tagReadLog << "##" << entry.key << ":" << entry.value << "\n";
            }
            
            writeBasicColumns();
            // printMessage();
        }

        virtual void writeTagReadLog(TagReadLog& data);
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
            insertByIndex("TumorSnpFile", params.tumorSnpFile, 1);
            insertByIndex("tumorBamFile", params.tumorBamFile, 4);
            insertByIndex("tagTumor", params.tagTumorSnp, 8);
            insertByIndex("somaticCallingThreshold", params.somaticCallingMpqThreshold, 9);
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
        // std::ofstream *tagResult;
        GermlineTagLog* tagResult;

        ReadStatistics localReadStats;
        
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
            GermlineTagLog *tagResult,
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
    protected:
        ReadStatistics& readStats;
        GermlineTagLog *tagResult;

        std::unique_ptr<ChromosomeProcessor> createProcessor(const std::string &chr) override{
            return std::unique_ptr<ChromosomeProcessor>(new GermlineHaplotagChrProcessor(writeOutputBam, mappingQualityFilter, readStats, tagResult));
        };

        // create tag read log
        virtual GermlineTagLog* createTagReadLog(const HaplotagParameters& params){
            std::cout << "create germline tag read log" << std::endl;
            return new GermlineTagLog(params);
        };

    public:
        GermlineHaplotagBamParser(
            ParsingBamMode mode,
            bool writeOutputBam,
            bool mappingQualityFilter,
            ReadStatistics& readStats
        );
        ~GermlineHaplotagBamParser() override;

        void createTagLog(const HaplotagParameters& params);
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
            std::cout << "create somatic tag read log" << std::endl;
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

class HaplotagProcess: public SomaticJudgeBase
{
    protected:
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

        ReadStatistics readStats;

        ReadHpDistriLog *hpBeforeInheritance;
        ReadHpDistriLog *hpAfterInheritance;
        
        SomaticReadBenchmark highConSomaticData;
        
        std::time_t processBegin;

        virtual void parseVariantFiles(VcfParser& vcfParser);
        // load SNP, SV, MOD vcf file
        void parseCommonVariantFiles(VcfParser& vcfParser);

        // decide which genome type chrVec and chrLength belong to
        void setChrVecAndChrLength(Genome tagGeneType);

        // update chromosome processing based on region
        void setProcessingChromRegion();

        // calculate SNP counts
        void calculateSnpCounts();

        virtual void prepareForHaplotag();

        virtual void tagRead(HaplotagParameters &params, const Genome& geneType);
        void printAlignmentStaristics();
    public:
        HaplotagProcess(HaplotagParameters &params);
        void taggingProcess();
        virtual ~HaplotagProcess();

};

#endif