#ifndef HAPLOTAGPROCESS_H
#define HAPLOTAGPROCESS_H

#include "HaplotagType.h"
#include "HaplotagParsingBam.h"
#include "HaplotagVcfParser.h"
#include "HaplotagLogging.h"

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

class HaplotagParamsMessage : public MessageManager{
    private:
        const HaplotagParameters& params;
    protected:

        void addCommonParamsMessage(){
            addEntry("phased SNP file", params.snpFile);
            addEntry("phased SV file", params.svFile);
            addEntry("phased MOD file", params.modFile);
            addEntry("input bam file", params.bamFile);
            addEntry("input ref file", params.fastaFile);
            addEntry("output bam file", params.resultPrefix + "." + params.outputFormat);
            addEntry("number of threads", params.numThreads);
            addEntry("write log file", params.writeReadLog);
            addEntry("log file", params.writeReadLog ? (params.resultPrefix+".out") : "");
            addEntry("#1", "-------------------------------------------");
            addEntry("tag region", !params.region.empty() ? params.region : "all");
            addEntry("filter mapping quality below", params.qualityThreshold);
            addEntry("percentage threshold", params.percentageThreshold);
            addEntry("tag supplementary", params.tagSupplementary);
            addEntry("#2", "-------------------------------------------");
        }

    public:
        HaplotagParamsMessage(const HaplotagParameters& params):params(params){}

        virtual void addParamsMessage(){
            addCommonParamsMessage();
        }

        void printParamsMessage(){
            // Print all entries with proper formatting
            for (const auto& entry : entries) {
                if (entry.key.at(0) == '#') {
                    std::cerr << entry.value << "\n";
                } else if(entry.key.at(0) == '[' || entry.key.at(0) == '\n') {
                    std::cerr << entry.key << "\n";
                } else {
                    std::cerr << std::left << std::setw(31) << entry.key << ": " << entry.value << "\n";
                }
            }
        }
};

// Normal header
class GermlineTagLog : public MessageManager {
    private:
        const HaplotagParameters& params;
    protected:

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
            sortEntries();

            // Write all entries
            for (const auto& entry : entries) {
                *tagReadLog << "##" << entry.key << ":" << entry.value << "\n";
            }
            
            writeBasicColumns();
            // printMessage();
        }

        virtual void writeTagReadLog(TagReadLog& data);
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
        const HaplotagParameters& params;
    protected:
        ReadStatistics& readStats;
        GermlineTagLog *tagResult;

        std::unique_ptr<ChromosomeProcessor> createProcessor(const std::string &chr) override{
            return std::unique_ptr<ChromosomeProcessor>(new GermlineHaplotagChrProcessor(writeOutputBam, mappingQualityFilter, readStats, tagResult));
        };

        // create tag read log
        virtual GermlineTagLog* createTagReadLog(){
            return new GermlineTagLog(params);
        };

    public:
        GermlineHaplotagBamParser(
            ParsingBamMode mode,
            bool writeOutputBam,
            bool mappingQualityFilter,
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
    protected:
        HaplotagParamsMessage paramsMessage;

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

        virtual void parseVariantFiles(VcfParser& vcfParser);
        // load SNP, SV, MOD vcf file
        void parseCommonVariantFiles(VcfParser& vcfParser);
        // decide which genome type chrVec and chrLength belong to
        virtual void setChrVecAndChrLength();
        // update chromosome processing based on region
        void setProcessingChromRegion();

        virtual void tagRead(HaplotagParameters &params, std::string& tagBamFile, const Genome& geneType);
        virtual void postprocessForHaplotag(){};
        void printAlignmentStaristics();

        virtual GermlineHaplotagBamParser* createHaplotagBamParser(
            ParsingBamMode mode,
            bool writeOutputBam,
            bool mappingQualityFilter,
            ReadStatistics& readStats
        ){
            return new GermlineHaplotagBamParser(mode, writeOutputBam, mappingQualityFilter, readStats, params);
        }

        virtual std::string getNormalSnpParsingMessage(){
            return "parsing SNP VCF ... ";
        };

        virtual std::string getTagReadStartMessage(){
            return "tag read start ...\n";
        };

    public:
        HaplotagProcess(HaplotagParameters &params);
        virtual void taggingProcess();
        virtual ~HaplotagProcess();

};

#endif