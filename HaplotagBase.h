#ifndef HAPLOTAG_BASE_H
#define HAPLOTAG_BASE_H

#include "Util.h"
#include "ParsingBam.h"
#include <cmath>
#include <iomanip>
#include <omp.h>
#include <climits>

class HaplotagBamParser;
class ChromosomeProcessor;
class CigarParser;


struct HaplotagParameters
{
    int numThreads;
    int qualityThreshold;
    
    double percentageThreshold;
    
    std::string snpFile;
    std::string tumorSnpFile;   
    std::string svFile;
    std::string modFile;
    std::string bamFile;
    std::string tumorBamFile;  
    std::string fastaFile;
    std::string resultPrefix;
    std::string region;
    std::string command;
    std::string version;
    std::string outputFormat;

    std::string benchmarkVcf;
    bool enableFilter;

    bool tagSupplementary;
    bool writeReadLog;
};

enum Genome
{
    NORMAL = 0,
    TUMOR = 1,
    HIGH_CON_SOMATIC = 2
};

enum Nitrogenous
{
    UNKOWN = 0,
    A = 1,
    C = 2,
    G = 3,
    T = 4,
};

enum VariantType
{
    NONE_VAR = 0,
    SNP = 1,
    INSERTION = 2,
    DELETION = 3,
    MNP = 4
};

enum SnpHP
{
    NONE_SNP = 0,
    GERMLINE_H1 = 1,
    GERMLINE_H2 = 2,
    SOMATIC_H3 = 3,
    SOMATIC_H4 = 4,
    SOMATIC_H5 = 5
};

enum ReadHP
{
    unTag = 0,
    H1 = 1,
    H2 = 2,
    H3 = 3,
    H4 = 4,
    H1_1 = 5,
    H1_2 = 6,
    H2_1 = 7,
    H2_2 = 8,
};

enum ParsingBamMode{
    SINGLE_THREAD = 0,
    MULTI_THREAD = 1
};

struct VarData{
    const static int NONE_PHASED_SET = -1;

    RefAlt allele;
    //phased set (-1 means no phase set)
    int PhasedSet;

    //haplotype
    std::string HP1;
    std::string HP2;

    VariantType variantType;

    bool is_phased_hetero;
    bool is_homozygous;
    bool is_unphased_hetero;

    bool isExistPhasedSet(){
        return PhasedSet != NONE_PHASED_SET;
    }

    void setVariantType(){
        if(allele.Ref.length() == 1 && allele.Alt.length() == 1){
            variantType = VariantType::SNP;
        }else if(allele.Ref.length() == 1 && allele.Alt.length() > 1){
            variantType = VariantType::INSERTION;
        }else if(allele.Ref.length() > 1 && allele.Alt.length() == 1){
            variantType = VariantType::DELETION;
        }else if((allele.Ref.length() > 1) && (allele.Ref.length() == allele.Alt.length())){
            variantType = VariantType::MNP;
        }else{
            throw std::runtime_error("(loadVariantType)Invalid allele: " + allele.Ref + " " + allele.Alt);
        }
    }

    VarData(): PhasedSet(NONE_PHASED_SET), HP1(""), HP2(""), variantType(VariantType::NONE_VAR)
             , is_phased_hetero(false), is_homozygous(false), is_unphased_hetero(false){}
};

struct MultiGenomeVar{
    // record the variants from the normal and tumor VCF files (normal:0, tumor:1, HighCon:2)
    std::map<Genome, VarData> Variant;
    bool isSomaticVariant;
    int somaticReadDeriveByHP;
    bool isExists(Genome type){
        return Variant.find(type) != Variant.end();
    }

    VarData& operator[](Genome type){
        return Variant[type];
    }
    
    MultiGenomeVar(): isSomaticVariant(false), somaticReadDeriveByHP(0){}
};

//record each type of base in specific position
struct PosBase{
    int A_count;
    int C_count;
    int G_count;
    int T_count;
    int unknow;
    int depth;
    int delCount;

    int MPQ_A_count;
    int MPQ_C_count;
    int MPQ_G_count;
    int MPQ_T_count;
    int MPQ_unknow;
    int filteredMpqDepth;

    float VAF;
    //Non-deletion Adjusted AF
    float nonDelVAF;
    float filteredMpqVAF;
    float lowMpqReadRatio;
    float delRatio;

    //germline haplotype imbalance ratio
    double germlineHaplotypeImbalanceRatio;
    double percentageOfGermlineHp;
    //snp position, read hp count
    std::map<int, int> ReadHpCount;
    
    PosBase(): A_count(0), C_count(0), G_count(0), T_count(0), unknow(0), depth(0), delCount(0)
             , MPQ_A_count(0), MPQ_C_count(0), MPQ_G_count(0), MPQ_T_count(0), MPQ_unknow(0), filteredMpqDepth(0) 
             , VAF(0.0), nonDelVAF(0.0), filteredMpqVAF(0.0), lowMpqReadRatio(0.0), delRatio(0.0)
             , germlineHaplotypeImbalanceRatio(0.0), percentageOfGermlineHp(0.0)
             , ReadHpCount(std::map<int, int>()){}

    int getBaseCount(const std::string& base) const {
        if (base == "A") return A_count;
        if (base == "T") return T_count;
        if (base == "C") return C_count;
        if (base == "G") return G_count;
        throw std::runtime_error("(getBaseCount)Invalid base: " + base);
    }

    int getMpqBaseCount(const std::string& base) const {
        if (base == "A") return MPQ_A_count;
        if (base == "T") return MPQ_T_count;
        if (base == "C") return MPQ_C_count;
        if (base == "G") return MPQ_G_count;
        throw std::runtime_error("(getMpqBaseCount)Invalid base: " + base);
    }
};

struct HP3_Info{
    //Case ratio 
    int totalCleanHP3Read;
    int pure_H1_1_read;
    int pure_H2_1_read;
    int pure_H3_read;
    int Mixed_HP_read;
    int unTag;

    int CaseReadCount;
    float pure_H1_1_readRatio;
    float pure_H2_1_readRatio;
    float pure_H3_readRatio;
    float Mixed_HP_readRatio;

    PosBase base;
    std::string GTtype;
    int somaticHp4Base;
    int somaticHp5Base;
    int somaticHp4BaseCount;
    int somaticHp5BaseCount;

    bool isHighConSomaticSNP;

    int somaticReadDeriveByHP;
    double shannonEntropy;
    int homopolymerLength;

    bool statisticPurity;

    // imbalance ratio for predict purity
    double allelicImbalanceRatio;
    double somaticHaplotypeImbalanceRatio;

    //interval snp filter information
    float MeanAltCountPerVarRead;
    float zScore;
    int intervalSnpCount;
    bool inDenseTumorInterval;
    
    // filter out by somatic feature filter
    bool isFilterOut;

    //readHp, count
    std::map<int, int> somaticReadHpCount;

    HP3_Info(): totalCleanHP3Read(0), pure_H1_1_read(0), pure_H2_1_read(0), pure_H3_read(0), Mixed_HP_read(0), unTag(0)
             , CaseReadCount(0), pure_H1_1_readRatio(0.0), pure_H2_1_readRatio(0.0), pure_H3_readRatio(0.0), Mixed_HP_readRatio(0.0)
             , base(), GTtype(""), somaticHp4Base(Nitrogenous::UNKOWN), somaticHp5Base(Nitrogenous::UNKOWN), somaticHp4BaseCount(0), somaticHp5BaseCount(0)
             , isHighConSomaticSNP(false), somaticReadDeriveByHP(0), shannonEntropy(0.0), homopolymerLength(0)
             , statisticPurity(false), allelicImbalanceRatio(0.0), somaticHaplotypeImbalanceRatio(0.0)
             , MeanAltCountPerVarRead(0.0), zScore(0.0), intervalSnpCount(0), inDenseTumorInterval(false)
             , isFilterOut(false), somaticReadHpCount(std::map<int, int>()){}
};

//record vcf information
struct VCF_Info 
{
    std::vector<std::string> chrVec;
    std::map<std::string, int> chrLength;
        
    // The number of SVs occurring on different haplotypes in a read
    // read id, sv haplotype, count
    std::map<std::string, std::map<int, int>> readSVHapCount;

    Genome gene_type;
};



struct ReadHpResult{
    std::map<int, int> readHpCounter;
    std::map<int, int> somaticBaseReadHpCounter;
    std::vector<float> deriveHPsimilarVec;
    int somaticSnpH3count;
    bool existDeriveByH1andH2;
    int deriveHP;
    int coverRegionStartPos;
    int coverRegionEndPos;
    ReadHpResult(): somaticSnpH3count(0), existDeriveByH1andH2(false), deriveHP(0),
                    coverRegionStartPos(INT_MAX), coverRegionEndPos(INT_MIN){}
};

struct chrReadHpResult{
    std::map<int, ReadHpResult> posReadHpResult;

    void recordReadHp(int &pos, int &hpResult, int &BaseHP);
    void recordDeriveHp(int &pos, int &deriveHP, float deriveHPsimilarity);
    void recordAlignCoverRegion(int& curVarPos, int &startPos, int &endPos);
};

class LogEntry {
    public:
        std::string key;
        std::string value;
        size_t order;

        LogEntry(const std::string& k, const std::string& v, size_t o) 
            : key(k), value(v), order(o) {}
};

class MessageManager{
    protected:
        std::vector<LogEntry> entries;
        size_t currentOrder;

        //Generally add entry (add at the end)
        template<typename T>
        void addEntry(const std::string& key, const T& value) {
            std::string valueStr = transformValueToString(value);

            entries.emplace_back(key, valueStr, currentOrder);
            // std::cerr << "Added entry: " << key << ":" << valueStr << " at index " << currentOrder << std::endl;
            currentOrder++;
        }

        // Insert entry after a specific order
        template<typename T>
        void insertByIndex(const std::string& key, const T& value, size_t insertIndex) {
            // move the order of the entries after the insertIndex
            for (auto& entry : entries) {
                if (entry.order >= insertIndex) {
                    entry.order ++;
                }
            }
            std::string valueStr = transformValueToString(value);
            entries.emplace_back(key, valueStr, insertIndex);
            // std::cerr << "Inserted entry: " << key << ":" << valueStr << " at index " << insertIndex << std::endl;
        }


        // Insert a new message after a specific key
        // If the target key is not found, append the message at the end
        template<typename T>
        void insertAfterKey(const std::string& newKey, const T& newValue, const std::string& targetKey) {
            // Find the position of the target key
            auto it = std::find_if(entries.begin(), entries.end(),
                [&targetKey](const LogEntry& entry) { return entry.key == targetKey; });
            
            if (it == entries.end()) {
                // If target key not found, append to the end
                std::cerr << "[ERROR](insertAfterKey) Target key not found: " << targetKey << std::endl;
                exit(1);
            }
            
            // Calculate insertion position (target key's order + 1)
            size_t insertIndex = it->order + 1;
            
            // Increment order for all entries after the insertion point
            for (auto& entry : entries) {
                if (entry.order >= insertIndex) {
                    entry.order++;
                }
            }
            
            // Insert the new entry
            std::string valueStr = transformValueToString(newValue);
            entries.emplace_back(newKey, valueStr, insertIndex);
        }

        void sortEntries(){
            std::sort(entries.begin(), entries.end(), [](const LogEntry& a, const LogEntry& b) {
                return a.order < b.order;
            });
        }

        void printMessage(){
            // Write all entries
            for (const auto& entry : entries) {
                std::cout << entry.key << ":" << entry.value << "\n";
            }
        }

        // basic type
        template<typename T>
        std::string transformValueToString(const T& value) {
            return std::to_string(value);
        }

        // Template Specialization for std::string
        std::string transformValueToString(const std::string& value) {
            return value;
        }

        // Template Specialization for bool
        std::string transformValueToString(const bool& value) {
            return value ? "true" : "false";
        }

        // Template Specialization for const char*
        std::string transformValueToString(const char* value) {
            return value;
        }

        // Template Specialization for char[]
        template<size_t N>
        std::string transformValueToString(const char (&value)[N]) {
            return std::string(value);
        }

        // Template Specialization for double
        std::string transformValueToString(const double& value) {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(2) << value;
            return oss.str();
        }

        // Template Specialization for float
        std::string transformValueToString(const float& value) {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(2) << value;
            return oss.str();
        }

    public:
        MessageManager() : currentOrder(0) {}

        virtual ~MessageManager() = default;

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

class ReadHpDistriLog{
    private :
        struct coverRegionInfo{
            int startPos;
            int endPos;
            int length;

            coverRegionInfo(): startPos(0), endPos(0), length(0){}
        };
        // chr, variant position (0-base), reads HP 
        // std::map<std::string, std::map<int, ReadHpResult>> chrVarReadHpResult;
        std::map<std::string, chrReadHpResult> chrVarReadHpResult;
    protected:

    public:
        ReadHpDistriLog();
        ~ReadHpDistriLog();

        // use in multi-thread scenario
        void loadChrKey(const std::string &chr);
        // Returns a pointer to ensure thread-safe access to chromosome results
        chrReadHpResult* getChrHpResultsPtr (const std::string &chr);

        // only use in single thread scenario
        void recordChrReadHp(const std::string &chr, int &pos, int &hpResult, int &BaseHP);
        void recordChrDeriveHp(const std::string &chr, int &pos, int &deriveHP, float deriveHPsimilarity);
        void recordChrAlignCoverRegion(const std::string &chr, int &pos, int &startPos, int &endPos);

        void writeReadHpDistriLog(const HaplotagParameters &params, std::string logPosfix, const std::vector<std::string> &chrVec);
        void writePosCoverRegionLog(const HaplotagParameters &params, std::string logPosfix, const std::vector<std::string> &chrVec);
        void writeTagReadCoverRegionLog(const HaplotagParameters &params, std::string logPosfix, const std::vector<std::string> &chrVec, std::map<std::string, int> &chrLength);
        void removeNotDeriveByH1andH2pos(const std::vector<std::string> &chrVec);
};


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
        
class VcfParser{
    private:
        Genome tagGeneType;
        bool parseSnpFile;
        bool parseSVFile;
        bool parseMODFile;
        bool integerPS;
        std::map<std::string, int> psIndex;

        void compressParser(std::string &variantFile, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
        void unCompressParser(std::string &variantFile, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
        virtual void parserProcess(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
    
    protected:

    public:
        VcfParser();
        VcfParser(Genome tagGeneType);
        virtual ~VcfParser();
        void setParseSnpFile(bool parseSnpFile);
        void setParseSVFile(bool parseSVFile);
        void setParseMODFile(bool parseMODFile);
        bool getParseSnpFile();
        void reset();
        void variantParser(std::string &variantFile, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat);
};

#endif
