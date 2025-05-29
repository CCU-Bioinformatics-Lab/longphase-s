#ifndef SOMATIC_HAPLOTAG_H
#define SOMATIC_HAPLOTAG_H

#include "Haplotag.h"
#include "ArgumentManager.h"
#include "SomaticHaplotagProcess.h"

template<>
struct ParamsHandler<SomaticHaplotagParameters>{
    static void initialize(SomaticHaplotagParameters& params) {

        ParamsHandler<HaplotagParameters>::initialize(params.basic);

        params.tumorPurity = 0.2;
        params.tumorPurity = 0.2;
        params.enableFilter = true;
        params.predictTumorPurity = true;
        params.onlyPredictTumorPurity = false;
    }

    static bool loadArgument(SomaticHaplotagParameters& params, char& opt, std::istringstream& arg) {
        // load base haplotag options
        bool isLoaded = ParamsHandler<HaplotagParameters>::loadArgument(params.basic, opt, arg);
        
        if(!isLoaded){
            //reset isLoaded
            isLoaded = true;
            //load somatic haplotag options
            switch (opt)
            {
                case HaplotagOption::TUM_SNP: arg >> params.tumorSnpFile; break;
                case HaplotagOption::TUM_BAM: arg >> params.tumorBamFile; break;
                case HaplotagOption::BENCHMARK_VCF: arg >> params.benchmarkVcf; break;
                case HaplotagOption::BENCHMARK_BED: arg >> params.benchmarkBedFile; break;
                case HaplotagOption::DISABLE_FILTER: params.enableFilter = false; break;
                case HaplotagOption::TUMOR_PURITY: 
                    arg >> params.tumorPurity; 
                    params.predictTumorPurity = false;
                    break;
                default: isLoaded = false; 
                break;
            }
        }
        return isLoaded;
    }

    static bool validateFiles(SomaticHaplotagParameters& params, const std::string& programName) {
        // validate base haplotag files
        bool isValid = ParamsHandler<HaplotagParameters>::validateFiles(params.basic, programName);
        // validate somatic haplotag files
        isValid &= FileValidator::validateRequiredFile(params.tumorSnpFile, "tumor SNP file", programName);
        isValid &= FileValidator::validateRequiredFile(params.tumorBamFile, "tumor BAM file", programName);
        isValid &= FileValidator::validateOptionalFile(params.benchmarkVcf, "benchmark VCF file", programName);
        isValid &= FileValidator::validateOptionalFile(params.benchmarkBedFile, "benchmark BED file", programName);
        return isValid;
    }

    static bool validateNumericParameter(SomaticHaplotagParameters& params, const std::string& programName) {
        bool isValid = ParamsHandler<HaplotagParameters>::validateNumericParameter(params.basic, programName);
    
        if (params.tumorPurity < 0.1 || params.tumorPurity > 1.0) {
            std::cerr << "[ERROR] " << programName << ": invalid tumor purity. value: " 
                    << params.tumorPurity 
                    << "\nthis value need: 0.1~1.0, --tumor-purity=Number\n";
            isValid = false;
        }
        
        return isValid;  
    }

    static void recordCommand(SomaticHaplotagParameters& params, int argc, char** argv) {
        for(int i = 0; i < argc; ++i){
            params.basic.command.append(argv[i]);
            params.basic.command.append(" ");
        }
    }
};

struct SomaticHaplotagParamHandler {
    static void initialize(SomaticHaplotagParameters& params) {

        HaplotagParamHandler::initialize(params.basic);

        params.tumorPurity = 0.2;
        params.predictTumorPurity = false;
        params.onlyPredictTumorPurity = false;
        params.enableFilter = false;
    }

    static bool loadArgument(SomaticHaplotagParameters& params, char& opt, std::istringstream& arg) {
        // load base haplotag options
        bool isLoaded = HaplotagParamHandler::loadArgument(params.basic, opt, arg);
        
        if(!isLoaded){
            //reset isLoaded
            isLoaded = true;
            //load somatic haplotag options
            switch (opt)
            {
                case HaplotagOption::TUM_SNP: arg >> params.tumorSnpFile; break;
                case HaplotagOption::TUM_BAM: arg >> params.tumorBamFile; break;
                case HaplotagOption::BENCHMARK_VCF: arg >> params.benchmarkVcf; break;
                case HaplotagOption::BENCHMARK_BED: arg >> params.benchmarkBedFile; break;
                case HaplotagOption::DISABLE_FILTER: params.enableFilter = false; break;
                case HaplotagOption::TUMOR_PURITY: 
                    arg >> params.tumorPurity; 
                    params.predictTumorPurity = false;
                    break;
                default: isLoaded = false; 
                break;
            }
        }
        return isLoaded;
    }

};

class SomaticHaplotagHelpManager : public HaplotagHelpManager {
    public:
        SomaticHaplotagHelpManager(const std::string& program) : HaplotagHelpManager(program) {}

        virtual void buildMessage() override;
    
};

class SomaticHaplotagArgumentManager : public HaplotagArgumentManager {
    SomaticHaplotagParameters ecParams;
    protected:
        virtual HelpMessageManager* createHelpManager(const std::string& program) override {
            return new SomaticHaplotagHelpManager(program);
        }

        virtual void initializeDefaultValues() override;

        virtual bool loadArgument(char& opt, std::istringstream& arg) override;

        virtual bool validateFiles() override;
        virtual bool validateNumericParameter() override;

    public:
        SomaticHaplotagArgumentManager(const std::string& program) : HaplotagArgumentManager(program) {}
        virtual void setOptions() override;
        virtual ~SomaticHaplotagArgumentManager() = default;
        SomaticHaplotagParameters& getParams() { return ecParams; }
};



// functions
int SomaticHaplotagMain(int argc, char** argv, std::string in_version);


#endif