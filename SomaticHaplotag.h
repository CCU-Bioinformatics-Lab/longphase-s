#ifndef SOMATIC_HAPLOTAG_H
#define SOMATIC_HAPLOTAG_H

#include "Haplotag.h"
#include "ArgumentManager.h"
#include "SomaticHaplotagProcess.h"

template<>
struct ParamsHandler<SomaticHaplotagParameters>{
    static void initialize(SomaticHaplotagParameters& params) {

        ParamsHandler<HaplotagParameters>::initialize(params.baseParams);

        params.tumorPurity = 0.2;
        params.predictTumorPurity = false;
        params.onlyPredictTumorPurity = false;
        params.enableFilter = false;
    }

    static bool loadArgument(SomaticHaplotagParameters& params, char& opt, std::istringstream& arg) {
        // load base haplotag options
        bool isLoaded = ParamsHandler<HaplotagParameters>::loadArgument(params.baseParams, opt, arg);
        
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

struct SomaticHaplotagParamHandler {
    static void initialize(SomaticHaplotagParameters& params) {

        HaplotagParamHandler::initialize(params.baseParams);

        params.tumorPurity = 0.2;
        params.predictTumorPurity = false;
        params.onlyPredictTumorPurity = false;
        params.enableFilter = false;
    }

    static bool loadArgument(SomaticHaplotagParameters& params, char& opt, std::istringstream& arg) {
        // load base haplotag options
        bool isLoaded = HaplotagParamHandler::loadArgument(params.baseParams, opt, arg);
        
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
};



// functions
int SomaticHaplotagMain(int argc, char** argv, std::string in_version);


#endif