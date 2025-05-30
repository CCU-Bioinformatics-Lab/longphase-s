#ifndef PURITY_PREDICTION_H
#define PURITY_PREDICTION_H

#include "Haplotag.h"
#include "ArgumentManager.h"
#include "PurityPredictionProcess.h"

template<>
struct ParamsHandler<PurityPredictionParameters>{
    static void initialize(PurityPredictionParameters& params, const std::string& version);

    static bool loadArgument(PurityPredictionParameters& params, char& opt, std::istringstream& arg);

    static bool validateFiles(PurityPredictionParameters& params, const std::string& programName);

    static bool validateNumericParameter(PurityPredictionParameters& params, const std::string& programName);

    static void recordCommand(PurityPredictionParameters& params, int argc, char** argv);

    static int getHelpEnumNum();
};


class PurityPredictionHelpManager : public HaplotagHelpManager {
    public:
        PurityPredictionHelpManager(const std::string& program) : HaplotagHelpManager(program) {}
        virtual void buildMessage() override;
        virtual ~PurityPredictionHelpManager() = default;
    
};

class PurityPredictionOptionDefiner : public HaplotagOptionDefiner {
    public:
        virtual void defineOptions(ArgumentManager& manager) override;
        virtual ~PurityPredictionOptionDefiner() = default;
};

class PurityPredictionArgumentManager : public ArgumentTemManager<PurityPredictionParameters> {
    protected:
        virtual HelpMessageManager* createHelpManager(const std::string& program) override {
            return new PurityPredictionHelpManager(program);
        }

        virtual OptionDefiner* createOptionDefiner() override {
            return new PurityPredictionOptionDefiner;
        }

    public:
        PurityPredictionArgumentManager(const std::string& program, const std::string& version)
         : ArgumentTemManager<PurityPredictionParameters>(program, version) {}
        virtual ~PurityPredictionArgumentManager() = default;
};



int PurityPredictionMain(int argc, char** argv, std::string in_version);

#endif