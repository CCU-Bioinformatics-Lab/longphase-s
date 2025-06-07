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


class PurityPredictionOptionDefiner : public HaplotagOptionDefiner {
    public:
        virtual void defineOptions(ArgumentManager& manager) override;
        virtual ~PurityPredictionOptionDefiner() = default;
};

class PurityPredictionArgumentManager : public ArgumentTemManager<PurityPredictionParameters> {
    protected:

        virtual OptionDefiner* createOptionDefiner() override {
            return new PurityPredictionOptionDefiner;
        }

    public:
        PurityPredictionArgumentManager(const std::string& program, const std::string& version, const char* HELP_MESSAGE)
         : ArgumentTemManager<PurityPredictionParameters>(program, version, HELP_MESSAGE) {}
        virtual ~PurityPredictionArgumentManager() = default;
};



int PurityPredictionMain(int argc, char** argv, std::string in_version);

#endif