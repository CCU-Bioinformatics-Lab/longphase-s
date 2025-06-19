#ifndef PURITY_ESTIMATION_H
#define PURITY_ESTIMATION_H

#include "Haplotag.h"
#include "ArgumentManager.h"
#include "PurityEstimationProcess.h"

template<>
struct ParamsHandler<PurityEstimParameters>{
    static void initialize(PurityEstimParameters& params, const std::string& version);

    static bool loadArgument(PurityEstimParameters& params, char& opt, std::istringstream& arg);

    static bool validateFiles(PurityEstimParameters& params, const std::string& programName);

    static bool validateNumericParams(PurityEstimParameters& params, const std::string& programName);

    static void recordCommand(PurityEstimParameters& params, int argc, char** argv);

    static int getHelpEnumNum();
};


class PurityEstimOptDefiner : public HaplotagOptionDefiner {
    public:
        virtual void defineOptions(ArgumentManager& manager) override;
        virtual ~PurityEstimOptDefiner() = default;
};

class PurityEstimArgsManager : public ArgumentTemManager<PurityEstimParameters> {
    protected:

        virtual OptionDefiner* createOptionDefiner() override {
            return new PurityEstimOptDefiner;
        }

    public:
        PurityEstimArgsManager(const std::string& program, const std::string& version, const char* HELP_MESSAGE)
         : ArgumentTemManager<PurityEstimParameters>(program, version, HELP_MESSAGE) {}
        virtual ~PurityEstimArgsManager() = default;
};



int PurityEstimMain(int argc, char** argv, std::string in_version);

#endif