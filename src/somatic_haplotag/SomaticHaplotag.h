#ifndef SOMATIC_HAPLOTAG_H
#define SOMATIC_HAPLOTAG_H

#include "../shared/ArgumentManager.h"
#include "../haplotag/Haplotag.h"
#include "SomaticHaplotagProcess.h"


template<>
struct ParamsHandler<SomaticHaplotagParameters>{
    static void initialize(SomaticHaplotagParameters& params, const std::string& version);

    static bool loadArgument(SomaticHaplotagParameters& params, char& opt, std::istringstream& arg);

    static bool validateFiles(SomaticHaplotagParameters& params, const std::string& programName);

    static bool validateNumericParams(SomaticHaplotagParameters& params, const std::string& programName);

    static void recordCommand(SomaticHaplotagParameters& params, int argc, char** argv);

    static int getHelpEnumNum();
};


class SomaticHaplotagOptionDefiner : public HaplotagOptionDefiner {
    public:
        virtual void defineOptions(ArgumentManager& manager) override;
        virtual ~SomaticHaplotagOptionDefiner() = default;
};

class SomaticHaplotagArgManager : public ArgumentTemManager<SomaticHaplotagParameters> {
    protected:

        virtual OptionDefiner* createOptionDefiner() override {
            return new SomaticHaplotagOptionDefiner;
        }

    public:
        SomaticHaplotagArgManager(const std::string& program, const std::string& version, const char* HELP_MESSAGE)
         : ArgumentTemManager<SomaticHaplotagParameters>(program, version, HELP_MESSAGE) {}
        virtual ~SomaticHaplotagArgManager() = default;
};



// functions
int SomaticHaplotagMain(int argc, char** argv, std::string in_version);


#endif