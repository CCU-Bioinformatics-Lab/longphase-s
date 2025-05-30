#ifndef SOMATIC_HAPLOTAG_H
#define SOMATIC_HAPLOTAG_H

#include "Haplotag.h"
#include "ArgumentManager.h"
#include "SomaticHaplotagProcess.h"

template<>
struct ParamsHandler<SomaticHaplotagParameters>{
    static void initialize(SomaticHaplotagParameters& params, const std::string& version);

    static bool loadArgument(SomaticHaplotagParameters& params, char& opt, std::istringstream& arg);

    static bool validateFiles(SomaticHaplotagParameters& params, const std::string& programName);

    static bool validateNumericParameter(SomaticHaplotagParameters& params, const std::string& programName);

    static void recordCommand(SomaticHaplotagParameters& params, int argc, char** argv);

    static int getHelpEnumNum();
};


class SomaticHaplotagHelpManager : public HaplotagHelpManager {
    public:
        SomaticHaplotagHelpManager(const std::string& program) : HaplotagHelpManager(program) {}
        virtual void buildMessage() override;
        virtual ~SomaticHaplotagHelpManager() = default;
    
};

class SomaticHaplotagOptionDefiner : public HaplotagOptionDefiner {
    public:
        virtual void defineOptions(ArgumentManager& manager) override;
        virtual ~SomaticHaplotagOptionDefiner() = default;
};

class SomaticHaplotagArgumentManager : public ArgumentTemManager<SomaticHaplotagParameters> {
    protected:
        virtual HelpMessageManager* createHelpManager(const std::string& program) override {
            return new SomaticHaplotagHelpManager(program);
        }

        virtual OptionDefiner* createOptionDefiner() override {
            return new SomaticHaplotagOptionDefiner;
        }

    public:
        SomaticHaplotagArgumentManager(const std::string& program, const std::string& version)
         : ArgumentTemManager<SomaticHaplotagParameters>(program, version) {}
        virtual ~SomaticHaplotagArgumentManager() = default;
};



// functions
int SomaticHaplotagMain(int argc, char** argv, std::string in_version);


#endif