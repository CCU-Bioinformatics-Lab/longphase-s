#ifndef HAPLOTAG_H
#define HAPLOTAG_H

#include "Util.h"
#include "HaplotagType.h"  // Include the header that defines HaplotagParameters
#include "ArgumentManager.h"
#include <getopt.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

// Option identifiers enum
enum HaplotagOption {
    OPT_HELP = 1,
    TAG_SUP,
    SV_FILE,
    REGION,
    LOG,
    MOD_FILE,
    CRAM,
    TUM_SNP,
    TUM_BAM,
    BENCHMARK_VCF,
    BENCHMARK_BED,
    DISABLE_FILTER,
    TUMOR_PURITY
};

class HaplotagHelpManager : public HelpMessageManager {
public:
    HaplotagHelpManager(const std::string& program) : HelpMessageManager(program) {}

    virtual void buildMessage() override;
};

// Haplotag-specific option definition manager
class HaplotagArgumentManager : public ArgumentManager {
    protected:
        HaplotagParameters ecParams;

        HelpMessageManager* helpManager;

        virtual void initializeDefaultValues() override;

        virtual bool loadOptions(char& opt, std::istringstream& arg);

        virtual void recordCommand(int argc, char** argv);

        // Validate all input files
        virtual bool validateFiles();

        virtual bool validateNumericParameter();

        virtual HelpMessageManager* createHelpManager(const std::string& program) override {
            return new HaplotagHelpManager(program);
        }


    public:
        HaplotagArgumentManager(const std::string& program);
        virtual ~HaplotagArgumentManager();

        virtual void setOptions() override;
        void setHelpMessage();
        void parseOptions(int argc, char** argv);

        void setVersion(std::string in_version) { ecParams.version = in_version; }
        HaplotagParameters getParams() const { return ecParams; }
};

// functions
int HaplotagMain(int argc, char** argv, std::string in_version);

#endif