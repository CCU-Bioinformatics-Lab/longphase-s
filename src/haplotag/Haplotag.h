#ifndef HAPLOTAG_H
#define HAPLOTAG_H

#include "../shared/Util.h"
#include "../shared/ArgumentManager.h"
#include "HaplotagType.h"
#include "HaplotagProcess.h"
#include <getopt.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

template<>
struct ParamsHandler<HaplotagParameters>{

    static void initialize(HaplotagParameters& params, const std::string& version);

    static bool loadArgument(HaplotagParameters& params, char& opt, std::istringstream& arg);

    static bool validateFiles(HaplotagParameters& params, const std::string& programName);

    static bool validateNumericParams(HaplotagParameters& params, const std::string& programName);

    static void recordCommand(HaplotagParameters& params, int argc, char** argv);

    static int getHelpEnumNum();
};

template<>
struct ParamsHandler<ParsingBamConfig>{

    static void initialize(ParsingBamConfig& params, const std::string& version);

    static bool loadArgument(ParsingBamConfig& params, char& opt, std::istringstream& arg);

    static bool validateNumericParams(ParsingBamConfig& params, const std::string& programName);

    static void recordCommand(ParsingBamConfig& params, int argc, char** argv);
};


class HaplotagOptionDefiner : public OptionDefiner {
    public:
        virtual void defineOptions(ArgumentManager& manager) override;
        virtual ~HaplotagOptionDefiner() = default;
};

// Haplotag-specific option definition manager
class HaplotagArgumentManager : public ArgumentTemManager<HaplotagParameters> {
    protected:

        virtual OptionDefiner* createOptionDefiner() override {
            return new HaplotagOptionDefiner();
        }

    public:
        HaplotagArgumentManager(const std::string& program, const std::string& version, const char* HELP_MESSAGE)
         : ArgumentTemManager<HaplotagParameters>(program, version, HELP_MESSAGE) {};
        ~HaplotagArgumentManager() = default;
};

// functions
int HaplotagMain(int argc, char** argv, std::string in_version);

#endif