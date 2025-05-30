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

template<>
struct ParamsHandler<HaplotagParameters>{

    static void initialize(HaplotagParameters& params, const std::string& version);

    static bool loadArgument(HaplotagParameters& params, char& opt, std::istringstream& arg);

    static bool validateFiles(HaplotagParameters& params, const std::string& programName);

    static bool validateNumericParameter(HaplotagParameters& params, const std::string& programName);

    static void recordCommand(HaplotagParameters& params, int argc, char** argv);

    static int getHelpEnumNum();
};


class HaplotagHelpManager : public HelpMessageManager {
public:
    HaplotagHelpManager(const std::string& program) : HelpMessageManager(program) {}

    virtual void buildMessage() override;
    virtual ~HaplotagHelpManager() = default;
};

class HaplotagOptionDefiner : public OptionDefiner {
    public:
        virtual void defineOptions(ArgumentManager& manager) override;
        virtual ~HaplotagOptionDefiner() = default;
};

// Haplotag-specific option definition manager
class HaplotagArgumentManager : public ArgumentTemManager<HaplotagParameters> {
    protected:

        virtual HelpMessageManager* createHelpManager(const std::string& program) override {
            return new HaplotagHelpManager(program);
        }

        virtual OptionDefiner* createOptionDefiner() override {
            return new HaplotagOptionDefiner();
        }

    public:
        HaplotagArgumentManager(const std::string& program, const std::string& version);
        virtual ~HaplotagArgumentManager();
};

// functions
int HaplotagMain(int argc, char** argv, std::string in_version);

#endif