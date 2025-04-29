#ifndef HAPLOTAG_H
#define HAPLOTAG_H
#include "Util.h"
#include "HaplotagBase.h"  // Include the header that defines HaplotagParameters
#include <getopt.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>

#define SUBPROGRAM "haplotag"

// 幫助訊息管理類
class HelpMessageManager {
    protected:
        struct HelpSection {
            std::string header;
            std::vector<std::string> items;
        };
        
        std::vector<HelpSection> sections;
        std::string programName;
        
        void addSection(const std::string& header) {
            sections.push_back({header, {}});
        }
        
        void addItem(const std::string& item) {
            if (!sections.empty()) {
                sections.back().items.push_back(item);
            }
        }

    public:
        HelpMessageManager(const std::string& program) : programName(program) {}
        
        virtual ~HelpMessageManager() = default;
        
        // 建立基本訊息結構
        virtual void buildBaseMessage() {}
        
        // 允許衍生類添加新的幫助訊息
        virtual void extendMessage() {}
        
        // 輸出幫助訊息
        virtual void printHelp() const {
            for (const auto& section : sections) {
                std::cout << section.header << std::endl;
                for (const auto& item : section.items) {
                    std::cout << item << std::endl;
                }
            }
        }
        
        // 獲取完整幫助訊息字符串
        virtual std::string getHelpMessage() const {
            std::ostringstream oss;
            for (const auto& section : sections) {
                oss << section.header << "\n";
                for (const auto& item : section.items) {
                    oss << item << "\n";
                }
            }
            return oss.str();
        }
};


class HaplotagHelpManager : public HelpMessageManager {
public:
    HaplotagHelpManager() : HelpMessageManager(SUBPROGRAM) {
        buildBaseMessage();
    }

    virtual void buildBaseMessage() override;
};


// Class to manage command line option definitions
class OptionManager {
    protected:
        std::vector<struct option> longOpts;
        std::string shortOpts;
        
        // Initialize the base set of command line options
        virtual void setOptions() {}
        
        // Add a new option to the list
        void addOption(const struct option& opt) {
            // If this is not the terminator, add it before the last element
            if (opt.name != NULL && !longOpts.empty()) {
                longOpts.insert(longOpts.end() - 1, opt);
            } else {
                longOpts.push_back(opt);
            }
        }

    public:
        // Getters for option definitions
        const struct option* getLongOpts() const { return longOpts.data(); }
        const char* getShortOpts() const { return shortOpts.c_str(); }
        
        // Allow derived classes to extend options
        virtual void extendOptions() {}

        OptionManager() {}
        
        virtual ~OptionManager() = default;
};

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
    SC_MPQ,
    TAG_TUM,
    HIGH_CON,
    DISABLE_FILTER
};

// Haplotag-specific option definition manager
class HaplotagOptionManager : public OptionManager {
    protected:
        HaplotagParameters ecParams;
    public:
        HaplotagOptionManager() : OptionManager() {
            setOptions();
        }
        void setOptions() override;
        // Parse command line options and return parameters
        void parseOptions(int argc, char** argv);
        // Set the version of the program
        void setVersion(std::string in_version) { ecParams.version = in_version; }
        // Get the parameters of the program
        HaplotagParameters getParams() const { return ecParams; }
};

// functions
int HaplotagMain(int argc, char** argv, std::string in_version);

#endif