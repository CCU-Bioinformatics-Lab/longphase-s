#ifndef HAPLOTAG_H
#define HAPLOTAG_H

#include "Util.h"
#include "HaplotagType.h"  // Include the header that defines HaplotagParameters
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

// Help message management class
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

        void addSectionItem(const std::string& sectionName, const std::string& newItem) {
            bool found = false;
            for (auto& section : sections) {
                if (section.header == sectionName) {
                    section.items.push_back(newItem);
                    found = true;
                    break;
                }
            }

            if (!found) {
                std::cerr << "Section '" << sectionName << "' not found." << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        void clearSectionItem(const std::string& sectionName) {
            bool found = false;
            for (auto& section : sections) {
                if (section.header == sectionName) {
                    section.items.clear();
                    found = true;
                    break;
                }
            }

            if (!found) {
                std::cerr << "Section '" << sectionName << "' not found." << std::endl;
                exit(EXIT_FAILURE);
            }
        }

    public:
        HelpMessageManager(const std::string& program) : programName(program) {}
        
        virtual ~HelpMessageManager() = default;
        
        // Build the basic message structure
        virtual void buildBaseMessage() {}
        
        // Allow derived classes to add new help messages
        virtual void modifyMessage() {}
        
        // Print the help message
        virtual void printHelp() const {
            for (const auto& section : sections) {
                std::cout << section.header << std::endl;
                for (const auto& item : section.items) {
                    std::cout << item << std::endl;
                }
            }
        }
        
        // Get the complete help message string
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

// Class to manage command line option definitions
class OptionManager {
    protected:
        std::string programName;
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

        // Validate if a required file exists
        bool validateRequiredFile(const std::string& filePath, const std::string& fileDescription);

        // Validate if an optional file exists (if specified)
        bool validateOptionalFile(const std::string& filePath, const std::string& fileDescription);

        virtual HelpMessageManager* createHelpManager(const std::string& program)= 0;
    public:
        // Getters for option definitions
        const struct option* getLongOpts() const { return longOpts.data(); }
        const char* getShortOpts() const { return shortOpts.c_str(); }
        
        // Allow derived classes to extend options
        virtual void extendOptions() {};


        OptionManager(const std::string& program) : programName(program) {};
        
        virtual ~OptionManager() = default;
};

class HaplotagHelpManager : public HelpMessageManager {
public:
    HaplotagHelpManager(const std::string& program) : HelpMessageManager(program) {
        buildBaseMessage();
    }

    virtual void buildBaseMessage() override;
};

// Haplotag-specific option definition manager
class HaplotagOptionManager : public OptionManager {
    protected:
        HaplotagParameters ecParams;

        HelpMessageManager* helpManager;

        virtual void initializeDefaultValues();

        bool loadHaplotagOptions(char& opt, std::istringstream& arg);

        // Validate all input files
        virtual bool validateFiles();

        virtual bool validateNumericParameter();

        virtual bool loadExtendOptions(char& opt, std::istringstream& arg);
        virtual bool validateExtendFiles();

        virtual HelpMessageManager* createHelpManager(const std::string& program) override {
            return new HaplotagHelpManager(program);
        }


    public:
        HaplotagOptionManager(const std::string& program);
        virtual ~HaplotagOptionManager();

        void setOptions() override;
        void setHelpMessage();
        void parseOptions(int argc, char** argv);

        void setVersion(std::string in_version) { ecParams.version = in_version; }
        HaplotagParameters getParams() const { return ecParams; }
};

// functions
int HaplotagMain(int argc, char** argv, std::string in_version);

#endif