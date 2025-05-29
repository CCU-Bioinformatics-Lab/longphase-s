#ifndef ARGUMENT_MANAGER_H
#define ARGUMENT_MANAGER_H

#include <getopt.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include "HaplotagType.h"

template<typename Params>
struct ParamsHandler{
    static void initialize(Params& params) {}

    static bool loadArgument(Params& params, char& opt, std::istringstream& arg) {return false;}

    static bool validateFiles(Params& params, const std::string& programName) {return false;}

    static bool validateNumericParameter(Params& params, const std::string& programName) {return false;}

    static void recordCommand(Params& params, int argc, char** argv) {}
};

namespace FileValidator{
    // Validate if a required file exists
    bool validateRequiredFile(const std::string& filePath, const std::string& fileDescription, const std::string& programName);

    // Validate if an optional file exists (if specified)
    bool validateOptionalFile(const std::string& filePath, const std::string& fileDescription, const std::string& programName);
}

// Help message management class
class HelpMessageManager {
    protected:
        struct HelpSection {
            std::string header;
            std::vector<std::string> items;
        };
        
        std::vector<HelpSection> sections;
        std::string programName;
        
        void addSection(const std::string& header);
        void addItem(const std::string& item);
        void addSectionItem(const std::string& sectionName, const std::string& newItem);
        void clearSectionItem(const std::string& sectionName);

    public:
        HelpMessageManager(const std::string& program) : programName(program) {}
        
        virtual ~HelpMessageManager() = default;
        
        // Build the basic message structure
        virtual void buildMessage() {}
        
        // Allow derived classes to add new help messages
        virtual void modifyMessage() {}
        
        // Print the help message
        virtual void printHelp() const;
        
        // Get the complete help message string
        virtual std::string getHelpMessage() const;
};

// Class to manage command line option definitions
class ArgumentManager {

    protected:
        std::string programName;
        std::vector<struct option> longOpts;
        std::string shortOpts;

        HelpMessageManager* helpManager;

        virtual void initializeDefaultValues(){};
        

        virtual int getHelpEnumNum()=0;
        
        // Getters for option definitions
        const struct option* getLongOpts() const { return longOpts.data(); }
        const char* getShortOpts() const { return shortOpts.c_str(); }

        // Add a new option to the list
        void addOption(const struct option& opt);

        // Remove an option from the list
        void removeOption(const std::string& name);

        // Validate if a required file exists
        bool validateRequiredFile(const std::string& filePath, const std::string& fileDescription);

        // Validate if an optional file exists (if specified)
        bool validateOptionalFile(const std::string& filePath, const std::string& fileDescription);

        virtual HelpMessageManager* createHelpManager(const std::string& program)= 0;



        virtual bool loadArgument(char& opt, std::istringstream& arg){return false;};

        virtual void recordCommand(int argc, char** argv){};

        // Validate all input files
        virtual bool validateFiles(){return false;};

        virtual bool validateNumericParameter(){return false;};
    public:
        // Initialize the base set of command line options
        virtual void setOptions() {}


        void setHelpMessage(){
            helpManager = createHelpManager(programName);
            helpManager->buildMessage();
        };

        void destroyHelpMessage(){
            if(helpManager) delete helpManager;
        };

        // Parse command line options
        virtual void parseOptions(int argc, char** argv);

        ArgumentManager(const std::string& program) : programName(program) {};
        
        virtual ~ArgumentManager(){
            destroyHelpMessage();
        };
};


#endif
