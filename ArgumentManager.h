#ifndef ARGUMENT_MANAGER_H
#define ARGUMENT_MANAGER_H

#include <getopt.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include "HaplotagType.h"


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

        virtual void initializeDefaultValues(){};
        
        // Initialize the base set of command line options
        virtual void setOptions() {}
        
        // Add a new option to the list
        void addOption(const struct option& opt);

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

        ArgumentManager(const std::string& program) : programName(program) {};
        
        virtual ~ArgumentManager() = default;
};


#endif
