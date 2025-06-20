#ifndef ARGUMENT_MANAGER_H
#define ARGUMENT_MANAGER_H

#include <getopt.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>


class ArgumentManager;
class OptionDefiner;

/**
 * @brief Template struct for handling parameter-specific operations
 * @tparam Params The parameter type to handle
 * 
 * This template struct provides a framework for parameter initialization,
 * argument loading, validation, and command recording. Derived classes
 * should specialize this template for their specific parameter types.
 */
template<typename Params>
struct ParamsHandler{
    // Initialize the parameters
    static void initialize(Params& params, const std::string& version) {}

    // Load the arguments
    static bool loadArgument(Params& params, char& opt, std::istringstream& arg) {return false;}

    // Validate the files
    static bool validateFiles(Params& params, const std::string& programName) {return false;}

    // Validate the numeric parameters
    static bool validateNumericParams(Params& params, const std::string& programName) {return false;}

    // Record the command
    static void recordCommand(Params& params, int argc, char** argv) {}

    // Get help enum number
    static int getHelpEnumNum() {return 0;}
};

/**
 * @brief Namespace containing file validation utilities
 * 
 * Provides functions to validate required and optional files
 */
namespace FileValidator{
    // Validate if a required file exists
    bool validateRequiredFile(const std::string& filePath, const std::string& fileDescription, const std::string& programName);

    // Validate if an optional file exists (if specified)
    bool validateOptionalFile(const std::string& filePath, const std::string& fileDescription, const std::string& programName);
}


/**
 * @brief Abstract base class for defining command line options
 * 
 * Derived classes must implement defineOptions() to specify their
 * command line argument structure
 */
class OptionDefiner {
    public:
        virtual ~OptionDefiner() = default;
        
        /**
         * @brief Define command line options
         * @param manager Reference to ArgumentManager for adding options
         * 
         * This pure virtual function must be implemented by derived classes to define
         * specific command line options. Derived classes can call manager.setShortOption()
         * and manager.addOption() to set short option strings and long option lists.
         */
        virtual void defineOptions(ArgumentManager& manager) = 0;
};

/**
 * @brief Base class for managing command line argument parsing
 * 
 * Provides functionality for parsing command line arguments, validating
 * files and parameters, and managing option definitions through OptionDefiner
 */
class ArgumentManager {

    protected:
        // program name and version
        std::string programName;
        std::string version;
        const char* HELP_MESSAGE;

        std::vector<struct option> longOpts;
        std::string shortOpts;

        OptionDefiner* optionDefiner;
        
        /**
         * @brief Factory method to create appropriate OptionDefiner
         * @return Pointer to OptionDefiner instance
         * 
         * This factory method should be overridden by derived classes to return
         * the appropriate OptionDefiner subclass that defines their specific
         * command-line options.
         */
        virtual OptionDefiner* createOptionDefiner() = 0;

        virtual void initializeDefaultValues(){};

        virtual int getHelpEnumNum()=0;
        
        // Getters for option definitions
        const struct option* getLongOpts() const { return longOpts.data(); }
        const char* getShortOpts() const { return shortOpts.c_str(); }


        // Validate if a required file exists
        bool validateRequiredFile(const std::string& filePath, const std::string& fileDescription);

        // Validate if an optional file exists (if specified)
        bool validateOptionalFile(const std::string& filePath, const std::string& fileDescription);


        virtual bool loadArgument(char& opt, std::istringstream& arg){return false;};

        virtual void recordCommand(int argc, char** argv){};

        // Validate all input files
        virtual bool validateFiles(){return false;};

        virtual bool validateNumericParams(){return false;};
    public:

        void setShortOption(const std::string& opt);

        void addOption(const struct option& opt);

        void setOptions() {
            optionDefiner = createOptionDefiner();
            optionDefiner->defineOptions(*this);
        }

        void destroy(){
            if(optionDefiner) delete optionDefiner;
        };

        // Parse command line options
        virtual void parseOptions(int argc, char** argv);

        ArgumentManager(const std::string& program, const std::string& version, const char* HELP_MESSAGE)
         : programName(program), version(version), HELP_MESSAGE(HELP_MESSAGE) {};
        
        virtual ~ArgumentManager(){
            destroy();
        };
};


/**
 * @brief Template class for managing typed parameters with argument parsing
 * @tparam Params The parameter type to manage
 * 
 * This template class extends ArgumentManager to provide type-safe parameter
 * management. It uses ParamsHandler<Params> to handle parameter-specific operations.
 */
template<typename Params>
class ArgumentTemManager : public ArgumentManager{
    private:
        Params params;
        ParamsHandler<Params> paramsHandler;
    
    protected:
        virtual void initializeDefaultValues() override {
            paramsHandler.initialize(params, version);
        };

        virtual bool loadArgument(char& opt, std::istringstream& arg) override {
            return paramsHandler.loadArgument(params, opt, arg);
        };
        
        virtual bool validateFiles() override {
            return paramsHandler.validateFiles(params, programName);
        };

        virtual bool validateNumericParams() override {
            return paramsHandler.validateNumericParams(params, programName);
        };

        virtual void recordCommand(int argc, char** argv) override {
            paramsHandler.recordCommand(params, argc, argv);
        };

        virtual int getHelpEnumNum() override {
            return ParamsHandler<Params>::getHelpEnumNum();
        };
        
    public:
        ArgumentTemManager(const std::string& program, const std::string& version, const char* HELP_MESSAGE)
         : ArgumentManager(program, version, HELP_MESSAGE) {};

        Params getParams() const {
            return params;
        }
        virtual ~ArgumentTemManager(){};  
};


#endif
