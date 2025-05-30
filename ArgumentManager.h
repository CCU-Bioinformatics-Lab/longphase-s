#ifndef ARGUMENT_MANAGER_H
#define ARGUMENT_MANAGER_H

#include <getopt.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>

/**
 * @file ArgumentManager.h
 * @brief Argument management system combining strategy pattern and template specialization
 * 
 * This file implements a argument management system that combines:
 * 1. Strategy Pattern for option definition (OptionDefiner)
 * 2. Template Specialization for parameter handling (ParamsHandler)
 * 3. Template inheritance for type-safe parameter management (ArgumentTemManager)
 * 
 * Architecture Overview:
 * - ArgumentManager: Base class for argument parsing and management
 * - OptionDefiner: Strategy pattern for defining command-line options (inheritance-based reuse)
 * - ParamsHandler: Template specialization for parameter processing (type-specific handling)
 * - ArgumentTemManager: Template class combining argument management with type-safe parameters
 * 
 * This design allows for:
 * - Reusable option definitions through inheritance (e.g., SomaticHaplotag inherits Haplotag options)
 * - Type-safe parameter handling through template specialization
 * - Clean separation of concerns between option definition and parameter processing
 * - Easy extension for new parameter types and option sets
 */

class ArgumentManager;
class OptionDefiner;
class HelpMessageManager;


/**
 * @brief Parameter handler template - Uses template specialization for type-specific parameter processing
 * @tparam Params The parameter structure type (e.g., HaplotagParameters, SomaticHaplotagParameters)
 * 
 * This template provides a unified interface for parameter processing operations.
 * Each parameter type should have its own template specialization that implements
 * the specific logic for that parameter type.
 * 
 * The default implementation provides empty/false returns, requiring specialization
 * for actual functionality.
 * 
 * Template specializations should implement:
 * - initialize(): Set default values for parameters
 * - loadArgument(): Parse and load command-line arguments into parameters
 * - validateFiles(): Validate input/output file paths
 * - validateNumericParameter(): Validate numeric parameter ranges
 * - recordCommand(): Record the command line for logging
 * - getHelpEnumNum(): Return the help option enum value
 * 
 * Example specialization:
 * ```cpp
 * template<>
 * struct ParamsHandler<MyParameters> {
 *     static void initialize(MyParameters& params) {
 *         params.value = 42;
 *     }
 *     // ... other methods
 * };
 * ```
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
    static bool validateNumericParameter(Params& params, const std::string& programName) {return false;}

    // Record the command
    static void recordCommand(Params& params, int argc, char** argv) {}

    // Get help enum number
    static int getHelpEnumNum() {return 0;}
};

namespace FileValidator{
    // Validate if a required file exists
    bool validateRequiredFile(const std::string& filePath, const std::string& fileDescription, const std::string& programName);

    // Validate if an optional file exists (if specified)
    bool validateOptionalFile(const std::string& filePath, const std::string& fileDescription, const std::string& programName);
}

/**
 * @brief Help message management class - Manages and formats help/usage information
 * 
 * This class provides a structured way to build and display help messages for command-line
 * applications. It organizes help content into sections with headers and items, allowing
 * for clean and consistent help message formatting.
 * 
 * Features:
 * - Section-based organization of help content
 * - Hierarchical structure with headers and items
 * - Dynamic content modification through inheritance
 * - Consistent formatting across different applications
 * 
 * The class uses a section-based approach where each section has:
 * - A header string (e.g., "Required Arguments:", "Optional Arguments:")
 * - A list of items describing individual options or usage information
 * 
 * Derived classes can override buildMessage() to create application-specific
 * help content and modifyMessage() to customize existing content.
 * 
 * Usage:
 * ```cpp
 * class MyHelpManager : public HelpMessageManager {
 * public:
 *     void buildMessage() override {
 *         addSection("Usage:");
 *         addItem("myprogram [OPTIONS] input.file");
 *         addSection("Options:");
 *         addItem("-h, --help    Show this help");
 *     }
 * };
 * ```
 */
// Help message management class
class HelpMessageManager {
    protected:
        /**
         * @brief Structure representing a help section with header and items
         */
        struct HelpSection {
            std::string header;
            std::vector<std::string> items;
        };
        
        std::vector<HelpSection> sections;
        std::string programName;
        
        /**
         * @brief Add a new section with the given header
         * @param header The section header text
         */
        void addSection(const std::string& header);
        
        /**
         * @brief Add an item to the current section
         * @param item The item text to add
         */
        void addItem(const std::string& item);
        
        /**
         * @brief Add an item to a specific section by name
         * @param sectionName The name of the target section
         * @param newItem The item text to add
         */
        void addSectionItem(const std::string& sectionName, const std::string& newItem);
        
        /**
         * @brief Clear all items from a specific section
         * @param sectionName The name of the section to clear
         */
        void clearSectionItem(const std::string& sectionName);

    public:
        HelpMessageManager(const std::string& program) : programName(program) {}
        
        virtual ~HelpMessageManager() = default;
        
        /**
         * @brief Build the basic message structure
         * 
         * This virtual method should be overridden by derived classes to build
         * the complete help message structure. The default implementation is empty.
         */
        virtual void buildMessage() {}
        
        /**
         * @brief Allow derived classes to modify existing help messages
         * 
         * This method can be used to customize or extend help messages after
         * the basic structure has been built.
         */
        virtual void modifyMessage() {}
        
        /**
         * @brief Print the help message to standard output
         */
        virtual void printHelp() const;
        
        /**
         * @brief Get the complete help message as a string
         * @return The formatted help message string
         */
        virtual std::string getHelpMessage() const;
};


/**
 * @brief Option definer base class - Uses strategy pattern to define command line options
 * 
 * This class is the base class for the strategy pattern, responsible for defining 
 * command line option strategies. By inheriting this class, different option sets 
 * can be defined for different applications.
 * 
 * Design Philosophy:
 * - Strategy Pattern: Separates option definition logic from parameter managers
 * - Inheritance Reuse: Derived classes can inherit basic options and add additional ones
 * - Separation of Concerns: Option definition and parameter processing are managed separately
 * 
 * Usage:
 * 1. Inherit from OptionDefiner class
 * 2. Implement defineOptions() method to define specific options
 * 3. Use strategy object in ArgumentManager to set options
 * 
 * Example:
 * ```cpp
 * class MyOptionDefiner : public OptionDefiner {
 * public:
 *     void defineOptions(ArgumentManager& manager) override {
 *         manager.setShortOption("abc:");
 *         manager.addOption({"help", no_argument, NULL, HELP});
 *         manager.addOption({"input", required_argument, NULL, 'a'});
 *     }
 * };
 * ```
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
 * @brief Base class for command-line argument management
 * 
 * This class serves as the foundation for command-line argument parsing and management.
 * It integrates three key components using composition and factory patterns:
 * 
 * 1. **OptionDefiner**: Strategy pattern for defining command-line options
 *    - Uses factory method createOptionDefiner() to create appropriate strategy
 *    - Calls optionDefiner->defineOptions() to set up options
 *    - Allows inheritance-based reuse of option definitions
 * 
 * 2. **HelpMessageManager**: Manages help/usage information display
 *    - Uses factory method createHelpManager() to create appropriate help manager
 *    - Calls helpManager->buildMessage() to construct help content
 *    - Supports customized help messages through inheritance
 * 
 * 3. **Parameter Processing**: Virtual methods for parameter handling
 *    - loadArgument(): Process individual command-line arguments
 *    - validateFiles(): Validate input/output file paths
 *    - validateNumericParameter(): Validate numeric parameter ranges
 * 
 * **Workflow**:
 * ```cpp
 * ArgumentManager* manager = new ConcreteArgumentManager("myapp");
 * manager->setOptions();        // Creates OptionDefiner and defines options
 * manager->setHelpMessage();    // Creates HelpMessageManager and builds help
 * manager->parseOptions(argc, argv); // Parses command line arguments
 * ```
 * 
 * **Factory Pattern Integration**:
 * - createOptionDefiner(): Returns strategy for option definition
 * - createHelpManager(): Returns manager for help message construction
 * 
 * **Template Integration**:
 * This class works seamlessly with ArgumentTemManager template to provide
 * type-safe parameter handling through ParamsHandler specializations.
 * 
 * Example implementation:
 * ```cpp
 * class MyArgumentManager : public ArgumentManager {
 * protected:
 *     HelpMessageManager* createHelpManager(const std::string& program) override {
 *         return new MyHelpManager(program);
 *     }
 *     OptionDefiner* createOptionDefiner() override {
 *         return new MyOptionDefiner();
 *     }
 * public:
 *     MyArgumentManager(const std::string& program) : ArgumentManager(program) {}
 * };
 * ```
 */
// Class to manage command line option definitions
class ArgumentManager {

    protected:
        std::string programName;
        std::string version;
        std::vector<struct option> longOpts;
        std::string shortOpts;

        HelpMessageManager* helpManager;
        OptionDefiner* optionDefiner;

        // Create the help message manager and option definer
        /**
         * @brief Factory method to create appropriate HelpMessageManager
         * @param program The program name for help message display
         * @return Pointer to HelpMessageManager instance
         * 
         * This factory method should be overridden by derived classes to return
         * the appropriate HelpMessageManager subclass for their specific needs.
         * The returned object will be used to build and display help messages.
         */
        virtual HelpMessageManager* createHelpManager(const std::string& program) = 0;
        
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

        virtual bool validateNumericParameter(){return false;};
    public:
        /**
         * @brief Set short option string for getopt
         * @param opt The short option string (e.g., "abc:" for -a, -b, -c with argument)
         */
        void setShortOption(const std::string& opt);

        /**
         * @brief Add a long option to the options list
         * @param opt The option structure for getopt_long
         */
        void addOption(const struct option& opt);

        /**
         * @brief Remove an option from the options list
         * @param name The name of the option to remove
         */
        void removeOption(const std::string& name);

        /**
         * @brief Initialize command-line options using OptionDefiner strategy
         * 
         * This method creates an OptionDefiner instance via createOptionDefiner()
         * and calls its defineOptions() method to set up all command-line options.
         * This implements the strategy pattern for option definition.
         */
        void setOptions() {
            optionDefiner = createOptionDefiner();
            optionDefiner->defineOptions(*this);
        }

        /**
         * @brief Initialize help message using HelpMessageManager
         * 
         * This method creates a HelpMessageManager instance via createHelpManager()
         * and calls its buildMessage() method to construct the help content.
         * The help manager can then be used to display help information.
         */
        void setHelpMessage(){
            helpManager = createHelpManager(programName);
            helpManager->buildMessage();
        };

        /**
         * @brief Clean up allocated resources
         * 
         * Safely deletes HelpMessageManager and OptionDefiner instances
         * if they have been allocated.
         */
        void destroy(){
            if(helpManager) delete helpManager;
            if(optionDefiner) delete optionDefiner;
        };

        // Parse command line options
        virtual void parseOptions(int argc, char** argv);

        ArgumentManager(const std::string& program, const std::string& version)
         : programName(program), version(version) {};
        
        virtual ~ArgumentManager(){
            destroy();
        };
};

/**
 * @brief Template-based argument manager - Combines ArgumentManager with type-safe parameter handling
 * @tparam Params The parameter structure type (e.g., HaplotagParameters, SomaticHaplotagParameters)
 * 
 * This template class extends ArgumentManager to provide type-safe parameter management
 * using template specialization. It automatically delegates parameter processing to the
 * appropriate ParamsHandler specialization for the given Params type.
 * 
 * Features:
 * - Type-safe parameter storage and access
 * - Automatic delegation to ParamsHandler specializations
 * - Unified interface for all parameter types
 * - Integration with strategy pattern for option definition
 * 
 * The class maintains a Params object and uses ParamsHandler<Params> to process
 * all parameter-related operations. This allows each parameter type to have its
 * own specialized behavior while maintaining a consistent interface.
 * 
 * Usage:
 * ```cpp
 * class MyArgumentManager : public ArgumentTemManager<MyParameters> {
 * protected:
 *     OptionDefiner* createOptionDefiner() override {
 *         return new MyOptionDefiner();
 *     }
 * public:
 *     MyArgumentManager(const std::string& program) 
 *         : ArgumentTemManager<MyParameters>(program) {}
 * };
 * ```
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

        virtual bool validateNumericParameter() override {
            return paramsHandler.validateNumericParameter(params, programName);
        };

        virtual void recordCommand(int argc, char** argv) override {
            paramsHandler.recordCommand(params, argc, argv);
        };

        virtual HelpMessageManager* createHelpManager(const std::string& program) override {
            return new HelpMessageManager(program);
        };

        virtual int getHelpEnumNum() override {
            return ParamsHandler<Params>::getHelpEnumNum();
        };
        
    public:
        ArgumentTemManager(const std::string& program, const std::string& version)
         : ArgumentManager(program, version) {};

        Params getParams() const {
            return params;
        }
        virtual ~ArgumentTemManager(){};  
};


#endif
