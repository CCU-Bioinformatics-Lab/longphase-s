#include "ArgumentManager.h"

/**
 * @brief Set the short options string for getopt
 * @param opt The short options string (e.g., "s:b:o:t:q:p:r:")
 */
void ArgumentManager::setShortOption(const std::string& opt) {
    shortOpts = opt;
}

/**
 * @brief Add a long option to the options list
 * @param opt The option structure to add
 * 
 * Inserts the option before the terminator if it's not the terminator itself
 */
void ArgumentManager::addOption(const struct option& opt) {
    // If this is not the terminator, add it before the last element
    if (opt.name != NULL && !longOpts.empty()) {
        longOpts.insert(longOpts.end() - 1, opt);
    } else {
        longOpts.push_back(opt);
    }
}

/**
 * @brief Validate that a required file exists and is accessible
 * @param filePath Path to the file to validate
 * @param fileDescription Description of the file for error messages
 * @return true if file exists and is accessible, false otherwise
 */
bool ArgumentManager::validateRequiredFile(const std::string& filePath, const std::string& fileDescription) {
    if(filePath.empty()) {
        std::cerr << "[ERROR] " << programName  << ": missing " << fileDescription << ".\n";

        return false;
    }
    
    std::ifstream openFile(filePath.c_str());
    if(!openFile.is_open()) {
        std::cerr << "[ERROR] " << programName  << ": " << fileDescription << ": " << filePath << " not exist.\n\n";
        return false;
    }
    return true;
}

/**
 * @brief Validate that an optional file exists and is accessible (if specified)
 * @param filePath Path to the file to validate
 * @param fileDescription Description of the file for error messages
 * @return true if file is not specified or exists and is accessible, false otherwise
 */
bool ArgumentManager::validateOptionalFile(const std::string& filePath, const std::string& fileDescription) {
    if(filePath.empty()) {
        // Optional file not specified, that's OK
        return true;  
    }
    
    std::ifstream openFile(filePath.c_str());
    if(!openFile.is_open()) {
        std::cerr << "[ERROR] " << programName << ": " << fileDescription << ": " << filePath << " not exist.\n\n";
        return false;
    }
    return true;
}

/**
 * @brief Parse command line arguments using getopt_long
 * @param argc Number of command line arguments
 * @param argv Array of command line argument strings
 * 
 * This function processes command line arguments, validates files and parameters,
 * and handles help requests. It calls virtual functions that should be overridden
 * by derived classes for specific parameter handling.
 */
void ArgumentManager::parseOptions(int argc, char** argv)
{

    if(!optionDefiner){
        std::cerr << "[ERROR] " << programName << ": option definer not set." << std::endl;
        exit(EXIT_FAILURE);
    }

    // Initialize default values
    initializeDefaultValues();

    optind = 1;    // Reset getopt

    bool die = false;
    for (char c; (c = getopt_long(argc, argv, getShortOpts(), getLongOpts(), NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");

        if(loadArgument(c, arg)){
            continue;
        }
        
        if(c == getHelpEnumNum()){
            std::cout << HELP_MESSAGE << std::endl;
            exit(EXIT_SUCCESS);
        }else{
            die = true;
        }
    }

    // Build command string
    recordCommand(argc, argv);

    // Validate arguments
    if (argc - optind < 0) {
        std::cerr << "[ERROR] " << programName << ": missing arguments\n";
        die = true;
    }
    
    // Validate all input files
    if (!validateFiles()) {
        die = true;
    }
    
    // Validate numeric parameters
    if (!validateNumericParams()) {
        die = true;
    }
    
    if (die)
    {
        std::cerr << "\n";
        std::cout << HELP_MESSAGE << std::endl;
        exit(EXIT_FAILURE);
    } 
}

/**
 * @brief Namespace containing file validation utilities
 * 
 * Provides standalone functions for file validation that can be used
 * without an ArgumentManager instance
 */
namespace FileValidator{
    /**
     * @brief Validate that a required file exists and is accessible
     * @param filePath Path to the file to validate
     * @param fileDescription Description of the file for error messages
     * @param programName Name of the program for error messages
     * @return true if file exists and is accessible, false otherwise
     */
    bool validateRequiredFile(const std::string& filePath, const std::string& fileDescription, const std::string& programName) {
        if(filePath.empty()) {
            std::cerr << "[ERROR] " << programName  << ": missing " << fileDescription << ".\n";

            return false;
        }
        
        std::ifstream openFile(filePath.c_str());
        if(!openFile.is_open()) {
            std::cerr << "[ERROR] " << programName  << ": " << fileDescription << ": " << filePath << " not exist.\n\n";
            return false;
        }
        return true;
    }

    /**
     * @brief Validate that an optional file exists and is accessible (if specified)
     * @param filePath Path to the file to validate
     * @param fileDescription Description of the file for error messages
     * @param programName Name of the program for error messages
     * @return true if file is not specified or exists and is accessible, false otherwise
     */
    bool validateOptionalFile(const std::string& filePath, const std::string& fileDescription, const std::string& programName) {
        if(filePath.empty()) {
            // Optional file not specified, that's OK
            return true;  
        }
        
        std::ifstream openFile(filePath.c_str());
        if(!openFile.is_open()) {
            std::cerr << "[ERROR] " << programName << ": " << fileDescription << ": " << filePath << " not exist.\n\n";
            return false;
        }
        return true;
    }
}
