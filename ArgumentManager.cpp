#include "ArgumentManager.h"


void ArgumentManager::setShortOption(const std::string& opt) {
    shortOpts = opt;
}

void ArgumentManager::addOption(const struct option& opt) {
    // If this is not the terminator, add it before the last element
    if (opt.name != NULL && !longOpts.empty()) {
        longOpts.insert(longOpts.end() - 1, opt);
    } else {
        longOpts.push_back(opt);
    }
}


// Validate if a required file exists
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

// Validate if an optional file exists (if specified)
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
    if (!validateNumericParameter()) {
        die = true;
    }
    
    if (die)
    {
        std::cerr << "\n";
        std::cout << HELP_MESSAGE << std::endl;
        exit(EXIT_FAILURE);
    } 
}

namespace FileValidator{
    // Validate if a required file exists
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

    // Validate if an optional file exists (if specified)
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
