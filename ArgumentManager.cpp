#include "ArgumentManager.h"

void HelpMessageManager::addSection(const std::string& header) {
    sections.push_back({header, {}});
}

void HelpMessageManager::addItem(const std::string& item) {
    if (!sections.empty()) {
        sections.back().items.push_back(item);
    }
}

void HelpMessageManager::addSectionItem(const std::string& sectionName, const std::string& newItem) {
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

void HelpMessageManager::clearSectionItem(const std::string& sectionName) {
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

void HelpMessageManager::printHelp() const {
    for (const auto& section : sections) {
        std::cout << section.header << std::endl;
        for (const auto& item : section.items) {
            std::cout << item << std::endl;
        }
    }
}

// Get the complete help message string
std::string HelpMessageManager::getHelpMessage() const {
    std::ostringstream oss;
    for (const auto& section : sections) {
        oss << section.header << "\n";
        for (const auto& item : section.items) {
            oss << item << "\n";
        }
    }
    return oss.str();
}


void ArgumentManager::addOption(const struct option& opt) {
    // If this is not the terminator, add it before the last element
    if (opt.name != NULL && !longOpts.empty()) {
        longOpts.insert(longOpts.end() - 1, opt);
    } else {
        longOpts.push_back(opt);
    }
}

void ArgumentManager::removeOption(const std::string& name) {
    auto it = std::remove_if(longOpts.begin(), longOpts.end() - 1, [&](const option& opt) {
        return opt.name != nullptr && name == opt.name;
    });
    if (it != longOpts.end() - 1) {
        longOpts.erase(it, longOpts.end() - 1);
    }else{
        std::cerr << "[ERROR](ArgumentManager) removeOption: Option " << name << " not found." << std::endl;
        exit(EXIT_FAILURE);
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
            if(helpManager) helpManager->printHelp();
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
        if(helpManager) helpManager->printHelp();
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
