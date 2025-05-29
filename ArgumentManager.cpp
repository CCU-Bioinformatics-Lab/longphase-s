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