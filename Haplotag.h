#ifndef HAPLOTAG_H
#define HAPLOTAG_H

#include "Util.h"
#include "HaplotagType.h"  // Include the header that defines HaplotagParameters
#include "ArgumentManager.h"
#include <getopt.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

template<>
struct ParamsHandler<HaplotagParameters>{
    static void initialize(HaplotagParameters& params) {
        // Initialize default values
        params.numThreads = 1;
        params.qualityThreshold = 1;
        params.percentageThreshold = 0.6;
        params.resultPrefix = "result";
        params.outputFormat = "bam";
        params.tagSupplementary = false;
        params.writeReadLog = false;
        params.command = "longphase ";
    }

    static bool loadArgument(HaplotagParameters& params, char& opt, std::istringstream& arg) {
        bool isLoaded = true;
        switch (opt)
        {
            case 's': arg >> params.snpFile; break;
            case 't': arg >> params.numThreads; break;
            case 'b': arg >> params.bamFile; break;
            case 'r': arg >> params.fastaFile; break; 
            case 'o': arg >> params.resultPrefix; break;
            case 'q': arg >> params.qualityThreshold; break;
            case 'p': arg >> params.percentageThreshold; break;
            case HaplotagOption::SV_FILE:  arg >> params.svFile; break;
            case HaplotagOption::MOD_FILE: arg >> params.modFile; break;     
            case HaplotagOption::TAG_SUP:  params.tagSupplementary = true; break;
            case HaplotagOption::REGION:   arg >> params.region; break;        
            case HaplotagOption::CRAM:     params.outputFormat = "cram"; break;
            case HaplotagOption::LOG:      params.writeReadLog = true; break;
            default: isLoaded = false; break;
        }
        return isLoaded;    
    }

    static bool validateFiles(HaplotagParameters& params, const std::string& programName) {
        bool isValid = true;
        
        // Required files
        isValid &= FileValidator::validateRequiredFile(params.snpFile, "SNP file", programName);
        isValid &= FileValidator::validateRequiredFile(params.bamFile, "BAM file", programName);
        isValid &= FileValidator::validateRequiredFile(params.fastaFile, "reference file", programName);
        
        // Optional files
        isValid &= FileValidator::validateOptionalFile(params.svFile, "SV file", programName);
        isValid &= FileValidator::validateOptionalFile(params.modFile, "MOD file", programName);
        
        return isValid;
    }

    static bool validateNumericParameter(HaplotagParameters& params, const std::string& programName) {
        bool isValid = true;
        
        if (params.numThreads < 1) {
            std::cerr << "[ERROR] " << programName << ": invalid threads. value: " 
                    << params.numThreads 
                    << "\nplease check -t, --threads=Num\n";
            isValid = false;
        }
        
        if (params.percentageThreshold > 1 || params.percentageThreshold < 0) {
            std::cerr << "[ERROR] " << programName << ": invalid percentage threshold. value: " 
                    << params.percentageThreshold
                    << "\nthis value need: 0~1, please check -p, --percentageThreshold=Num\n";
            isValid = false;
        }

        return isValid;   
    }

    static void recordCommand(HaplotagParameters& params, int argc, char** argv) {
        for(int i = 0; i < argc; ++i){
            params.command.append(argv[i]);
            params.command.append(" ");
        }
    }

};

struct HaplotagParamHandler{
    static void initialize(HaplotagParameters& params) {
        // Initialize default values
        params.numThreads = 1;
        params.qualityThreshold = 1;
        params.percentageThreshold = 0.6;
        params.resultPrefix = "result";
        params.outputFormat = "bam";
        params.tagSupplementary = false;
        params.writeReadLog = false;
        params.command = "longphase ";
    }

    static bool loadArgument(HaplotagParameters& params, char& opt, std::istringstream& arg) {
        bool isLoaded = true;
        switch (opt)
        {
            case 's': arg >> params.snpFile; break;
            case 't': arg >> params.numThreads; break;
            case 'b': arg >> params.bamFile; break;
            case 'r': arg >> params.fastaFile; break; 
            case 'o': arg >> params.resultPrefix; break;
            case 'q': arg >> params.qualityThreshold; break;
            case 'p': arg >> params.percentageThreshold; break;
            case HaplotagOption::SV_FILE:  arg >> params.svFile; break;
            case HaplotagOption::MOD_FILE: arg >> params.modFile; break;     
            case HaplotagOption::TAG_SUP:  params.tagSupplementary = true; break;
            case HaplotagOption::REGION:   arg >> params.region; break;        
            case HaplotagOption::CRAM:     params.outputFormat = "cram"; break;
            case HaplotagOption::LOG:      params.writeReadLog = true; break;
            default: isLoaded = false; break;
        }
        return isLoaded;    
    }

    static bool validateFiles(HaplotagParameters& params, const std::string& programName) {
        bool isValid = true;
        
        // Required files
        isValid &= FileValidator::validateRequiredFile(params.snpFile, "SNP file", programName);
        isValid &= FileValidator::validateRequiredFile(params.bamFile, "BAM file", programName);
        isValid &= FileValidator::validateRequiredFile(params.fastaFile, "reference file", programName);
        
        // Optional files
        isValid &= FileValidator::validateOptionalFile(params.svFile, "SV file", programName);
        isValid &= FileValidator::validateOptionalFile(params.modFile, "MOD file", programName);
        
        return isValid;
    }

    static bool validateNumericParameter(HaplotagParameters& params, const std::string& programName) {
        bool isValid = true;
        
        if (params.numThreads < 1) {
            std::cerr << "[ERROR] " << programName << ": invalid threads. value: " 
                    << params.numThreads 
                    << "\nplease check -t, --threads=Num\n";
            isValid = false;
        }
        
        if (params.percentageThreshold > 1 || params.percentageThreshold < 0) {
            std::cerr << "[ERROR] " << programName << ": invalid percentage threshold. value: " 
                    << params.percentageThreshold
                    << "\nthis value need: 0~1, please check -p, --percentageThreshold=Num\n";
            isValid = false;
        }

        return isValid;   
    }

    static void recordCommand(HaplotagParameters& params, int argc, char** argv) {
        for(int i = 0; i < argc; ++i){
            params.command.append(argv[i]);
            params.command.append(" ");
        }
    }

};

class HaplotagHelpManager : public HelpMessageManager {
public:
    HaplotagHelpManager(const std::string& program) : HelpMessageManager(program) {}

    virtual void buildMessage() override;
};

// Haplotag-specific option definition manager
class HaplotagArgumentManager : public ArgumentManager {
    protected:
        HaplotagParameters ecParams;

        HelpMessageManager* helpManager;

        virtual void initializeDefaultValues();

        virtual bool loadArgument(char& opt, std::istringstream& arg);

        virtual void recordCommand(int argc, char** argv);

        // Validate all input files
        virtual bool validateFiles();

        virtual bool validateNumericParameter();

        virtual HelpMessageManager* createHelpManager(const std::string& program) override {
            return new HaplotagHelpManager(program);
        }

        virtual int getHelpEnumNum() override {
            return HaplotagOption::OPT_HELP;
        }

    public:
        HaplotagArgumentManager(const std::string& program);
        virtual ~HaplotagArgumentManager();

        virtual void setOptions() override;
        // void setHelpMessage();
        // void parseOptions(int argc, char** argv);

        void setVersion(std::string in_version) { ecParams.version = in_version; }
        HaplotagParameters getParams() const { return ecParams; }
};

// functions
int HaplotagMain(int argc, char** argv, std::string in_version);

#endif