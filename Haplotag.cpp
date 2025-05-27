#include "Haplotag.h"
#include "HaplotagProcess.h"
#include "Util.h"
#include <getopt.h>

#define SUBPROGRAM "haplotag"

void HaplotagHelpManager::buildBaseMessage() {
    // Usage
    addSection("Usage: " + programName + " [OPTION] ... READSFILE");
    addItem("      --help                          display this help and exit.\n");
    
    // Required arguments - General mode
    addSection("required arguments:");
    addItem("      -s, --snp-file=NAME             input SNP vcf file.");
    addItem("      -b, --bam-file=NAME             input bam file.");
    addItem("      -r, --reference=NAME            reference FASTA.\n");
        
    // Optional arguments
    addSection("optional arguments:");
    addItem("      --tagSupplementary              tag supplementary alignment. default:false");
    addItem("      --sv-file=NAME                  input phased SV vcf file.");
    addItem("      --mod-file=NAME                 input a modified VCF file (produced by longphase modcall and processed by longphase phase).");
    addItem("      -q, --qualityThreshold=Num      not tag alignment if the mapping quality less than threshold. default:1");
    addItem("      -p, --percentageThreshold=Num   the alignment will be tagged according to the haplotype corresponding to most alleles.");
    addItem("                                      if the alignment has no obvious corresponding haplotype, it will not be tagged. default:0.6");
    addItem("      -t, --threads=Num               number of thread. default:1");
    addItem("      -o, --out-prefix=NAME           prefix of phasing result. default:result");
    addItem("      --cram                          the output file will be in the cram format. default:bam");
    addItem("      --region=REGION                 tagging include only reads/variants overlapping those regions. default:\"\"(all regions)");
    addItem("                                      input format:chrom (consider entire chromosome)");
    addItem("                                                   chrom:start (consider region from this start to end of chromosome)");
    addItem("                                                   chrom:start-end");
    addItem("      --log                           an additional log file records the result of each read. default:false");
}

HaplotagOptionManager::HaplotagOptionManager(const std::string& program) : OptionManager(program) {
}

void HaplotagOptionManager::setHelpMessage() {
    helpManager = createHelpManager(programName);
    helpManager->modifyMessage();
}

HaplotagOptionManager::~HaplotagOptionManager() {
   if(helpManager) delete helpManager;
}


void HaplotagOptionManager::setOptions() {
    // Initialize short options string
    shortOpts = "s:b:o:t:q:p:r:";

    // Help option
    addOption({"help", no_argument, NULL, OPT_HELP});
    
    // Input/Output files
    addOption({"snp-file", required_argument, NULL, 's'});
    addOption({"bam-file", required_argument, NULL, 'b'});
    addOption({"reference", required_argument, NULL, 'r'});
    addOption({"sv-file", required_argument, NULL, SV_FILE});
    addOption({"mod-file", required_argument, NULL, MOD_FILE});
        
    // Processing options
    addOption({"threads", required_argument, NULL, 't'});
    addOption({"qualityThreshold", required_argument, NULL, 'q'});
    addOption({"percentageThreshold", required_argument, NULL, 'p'});
    addOption({"tagSupplementary", no_argument, NULL, TAG_SUP});
    addOption({"out-prefix", required_argument, NULL, 'o'});
    addOption({"region", required_argument, NULL, REGION});
    addOption({"cram", no_argument, NULL, CRAM});
    addOption({"log", no_argument, NULL, LOG});
    
    // Add terminator
    addOption({NULL, 0, NULL, 0});

    //extend haplotag options
    extendOptions();
}

void HaplotagOptionManager::initializeDefaultValues() {
    // Initialize default values
    ecParams.numThreads = 1;
    ecParams.qualityThreshold = 1;
    ecParams.percentageThreshold = 0.6;
    ecParams.resultPrefix = "result";
    ecParams.outputFormat = "bam";
    ecParams.tagSupplementary = false;
    ecParams.writeReadLog = false;
    ecParams.command = "longphase ";
}

void HaplotagOptionManager::parseOptions(int argc, char** argv)
{

    // Initialize default values
    initializeDefaultValues();

    optind = 1;    // Reset getopt

    bool die = false;
    for (char c; (c = getopt_long(argc, argv, getShortOpts(), getLongOpts(), NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");

        if(loadHaplotagOptions(c, arg)){
            continue;
        }

        if(loadExtendOptions(c, arg)){
            continue;
        }
        
        if(c == HaplotagOption::OPT_HELP){
            helpManager->printHelp();
            exit(EXIT_SUCCESS);
        }else{
            die = true;
        }
    }

    // Build command string
    for(int i = 0; i < argc; ++i){
        ecParams.command.append(argv[i]);
        ecParams.command.append(" ");
    }

    // Validate arguments
    if (argc - optind < 0) {
        std::cerr << "[ERROR] " << programName << ": missing arguments\n";
        die = true;
    }
    
    // Validate all input files
    if (!validateFiles()) {
        die = true;
    }

    if (!validateExtendFiles()) {
        die = true;
    }
    
    // Validate numeric parameters
    if (!validateNumericParameter()) {
        die = true;
    }
    
    if (die)
    {
        std::cerr << "\n";
        helpManager->printHelp();
        exit(EXIT_FAILURE);
    } 
}

bool HaplotagOptionManager::loadHaplotagOptions(char& opt, std::istringstream& arg) {
    bool isLoaded = true;
    switch (opt)
    {
        case 's': arg >> ecParams.snpFile; break;
        case 't': arg >> ecParams.numThreads; break;
        case 'b': arg >> ecParams.bamFile; break;
        case 'r': arg >> ecParams.fastaFile; break; 
        case 'o': arg >> ecParams.resultPrefix; break;
        case 'q': arg >> ecParams.qualityThreshold; break;
        case 'p': arg >> ecParams.percentageThreshold; break;
        case HaplotagOption::SV_FILE:  arg >> ecParams.svFile; break;
        case HaplotagOption::MOD_FILE: arg >> ecParams.modFile; break;     
        case HaplotagOption::TAG_SUP:  ecParams.tagSupplementary = true; break;
        case HaplotagOption::REGION:   arg >> ecParams.region; break;        
        case HaplotagOption::CRAM:     ecParams.outputFormat = "cram"; break;
        case HaplotagOption::LOG:      ecParams.writeReadLog = true; break;
        default: isLoaded = false; break;
    }
    return isLoaded;
}


bool HaplotagOptionManager::validateFiles() {
    bool isValid = true;
    
    // Required files
    isValid &= validateRequiredFile(ecParams.snpFile, "SNP file");
    isValid &= validateRequiredFile(ecParams.bamFile, "BAM file");
    isValid &= validateRequiredFile(ecParams.fastaFile, "reference file");
    
    // Optional files
    isValid &= validateOptionalFile(ecParams.svFile, "SV file");
    isValid &= validateOptionalFile(ecParams.modFile, "MOD file");
    
    return isValid;
}

bool HaplotagOptionManager::loadExtendOptions(char& opt, std::istringstream& arg) {
    return false;
}

bool HaplotagOptionManager::validateExtendFiles() {
    return true;
}

bool HaplotagOptionManager::validateNumericParameter() {
    bool isValid = true;
    
    if (ecParams.numThreads < 1) {
        std::cerr << "[ERROR] " << programName << ": invalid threads. value: " 
                << ecParams.numThreads 
                << "\nplease check -t, --threads=Num\n";
        isValid = false;
    }
    
    if (ecParams.percentageThreshold > 1 || ecParams.percentageThreshold < 0) {
        std::cerr << "[ERROR] " << programName << ": invalid percentage threshold. value: " 
                << ecParams.percentageThreshold
                << "\nthis value need: 0~1, please check -p, --percentageThreshold=Num\n";
        isValid = false;
    }

    return isValid;
}


// Validate if a required file exists
bool OptionManager::validateRequiredFile(const std::string& filePath, const std::string& fileDescription) {
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
bool OptionManager::validateOptionalFile(const std::string& filePath, const std::string& fileDescription) {
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


int HaplotagMain(int argc, char** argv, std::string in_version)
{
    HaplotagOptionManager optionManager(SUBPROGRAM);

    optionManager.setOptions();
    optionManager.setHelpMessage();

    optionManager.parseOptions(argc, argv);
    optionManager.setVersion(in_version);

    HaplotagParameters ecParams = optionManager.getParams();
    
    HaplotagProcess processor(ecParams);
    processor.taggingProcess();

    return 0;
}
