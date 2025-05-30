#include "Haplotag.h"
#include "HaplotagProcess.h"
#include "Util.h"
#include <getopt.h>

#define SUBPROGRAM "haplotag"

void HaplotagHelpManager::buildMessage() {
    // Usage
    addSection("Usage: " + programName + " [OPTION] ... READSFILE");
    addItem("      --help                          display this help and exit.\n");
    
    // Required arguments - General mode
    addSection(REQUIRED_SECTION);
    addItem("      -s, --snp-file=NAME             input SNP vcf file.");
    addItem("      -b, --bam-file=NAME             input bam file.");
    addItem("      -r, --reference=NAME            reference FASTA.\n");
        
    // Optional arguments
    addSection(OPTIONAL_SECTION);
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

void HaplotagOptionDefiner::defineOptions(ArgumentManager& manager) {
    // Initialize short options string
    manager.setShortOption("s:b:o:t:q:p:r:");

    // Help option
    manager.addOption({"help", no_argument, NULL, OPT_HELP});
    
    // Input/Output files
    manager.addOption({"snp-file", required_argument, NULL, 's'});
    manager.addOption({"bam-file", required_argument, NULL, 'b'});
    manager.addOption({"reference", required_argument, NULL, 'r'});
    manager.addOption({"sv-file", required_argument, NULL, SV_FILE});
    manager.addOption({"mod-file", required_argument, NULL, MOD_FILE});
        
    // Processing options
    manager.addOption({"threads", required_argument, NULL, 't'});
    manager.addOption({"qualityThreshold", required_argument, NULL, 'q'});
    manager.addOption({"percentageThreshold", required_argument, NULL, 'p'});
    manager.addOption({"tagSupplementary", no_argument, NULL, TAG_SUP});
    manager.addOption({"out-prefix", required_argument, NULL, 'o'});
    manager.addOption({"region", required_argument, NULL, REGION});
    manager.addOption({"cram", no_argument, NULL, CRAM});
    manager.addOption({"log", no_argument, NULL, LOG});
    
    // Add terminator
    manager.addOption({NULL, 0, NULL, 0});
}

HaplotagArgumentManager::HaplotagArgumentManager(const std::string& program, const std::string& version)
 : ArgumentTemManager<HaplotagParameters>(program, version) {
}


HaplotagArgumentManager::~HaplotagArgumentManager() {
}

void ParamsHandler<HaplotagParameters>::initialize(HaplotagParameters& params, const std::string& version) {
    // Initialize default values
    params.numThreads = 1;
    params.qualityThreshold = 1;
    params.percentageThreshold = 0.6;
    params.resultPrefix = "result";
    params.outputFormat = "bam";
    params.tagSupplementary = false;
    params.writeReadLog = false;
    params.command = "longphase-s ";
    params.version = version;
}

bool ParamsHandler<HaplotagParameters>::loadArgument(HaplotagParameters& params, char& opt, std::istringstream& arg) {
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

bool ParamsHandler<HaplotagParameters>::validateFiles(HaplotagParameters& params, const std::string& programName) {
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

bool ParamsHandler<HaplotagParameters>::validateNumericParameter(HaplotagParameters& params, const std::string& programName) {
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

void ParamsHandler<HaplotagParameters>::recordCommand(HaplotagParameters& params, int argc, char** argv) {
    for(int i = 0; i < argc; ++i){
        params.command.append(argv[i]);
        params.command.append(" ");
    }
}

int ParamsHandler<HaplotagParameters>::getHelpEnumNum() {
    return HaplotagOption::OPT_HELP;
}

int HaplotagMain(int argc, char** argv, std::string in_version)
{
    HaplotagArgumentManager optionManager(SUBPROGRAM, in_version);

    optionManager.setOptions();
    optionManager.setHelpMessage();

    optionManager.parseOptions(argc, argv);
    // optionManager.setVersion(in_version);

    HaplotagParameters ecParams = optionManager.getParams();
    
    HaplotagProcess processor(ecParams);
    processor.taggingProcess();

    return 0;
}
