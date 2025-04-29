#include "Haplotag.h"
#include "HaplotagProcess.h"
#include "SomaticHaplotagProcess.h"
#include "Util.h"
#include <getopt.h>

void HaplotagHelpManager::buildBaseMessage() {
    // Usage
    addSection("Usage: " + programName + " [OPTION] ... READSFILE");
    addItem("      --help                          display this help and exit.\n");
    
    // Required arguments - General mode
    addSection("required arguments:");
    addItem("   [General mode]");
    addItem("      -s, --snp-file=NAME             input SNP vcf file.");
    addItem("      -b, --bam-file=NAME             input bam file.");
    addItem("      -r, --reference=NAME            reference FASTA.\n");
    
    // Required arguments - Somatic mode
    addItem("   [Somatic mode] (tumor/normal pair data):");
    addItem("      --somaticMode                   enable somatic mutation tagging. default: false (disabled)");
    addItem("      -s, --snp-file=NAME             input normal sample SNP VCF file.");
    addItem("      -b, --bam-file=NAME             input normal sample BAM file (used as a reference for comparison).");
    addItem("      --tumor-snp-file=NAME           input tumor sample SNP VCF file.");
    addItem("      --tumor-bam-file=NAME           input tumor sample BAM file for mutation tagging.");
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

void HaplotagOptionManager::setOptions() {
    // Initialize short options string
    shortOpts = "s:b:o:t:q:p:r:";

    // Help option
    addOption({"help", no_argument, NULL, OPT_HELP});
    
    // Input/Output files
    addOption({"snp-file", required_argument, NULL, 's'});
    addOption({"bam-file", required_argument, NULL, 'b'});
    addOption({"reference", required_argument, NULL, 'r'});
    addOption({"out-prefix", required_argument, NULL, 'o'});
    
    // Tumor-specific options
    addOption({"tumor-snp-file", required_argument, NULL, TUM_SNP});
    addOption({"tumor-bam-file", required_argument, NULL, TUM_BAM});
    addOption({"somaticMode", no_argument, NULL, TAG_TUM});
    
    // Processing options
    addOption({"tagSupplementary", no_argument, NULL, TAG_SUP});
    addOption({"disableFilter", no_argument, NULL, DISABLE_FILTER});
    addOption({"sv-file", required_argument, NULL, SV_FILE});
    addOption({"mod-file", required_argument, NULL, MOD_FILE});
    addOption({"cram", no_argument, NULL, CRAM});
    
    // Configuration options
    addOption({"threads", required_argument, NULL, 't'});
    addOption({"qualityThreshold", required_argument, NULL, 'q'});
    addOption({"percentageThreshold", required_argument, NULL, 'p'});
    addOption({"somaticCallingMPQ", required_argument, NULL, SC_MPQ});
    addOption({"region", required_argument, NULL, REGION});
    addOption({"log", no_argument, NULL, LOG});
    addOption({"highCon-snp", required_argument, NULL, HIGH_CON});
    
    // Add terminator
    addOption({NULL, 0, NULL, 0});
}


void HaplotagOptionManager::parseOptions(int argc, char** argv)
{
    static HaplotagHelpManager helpManager;

    // Initialize default values
    ecParams.numThreads = 1;
    ecParams.qualityThreshold = 1;
    ecParams.somaticCallingMpqThreshold = 1;
    ecParams.percentageThreshold = 0.6;
    ecParams.resultPrefix = "result";
    ecParams.outputFormat = "bam";
    ecParams.tagSupplementary = false;
    ecParams.writeReadLog = false;
    ecParams.tagTumorSnp = false;
    ecParams.command = "longphase ";
    ecParams.enableFilter = true;

    optind = 1;    // Reset getopt
    
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, getShortOpts(), getLongOpts(), NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
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
            case HaplotagOption::TUM_SNP: arg >> ecParams.tumorSnpFile; break;
            case HaplotagOption::TUM_BAM: arg >> ecParams.tumorBamFile; break;
            case HaplotagOption::REGION:   arg >> ecParams.region; break;        
            case HaplotagOption::TAG_TUM: ecParams.tagTumorSnp = true; break;
            case HaplotagOption::TAG_SUP:  ecParams.tagSupplementary = true; break;
            case HaplotagOption::SC_MPQ: arg >> ecParams.somaticCallingMpqThreshold; break;
            case HaplotagOption::HIGH_CON: arg >> ecParams.benchmarkVcf; break;
            case HaplotagOption::DISABLE_FILTER: ecParams.enableFilter = false; break;
            case HaplotagOption::CRAM:     ecParams.outputFormat = "cram"; break;
            case HaplotagOption::LOG:      ecParams.writeReadLog = true; break;
            case HaplotagOption::OPT_HELP:
                helpManager.printHelp();
                exit(EXIT_SUCCESS);
            default: die = true; break;
        }
    }

    ecParams.somaticCallingMpqThreshold = ecParams.qualityThreshold;

    // Build command string
    for(int i = 0; i < argc; ++i){
        ecParams.command.append(argv[i]);
        ecParams.command.append(" ");
    }

    // Validate arguments
    if (argc - optind < 0) {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    }
    
    // Validate input files
    if(ecParams.snpFile != "") {
        std::ifstream openFile(ecParams.snpFile.c_str());
        if(!openFile.is_open()) {
            std::cerr << "File " << ecParams.snpFile << " not exist.\n\n";
            die = true;
        }
    } else {
        std::cerr << SUBPROGRAM ": missing SNP file.\n";
        die = true;
    }
    
    if(ecParams.svFile != "")
    {
        std::ifstream openFile(ecParams.svFile.c_str());
        if(!openFile.is_open())
        {
            std::cerr<< "File " << ecParams.svFile << " not exist.\n\n";
            die = true;
        }
    }
    
    if(ecParams.modFile != "")
    {
        std::ifstream openFile(ecParams.modFile.c_str());
        if(!openFile.is_open())
        {
            std::cerr<< "File " << ecParams.modFile << " not exist.\n\n";
            die = true;
        }
    }
    
    if(ecParams.bamFile != "")
    {
        std::ifstream openFile(ecParams.bamFile.c_str());
        if(!openFile.is_open())
        {
            std::cerr<< "File " << ecParams.bamFile << " not exist.\n\n";
            die = true;
        }
    }
    else{
        std::cerr << SUBPROGRAM ": missing bam file.\n";
        die = true;
    }
    
    if(ecParams.fastaFile != "")
    {
        std::ifstream openFile(ecParams.snpFile.c_str());
        if(!openFile.is_open())
        {
            std::cerr<< "File " << ecParams.fastaFile << " not exist.\n\n";
            die = true;
        }
    }
    else{
        std::cerr << SUBPROGRAM ": missing reference.\n";
        die = true;
    }

    //somatic mode  
    if(ecParams.tagTumorSnp == true){
        // tumor vcf
        if(ecParams.tumorSnpFile != "")
        {
            std::ifstream openFile(ecParams.tumorSnpFile.c_str());
            if(!openFile.is_open())
            {
                std::cerr<< "File " << ecParams.tumorSnpFile << " not exist.\n\n";
                die = true;
            }
        }else{
            std::cerr << SUBPROGRAM ": missing tumor genome SNP file.\n";
            die = true;
        }

        //tumor bam
        if(ecParams.tumorBamFile != "")
        {
            std::ifstream openFile(ecParams.tumorBamFile.c_str());
            if(!openFile.is_open())
            {
                std::cerr<< "File " << ecParams.tumorBamFile << " not exist.\n\n";
                die = true;
            }
        }else{
            std::cerr << SUBPROGRAM ": missing tumor bam file.\n";
            die = true;
        }

        //benchmarking
        if(ecParams.benchmarkVcf != ""){
            std::ifstream openFile(ecParams.benchmarkVcf.c_str());
            if(!openFile.is_open())
            {
                std::cerr<< "File " << ecParams.benchmarkVcf << " not exist.\n\n";
                die = true;
            }
        }
    }  
    
    if (ecParams.numThreads < 1){
        std::cerr << SUBPROGRAM " invalid threads. value: " 
                  << ecParams.numThreads 
                  << "\nplease check -t, --threads=Num\n";
        die = true;
    }
    
    if (ecParams.percentageThreshold > 1 || ecParams.percentageThreshold < 0){
        std::cerr << SUBPROGRAM " invalid percentage threshold. value: " 
                  << ecParams.percentageThreshold
                  << "\nthis value need: 0~1, please check -p, --percentageThreshold=Num\n";
        die = true;
    }
    
    if (die)
    {
        std::cerr << "\n";
        helpManager.printHelp();
        exit(EXIT_FAILURE);
    }

    return;
}

int HaplotagMain(int argc, char** argv, std::string in_version)
{
    HaplotagOptionManager optionManager;

    optionManager.parseOptions(argc, argv);
    optionManager.setVersion(in_version);

    HaplotagParameters ecParams = optionManager.getParams();
    
    if(ecParams.tagTumorSnp){
        SomaticHaplotagProcess processor(ecParams);
        processor.taggingProcess();
    }else{
        HaplotagProcess processor(ecParams);
        processor.taggingProcess();
    }

    return 0;
}
