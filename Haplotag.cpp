#include "Haplotag.h"
#include "HaplotagProcess.h"
#include "Util.h"
#include <getopt.h>


#define SUBPROGRAM "haplotag"

static const char *CORRECT_USAGE_MESSAGE =
"Usage: "  " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"      --help                          display this help and exit.\n\n"
"required arguments:\n"
"   [General mode]\n"
"      -s, --snp-file=NAME             input SNP vcf file.\n"
"      -b, --bam-file=NAME             input bam file.\n"
"      -r, --reference=NAME            reference FASTA.\n\n"
//"required arguments for somatic mode (tumor/normal pair data):\n"
"   [Somatic mode] (tumor/normal pair data):\n"
"      --somaticMode                   enable somatic mutation tagging. default: false (disabled)\n"
"      -s, --snp-file=NAME             input normal sample SNP VCF file.\n"
"      -b, --bam-file=NAME             input normal sample BAM file (used as a reference for comparison).\n"
"      --tumor-snp-file=NAME           input tumor sample SNP VCF file.\n"
"      --tumor-bam-file=NAME           input tumor sample BAM file for mutation tagging.\n"
"      -r, --reference=NAME            reference FASTA.\n\n"
//"      --somaticCallingMPQ=Num         mapping quality threshold for calling somatic SNPs. default: 40\n\n"
"optional arguments:\n"
"      --tagSupplementary              tag supplementary alignment. default:false\n"
"      --sv-file=NAME                  input phased SV vcf file.\n"
"      --mod-file=NAME                 input a modified VCF file (produced by longphase modcall and processed by longphase phase).\n"
"      -q, --qualityThreshold=Num      not tag alignment if the mapping quality less than threshold. default:1\n"
"      -p, --percentageThreshold=Num   the alignment will be tagged according to the haplotype corresponding to most alleles.\n"
"                                      if the alignment has no obvious corresponding haplotype, it will not be tagged. default:0.6\n"
"      -t, --threads=Num               number of thread. default:1\n"
"      -o, --out-prefix=NAME           prefix of phasing result. default:result\n"
"      --cram                          the output file will be in the cram format. default:bam\n"
"      --region=REGION                 tagging include only reads/variants overlapping those regions. default:\"\"(all regions)\n"
"                                      input format:chrom (consider entire chromosome)\n"
"                                                   chrom:start (consider region from this start to end of chromosome)\n"
"                                                   chrom:start-end\n"
"      --log                           an additional log file records the result of each read. default:false\n";


static const char* shortopts = "s:b:o:t:q:p:r:";

enum { OPT_HELP = 1, TAG_SUP, SV_FILE, REGION, LOG, MOD_FILE, CRAM, TUM_SNP, TUM_BAM, SC_MPQ, TAG_TUM, HIGH_CON, DISABLE_FILTER};

static const struct option longopts[] = { 
    { "help",                 no_argument,        NULL, OPT_HELP },
    { "snp-file",             required_argument,  NULL, 's' },
    { "bam-file",             required_argument,  NULL, 'b' },
    { "tumor-snp-file",       required_argument,  NULL, TUM_SNP },   //new
    { "tumor-bam-file",       required_argument,  NULL, TUM_BAM },   //new
    { "reference",            required_argument,  NULL, 'r' },
    { "tagSupplementary",     no_argument,        NULL, TAG_SUP },
    { "somaticMode",          no_argument,        NULL, TAG_TUM },  //new
    { "disableFilter",        no_argument,        NULL, DISABLE_FILTER },  //new
    { "sv-file",              required_argument,  NULL, SV_FILE },
    { "mod-file",             required_argument,  NULL, MOD_FILE },
    { "out-prefix",           required_argument,  NULL, 'o' },
    { "cram",                 no_argument,        NULL, CRAM },
    { "threads",              required_argument,  NULL, 't' },
    { "qualityThreshold",     required_argument,  NULL, 'q' },
    { "percentageThreshold",  required_argument,  NULL, 'p' },
    { "somaticCallingMPQ",    required_argument,  NULL, SC_MPQ },   //new
    { "region",               required_argument,  NULL, REGION },
    { "log",                  no_argument,        NULL, LOG },
    { "highCon-snp",          required_argument,  NULL, HIGH_CON},
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static int numThreads = 1;
    static int qualityThreshold = 1;
    //static int somaticCallingMpqThreshold = 1;  //new
    static double percentageThreshold = 0.6;
    static std::string snpFile="";
    static std::string svFile="";
    static std::string modFile="";
    static std::string bamFile="";
    static std::string tumorSnpFile="";  //new
    static std::string tumorBamFile="";  //new
    static std::string fastaFile="";
    static std::string resultPrefix="result";
    static std::string region="";
    static std::string outputFormat="bam";
    static bool tagSupplementary = false;
    static bool writeReadLog = false;
    static bool tumorMode = false;  //new
    static std::string command="longphase ";
    static std::string highConSnp="";
    static bool enableFilter = true;
}

void HaplotagOptions(int argc, char** argv)
{
    optind=1;    //reset getopt
    
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 's': arg >> opt::snpFile; break;
            case 't': arg >> opt::numThreads; break;
            case 'b': arg >> opt::bamFile; break;
            case 'r': arg >> opt::fastaFile; break; 
            case 'o': arg >> opt::resultPrefix; break;
            case 'q': arg >> opt::qualityThreshold; break;
            case 'p': arg >> opt::percentageThreshold; break;
            case SV_FILE:  arg >> opt::svFile; break;
            case MOD_FILE: arg >> opt::modFile; break;     
            case TUM_SNP: arg >> opt::tumorSnpFile; break;  //new
            case TUM_BAM: arg >> opt::tumorBamFile; break;  //new
            case REGION:   arg >> opt::region; break;        
            case TAG_TUM: opt::tumorMode = true; break;  //new
            case TAG_SUP:  opt::tagSupplementary = true; break;
            case SC_MPQ: arg >> opt::qualityThreshold; break; //new
            case HIGH_CON: arg >> opt::highConSnp; break;
            case DISABLE_FILTER: opt::enableFilter = false; break;
            case CRAM:     opt::outputFormat = "cram"; break;
            case LOG:      opt::writeReadLog = true; break;
            case OPT_HELP:
                std::cout << CORRECT_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            default: die = true; break;
        }
    }

    for(int i = 0; i < argc; ++i){
        opt::command.append(argv[i]);
        opt::command.append(" ");
    }

    if (argc - optind < 0 )
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    }
    
    if( opt::snpFile != "")
    {
        std::ifstream openFile( opt::snpFile.c_str() );
        if( !openFile.is_open() )
        {
            std::cerr<< "File " << opt::snpFile << " not exist.\n\n";
            die = true;
        }
    }
    else{
        std::cerr << SUBPROGRAM ": missing SNP file.\n";
        die = true;
    }
    
    if( opt::svFile != "")
    {
        std::ifstream openFile( opt::svFile.c_str() );
        if( !openFile.is_open() )
        {
            std::cerr<< "File " << opt::svFile << " not exist.\n\n";
            die = true;
        }
    }
    
    if( opt::modFile != "")
    {
        std::ifstream openFile( opt::modFile.c_str() );
        if( !openFile.is_open() )
        {
            std::cerr<< "File " << opt::modFile << " not exist.\n\n";
            die = true;
        }
    }
    
    if( opt::bamFile != "")
    {
        std::ifstream openFile( opt::bamFile.c_str() );
        if( !openFile.is_open() )
        {
            std::cerr<< "File " << opt::bamFile << " not exist.\n\n";
            die = true;
        }
    }
    else{
        std::cerr << SUBPROGRAM ": missing bam file.\n";
        die = true;
    }
    
    if( opt::fastaFile != "")
    {
        std::ifstream openFile( opt::snpFile.c_str() );
        if( !openFile.is_open() )
        {
            std::cerr<< "File " << opt::fastaFile << " not exist.\n\n";
            die = true;
        }
    }
    else{
        std::cerr << SUBPROGRAM ": missing reference.\n";
        die = true;
    }

    //somatic mode  
    if(opt::tumorMode == true){
        // tumor vcf
        if( opt::tumorSnpFile != "")
        {
            std::ifstream openFile( opt::tumorSnpFile.c_str() );
            if( !openFile.is_open() )
            {
                std::cerr<< "File " << opt::tumorSnpFile << " not exist.\n\n";
                die = true;
            }
        }else{
            std::cerr << SUBPROGRAM ": missing tumor genome SNP file.\n";
            die = true;
        }

        //tumor bam
        if( opt::tumorBamFile != "")
        {
            std::ifstream openFile( opt::tumorBamFile.c_str() );
            if( !openFile.is_open() )
            {
                std::cerr<< "File " << opt::tumorBamFile << " not exist.\n\n";
                die = true;
            }
        }else{
            std::cerr << SUBPROGRAM ": missing tumor bam file.\n";
            die = true;
        }

        //benchmarking
        if( opt::highConSnp != ""){
            std::ifstream openFile( opt::highConSnp.c_str() );
            if( !openFile.is_open() )
            {
                std::cerr<< "File " << opt::highConSnp << " not exist.\n\n";
                die = true;
            }
        }
    }  
    
    if ( opt::numThreads < 1 ){
        std::cerr << SUBPROGRAM " invalid threads. value: " 
                  << opt::numThreads 
                  << "\nplease check -t, --threads=Num\n";
        die = true;
    }
    
    if ( opt::percentageThreshold > 1 || opt::percentageThreshold < 0 ){
        std::cerr << SUBPROGRAM " invalid percentage threshold. value: " 
                  << opt::percentageThreshold
                  << "\nthis value need: 0~1, please check -p, --percentageThreshold=Num\n";
        die = true;
    }
    
    if (die)
    {
        std::cerr << "\n" << CORRECT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

}

int HaplotagMain(int argc, char** argv, std::string in_version)
{
    HaplotagParameters ecParams;
    // set parameters
    HaplotagOptions(argc, argv);

    ecParams.numThreads=opt::numThreads;
    ecParams.qualityThreshold=opt::qualityThreshold;
    ecParams.somaticCallingMpqThreshold=opt::qualityThreshold;  //new
    ecParams.snpFile=opt::snpFile;
    ecParams.svFile=opt::svFile;
    ecParams.modFile=opt::modFile;
    ecParams.bamFile=opt::bamFile;
    ecParams.fastaFile=opt::fastaFile;
    ecParams.tumorSnpFile=opt::tumorSnpFile;  //new
    ecParams.tumorBamFile=opt::tumorBamFile;  //new
    ecParams.tagTumorSnp=opt::tumorMode;    //new
    ecParams.resultPrefix=opt::resultPrefix;
    ecParams.tagSupplementary=opt::tagSupplementary;
    ecParams.percentageThreshold=opt::percentageThreshold;
    ecParams.region=opt::region;
    ecParams.writeReadLog=opt::writeReadLog;
    ecParams.version=in_version;
    ecParams.command=opt::command;
    ecParams.outputFormat=opt::outputFormat;
    ecParams.benchmarkVcf=opt::highConSnp;
    ecParams.enableFilter=opt::enableFilter;
    HaplotagProcess processor;
    processor.taggingProcess(ecParams);

    return 0;
}
