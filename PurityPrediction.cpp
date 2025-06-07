#include "PurityPrediction.h"

#define SUBPROGRAM "purity prediction"

constexpr const char *CORRECT_USAGE_MESSAGE =
"Usage: "  " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"      --help                          display this help and exit.\n\n"
"required arguments:\n"
"      -s, --snp-file=NAME             input normal sample SNP VCF file.\n"
"      -b, --bam-file=NAME             input normal sample BAM file.\n"
"      --tumor-snp-file=NAME           input tumor sample SNP VCF file.\n"
"      --tumor-bam-file=NAME           input tumor sample BAM file.\n"
"      -r, --reference=NAME            reference FASTA.\n\n"
"optional arguments:\n"
"      --tagSupplementary              include supplementary alignments in haplotype statistics. default:true\n"
"      -q, --qualityThreshold=Num      exclude reads with mapping quality below threshold. default:20\n"
"      -p, --percentageThreshold=Num   include alignments in statistics based on haplotype allele percentage.\n"
"                                      alignments without clear haplotype assignment are excluded. default:0.6\n"
"      -t, --threads=Num               number of thread. default:1\n"
"      -o, --out-prefix=NAME           prefix of tumor purity prediction result. default:result\n"
"      --region=REGION                 tumor purity prediction include only reads/variants overlapping those regions. default:\"\"(all regions)\n"
"                                      input format:chrom (consider entire chromosome)\n"
"                                                   chrom:start (consider region from this start to end of chromosome)\n"
"                                                   chrom:start-end\n";

void PurityPredictionOptionDefiner::defineOptions(ArgumentManager& manager) {
    // base haplotag options
    HaplotagOptionDefiner::defineOptions(manager);

    // somatic haplotag-specific options
    manager.addOption({"tumor-snp-file", required_argument, NULL, TUM_SNP});
    manager.addOption({"tumor-bam-file", required_argument, NULL, TUM_BAM});
}

void ParamsHandler<PurityPredictionParameters>::initialize(PurityPredictionParameters& params, const std::string& version) {

    ParamsHandler<HaplotagParameters>::initialize(params.basic, version);
    params.basic.config.qualityThreshold = 20;
    params.basic.config.tagSupplementary = true;
}

bool ParamsHandler<PurityPredictionParameters>::loadArgument(PurityPredictionParameters& params, char& opt, std::istringstream& arg) {
    // load base haplotag options
    bool isLoaded = ParamsHandler<HaplotagParameters>::loadArgument(params.basic, opt, arg);
    
    if(!isLoaded){
        //reset isLoaded
        isLoaded = true;
        //load somatic haplotag options
        switch (opt)
        {
            case SomaticHaplotagOption::TUM_SNP: arg >> params.tumorSnpFile; break;
            case SomaticHaplotagOption::TUM_BAM: arg >> params.tumorBamFile; break;
            default: isLoaded = false; 
            break;
        }
    }
    return isLoaded;
}

bool ParamsHandler<PurityPredictionParameters>::validateFiles(PurityPredictionParameters& params, const std::string& programName) {
    // validate base haplotag files
    bool isValid = ParamsHandler<HaplotagParameters>::validateFiles(params.basic, programName);
    // validate somatic haplotag files
    isValid &= FileValidator::validateRequiredFile(params.tumorSnpFile, "tumor SNP file", programName);
    isValid &= FileValidator::validateRequiredFile(params.tumorBamFile, "tumor BAM file", programName);
    return isValid;
}

bool ParamsHandler<PurityPredictionParameters>::validateNumericParameter(PurityPredictionParameters& params, const std::string& programName) {
    bool isValid = ParamsHandler<HaplotagParameters>::validateNumericParameter(params.basic, programName);
    
    return isValid;  
}

void ParamsHandler<PurityPredictionParameters>::recordCommand(PurityPredictionParameters& params, int argc, char** argv) {
    for(int i = 0; i < argc; ++i){
        params.basic.config.command.append(argv[i]);
        params.basic.config.command.append(" ");
    }
}

int ParamsHandler<PurityPredictionParameters>::getHelpEnumNum(){
    return HaplotagOption::OPT_HELP;
}

int PurityPredictionMain(int argc, char** argv, std::string in_version){

    PurityPredictionArgumentManager optionManager(SUBPROGRAM, in_version, CORRECT_USAGE_MESSAGE);

    optionManager.setOptions();

    optionManager.parseOptions(argc, argv);

    PurityPredictionParameters ecParams = optionManager.getParams();

    PurityPredictionProcess processor(ecParams);

    processor.predictPurity();

    return 0;
}