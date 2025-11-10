#include "PurityEstimation.h"

#define SUBPROGRAM "estimate_purity"

constexpr const char *CORRECT_USAGE_MESSAGE =
"Usage: "  " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"      --help                          display this help and exit.\n\n"
"required arguments:\n"
"      -s, --snp-file=NAME             input phased normal sample SNP VCF file.\n"
"      -b, --bam-file=NAME             input normal sample BAM file.\n"
"      --tumor-snv-file=NAME           input tumor sample SNV VCF file.\n"
"      --tumor-bam-file=NAME           input tumor sample BAM file.\n"
"      -r, --reference=NAME            reference FASTA.\n\n"
"optional arguments:\n"
"      --tagSupplementary              include supplementary alignments in haplotype statistics. default:true\n"
"      -q, --qualityThreshold=Num      exclude reads with mapping quality below threshold. default:20\n"
"      -p, --percentageThreshold=Num   include alignments in statistics based on haplotype allele percentage.\n"
"                                      alignments without clear haplotype assignment are excluded. default:0.6\n"
"      -t, --threads=Num               number of thread. default:1\n"
"      -o, --out-prefix=NAME           prefix of tumor purity estimation result. default:result\n"
"      --region=REGION                 tumor purity estimation include only reads/variants overlapping those regions. default:\"\"(all regions)\n"
"                                      input format:chrom (consider entire chromosome)\n"
"                                                   chrom:start (consider region from this start to end of chromosome)\n"
"                                                   chrom:start-end\n";

void PurityEstimOptDefiner::defineOptions(ArgumentManager& manager) {
    // base haplotag options
    HaplotagOptionDefiner::defineOptions(manager);

    // somatic haplotag-specific options
    manager.addOption({"tumor-snv-file", required_argument, NULL, TUM_SNP});
    manager.addOption({"tumor-bam-file", required_argument, NULL, TUM_BAM});
}

void ParamsHandler<PurityEstimParameters>::initialize(PurityEstimParameters& params, const std::string& version) {

    ParamsHandler<HaplotagParameters>::initialize(params.basic, version);
    params.basic.bamCfg.qualityThreshold = 20;
    params.basic.bamCfg.tagSupplementary = true;
}

bool ParamsHandler<PurityEstimParameters>::loadArgument(PurityEstimParameters& params, char& opt, std::istringstream& arg) {
    // load base haplotag options
    bool isLoaded = ParamsHandler<HaplotagParameters>::loadArgument(params.basic, opt, arg);
    
    if(!isLoaded){
        //reset isLoaded
        isLoaded = true;
        //load somatic haplotag options
        switch (opt)
        {
            case SomaticHaplotagOption::TUM_SNP: arg >> params.tumorSnvFile; break;
            case SomaticHaplotagOption::TUM_BAM: arg >> params.tumorBamFile; break;
            default: isLoaded = false; 
            break;
        }
    }
    return isLoaded;
}

bool ParamsHandler<PurityEstimParameters>::validateFiles(PurityEstimParameters& params, const std::string& programName) {
    // validate base haplotag files
    bool isValid = ParamsHandler<HaplotagParameters>::validateFiles(params.basic, programName);
    // validate somatic haplotag files
    isValid &= FileValidator::validateRequiredFile(params.tumorSnvFile, "tumor SNV file", programName);
    isValid &= FileValidator::validateRequiredFile(params.tumorBamFile, "tumor BAM file", programName);
    return isValid;
}

bool ParamsHandler<PurityEstimParameters>::validateNumericParams(PurityEstimParameters& params, const std::string& programName) {
    bool isValid = ParamsHandler<HaplotagParameters>::validateNumericParams(params.basic, programName);
    
    return isValid;  
}

void ParamsHandler<PurityEstimParameters>::recordCommand(PurityEstimParameters& params, int argc, char** argv) {
    ParamsHandler<ParsingBamConfig>::recordCommand(params.basic.bamCfg, argc, argv);
}

int ParamsHandler<PurityEstimParameters>::getHelpEnumNum(){
    return HaplotagOption::OPT_HELP;
}

int PurityEstimMain(int argc, char** argv, std::string in_version){

    PurityEstimArgsManager optionManager(SUBPROGRAM, in_version, CORRECT_USAGE_MESSAGE);

    optionManager.setOptions();

    optionManager.parseOptions(argc, argv);

    PurityEstimParameters ecParams = optionManager.getParams();

    PurityEstimProcess processor(ecParams);

    processor.estimatePurity();

    return 0;
}