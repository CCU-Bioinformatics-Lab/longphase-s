#include "SomaticHaplotag.h"

#define SUBPROGRAM "somatic_haplotag"

constexpr const char *CORRECT_USAGE_MESSAGE =
"Usage: "  " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"      --help                          display this help and exit.\n\n"
"required arguments:\n"
"      -s, --snp-file=NAME             input normal sample SNP VCF file.\n"
"      -b, --bam-file=NAME             input normal sample BAM file.\n"
"      --tumor-snp-file=NAME           input tumor sample SNP VCF file.\n"
"      --tumor-bam-file=NAME           input tumor sample BAM file for somatic haplotag.\n"
"      -r, --reference=NAME            reference FASTA.\n\n"
"optional arguments:\n"
"      --tagSupplementary              tag supplementary alignment. default:false\n"
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
"      --log                           an additional log file records the result of each read. default:false\n\n"
"somatic variant calling arguments:\n"
"      --tumor-purity=Num              tumor purity (0.1~1.0) for adjusting somatic variant filtering.\n"
"                                      lower values apply more conservative filters. default: auto-prediction.\n"
"      --disableFilter                 disable somatic variant filtering and accept all tumor VCF variants. default: false.\n\n"
"somatic haplotag benchmark arguments:\n"
"      --truth-vcf=NAME                truth somatic variants VCF file for somatic haplotag evaluation.\n"
"      --truth-bed=NAME                confident regions BED file for somatic haplotag evaluation.\n";

void SomaticHaplotagOptionDefiner::defineOptions(ArgumentManager& manager) {
    // base haplotag options
    HaplotagOptionDefiner::defineOptions(manager);

    // somatic haplotag-specific options
    manager.addOption({"tumor-snp-file", required_argument, NULL, TUM_SNP});
    manager.addOption({"tumor-bam-file", required_argument, NULL, TUM_BAM});
    
    manager.addOption({"disableFilter", no_argument, NULL, DISABLE_FILTER});
    manager.addOption({"tumor-purity", required_argument, NULL, TUMOR_PURITY});

    manager.addOption({"truth-vcf", required_argument, NULL, BENCHMARK_VCF});
    manager.addOption({"truth-bed", required_argument, NULL, BENCHMARK_BED});
}

void ParamsHandler<SomaticHaplotagParameters>::initialize(SomaticHaplotagParameters& params, const std::string& version) {

    ParamsHandler<HaplotagParameters>::initialize(params.basic, version);

    params.callerCfg.tumorPurity = 0.2;
    params.callerCfg.enableFilter = true;
    params.callerCfg.predictTumorPurity = true;
    params.metricsSuffix = "_somatic_haplotag.metrics";
}

bool ParamsHandler<SomaticHaplotagParameters>::loadArgument(SomaticHaplotagParameters& params, char& opt, std::istringstream& arg) {
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
            case SomaticHaplotagOption::BENCHMARK_VCF: arg >> params.benchmarkVcf; break;
            case SomaticHaplotagOption::BENCHMARK_BED: arg >> params.benchmarkBedFile; break;
            case SomaticHaplotagOption::DISABLE_FILTER: params.callerCfg.enableFilter = false; break;
            case SomaticHaplotagOption::TUMOR_PURITY: 
                arg >> params.callerCfg.tumorPurity; 
                params.callerCfg.predictTumorPurity = false;
                break;
            default: isLoaded = false; 
            break;
        }
    }
    return isLoaded;
}

bool ParamsHandler<SomaticHaplotagParameters>::validateFiles(SomaticHaplotagParameters& params, const std::string& programName) {
    // validate base haplotag files
    bool isValid = ParamsHandler<HaplotagParameters>::validateFiles(params.basic, programName);
    // validate somatic haplotag files
    isValid &= FileValidator::validateRequiredFile(params.tumorSnpFile, "tumor SNP file", programName);
    isValid &= FileValidator::validateRequiredFile(params.tumorBamFile, "tumor BAM file", programName);
    isValid &= FileValidator::validateOptionalFile(params.benchmarkVcf, "benchmark VCF file", programName);
    isValid &= FileValidator::validateOptionalFile(params.benchmarkBedFile, "benchmark BED file", programName);
    return isValid;
}

bool ParamsHandler<SomaticHaplotagParameters>::validateNumericParameter(SomaticHaplotagParameters& params, const std::string& programName) {
    bool isValid = ParamsHandler<HaplotagParameters>::validateNumericParameter(params.basic, programName);

    if (params.callerCfg.tumorPurity < 0.1 || params.callerCfg.tumorPurity > 1.0) {
        std::cerr << "[ERROR] " << programName << ": invalid tumor purity. value: " 
                << params.callerCfg.tumorPurity 
                << "\nthis value need: 0.1~1.0, --tumor-purity=Number\n";
        isValid = false;
    }
    
    return isValid;  
}

void ParamsHandler<SomaticHaplotagParameters>::recordCommand(SomaticHaplotagParameters& params, int argc, char** argv) {
    for(int i = 0; i < argc; ++i){
        params.basic.config.command.append(argv[i]);
        params.basic.config.command.append(" ");
    }
}

int ParamsHandler<SomaticHaplotagParameters>::getHelpEnumNum() {
    return HaplotagOption::OPT_HELP;
}

int SomaticHaplotagMain(int argc, char** argv, std::string in_version){
    
    SomaticHaplotagArgumentManager optionManager(SUBPROGRAM, in_version, CORRECT_USAGE_MESSAGE);

    optionManager.setOptions();

    optionManager.parseOptions(argc, argv);

    SomaticHaplotagParameters ecParams = optionManager.getParams();
    
    SomaticHaplotagProcess processor(ecParams);
    processor.taggingProcess();

    return 0;
}