#include "SomaticHaplotag.h"

#define SUBPROGRAM "somatic haplotag"

void SomaticHaplotagHelpManager::buildMessage() {
    // build base message
    HaplotagHelpManager::buildMessage();

    const std::string requireSection = "required arguments:";
    const std::string optionalSection = "optional arguments:";
    //clear the section item
    clearSectionItem(requireSection);
    // Required arguments - Somatic mode
    addSectionItem(requireSection, "      -s, --snp-file=NAME             input normal sample SNP VCF file.");
    addSectionItem(requireSection, "      -b, --bam-file=NAME             input normal sample BAM file (used as a reference for comparison).");
    addSectionItem(requireSection, "      --tumor-snp-file=NAME           input tumor sample SNP VCF file.");
    addSectionItem(requireSection, "      --tumor-bam-file=NAME           input tumor sample BAM file for mutation tagging.");
    addSectionItem(requireSection, "      -r, --reference=NAME            reference FASTA.\n");

    addSectionItem(optionalSection, " ");
    
    addSection("somatic variant calling arguments:");
    addItem("      --tumor-purity=Num              tumor purity value (0.1~1.0) used to adjust somatic variant calling sensitivity and specificity");
}

void SomaticHaplotagOptionDefiner::defineOptions(ArgumentManager& manager) {
    // base haplotag options
    HaplotagOptionDefiner::defineOptions(manager);

    // somatic haplotag-specific options
    manager.addOption({"tumor-snp-file", required_argument, NULL, TUM_SNP});
    manager.addOption({"tumor-bam-file", required_argument, NULL, TUM_BAM});
    manager.addOption({"benchmark-snp", required_argument, NULL, BENCHMARK_VCF});
    manager.addOption({"benchmark-bed", required_argument, NULL, BENCHMARK_BED});
    manager.addOption({"disableFilter", no_argument, NULL, DISABLE_FILTER});
    manager.addOption({"tumor-purity", required_argument, NULL, TUMOR_PURITY});
}

void ParamsHandler<SomaticHaplotagParameters>::initialize(SomaticHaplotagParameters& params, const std::string& version) {

    ParamsHandler<HaplotagParameters>::initialize(params.basic, version);

    params.tumorPurity = 0.2;
    params.tumorPurity = 0.2;
    params.enableFilter = true;
    params.predictTumorPurity = true;
    params.onlyPredictTumorPurity = false;
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
            case HaplotagOption::TUM_SNP: arg >> params.tumorSnpFile; break;
            case HaplotagOption::TUM_BAM: arg >> params.tumorBamFile; break;
            case HaplotagOption::BENCHMARK_VCF: arg >> params.benchmarkVcf; break;
            case HaplotagOption::BENCHMARK_BED: arg >> params.benchmarkBedFile; break;
            case HaplotagOption::DISABLE_FILTER: params.enableFilter = false; break;
            case HaplotagOption::TUMOR_PURITY: 
                arg >> params.tumorPurity; 
                params.predictTumorPurity = false;
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

    if (params.tumorPurity < 0.1 || params.tumorPurity > 1.0) {
        std::cerr << "[ERROR] " << programName << ": invalid tumor purity. value: " 
                << params.tumorPurity 
                << "\nthis value need: 0.1~1.0, --tumor-purity=Number\n";
        isValid = false;
    }
    
    return isValid;  
}

void ParamsHandler<SomaticHaplotagParameters>::recordCommand(SomaticHaplotagParameters& params, int argc, char** argv) {
    for(int i = 0; i < argc; ++i){
        params.basic.command.append(argv[i]);
        params.basic.command.append(" ");
    }
}

int ParamsHandler<SomaticHaplotagParameters>::getHelpEnumNum() {
    return HaplotagOption::OPT_HELP;
}

int SomaticHaplotagMain(int argc, char** argv, std::string in_version){
    
    SomaticHaplotagArgumentManager optionManager(SUBPROGRAM, in_version);

    optionManager.setOptions();
    optionManager.setHelpMessage();

    optionManager.parseOptions(argc, argv);
    // optionManager.setVersion(in_version);

    SomaticHaplotagParameters ecParams = optionManager.getParams();
    
    SomaticHaplotagProcess processor(ecParams);
    processor.taggingProcess();

    return 0;
}