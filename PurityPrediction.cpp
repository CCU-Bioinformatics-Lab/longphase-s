#include "PurityPrediction.h"

#define SUBPROGRAM "purity prediction"

void PurityPredictionHelpManager::buildMessage() {
    // build base message
    HaplotagHelpManager::buildMessage();

    //clear the section item
    clearSectionItem(REQUIRED_SECTION);

    addSectionItem(REQUIRED_SECTION, "      -s, --snp-file=NAME             input normal sample SNP VCF file.");
    addSectionItem(REQUIRED_SECTION, "      -b, --bam-file=NAME             input normal sample BAM file");
    addSectionItem(REQUIRED_SECTION, "      --tumor-snp-file=NAME           input tumor sample SNP VCF file.");
    addSectionItem(REQUIRED_SECTION, "      --tumor-bam-file=NAME           input tumor sample BAM file.");
    addSectionItem(REQUIRED_SECTION, "      -r, --reference=NAME            reference FASTA.\n");

    clearSectionItem(OPTIONAL_SECTION);
    addSectionItem(OPTIONAL_SECTION, "      -t, --threads=Num               number of thread. default:1");
    addSectionItem(OPTIONAL_SECTION, "      -o, --out-prefix=NAME           prefix of tumor purity prediction result. default:result");
}

void PurityPredictionOptionDefiner::defineOptions(ArgumentManager& manager) {
    // base haplotag options
    HaplotagOptionDefiner::defineOptions(manager);

    // somatic haplotag-specific options
    manager.addOption({"tumor-snp-file", required_argument, NULL, TUM_SNP});
    manager.addOption({"tumor-bam-file", required_argument, NULL, TUM_BAM});
}

void ParamsHandler<PurityPredictionParameters>::initialize(PurityPredictionParameters& params, const std::string& version) {

    ParamsHandler<HaplotagParameters>::initialize(params.basic, version);
    params.basic.qualityThreshold = 20;
    params.basic.tagSupplementary = true;
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
        params.basic.command.append(argv[i]);
        params.basic.command.append(" ");
    }
}

int ParamsHandler<PurityPredictionParameters>::getHelpEnumNum(){
    return HaplotagOption::OPT_HELP;
}

int PurityPredictionMain(int argc, char** argv, std::string in_version){

    PurityPredictionArgumentManager optionManager(SUBPROGRAM, in_version);

    optionManager.setOptions();
    optionManager.setHelpMessage();

    optionManager.parseOptions(argc, argv);

    PurityPredictionParameters ecParams = optionManager.getParams();

    PurityPredictionProcess processor(ecParams);

    processor.predictPurity();

    return 0;
}