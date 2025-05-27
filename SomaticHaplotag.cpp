#include "SomaticHaplotag.h"

#define SUBPROGRAM "somatic haplotag"

void SomaticHaplotagHelpManager::modifyMessage() {
    const std::string requireSection = "required arguments:";
    //clear the section item
    clearSectionItem(requireSection);
    // Required arguments - Somatic mode
    addSectionItem(requireSection, "      -s, --snp-file=NAME             input normal sample SNP VCF file.");
    addSectionItem(requireSection, "      -b, --bam-file=NAME             input normal sample BAM file (used as a reference for comparison).");
    addSectionItem(requireSection, "      --tumor-snp-file=NAME           input tumor sample SNP VCF file.");
    addSectionItem(requireSection, "      --tumor-bam-file=NAME           input tumor sample BAM file for mutation tagging.");
    addSectionItem(requireSection, "      -r, --reference=NAME            reference FASTA.\n");

    const std::string optionalSection = "optional arguments:";
    addSectionItem(optionalSection, "      --tumor-purity=Num              tumor purity value (0.1~1.0) used to adjust somatic variant calling sensitivity and specificity");
}

void SomaticHaplotagOptionManager::extendOptions() {
    // Tumor-specific options
    addOption({"tumor-snp-file", required_argument, NULL, TUM_SNP});
    addOption({"tumor-bam-file", required_argument, NULL, TUM_BAM});
    addOption({"benchmark-snp", required_argument, NULL, BENCHMARK_VCF});
    addOption({"benchmark-bed", required_argument, NULL, BENCHMARK_BED});
    addOption({"disableFilter", no_argument, NULL, DISABLE_FILTER});
    addOption({"tumor-purity", required_argument, NULL, TUMOR_PURITY});
}

bool SomaticHaplotagOptionManager::validateExtendFiles() {
    bool isValid = true;
    isValid &= validateRequiredFile(ecParams.tumorSnpFile, "tumor SNP file");
    isValid &= validateRequiredFile(ecParams.tumorBamFile, "tumor BAM file");
    isValid &= validateOptionalFile(ecParams.benchmarkVcf, "benchmark VCF file");
    isValid &= validateOptionalFile(ecParams.benchmarkBedFile, "benchmark BED file");
    return isValid;
}

void SomaticHaplotagOptionManager::initializeDefaultValues() {
    HaplotagOptionManager::initializeDefaultValues();

    ecParams.tumorPurity = 0.2;
    ecParams.enableFilter = true;
    ecParams.predictTumorPurity = true;
}

bool SomaticHaplotagOptionManager::loadExtendOptions(char& opt, std::istringstream& arg) {
    bool isLoaded = true;
    switch (opt)
    {
        case HaplotagOption::TUM_SNP: arg >> ecParams.tumorSnpFile; break;
        case HaplotagOption::TUM_BAM: arg >> ecParams.tumorBamFile; break;
        case HaplotagOption::BENCHMARK_VCF: arg >> ecParams.benchmarkVcf; break;
        case HaplotagOption::BENCHMARK_BED: arg >> ecParams.benchmarkBedFile; break;
        case HaplotagOption::DISABLE_FILTER: ecParams.enableFilter = false; break;
        case HaplotagOption::TUMOR_PURITY: 
            arg >> ecParams.tumorPurity; 
            ecParams.predictTumorPurity = false;
            break;
        default: isLoaded = false; break;
    }
    return isLoaded;
}

bool SomaticHaplotagOptionManager::validateNumericParameter() {
    
    bool isValid = HaplotagOptionManager::validateNumericParameter();
    
    if (ecParams.tumorPurity < 0.1 || ecParams.tumorPurity > 1.0) {
        std::cerr << "[ERROR] " << programName << ": invalid tumor purity. value: " 
                << ecParams.tumorPurity 
                << "\nthis value need: 0.1~1.0, --tumor-purity=Number\n";
        isValid = false;
    }
    
    return isValid;
}
int SomaticHaplotagMain(int argc, char** argv, std::string in_version){
    
    SomaticHaplotagOptionManager optionManager(SUBPROGRAM);

    optionManager.setOptions();
    optionManager.setHelpMessage();

    optionManager.parseOptions(argc, argv);
    optionManager.setVersion(in_version);

    HaplotagParameters ecParams = optionManager.getParams();
    
    SomaticHaplotagProcess processor(ecParams);
    processor.taggingProcess();

    return 0;
}