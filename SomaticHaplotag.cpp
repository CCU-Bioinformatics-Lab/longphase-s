#include "SomaticHaplotag.h"

#define SUBPROGRAM "somatic haplotag"

void SomaticHaplotagHelpManager::modifyMessage() {
    const std::string sectionName = "required arguments:";
    //clear the section item
    clearSectionItem(sectionName);

    // Required arguments - Somatic mode
    addSectionItem(sectionName, "   [Somatic mode] (tumor/normal pair data):");
    addSectionItem(sectionName, "      -s, --snp-file=NAME             input normal sample SNP VCF file.");
    addSectionItem(sectionName, "      -b, --bam-file=NAME             input normal sample BAM file (used as a reference for comparison).");
    addSectionItem(sectionName, "      --tumor-snp-file=NAME           input tumor sample SNP VCF file.");
    addSectionItem(sectionName, "      --tumor-bam-file=NAME           input tumor sample BAM file for mutation tagging.");
    addSectionItem(sectionName, "      -r, --reference=NAME            reference FASTA.\n");
}

void SomaticHaplotagOptionManager::extendOptions() {
    // Tumor-specific options
    addOption({"tumor-snp-file", required_argument, NULL, TUM_SNP});
    addOption({"tumor-bam-file", required_argument, NULL, TUM_BAM});
    addOption({"highCon-snp", required_argument, NULL, BENCHMARK_VCF});
    addOption({"benchmark-bed", required_argument, NULL, BENCHMARK_BED});
    addOption({"disableFilter", no_argument, NULL, DISABLE_FILTER});
}

bool SomaticHaplotagOptionManager::validateExtendFiles() {
    bool isValid = true;
    isValid &= validateRequiredFile(ecParams.tumorSnpFile, "tumor SNP file");
    isValid &= validateRequiredFile(ecParams.tumorBamFile, "tumor BAM file");
    isValid &= validateOptionalFile(ecParams.benchmarkVcf, "benchmark VCF file");
    isValid &= validateOptionalFile(ecParams.benchmarkBedFile, "benchmark BED file");
    return isValid;
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
        default: isLoaded = false; break;
    }
    return isLoaded;
}


int SomaticHaplotagMain(int argc, char** argv, std::string in_version){
    
    SomaticHaplotagOptionManager optionManager;

    optionManager.setOptions();
    optionManager.setHelpMessage();

    optionManager.parseOptions(argc, argv);
    optionManager.setVersion(in_version);

    HaplotagParameters ecParams = optionManager.getParams();
    
    SomaticHaplotagProcess processor(ecParams);
    processor.taggingProcess();

    return 0;
}