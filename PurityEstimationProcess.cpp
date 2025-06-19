#include "PurityEstimationProcess.h"

PurityEstimProcess::PurityEstimProcess(PurityEstimParameters& params)
    : HaplotagProcess(params.basic), params(params) {
    tumorPurity = 0.0;
}

PurityEstimProcess::~PurityEstimProcess() {

}

void PurityEstimProcess::printParamsMessage() {
    std::cerr<< "LongPhase-S v" << params.basic.bamCfg.version << " - Estimate Tumor Purity \n";
    std::cerr<< "\n";
    std::cerr<< "[Input Files]\n";
    std::cerr<< "phased normal SNP file       : " << params.basic.snpFile << "\n";
    std::cerr<< "tumor SNP file               : " << params.tumorSnpFile << "\n";
    std::cerr<< "normal BAM file              : " << params.basic.bamFile << "\n";
    std::cerr<< "tumor BAM file               : " << params.tumorBamFile << "\n";
    std::cerr<< "reference file               : " << params.basic.fastaFile << "\n\n";
    std::cerr<< "[Output Files]\n";
    std::cerr<< "purity estimation file       :" << params.basic.bamCfg.resultPrefix + "_purity.out" << "\n";
    std::cerr<< "-------------------------------------------\n";
    std::cerr<< "[Purity Estimation Params] " << "\n";
    std::cerr<< "number of threads            : " << params.basic.bamCfg.numThreads << "\n";
    std::cerr<< "estimation region            : " << (!params.basic.bamCfg.region.empty() ? params.basic.bamCfg.region : "all") << "\n";
    std::cerr<< "filter mapping quality below : " << params.basic.bamCfg.qualityThreshold << "\n";
    std::cerr<< "percentage threshold         : " << params.basic.bamCfg.percentageThreshold << "\n";
    std::cerr<< "include supplementary reads  : " << (params.basic.bamCfg.tagSupplementary ? "enabled" : "disabled") << "\n";
    std::cerr<< "-------------------------------------------\n";

}

void PurityEstimProcess::parseVariantFiles(VcfParser& vcfParser) {
    // parse normal SNP, SV, MOD vcf file
    HaplotagProcess::parseVariantFiles(vcfParser);

    //load tumor snp vcf
    if(params.tumorSnpFile != ""){
        std::time_t begin = time(NULL);
        std::cerr<< "parsing tumor SNP VCF ... ";
        vcfParser.setParseSnpFile(true);
        vcfParser.parsingVCF(params.tumorSnpFile, vcfSet[Genome::TUMOR], *mergedChrVarinat);
        vcfParser.reset();
        std::cerr<< difftime(time(NULL), begin) << "s\n";
    }
}

void PurityEstimProcess::estimatePurity() {
    // set the tumor genome sample for parsing all 
    printParamsMessage();

    Genome geneSample = TUMOR;

    VcfParser vcfParser(geneSample);
    // load SNP vcf file
    parseVariantFiles(vcfParser);
    //decide which genome sample chrVec and chrLength belong to
    setChrVecAndChrLength();
    // update chromosome processing based on region
    setProcessingChromRegion();

    //estimate tumor purity 
    CallerConfig callerCfg;
    SomaticVarCaller *somaticVarCaller = new SomaticVarCaller(callerCfg, params.basic.bamCfg, *chrVec);
    somaticVarCaller->extractSomaticData(params.basic.bamFile, params.tumorBamFile, params.basic.fastaFile, params.basic.bamCfg, *chrVec, *chrLength, *mergedChrVarinat, vcfSet);
    tumorPurity = somaticVarCaller->runTumorPurityEstimator(params.basic.bamCfg.writeReadLog, params.basic.bamCfg.resultPrefix);
    delete somaticVarCaller;

    printExecutionReport();
}

void PurityEstimProcess::printExecutionReport(){
    std::cerr<< "-------------------------------------------\n";
    std::cerr<< "total process time:    " << difftime(time(NULL), processBegin) << "s\n";
    std::cerr<< "estimated tumor purity: " << tumorPurity << "\n";
    std::cerr<< "-------------------------------------------\n";
}