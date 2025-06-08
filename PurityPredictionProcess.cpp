#include "PurityPredictionProcess.h"

PurityPredictionProcess::PurityPredictionProcess(PurityPredictionParameters& params)
    : HaplotagProcess(params.basic), params(params) {
    tumorPurity = 0.0;
}

PurityPredictionProcess::~PurityPredictionProcess() {

}

void PurityPredictionProcess::printParamsMessage() {
    std::cerr<< "LongPhase-S v" << params.basic.config.version << " - Predict Tumor Purity\n";
    std::cerr<< "\n";
    std::cerr<< "[Input Files]\n";
    std::cerr<< "phased normal SNP file       : " << params.basic.snpFile << "\n";
    std::cerr<< "tumor SNP file               : " << params.tumorSnpFile << "\n";
    std::cerr<< "normal BAM file              : " << params.basic.bamFile << "\n";
    std::cerr<< "tumor BAM file               : " << params.tumorBamFile << "\n";
    std::cerr<< "reference file               : " << params.basic.fastaFile << "\n\n";
    std::cerr<< "[Output Files]\n";
    std::cerr<< "purity prediction file       :" << params.basic.config.resultPrefix + "_purity.out" << "\n";
    std::cerr<< "-------------------------------------------\n";
    std::cerr<< "[Purity Prediction Params] " << "\n";
    std::cerr<< "number of threads            : " << params.basic.config.numThreads << "\n";
    std::cerr<< "prediction region            : " << (!params.basic.config.region.empty() ? params.basic.config.region : "all") << "\n";
    std::cerr<< "filter mapping quality below : " << params.basic.config.qualityThreshold << "\n";
    std::cerr<< "percentage threshold         : " << params.basic.config.percentageThreshold << "\n";
    std::cerr<< "tag supplementary            : " << (params.basic.config.tagSupplementary ? "enabled" : "disabled") << "\n";
    std::cerr<< "-------------------------------------------\n";

}

void PurityPredictionProcess::parseVariantFiles(VcfParser& vcfParser) {
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

void PurityPredictionProcess::predictPurity() {
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

    //predict tumor purity 
    CallerConfig callerCfg;
    SomaticVarCaller *somaticVarCaller = new SomaticVarCaller(callerCfg, params.basic.config, *chrVec);
    somaticVarCaller->extractSomaticData(params.basic.bamFile, params.tumorBamFile, params.basic.fastaFile, params.basic.config, *chrVec, *chrLength, *mergedChrVarinat, vcfSet);
    tumorPurity = somaticVarCaller->runTumorPurityPredictor(params.basic.config.writeReadLog, params.basic.config.resultPrefix);
    delete somaticVarCaller;

    printExecutionReport();
}

void PurityPredictionProcess::printExecutionReport(){
    std::cerr<< "-------------------------------------------\n";
    std::cerr<< "total process time:    " << difftime(time(NULL), processBegin) << "s\n";
    std::cerr<< "predicted tumor purity: " << tumorPurity << "\n";
    std::cerr<< "-------------------------------------------\n";
}