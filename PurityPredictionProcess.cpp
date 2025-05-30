#include "PurityPredictionProcess.h"

PurityPredictionProcess::PurityPredictionProcess(PurityPredictionParameters& params)
    : HaplotagProcess(params.basic), params(params) {
}

PurityPredictionProcess::~PurityPredictionProcess() {

}

void PurityPredictionProcess::parseVariantFiles(VcfParser& vcfParser) {
    // parse normal SNP, SV, MOD vcf file
    HaplotagProcess::parseVariantFiles(vcfParser);

    //load tumor snp vcf
    if(params.tumorSnpFile != ""){
        std::time_t begin = time(NULL);
        std::cerr<< "parsing tumor SNP VCF ... ";
        vcfParser.setParseSnpFile(true);
        vcfParser.variantParser(params.tumorSnpFile, vcfSet[Genome::TUMOR], *mergedChrVarinat);
        vcfParser.reset();
        std::cerr<< difftime(time(NULL), begin) << "s\n";
    }
}

void PurityPredictionProcess::predictPurity() {
    // set the tumor genome type for parsing all Genotype

    paramsMessage.addParamsMessage();
    paramsMessage.printParamsMessage();

    Genome geneType = TUMOR;

    VcfParser vcfParser(geneType);
    // load SNP, SV, MOD vcf file
    parseVariantFiles(vcfParser);
    //decide which genome type chrVec and chrLength belong to
    setChrVecAndChrLength();
    // update chromosome processing based on region
    setProcessingChromRegion();

    //predict tumor purity 
    SomaticVarCaller *somaticVarCaller = new SomaticVarCaller(*chrVec);
    somaticVarCaller->extractSomaticData(params.basic.bamFile, params.tumorBamFile, params.basic, *chrVec, *chrLength, *mergedChrVarinat, vcfSet);
    somaticVarCaller->runTumorPurityPredictor(params.basic.writeReadLog, params.basic.resultPrefix+"_test");

    delete somaticVarCaller;
}