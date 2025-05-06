#include "SomaticBenchmark.h"

SomaticReadVerifier::SomaticReadVerifier(bool openTestingFunc, SomaticReadMetrics *metrics): 
openTestingFunc(openTestingFunc),
metrics(metrics)
{}

SomaticReadVerifier::~SomaticReadVerifier(){

}

void SomaticReadVerifier::recordDelReadCount(const std::string &chr, std::map<int, MultiGenomeVar>::iterator &currentVariantIter){
    if(!openTestingFunc) return;
    
    if(currentVariantIter->second.isExists(HIGH_CON_SOMATIC)){
        int pos = (*currentVariantIter).first;
        metrics->posAltRefDelCount[pos].delCount++;

        //record somatic position for record crossing high con snp read
        metrics->highConSomaticPos.push_back(std::make_pair(pos, SnpHP::NONE_SNP));
    }
}

void SomaticReadVerifier::recordRefAltAlleleCount(const std::string &chr, std::string &base, std::map<int, MultiGenomeVar>::iterator &currentVariantIter){
    if(!openTestingFunc) return;

    if(currentVariantIter->second.isExists(HIGH_CON_SOMATIC)){
        int pos = currentVariantIter->first;
        std::string& refAllele = currentVariantIter->second.Variant[HIGH_CON_SOMATIC].allele.Ref;
        std::string& altAllele = currentVariantIter->second.Variant[HIGH_CON_SOMATIC].allele.Alt;

        int baseHP = SnpHP::NONE_SNP;

        if(base == refAllele){
            metrics->posAltRefDelCount[pos].refCount++;
        }else if(base == altAllele){
            metrics->posAltRefDelCount[pos].altCount++;
            baseHP = SnpHP::SOMATIC_H3;
        }
        //record somatic position for record crossing high con snp read
        metrics->highConSomaticPos.push_back(std::make_pair(pos, baseHP));
    }
}

SomaticReadLog SomaticReadVerifier::createBasicSomaticReadLog(const std::string &chr, std::string &readID, std::string &hpResult, double &norHPsimilarity, float &deriveByHpSimilarity, std::map<int, int> &hpCount){
    SomaticReadLog tmp;
    tmp.chr = chr;
    tmp.readID = readID;
    tmp.hpResult = hpResult;
    tmp.germlineVarSimilarity = norHPsimilarity;
    tmp.deriveByHpSimilarity = deriveByHpSimilarity;
    tmp.germlineSnpCount = hpCount[SnpHP::GERMLINE_H1] + hpCount[SnpHP::GERMLINE_H2];
    tmp.tumorSnpCount = hpCount[SnpHP::SOMATIC_H3];

    return tmp;
}

void SomaticReadVerifier::recordCrossingHighConSnpRead(const std::string &chr, std::string &readID, std::string &hpResult, std::map<int, int> &variantsHP, std::map<int, int> &hpCount, double &norHPsimilarity, float &deriveByHpSimilarity, std::map<int, MultiGenomeVar> &currentChrVariants){
    // if not open testing function, return
    if(!openTestingFunc) return;

    SomaticReadLog tmp = createBasicSomaticReadLog(chr, readID, hpResult, norHPsimilarity, deriveByHpSimilarity, hpCount);

    bool isCrossHighConSomatic = false;
    bool existHighConVaraints = false;

    for(auto varIter : metrics->highConSomaticPos){
        int pos = varIter.first;
        int baseHP = varIter.second;

        //transfer hpResult to H3
        if(baseHP == SnpHP::SOMATIC_H3){
            existHighConVaraints = true;
        }

        tmp.somaticSnpHp[pos] = baseHP;

        isCrossHighConSomatic = true;
    }

    if(isCrossHighConSomatic){
        //exist high con variants alt allele
        if(existHighConVaraints){
            //correction hpResult that exist high con variants
            if(hpResult == "1"){
                tmp.hpResult = "1-1";
            }else if(hpResult == "2"){
                tmp.hpResult = "2-1";
            }else if(hpResult == "."){
                tmp.hpResult = "3";
            }
        }else{
            //correction hpResult that not exist high con variants
            if(hpResult == "2-1"){
                tmp.hpResult = "2";
            }else if(hpResult == "1-1"){
                tmp.hpResult = "1";
            }else if(hpResult == "3"){
                tmp.hpResult = ".";
            }
        }
    }

    if(isCrossHighConSomatic){
        metrics->readsCrossingHighConSnpVec.push_back(tmp);
    }
    
    //clear high con somatic position in current read for next read
    if(!metrics->highConSomaticPos.empty()){
        metrics->highConSomaticPos.clear();
    }
}

void SomaticReadVerifier::recordTaggedRead(const std::string &chr, std::string &readID, std::string &hpResult, std::map<int, int> &variantsHP, std::map<int, int> &hpCount, double &norHPsimilarity, float &deriveByHpSimilarity, std::map<int, MultiGenomeVar> &currentChrVariants){
    // if not open testing function, return
    if(!openTestingFunc || hpResult == ".") return;

    SomaticReadLog tmp = createBasicSomaticReadLog(chr, readID, hpResult, norHPsimilarity, deriveByHpSimilarity, hpCount);

    bool readExistHighConSomatic = false;

    auto varIter = variantsHP.begin();
    while(varIter != variantsHP.end()){
        int pos = varIter->first;
        int snpHP = varIter->second;
        if(currentChrVariants.find(pos) != currentChrVariants.end()){
            if(currentChrVariants[pos].isExists(HIGH_CON_SOMATIC) && (snpHP == SnpHP::SOMATIC_H3)){
                tmp.somaticSnpHp[pos] = snpHP;
                readExistHighConSomatic = true;
            }
        }
        varIter++;
    }

    if(readExistHighConSomatic){
        metrics->taggedSomaticReadVec.push_back(tmp);
    }

    metrics->totalReadVec.push_back(tmp);
}

SomaticReadBenchmark::SomaticReadBenchmark(){
    setParseSnpFile(true);
    openTestingFunc = false;
}
SomaticReadBenchmark::~SomaticReadBenchmark(){

}

void SomaticReadBenchmark::setTestingFunc(bool openTestingFunc){
    this->openTestingFunc = openTestingFunc;
}

bool SomaticReadBenchmark::getTestingFunc(){
    return openTestingFunc;
}

void SomaticReadBenchmark::loadChrKey(const std::string &chr){
    chrMetrics[chr] = SomaticReadMetrics();
}

void SomaticReadBenchmark::loadHighConSomatic(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat){
    if(!openTestingFunc) return;
    variantParser(input, Info, mergedChrVarinat);
}

void SomaticReadBenchmark::parserProcess(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat){
    if(!openTestingFunc) return;
    
    if( input.substr(0, 2) == "##" && getParseSnpFile()){
        if( input.find("contig=")!= std::string::npos ){
            int id_start  = input.find("ID=")+3;
            int id_end    = input.find(",length=");
            int len_start = id_end+8;
            int len_end   = input.find(">");
            
            std::string chr = input.substr(id_start,id_end-id_start);
            int chrLen = std::stoi( input.substr(len_start,len_end-len_start) );

            Info.chrVec.push_back(chr);
            Info.chrLength[chr]=chrLen;                
        }
    }
    else if ( input.substr(0, 1) == "#" ){
        
    }
    else{
        std::istringstream iss(input);
        std::vector<std::string> fields((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());

        if( fields.size() == 0 )
            return;
            
        // trans to 0-base
        int pos = std::stoi( fields[1] ) - 1;
        std::string chr = fields[0];

        VarData varData;
        varData.allele.Ref = fields[3];
        varData.allele.Alt = fields[4];
        mergedChrVarinat[chr][pos].Variant[Genome::HIGH_CON_SOMATIC] = varData;
    }
}

void SomaticReadBenchmark::parseBedFile(const std::string& bedFile) {
    std::ifstream file(bedFile);
    if (!file.is_open()) {
        std::cerr << "Failed to open BED file: " << bedFile << std::endl;
        return;
    }
    
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        processBedLine(line);
    }
}

void SomaticReadBenchmark::processBedLine(const std::string& line) {
    std::istringstream iss(line);
    std::string chr;
    int start, end;
    iss >> chr >> start >> end;
    
    bedRegions[chr].push_back(BedRegion{start, end - 1});
}

void SomaticReadBenchmark::markVariantsInBedRegions(std::vector<std::string> &chrVec, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat) {
    int tumorInBedRegionCount = 0;
    int tumorOutBedRegionCount = 0;
    int normalInBedRegionCount = 0;
    int normalOutBedRegionCount = 0;
    int benchmarkInBedRegionCount = 0;
    int benchmarkOutBedRegionCount = 0;

    for (auto& chr : chrVec) {
        auto& chrPosVariants = mergedChrVarinat[chr];
        
        // check if this chromosome has bed regions
        auto bedIt = bedRegions.find(chr);
        if (bedIt == bedRegions.end()) {
            // if this chromosome has no bed regions, mark all points as false
            for (auto& curVar : chrPosVariants) {
                curVar.second.isInBedRegion = false;
                if(curVar.second.isExists(Genome::TUMOR)) {
                    tumorOutBedRegionCount++;
                }
                if(curVar.second.isExists(Genome::NORMAL)) {
                    normalOutBedRegionCount++;
                }
                if(curVar.second.isExists(Genome::HIGH_CON_SOMATIC)){
                    benchmarkOutBedRegionCount++;
                }
            }
            continue;
        }

        const auto& regions = bedIt->second;
        if (regions.empty()) {
            // If there are no bed regions, mark all points as false
            for (auto& curVar : chrPosVariants) {
                curVar.second.isInBedRegion = false;
                if(curVar.second.isExists(Genome::TUMOR)) {
                    tumorOutBedRegionCount++;
                }
                if(curVar.second.isExists(Genome::NORMAL)) {
                    normalOutBedRegionCount++;
                }
                if(curVar.second.isExists(Genome::HIGH_CON_SOMATIC)){
                    benchmarkOutBedRegionCount++;
                }
            }
            continue;
        }

        auto regionIter = regions.begin();
        auto varPosIter = chrPosVariants.begin();
        
        // traverse bed regions and variants at the same time
        while (varPosIter != chrPosVariants.end()) {
            int variantPos = varPosIter->first;
            
            // if variant position is out of current bed region range
            while (regionIter != regions.end() && variantPos > regionIter->end) {
                ++regionIter;
            }
            
            // check if variant is in current bed region
            if (regionIter != regions.end() && 
                variantPos >= regionIter->start && 
                variantPos <= regionIter->end) {
                varPosIter->second.isInBedRegion = true;
                if(varPosIter->second.isExists(Genome::TUMOR)) {
                    tumorInBedRegionCount++;
                }
                if(varPosIter->second.isExists(Genome::NORMAL)) {
                    normalInBedRegionCount++;
                }
                if(varPosIter->second.isExists(Genome::HIGH_CON_SOMATIC)){
                    benchmarkInBedRegionCount++;
                }
            } else {
                varPosIter->second.isInBedRegion = false;
                if(varPosIter->second.isExists(Genome::TUMOR)) {
                    tumorOutBedRegionCount++;
                }
                if(varPosIter->second.isExists(Genome::NORMAL)) {
                    normalOutBedRegionCount++;
                }
                if(varPosIter->second.isExists(Genome::HIGH_CON_SOMATIC)){
                    benchmarkOutBedRegionCount++;
                }
            }
            
            ++varPosIter;
        }
    }
    
    std::cout << "Tumor in bed region count: " << tumorInBedRegionCount << "\n";
    std::cout << "Tumor out bed region count: " << tumorOutBedRegionCount << "\n";
    std::cout << "Normal in bed region count: " << normalInBedRegionCount << "\n";
    std::cout << "Normal out bed region count: " << normalOutBedRegionCount << "\n";
    std::cout << "Benchmark in bed region count: " << benchmarkInBedRegionCount << "\n";
    std::cout << "Benchmark out bed region count: " << benchmarkOutBedRegionCount << "\n";
}

void SomaticReadBenchmark::removeVariantsOutBedRegion(std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat){
    //remove variants out bed region
    for (auto& chrPair : mergedChrVarinat) {
        auto& chrPosVariants = chrPair.second;
        
        auto varPosIter = chrPosVariants.begin();
        while (varPosIter != chrPosVariants.end()) {
            if (!varPosIter->second.isInBedRegion) {
                bool hasTumor = varPosIter->second.isExists(Genome::TUMOR);
                bool hasHighCon = varPosIter->second.isExists(Genome::HIGH_CON_SOMATIC);
                
                if (hasTumor || hasHighCon) {
                    // if no exist NORMAL data, remove the whole position
                    if (!varPosIter->second.isExists(Genome::NORMAL)) {
                        varPosIter = chrPosVariants.erase(varPosIter);
                        continue;
                    } else {
                        // if there is NORMAL data, only remove the data to be cleared
                        if (hasTumor) {
                            varPosIter->second.Variant.erase(Genome::TUMOR);
                        }
                        if (hasHighCon) {
                            varPosIter->second.Variant.erase(Genome::HIGH_CON_SOMATIC);
                        }
                        ++varPosIter;
                    }
                } else {
                    ++varPosIter;
                }
            } else {
                ++varPosIter;
            }
        }
    }
}

void SomaticReadBenchmark::writeBedRegionLog(const std::vector<std::string>& chrVec,
                                           const std::map<std::string, std::map<int, MultiGenomeVar>>& mergedChrVarinat,
                                           const std::string& outPrefix) {
    std::ofstream inBedLog(outPrefix + "_in_bed.out");
    std::ofstream outBedLog(outPrefix + "_out_bed.out");
    
    std::string header = "#Chr\tPosition\tRef\tAlt\tVariant_Type\n";
    inBedLog << header;
    outBedLog << header;

    for (const auto& chr : chrVec) {
        auto chrPosVariants = mergedChrVarinat.at(chr);

        auto varPosIter = chrPosVariants.begin();

        while(varPosIter != chrPosVariants.end()){
            int pos = varPosIter->first + 1;
            
            std::string baseInfo = chr + "\t" + std::to_string(pos) + "\t";
            
            if(varPosIter->second.isExists(Genome::TUMOR)){
                auto& tumorVar = varPosIter->second.Variant.at(TUMOR);
                std::string varInfo = baseInfo + 
                                    tumorVar.allele.Ref + "\t" + 
                                    tumorVar.allele.Alt + "\t" +
                                    "TUMOR\n";

                if (varPosIter->second.isInBedRegion) {
                    inBedLog << varInfo;
                } else {
                    outBedLog << varInfo;
                }
            }

            varPosIter++;
        }
    }
    
    inBedLog.close();
    outBedLog.close();

}

SomaticReadMetrics* SomaticReadBenchmark::getMetricsPtr(const std::string &chr){
    return &(chrMetrics[chr]);
}

void SomaticReadBenchmark::writePosAlleleCountLog(
    std::vector<std::string>& chrVec,
    HaplotagParameters &params,
    std::string logPosfix,
    std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat
){
    // if not open testing function, return
    if(!openTestingFunc) return;

    std::ofstream *refAltCountLog=NULL;
    refAltCountLog=new std::ofstream(params.resultPrefix + logPosfix);
    int totalVariantCount = 0;

    for(auto chr: chrVec){
        totalVariantCount += chrMetrics[chr].posAltRefDelCount.size();
    }

    if(!refAltCountLog->is_open()){
        std::cerr<< "Fail to open write file: " << params.resultPrefix + logPosfix << "\n";
        exit(1);
    }else{
        (*refAltCountLog) << "#############################\n";
        (*refAltCountLog) << "# Somatic SNP allele count #\n";
        (*refAltCountLog) << "#############################\n";
        (*refAltCountLog) << "##High confidence VCF:"  << params.benchmarkVcf << "\n";
        (*refAltCountLog) << "##MappingQualityThreshold:"  << params.qualityThreshold << "\n";
        (*refAltCountLog) << "##Tatal variants:"  << totalVariantCount << "\n";
        (*refAltCountLog) << "#CHROM\t"
                          << "POS\t"
                          << "REF\t"
                          << "ALT\t"
                          << "REF_COUNT\t"
                          << "ALT_COUNT\t"
                          << "DEL_COUNT\n";
    }

    for(auto chr: chrVec){
        for(auto &posIter : chrMetrics[chr].posAltRefDelCount){
            (*refAltCountLog) << chr << "\t"
                              << posIter.first << "\t"
                              << mergedChrVarinat[chr][posIter.first].Variant[HIGH_CON_SOMATIC].allele.Ref << "\t"
                              << mergedChrVarinat[chr][posIter.first].Variant[HIGH_CON_SOMATIC].allele.Alt << "\t"
                              << posIter.second.refCount << "\t"
                              << posIter.second.altCount << "\t"
                              << posIter.second.delCount << "\n";
        }
    }

    refAltCountLog->close();
    delete refAltCountLog;
    refAltCountLog = nullptr;
}

void SomaticReadBenchmark::writeTaggedReadLog(
    const std::vector<std::string>& chrVec,
    HaplotagParameters &params,
    std::string logPosfix
){
    // if not open testing function, return
    if(!openTestingFunc) return;

    std::map<std::string, std::vector<SomaticReadLog>*> somaticReadVecMap;

    auto metricsIter = chrMetrics.begin();
    while(metricsIter != chrMetrics.end()){
        setChrSomaticReadVecPtr(metricsIter->first, somaticReadVecMap, metricsIter->second.totalReadVec);
        metricsIter++;
    }

    writeReadLog(chrVec, params, logPosfix, somaticReadVecMap);
}
void SomaticReadBenchmark::writeTaggedSomaticReadLog(
    const std::vector<std::string>& chrVec,
    HaplotagParameters &params,
    std::string logPosfix
){
    // if not open testing function, return
    if(!openTestingFunc) return;

    std::map<std::string, std::vector<SomaticReadLog>*> somaticReadVecMap;

    auto metricsIter = chrMetrics.begin();
    while(metricsIter != chrMetrics.end()){
        setChrSomaticReadVecPtr(metricsIter->first, somaticReadVecMap, metricsIter->second.taggedSomaticReadVec);
        metricsIter++;
    }

    writeReadLog(chrVec, params, logPosfix, somaticReadVecMap);
}

void SomaticReadBenchmark::writeCrossHighConSnpReadLog(
    const std::vector<std::string>& chrVec,
    HaplotagParameters &params,
    std::string logPosfix
){
    // if not open testing function, return
    if(!openTestingFunc) return;

    std::map<std::string, std::vector<SomaticReadLog>*> somaticReadVecMap;

    auto metricsIter = chrMetrics.begin();
    while(metricsIter != chrMetrics.end()){
        setChrSomaticReadVecPtr(metricsIter->first, somaticReadVecMap, metricsIter->second.readsCrossingHighConSnpVec);
        metricsIter++;
    }

    writeReadLog(chrVec, params, logPosfix, somaticReadVecMap);
}

void SomaticReadBenchmark::setChrSomaticReadVecPtr(
    const std::string &chr,
    std::map<std::string, std::vector<SomaticReadLog>*> &somaticReadVecMap,
    std::vector<SomaticReadLog> &somaticReadVec
){
    somaticReadVecMap[chr] = &somaticReadVec;
}

void SomaticReadBenchmark::writeReadLog(
    const std::vector<std::string>& chrVec,
    HaplotagParameters &params,
    std::string logPosfix,
    std::map<std::string, std::vector<SomaticReadLog>*> &somaticReadVecMap
){
    // if not open testing function, return
    if(!openTestingFunc) return;

    std::ofstream *somaticReadLog=NULL;
    somaticReadLog=new std::ofstream(params.resultPrefix + logPosfix);

    int totalTruthSomaticReads = 0;
    int totalTaggedSomaticReads = 0;
    int totalReads = 0;
    for(auto chr: chrVec){
        for(auto readIter :chrMetrics[chr].readsCrossingHighConSnpVec){
            if(readIter.hpResult == "1-1" || readIter.hpResult == "2-1" || readIter.hpResult == "3"){
                totalTruthSomaticReads++;
            }
        }

        for(auto readIter :chrMetrics[chr].taggedSomaticReadVec){
            if(readIter.hpResult == "1-1" || readIter.hpResult == "2-1" || readIter.hpResult == "3"){
                totalTaggedSomaticReads++;
            }
        }
        totalReads += chrMetrics[chr].totalReadVec.size();
    }


    if(!somaticReadLog->is_open()){
        std::cerr<< "Fail to open write file: " << params.resultPrefix + logPosfix << "\n";
        exit(1);
    }else{
        (*somaticReadLog) << "#####################\n";
        (*somaticReadLog) << "# Somatic Reads Log #\n";
        (*somaticReadLog) << "#####################\n";
        (*somaticReadLog) << "##High confidence VCF: "  << params.benchmarkVcf << "\n";
        (*somaticReadLog) << "##MappingQualityThreshold: "  << params.qualityThreshold << "\n";
        (*somaticReadLog) << "##Tatal reads: "  << totalReads << "\n";
        (*somaticReadLog) << "##Tatal tagged somatic reads: "  << totalTaggedSomaticReads << "\n";
        (*somaticReadLog) << "##Tatal truth somatic reads: "  << totalTruthSomaticReads << "\n";
        (*somaticReadLog) << "##Truth somatic read ratio: "  << (float)totalTaggedSomaticReads / (float)totalTruthSomaticReads << "\n";
        (*somaticReadLog) << "#CHROM\t"
                          << "ReadID\t"
                          << "germlineVarSimilarity\t"
                          << "deriveByHpSimilarity\t"
                          << "germlineSnpCount\t"
                          << "tumorSnpCount\t"
                          << "Haplotype\t"
                          << "somaticVariant,HP\n";
    }

    for(auto chr: chrVec){
        for(auto somaticRead: *somaticReadVecMap[chr]){
            (*somaticReadLog) << somaticRead.chr << "\t"
                            << somaticRead.readID << "\t"
                            << somaticRead.germlineVarSimilarity << "\t"
                            << somaticRead.deriveByHpSimilarity << "\t"
                            << somaticRead.germlineSnpCount << "\t"
                            << somaticRead.tumorSnpCount << "\t"
                            << "H" << somaticRead.hpResult << "\t";

            for(auto SnpHp: somaticRead.somaticSnpHp){
                (*somaticReadLog) << SnpHp.first+1 << "," << SnpHp.second << "\t";
            }
                (*somaticReadLog) << "\n";
        }
    }
    (*somaticReadLog).close();
    delete somaticReadLog;
    somaticReadLog = nullptr;
}

void SomaticReadBenchmark::displaySomaticVarCount(std::vector<std::string> &chrVec, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat){
    if(!openTestingFunc) return;
    int totalVariantCount = 0;
    
    for(auto &chr : chrVec){
        std::map<int, MultiGenomeVar>::iterator chrVariantIter = mergedChrVarinat[chr].begin();
        while(chrVariantIter != mergedChrVarinat[chr].end()){
            if(chrVariantIter->second.isExists(HIGH_CON_SOMATIC)){
                totalVariantCount++;
            }
            chrVariantIter++;
        }
    }

    int totalBedRegionCount = 0;
    for(auto &chr : chrVec){
        totalBedRegionCount += bedRegions[chr].size();
    }
    std::cout << "Total somatic variants: " << totalVariantCount << "\n";
    std::cout << "Total bed regions: " << totalBedRegionCount << "\n";
}