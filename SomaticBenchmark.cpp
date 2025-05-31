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

SomaticReadBenchmark::SomaticReadBenchmark(std::string benchmarkVcf, std::string benchmarkBed, int mappingQualityThreshold){
    setParseSnpFile(true);
    openTestingFunc = false;
    this->benchmarkVcf = benchmarkVcf;
    this->benchmarkBed = benchmarkBed;
    this->mappingQualityThreshold = mappingQualityThreshold;
}
SomaticReadBenchmark::~SomaticReadBenchmark(){

}

void SomaticReadBenchmark::setEnabled(bool openTestingFunc){
    this->openTestingFunc = openTestingFunc;
}

bool SomaticReadBenchmark::isEnabled(){
    return openTestingFunc;
}

bool SomaticReadBenchmark::isLoadBedFile(){
    return loadedBedFile;
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
        if( fields.size() == 0 ){
            return;
        }
            
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
    if(!openTestingFunc) return;

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
    loadedBedFile = true;
}

void SomaticReadBenchmark::processBedLine(const std::string& line) {
    std::istringstream iss(line);
    std::string chr;
    int start, end;
    iss >> chr >> start >> end;
    
    bedRegions[chr].push_back(BedRegion{start, end - 1});
}

void SomaticReadBenchmark::markVariantsInBedRegions(std::vector<std::string> &chrVec, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat) {
    if(!openTestingFunc || !loadedBedFile) return;
    for (auto& chr : chrVec) {
        auto& chrPosVariants = mergedChrVarinat[chr];
        
        // check if this chromosome has bed regions
        auto bedIt = bedRegions.find(chr);
        if (bedIt == bedRegions.end()) {
            // if this chromosome has no bed regions, mark all points as false
            for (auto& curVar : chrPosVariants) {
                curVar.second.isInBedRegion = false;
                if(curVar.second.isExists(Genome::TUMOR)) {
                    variantOutBedRegionCount[Genome::TUMOR]++;
                }
                if(curVar.second.isExists(Genome::NORMAL)) {
                    variantOutBedRegionCount[Genome::NORMAL]++;
                }
                if(curVar.second.isExists(Genome::HIGH_CON_SOMATIC)){
                    variantOutBedRegionCount[Genome::HIGH_CON_SOMATIC]++;
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
                    variantOutBedRegionCount[Genome::TUMOR]++;
                }
                if(curVar.second.isExists(Genome::NORMAL)) {
                    variantOutBedRegionCount[Genome::NORMAL]++;
                }
                if(curVar.second.isExists(Genome::HIGH_CON_SOMATIC)){
                    variantOutBedRegionCount[Genome::HIGH_CON_SOMATIC]++;
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
                    variantInBedRegionCount[Genome::TUMOR]++;
                }
                if(varPosIter->second.isExists(Genome::NORMAL)) {
                    variantInBedRegionCount[Genome::NORMAL]++;
                }
                if(varPosIter->second.isExists(Genome::HIGH_CON_SOMATIC)){
                    variantInBedRegionCount[Genome::HIGH_CON_SOMATIC]++;
                }
            } else {
                varPosIter->second.isInBedRegion = false;
                if(varPosIter->second.isExists(Genome::TUMOR)) {
                    variantOutBedRegionCount[Genome::TUMOR]++;
                }
                if(varPosIter->second.isExists(Genome::NORMAL)) {
                    variantOutBedRegionCount[Genome::NORMAL]++;
                }
                if(varPosIter->second.isExists(Genome::HIGH_CON_SOMATIC)){
                    variantOutBedRegionCount[Genome::HIGH_CON_SOMATIC]++;
                }
            }
            
            ++varPosIter;
        }
    }

}

void SomaticReadBenchmark::removeVariantsOutBedRegion(std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat){
    if(!openTestingFunc || !loadedBedFile) return;

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
    if(!openTestingFunc || !loadedBedFile) return;
    std::ofstream inBedLog(outPrefix + "_var_in_bed.out");
    std::ofstream outBedLog(outPrefix + "_var_out_bed.out");
    
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
    std::string outputFileName,
    std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat
){
    // if not open testing function, return
    if(!openTestingFunc) return;

    std::ofstream *refAltCountLog=NULL;
    refAltCountLog=new std::ofstream(outputFileName);
    int totalVariantCount = 0;

    for(auto chr: chrVec){
        totalVariantCount += chrMetrics[chr].posAltRefDelCount.size();
    }

    if(!refAltCountLog->is_open()){
        std::cerr<< "Fail to open write file: " << outputFileName << "\n";
        exit(1);
    }else{
        (*refAltCountLog) << "#############################\n";
        (*refAltCountLog) << "# Somatic SNP allele count #\n";
        (*refAltCountLog) << "#############################\n";
        (*refAltCountLog) << "##Benchmark VCF:"  << benchmarkVcf << "\n";
        (*refAltCountLog) << "##MappingQualityThreshold:"  << mappingQualityThreshold << "\n";
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

void SomaticReadBenchmark::writeTaggedReadReport(
    const std::vector<std::string>& chrVec,
    std::string outputFileName
){
    // if not open testing function, return
    if(!openTestingFunc) return;

    std::map<std::string, std::vector<SomaticReadLog>*> somaticReadVecMap;

    auto metricsIter = chrMetrics.begin();
    while(metricsIter != chrMetrics.end()){
        setChrSomaticReadVecPtr(metricsIter->first, somaticReadVecMap, metricsIter->second.totalReadVec);
        metricsIter++;
    }

    writeReadLog(chrVec, outputFileName, somaticReadVecMap);
}
void SomaticReadBenchmark::writeTaggedSomaticReadReport(
    const std::vector<std::string>& chrVec,
    std::string outputFileName
){
    // if not open testing function, return
    if(!openTestingFunc) return;

    time_t begin = time(NULL);
    std::cerr << "write somatic haplotag metrics report... ";

    std::map<std::string, std::vector<SomaticReadLog>*> somaticReadVecMap;

    auto metricsIter = chrMetrics.begin();
    while(metricsIter != chrMetrics.end()){
        setChrSomaticReadVecPtr(metricsIter->first, somaticReadVecMap, metricsIter->second.taggedSomaticReadVec);
        metricsIter++;
    }

    writeReadLog(chrVec, outputFileName, somaticReadVecMap);
    std::cerr<< difftime(time(NULL), begin) << "s\n";
}

void SomaticReadBenchmark::writeTotalTruthSomaticReadReport(
    const std::vector<std::string>& chrVec,
    std::string outputFileName
){
    // if not open testing function, return
    if(!openTestingFunc) return;


    std::map<std::string, std::vector<SomaticReadLog>*> somaticReadVecMap;

    auto metricsIter = chrMetrics.begin();
    while(metricsIter != chrMetrics.end()){
        setChrSomaticReadVecPtr(metricsIter->first, somaticReadVecMap, metricsIter->second.readsCrossingHighConSnpVec);
        metricsIter++;
    }

    writeReadLog(chrVec, outputFileName, somaticReadVecMap);

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
    std::string outputFileName,
    std::map<std::string, std::vector<SomaticReadLog>*> &somaticReadVecMap
){
    // if not open testing function, return
    if(!openTestingFunc) return;

    std::ofstream *somaticReadLog=NULL;
    somaticReadLog=new std::ofstream(outputFileName);

    int totalReads = 0;

    // all truth somatic reads
    int totalTruthSomaticReads = 0;
    std::map<std::string, int> truthSomaticReadsMap;

    // all tagged truth somatic reads(TP)
    int totalTaggedTruthSomaticReads = 0;
    std::map<std::string, int> taggedTruthSomaticReadsMap;

    // all tagged somatic reads(TP+FP)
    int totalTaggedSomaticReads = 0;
    std::map<std::string, int> toatlTaggedSomaticReadsMap;


    for(auto chr: chrVec){
        for(auto readIter :chrMetrics[chr].readsCrossingHighConSnpVec){
            if(readIter.hpResult == "1-1" || readIter.hpResult == "2-1" || readIter.hpResult == "3"){
                truthSomaticReadsMap[readIter.hpResult]++;
                totalTruthSomaticReads++;
            }
        }

        for(auto readIter :chrMetrics[chr].taggedSomaticReadVec){
            if(readIter.hpResult == "1-1" || readIter.hpResult == "2-1" || readIter.hpResult == "3"){
                totalTaggedTruthSomaticReads++;
                taggedTruthSomaticReadsMap[readIter.hpResult]++;
            }
        }
        for(auto readIter :chrMetrics[chr].totalReadVec){
            if(readIter.hpResult == "1-1" || readIter.hpResult == "2-1" || readIter.hpResult == "3"){
                toatlTaggedSomaticReadsMap[readIter.hpResult]++;
                totalTaggedSomaticReads++;
            }
            totalReads++;
        }
    }

    if(!somaticReadLog->is_open()){
        std::cerr<< "Fail to open write file: " << outputFileName << "\n";
        exit(1);
    }
    float recall = calculateRecall(totalTaggedTruthSomaticReads, totalTruthSomaticReads);;
    float precision = calculatePrecision(totalTaggedTruthSomaticReads, totalTaggedSomaticReads);
    float f1_score = calculateF1Score(recall, precision);


    (*somaticReadLog) << "############################\n";
    (*somaticReadLog) << "# Somatic Haplotag Metrics #\n";
    (*somaticReadLog) << "############################\n";
    (*somaticReadLog) << "##Benchmark VCF File: "  << benchmarkVcf << "\n";
    (*somaticReadLog) << "##Benchmark Bed File: "  << benchmarkBed << "\n";
    (*somaticReadLog) << "##MappingQualityThreshold: "  << mappingQualityThreshold << "\n";
    (*somaticReadLog) << "##Tatal reads: "  << totalReads << "\n";
    (*somaticReadLog) << "##Tatal truth somatic reads: "  << totalTruthSomaticReads << "\n";
    (*somaticReadLog) << "##Tatal truth H1-1: "  << truthSomaticReadsMap["1-1"] << "\n";
    (*somaticReadLog) << "##Tatal truth H2-1: "  << truthSomaticReadsMap["2-1"] << "\n";
    (*somaticReadLog) << "##Tatal truth H3: "  << truthSomaticReadsMap["3"] << "\n";
    
    int separatorLength = 95;
    int columnWidth = 15;

    (*somaticReadLog) << std::left << std::setw(columnWidth) << "## Haplotype" 
                    << std::setw(columnWidth) << "Precision" 
                    << std::setw(columnWidth) << "Recall" 
                    << std::setw(columnWidth) << "F1-Score"
                    << std::setw(columnWidth) << "TP" 
                    << std::setw(columnWidth) << "FP" 
                    << std::setw(columnWidth) << "FN" << "\n";
    (*somaticReadLog) << "##" << std::string(separatorLength, '-') << "\n";
    
    std::vector<std::string> haplotypes = {"1-1", "2-1", "3"};
    for(const auto& hp : haplotypes) {
        int tp = taggedTruthSomaticReadsMap[hp];
        int fp = toatlTaggedSomaticReadsMap[hp] - taggedTruthSomaticReadsMap[hp];
        int fn = truthSomaticReadsMap[hp] - taggedTruthSomaticReadsMap[hp];
        
        float precision = calculatePrecision(tp, tp + fp);
        float recall = calculateRecall(tp, tp + fn);
        float f1 = calculateF1Score(recall, precision);

        
        (*somaticReadLog) << std::left << std::setw(columnWidth) << ("## H" + hp)
                        << std::setw(columnWidth) << std::fixed << std::setprecision(4) << precision
                        << std::setw(columnWidth) << std::fixed << std::setprecision(4) << recall
                        << std::setw(columnWidth) << std::fixed << std::setprecision(4) << f1 
                        << std::setw(columnWidth) << tp
                        << std::setw(columnWidth) << fp  
                        << std::setw(columnWidth) << fn << "\n";
    }

    // overall statistics
    (*somaticReadLog) << "##" << std::string(separatorLength, '-') << "\n";
    (*somaticReadLog) << std::left << std::setw(columnWidth) << "## Overall"
                    << std::setw(columnWidth) << std::fixed << std::setprecision(4) << precision
                    << std::setw(columnWidth) << std::fixed << std::setprecision(4) << recall
                    << std::setw(columnWidth) << std::fixed << std::setprecision(4) << f1_score 
                    << std::setw(columnWidth) << totalTaggedTruthSomaticReads
                    << std::setw(columnWidth) << totalTaggedSomaticReads - totalTaggedTruthSomaticReads
                    << std::setw(columnWidth) << totalTruthSomaticReads - totalTaggedTruthSomaticReads << "\n";
    
    (*somaticReadLog) << "##\n";

    (*somaticReadLog) << "#CHROM\t"
                        << "ReadID\t"
                        << "germlineVarSimilarity\t"
                        << "deriveByHpSimilarity\t"
                        << "germlineSnpCount\t"
                        << "tumorSnpCount\t"
                        << "Haplotype\t"
                        << "somaticVariant,HP\n";

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
            if(chrVariantIter->second.isExists(Genome::HIGH_CON_SOMATIC)){
                totalVariantCount++;
            }
            chrVariantIter++;
        }
    }

    std::cout << "Total somatic variants: " << totalVariantCount << "\n";
}

void SomaticReadBenchmark::displayBedRegionCount(std::vector<std::string> &chrVec){
    if(!openTestingFunc || !loadedBedFile) return;
    int totalBedRegionCount = 0;
    for(auto &chr : chrVec){
        totalBedRegionCount += bedRegions[chr].size();
    }
    std::cout << "Total bed regions: " << totalBedRegionCount << "\n";
    std::cout << "--------------------------------" << "\n";
    std::cout << "Variant in bed region count: " << "\n";
    std::cout << "  -Tumor: " << variantInBedRegionCount[Genome::TUMOR] << "\n";
    std::cout << "  -Normal: " << variantInBedRegionCount[Genome::NORMAL] << "\n";
    std::cout << "  -Benchmark: " << variantInBedRegionCount[Genome::HIGH_CON_SOMATIC] << "\n";
    std::cout << "\n";
    std::cout << "Variant out bed region count: " << "\n";
    std::cout << "  -Tumor: " << variantOutBedRegionCount[Genome::TUMOR] << "\n";
    std::cout << "  -Normal: " << variantOutBedRegionCount[Genome::NORMAL] << "\n";
    std::cout << "  -Benchmark: " << variantOutBedRegionCount[Genome::HIGH_CON_SOMATIC] << "\n";
    std::cout << "--------------------------------" << "\n";
}