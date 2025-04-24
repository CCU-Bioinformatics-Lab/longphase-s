#include "SomaticBenchmark.h"

SomaticReadVerifier::SomaticReadVerifier(){
    setParseSnpFile(true);
    openTestingFunc = false;
}
SomaticReadVerifier::~SomaticReadVerifier(){

}

void SomaticReadVerifier::setTestingFunc(bool openTestingFunc){
    this->openTestingFunc = openTestingFunc;
}

void SomaticReadVerifier::loadHighConSomatic(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat){
    if(!openTestingFunc) return;
    variantParser(input, Info, mergedChrVarinat);
}

void SomaticReadVerifier::parserProcess(std::string &input, VCF_Info &Info, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat){
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

void SomaticReadVerifier::recordDelReadCount(const std::string &chr, std::map<int, MultiGenomeVar>::iterator &currentVariantIter){
    if(!openTestingFunc) return;
    
    if(currentVariantIter->second.isExists(HIGH_CON_SOMATIC)){
        int pos = (*currentVariantIter).first;
        posAltRefDelCount[chr][pos].delCount++;

        //record somatic position for record crossing high con snp read
        highConSomaticPos.push_back(std::make_pair(pos, SnpHP::NONE_SNP));
    }
}

void SomaticReadVerifier::recordRefAltAlleleCount(const std::string &chr, std::string &base, std::map<int, MultiGenomeVar>::iterator &currentVariantIter){
    if(!openTestingFunc) return;

    if(currentVariantIter->second.isExists(HIGH_CON_SOMATIC)){
        int pos = currentVariantIter->first;
        std::string refAllele = currentVariantIter->second.Variant[HIGH_CON_SOMATIC].allele.Ref;
        std::string altAllele = currentVariantIter->second.Variant[HIGH_CON_SOMATIC].allele.Alt;

        int baseHP = SnpHP::NONE_SNP;

        if(base == refAllele){
            posAltRefDelCount[chr][pos].refCount++;
        }else if(base == altAllele){
            posAltRefDelCount[chr][pos].altCount++;
            baseHP = SnpHP::SOMATIC_H3;
        }
        //record somatic position for record crossing high con snp read
        highConSomaticPos.push_back(std::make_pair(pos, baseHP));
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

    for(auto varIter : highConSomaticPos){
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
        readsCrossingHighConSnpVec.push_back(tmp);
    }
    
    //clear high con somatic position in current read for next read
    if(!highConSomaticPos.empty()){
        highConSomaticPos.clear();
    }
}

void SomaticReadVerifier::recordTaggedSomaticRead(const std::string &chr, std::string &readID, std::string &hpResult, std::map<int, int> &variantsHP, std::map<int, int> &hpCount, double &norHPsimilarity, float &deriveByHpSimilarity, std::map<int, MultiGenomeVar> &currentChrVariants){
    // if not open testing function, return
    if(!openTestingFunc) return;

    SomaticReadLog tmp = createBasicSomaticReadLog(chr, readID, hpResult, norHPsimilarity, deriveByHpSimilarity, hpCount);

    bool readExistHighConSomatic = false;

    auto varIter = variantsHP.begin();
    while(varIter != variantsHP.end()){
        int pos = varIter->first;
        int snpHP = varIter->second;
        if(currentChrVariants.find(pos) != currentChrVariants.end()){
            if(currentChrVariants[pos].isExists(HIGH_CON_SOMATIC) && (snpHP == SnpHP::SOMATIC_H3 || snpHP == SnpHP::SOMATIC_H4)){
                tmp.somaticSnpHp[pos] = snpHP;
                readExistHighConSomatic = true;
            }
        }
        varIter++;
    }

    if(readExistHighConSomatic){
        taggedSomaticReadVec.push_back(tmp);
    }
}

void SomaticReadVerifier::recordTaggedRead(const std::string &chr, std::string &readID, std::string &hpResult, std::map<int, int> &variantsHP, std::map<int, int> &hpCount, double &norHPsimilarity, float &deriveByHpSimilarity, std::map<int, MultiGenomeVar> &currentChrVariants){
    // if not open testing function, return
    if(!openTestingFunc) return;
    
    if(hpResult != "."){
        SomaticReadLog tmp = createBasicSomaticReadLog(chr, readID, hpResult, norHPsimilarity, deriveByHpSimilarity, hpCount);

        auto varIter = variantsHP.begin();
        while(varIter != variantsHP.end()){
            int pos = varIter->first;
            int snpHP = varIter->second;
            if(currentChrVariants.find(pos) != currentChrVariants.end()){
                if(currentChrVariants[pos].isExists(HIGH_CON_SOMATIC) && (snpHP == SnpHP::SOMATIC_H3 || snpHP == SnpHP::SOMATIC_H4)){
                    tmp.somaticSnpHp[pos] = snpHP;
                }
            }
            varIter++;
        }

        totalReadVec.push_back(tmp);
    }
}

void SomaticReadVerifier::writePosAlleleCountLog(std::vector<std::string> &chrVec, HaplotagParameters &params, std::string logPosfix, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat){
    // if not open testing function, return
    if(!openTestingFunc) return;

    std::ofstream *refAltCountLog=NULL;
    refAltCountLog=new std::ofstream(params.resultPrefix + logPosfix);
    int totalVariantCount = 0;

    for(auto &chr : chrVec){
        totalVariantCount += posAltRefDelCount[chr].size();
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

    for(auto &chr : chrVec){
        for(auto &posIter : posAltRefDelCount[chr]){
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

void SomaticReadVerifier::writeTaggedReadLog(HaplotagParameters &params, std::string logPosfix){
    // if not open testing function, return
    if(!openTestingFunc) return;
    writeReadLog(params, logPosfix, totalReadVec);
}
void SomaticReadVerifier::writeTaggedSomaticReadLog(HaplotagParameters &params, std::string logPosfix){
    // if not open testing function, return
    if(!openTestingFunc) return;
    writeReadLog(params, logPosfix, taggedSomaticReadVec);
}

void SomaticReadVerifier::writeCrossHighConSnpReadLog(HaplotagParameters &params, std::string logPosfix){
    // if not open testing function, return
    if(!openTestingFunc) return;
    writeReadLog(params, logPosfix, readsCrossingHighConSnpVec);
}

void SomaticReadVerifier::writeReadLog(HaplotagParameters &params, std::string logPosfix, std::vector<SomaticReadLog> &somaticReadVec){
    // if not open testing function, return
    if(!openTestingFunc) return;

    std::ofstream *somaticReadLog=NULL;
    somaticReadLog=new std::ofstream(params.resultPrefix + logPosfix);

    int totalTruthSomaticReads = 0;
    for(auto readIter :readsCrossingHighConSnpVec){
        if(readIter.hpResult == "1-1" || readIter.hpResult == "2-1" || readIter.hpResult == "3"){
            totalTruthSomaticReads++;
        }
    }

    int totalTaggedSomaticReads = 0;
    for(auto readIter :taggedSomaticReadVec){
        if(readIter.hpResult == "1-1" || readIter.hpResult == "2-1" || readIter.hpResult == "3"){
            totalTaggedSomaticReads++;
        }
    }

    int totalReads = totalReadVec.size();

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

    for(auto somaticRead: somaticReadVec){

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
    (*somaticReadLog).close();
    delete somaticReadLog;
    somaticReadLog = nullptr;
}

void SomaticReadVerifier::displaySomaticVarCount(std::vector<std::string> &chrVec, std::map<std::string, std::map<int, MultiGenomeVar>> &mergedChrVarinat){
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
    std::cout << "Total somatic variants: " << totalVariantCount << "\n";
}