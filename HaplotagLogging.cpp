#include "HaplotagLogging.h"

ReadHpDistriLog::ReadHpDistriLog(){

}

ReadHpDistriLog::~ReadHpDistriLog(){

}

void ReadHpDistriLog::loadChrKey(const std::string &chr){
    chrVarReadHpResult[chr] = chrReadHpResult();
}

chrReadHpResult* ReadHpDistriLog::getChrHpResultsPtr (const std::string &chr){
    return &(chrVarReadHpResult[chr]);
}

void ReadHpDistriLog::recordChrReadHp(const std::string &chr, int &pos, int &hpResult, int &BaseHP){
    //only use in single thread scenario
    chrVarReadHpResult[chr].recordReadHp(pos, hpResult, BaseHP);
}

void ReadHpDistriLog::recordChrDeriveHp(const std::string &chr, int &pos, int &deriveHP, float deriveHPsimilarity){
    //only use in single thread scenario
    chrVarReadHpResult[chr].recordDeriveHp(pos, deriveHP, deriveHPsimilarity);
}

void ReadHpDistriLog::recordChrAlignCoverRegion(const std::string &chr, int &pos, int &startPos, int &endPos){
    //only use in single thread scenario
    chrVarReadHpResult[chr].recordAlignCoverRegion(pos, startPos, endPos);
}


void ReadHpDistriLog::writeReadHpDistriLog(const HaplotagParameters &params, std::string logPosfix, const std::vector<std::string> &chrVec){
    std::ofstream *readHpDistriLog=NULL;
    readHpDistriLog=new std::ofstream(params.resultPrefix + logPosfix);

    int somaticSnpCount = 0;
    for(auto chr: chrVec){
        if(!chrVarReadHpResult[chr].posReadHpResult.empty()){
            somaticSnpCount += chrVarReadHpResult[chr].posReadHpResult.size();
        }
    }

    if(!readHpDistriLog->is_open()){
        std::cerr<< "Fail to open write file: " << params.resultPrefix + logPosfix << "\n";
        exit(1);
    }else{
        (*readHpDistriLog) << "####################################\n";
        (*readHpDistriLog) << "# Somatic SNP read HP distribution #\n";
        (*readHpDistriLog) << "####################################\n";
        (*readHpDistriLog) << "##MappingQualityThreshold:"        << params.qualityThreshold << "\n";
        (*readHpDistriLog) << "##SomaticSNP: " << somaticSnpCount << "\n";
        (*readHpDistriLog) << "#Chr\t"
                            << "Pos\t"
                            << "DeriveHP\t"
                            << "DeriveHPsimilarity\t\t"
                            << "AltCount\t"
                            << "somaticBase_H1-1\t"
                            << "somaticBase_H2-1\t"
                            << "somaticBase_H3\t\t"
                            << "HP1read\t"
                            << "HP2read\t"
                            << "HP1-1read\t"
                            << "HP2-1read\t"
                            << "HP3read\t"
                            << "untagRead\t"
                            << "HP1ratio\t"
                            << "HP2ratio\t"
                            << "HP1-1ratio\t"
                            << "HP2-1ratio\t"
                            << "HP3ratio\n";
    }

    for(auto chr: chrVec){
        std::map<int, ReadHpResult>::iterator curVarReadHpIter = chrVarReadHpResult[chr].posReadHpResult.begin();
        while(curVarReadHpIter != chrVarReadHpResult[chr].posReadHpResult.end()){
            int pos = (*curVarReadHpIter).first + 1;
            int HP1readCount = (*curVarReadHpIter).second.readHpCounter[ReadHP::H1];
            int HP1_1readCount = (*curVarReadHpIter).second.readHpCounter[ReadHP::H1_1];

            int HP2readCount = (*curVarReadHpIter).second.readHpCounter[ReadHP::H2];
            int HP2_1readCount = (*curVarReadHpIter).second.readHpCounter[ReadHP::H2_1];

            int HP3readCount = (*curVarReadHpIter).second.readHpCounter[ReadHP::H3];

            int totaltagRead = HP1readCount + HP2readCount + HP3readCount + HP1_1readCount + HP2_1readCount;

            float HP1readRatio = (float)HP1readCount / (float)totaltagRead; 
            float HP1_1readRatio = (float)HP1_1readCount / (float)totaltagRead; 

            float HP2readRatio = (float)HP2readCount / (float)totaltagRead; 
            float HP2_1readRatio = (float)HP2_1readCount / (float)totaltagRead; 

            float HP3readRatio = (float)HP3readCount / (float)totaltagRead; 

            float meanDeriveHPsimilarity = 0.0;
            if((*curVarReadHpIter).second.deriveHPsimilarVec.size() != 0){
                float size = (*curVarReadHpIter).second.deriveHPsimilarVec.size();
                for(float deriveHPsimilarity: (*curVarReadHpIter).second.deriveHPsimilarVec){
                    meanDeriveHPsimilarity += deriveHPsimilarity;
                    // std::cout << deriveHPsimilarity << "\t";
                }
                meanDeriveHPsimilarity /= size;
            }
            (*readHpDistriLog) << std::fixed << std::setprecision(3) 
                                << chr << "\t"
                                << pos << "\t"
                                << "H" << (*curVarReadHpIter).second.deriveHP << "\t"
                                << meanDeriveHPsimilarity << "\t\t"
                                << (*curVarReadHpIter).second.somaticSnpH3count << "\t"
                                << (*curVarReadHpIter).second.somaticBaseReadHpCounter[ReadHP::H1_1] << "\t"
                                << (*curVarReadHpIter).second.somaticBaseReadHpCounter[ReadHP::H2_1] << "\t"
                                << (*curVarReadHpIter).second.somaticBaseReadHpCounter[ReadHP::H3] << "\t\t"
                                << HP1readCount << "\t"
                                << HP2readCount << "\t\t"
                                << HP1_1readCount << "\t"
                                << HP2_1readCount << "\t"
                                << HP3readCount << "\t"
                                << (*curVarReadHpIter).second.readHpCounter[ReadHP::unTag] << "\t"
                                << HP1readRatio << "\t"
                                << HP2readRatio << "\t"
                                << HP1_1readRatio << "\t"
                                << HP2_1readRatio << "\t"
                                << HP3readRatio << "\n";
            curVarReadHpIter++;
        }
    }
    (*readHpDistriLog).close();
    delete readHpDistriLog;
    readHpDistriLog = nullptr;
}

void ReadHpDistriLog::writePosCoverRegionLog(const HaplotagParameters &params, std::string logPosfix, const std::vector<std::string> &chrVec){
    std::ofstream *posCoverRegionLog=NULL;
    posCoverRegionLog=new std::ofstream(params.resultPrefix + logPosfix);

    int somaticSnpCount = 0;
    for(auto chr: chrVec){
        if(!chrVarReadHpResult[chr].posReadHpResult.empty()){
            somaticSnpCount += chrVarReadHpResult[chr].posReadHpResult.size();
        }
    }

    if(!posCoverRegionLog->is_open()){
        std::cerr<< "Fail to open write file: " << params.resultPrefix + logPosfix << "\n";
        exit(1);
    }else{
        (*posCoverRegionLog) << "############################\n";
        (*posCoverRegionLog) << "# Somatic SNP cover region #\n";
        (*posCoverRegionLog) << "############################\n";
        (*posCoverRegionLog) << "##MappingQualityThreshold:"        << params.qualityThreshold << "\n";
        (*posCoverRegionLog) << "##SomaticSNP: " << somaticSnpCount << "\n";
        (*posCoverRegionLog) << "#Chr\t"
                              << "Pos\t"
                              << "Type\t"
                              << "StartPos\t"
                              << "EndPos\n";
    }

    for(auto chr: chrVec){
        std::map<int, ReadHpResult>::iterator curVarReadHpIter = chrVarReadHpResult[chr].posReadHpResult.begin();
        while(curVarReadHpIter != chrVarReadHpResult[chr].posReadHpResult.end()){
            int pos = (*curVarReadHpIter).first + 1;

            (*posCoverRegionLog) << std::fixed << std::setprecision(3) 
                                << chr << "\t"
                                << pos << "\t"
                                << "somatic" << "\t"
                                << (*curVarReadHpIter).second.coverRegionStartPos << "\t"
                                << (*curVarReadHpIter).second.coverRegionEndPos << "\n";
            curVarReadHpIter++;
        }
    }
    (*posCoverRegionLog).close();
    delete posCoverRegionLog;
    posCoverRegionLog = nullptr;
}

void ReadHpDistriLog::writeTagReadCoverRegionLog(const HaplotagParameters &params, std::string logPosfix, const std::vector<std::string> &chrVec, std::map<std::string, int> &chrLength){
    std::ofstream *tagReadCoverRegionLog=NULL;
    tagReadCoverRegionLog=new std::ofstream(params.resultPrefix + logPosfix);

    std::map<std::string, std::vector<coverRegionInfo>> coverRegion;

    // merge cover region from different SNPs
    for (const auto& chr : chrVec){
        auto curVarReadHpIter = chrVarReadHpResult[chr].posReadHpResult.begin();

        if (curVarReadHpIter == chrVarReadHpResult[chr].posReadHpResult.end()) {
            continue;
        }

        int curStartPos = curVarReadHpIter->second.coverRegionStartPos;
        int curEndPos = curVarReadHpIter->second.coverRegionEndPos;

        while (curVarReadHpIter != chrVarReadHpResult[chr].posReadHpResult.end()){

            auto nextVarReadHpIter = std::next(curVarReadHpIter);
            if(nextVarReadHpIter != chrVarReadHpResult[chr].posReadHpResult.end()){
                int nextStartPos = nextVarReadHpIter->second.coverRegionStartPos;
                int nextEndPos = nextVarReadHpIter->second.coverRegionEndPos;

                // region not overlap
                if (curEndPos < nextStartPos) {
                    coverRegionInfo regionInfo;
                    regionInfo.startPos = curStartPos;
                    regionInfo.endPos = curEndPos;
                    regionInfo.length = curEndPos - curStartPos + 1;
                    coverRegion[chr].emplace_back(regionInfo);
                    curStartPos = nextStartPos;
                    curEndPos = nextEndPos;
                }else{
                    // regions overlap
                    curStartPos = std::min(curStartPos, nextStartPos);
                    curEndPos = std::max(curEndPos, nextEndPos);
                }
            // last region
            }else{
                coverRegionInfo regionInfo;
                regionInfo.startPos = curStartPos;
                regionInfo.endPos = curEndPos;
                regionInfo.length = curEndPos - curStartPos + 1;
                coverRegion[chr].emplace_back(regionInfo);
                break;
            }

            curVarReadHpIter++;
        }
    }

    std::map<std::string, float> coverRegionRatio;
    long long totalChrLength = 0;
    long long totalChrCoverLength = 0;
    double totalChrCoverageRatio = 0.0;

    // calculate cover region ratio
    for(auto chr: chrVec){
        int totalCoverLength = 0;
        auto coverRegionIter = coverRegion[chr].begin();
        while(coverRegionIter != coverRegion[chr].end()){
            totalCoverLength += (*coverRegionIter).length;
            coverRegionIter++;
        }
        coverRegionRatio[chr] = (float)totalCoverLength / (float)chrLength[chr];
        totalChrLength += chrLength[chr];
        totalChrCoverLength += totalCoverLength;
    }
    totalChrCoverageRatio = (double)totalChrCoverLength / (double)totalChrLength;

    if(!tagReadCoverRegionLog->is_open()){
        std::cerr<< "Fail to open write file: " << params.resultPrefix + logPosfix << "\n";
        exit(1);
    }else{
        (*tagReadCoverRegionLog) << "##################################\n";
        (*tagReadCoverRegionLog) << "# Somatic reads cover region bed #\n";
        (*tagReadCoverRegionLog) << "##################################\n";
        (*tagReadCoverRegionLog) << "##MappingQualityThreshold: "   << params.qualityThreshold << "\n";
        (*tagReadCoverRegionLog) << "##----Chr coverage ratio----\n";
        (*tagReadCoverRegionLog) << "##Total chr coverage ratio: " << totalChrCoverageRatio << "\n";
        for(auto chr: chrVec){
            (*tagReadCoverRegionLog) <<"##" << chr << ":" << coverRegionRatio[chr] << "\n";
        }
        (*tagReadCoverRegionLog) << "#Chr\t"
                                 << "StartPos\t"
                                 << "EndPos\n";
                                // << "length\n";
    }

    for(auto chr: chrVec){
        auto coverRegionIter = coverRegion[chr].begin();
        while(coverRegionIter != coverRegion[chr].end()){

            (*tagReadCoverRegionLog) << std::fixed << std::setprecision(3) 
                                     << chr << "\t"
                                     << (*coverRegionIter).startPos << "\t"
                                     << (*coverRegionIter).endPos << "\n";
                                     //<< (*coverRegionIter).length << "\n";
            coverRegionIter++;
        }
    }
    (*tagReadCoverRegionLog).close();
    delete tagReadCoverRegionLog;
    tagReadCoverRegionLog = nullptr;
}


void ReadHpDistriLog::removeNotDeriveByH1andH2pos(const std::vector<std::string> &chrVec){
    for (const auto& chr : chrVec) {
        auto& chrResult = chrVarReadHpResult[chr];
        for (auto it = chrResult.posReadHpResult.begin(); it != chrResult.posReadHpResult.end(); ) {
            if (!it->second.existDeriveByH1andH2) {
                //std::cerr << "Removed position not derived by H1 and H2: " << chr << " " << it->first << std::endl;
                it = chrResult.posReadHpResult.erase(it);
            } else {
                ++it;
            }
        }
    }
}
