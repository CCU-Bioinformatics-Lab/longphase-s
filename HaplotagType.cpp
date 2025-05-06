#include "HaplotagType.h"

void chrReadHpResult::recordReadHp(int &pos, int &hpResult, int &BaseHP){
    posReadHpResult[pos].readHpCounter[hpResult]++;
    
    if(hpResult != ReadHP::unTag){
        if(BaseHP == SnpHP::SOMATIC_H3){
            if(hpResult != ReadHP::H1_1 && hpResult != ReadHP::H2_1 && hpResult != ReadHP::H3){
                std::cerr << "Error(recordReadHp) => error read hp : BaseHP: " <<BaseHP << " readHP: " << hpResult << " pos: " << pos+1 << std::endl; 
                exit(1);
            }
            posReadHpResult[pos].somaticSnpH3count++;
            posReadHpResult[pos].somaticBaseReadHpCounter[hpResult]++;
        }
    }
}

void chrReadHpResult::recordDeriveHp(int &pos, int &deriveHP, float deriveHPsimilarity){
    if(deriveHP != SnpHP::GERMLINE_H1 && deriveHP != SnpHP::GERMLINE_H2 && deriveHP != SnpHP::NONE_SNP){
        std::cerr << "Error(recordDeriveHp) => error derive hp : pos: " <<pos+1 << " deriveHP: " << deriveHP << std::endl; 
        exit(1);        
    }
    posReadHpResult[pos].deriveHP = deriveHP;
    if(deriveHPsimilarity != 0.0){
        posReadHpResult[pos].deriveHPsimilarVec.emplace_back(deriveHPsimilarity);
        if(deriveHPsimilarity != 1.0){
            // std::cout << "deriveHPsimilarity: " << deriveHPsimilarity << "\n";
            // std::cout << "deriveHPsimilarityVec: " << varReadHpResult[pos].deriveHPsimilarVec.back() << "\n";
        }
    }
}

void chrReadHpResult::recordAlignCoverRegion(int& curVarPos, int &startPos, int &endPos){
    if(posReadHpResult[curVarPos].coverRegionStartPos > startPos){
        posReadHpResult[curVarPos].coverRegionStartPos = startPos;
    }
    if(posReadHpResult[curVarPos].coverRegionEndPos < endPos){
        posReadHpResult[curVarPos].coverRegionEndPos = endPos;
    }
}


