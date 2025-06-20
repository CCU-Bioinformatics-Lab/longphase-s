#include "HaplotagLogging.h"

/**
 * @brief Record read haplotype assignment for a specific position
 * 
 * Updates haplotype counters and performs validation for somatic variants.
 * For somatic H3 variants, validates that read haplotypes are valid (H1_1, H2_1, or H3).
 * 
 * @param pos Variant position (0-based)
 * @param hpResult Assigned haplotype for the read
 * @param BaseHP Base haplotype type (germline or somatic)
 */
void chrReadHpResult::recordReadHp(int &pos, int &hpResult, int &BaseHP){
    posReadHpResult[pos].readHpCounter[hpResult]++;
    
    if(hpResult != ReadHP::unTag){
        if(BaseHP == SnpHP::SOMATIC_H3){
            // Validate read haplotype for somatic H3 variants
            if(hpResult != ReadHP::H1_1 && hpResult != ReadHP::H2_1 && hpResult != ReadHP::H3){
                std::cerr << "[ERROR](recordReadHp) => error read hp : BaseHP: " <<BaseHP << " readHP: " << hpResult << " pos: " << pos+1 << std::endl; 
                exit(1);
            }
            posReadHpResult[pos].somaticSnpH3count++;
            posReadHpResult[pos].somaticBaseReadHpCounter[hpResult]++;
        }
    }
}

/**
 * @brief Record derived haplotype information for a position
 * 
 * Validates derived haplotype values and stores similarity scores.
 * Only accepts valid derived haplotypes: GERMLINE_H1, GERMLINE_H2, or NONE_SNP.
 * 
 * @param pos Variant position (0-based)
 * @param deriveHP Derived haplotype (H1, H2, or NONE)
 * @param deriveHPsimilarity Similarity score for derived haplotype
 */
void chrReadHpResult::recordDeriveHp(int &pos, int &deriveHP, float deriveHPsimilarity){
    // Validate derived haplotype values
    if(deriveHP != SnpHP::GERMLINE_H1 && deriveHP != SnpHP::GERMLINE_H2 && deriveHP != SnpHP::NONE_SNP){
        std::cerr << "[ERROR](recordDeriveHp) => error derive hp : pos: " <<pos+1 << " deriveHP: " << deriveHP << std::endl; 
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

/**
 * @brief Record alignment coverage region for a variant position
 * 
 * Updates the coverage region boundaries for a variant position.
 * Expands the region if the new start/end positions extend beyond current boundaries.
 * 
 * @param curVarPos Current variant position
 * @param startPos Start position of coverage region
 * @param endPos End position of coverage region
 */
void chrReadHpResult::recordAlignCoverRegion(int& curVarPos, int &startPos, int &endPos){
    if(posReadHpResult[curVarPos].coverRegionStartPos > startPos){
        posReadHpResult[curVarPos].coverRegionStartPos = startPos;
    }
    if(posReadHpResult[curVarPos].coverRegionEndPos < endPos){
        posReadHpResult[curVarPos].coverRegionEndPos = endPos;
    }
}

ReadHpDistriLog::ReadHpDistriLog(){

}

ReadHpDistriLog::~ReadHpDistriLog(){

}

/**
 * @brief Initialize chromosome key for multi-threaded access
 * 
 * Creates an empty chrReadHpResult entry for the specified chromosome.
 * This method should be called before multi-threaded processing to ensure
 * thread-safe access to chromosome data.
 * 
 * @param chr Chromosome name
 */
void ReadHpDistriLog::loadChrKey(const std::string &chr){
    chrVarReadHpResult[chr] = chrReadHpResult();
}

/**
 * @brief Get thread-safe pointer to chromosome haplotype results
 * 
 * Returns a pointer to the chromosome's haplotype results for safe
 * multi-threaded access. The pointer remains valid throughout the
 * chromosome's processing.
 * 
 * @param chr Chromosome name
 * @return Pointer to chromosome haplotype results
 */
chrReadHpResult* ReadHpDistriLog::getChrHpResultsPtr (const std::string &chr){
    return &(chrVarReadHpResult[chr]);
}

/**
 * @brief Record read haplotype for a chromosome (single-thread use only)
 * 
 * Delegates to chrReadHpResult::recordReadHp for the specified chromosome.
 * This method is not thread-safe and should only be used in single-threaded scenarios.
 * 
 * @param chr Chromosome name
 * @param pos Variant position (0-based)
 * @param hpResult Assigned haplotype for the read
 * @param BaseHP Base haplotype type
 */
void ReadHpDistriLog::recordChrReadHp(const std::string &chr, int &pos, int &hpResult, int &BaseHP){
    //only use in single thread scenario
    chrVarReadHpResult[chr].recordReadHp(pos, hpResult, BaseHP);
}

/**
 * @brief Record derived haplotype for a chromosome (single-thread use only)
 * 
 * Delegates to chrReadHpResult::recordDeriveHp for the specified chromosome.
 * This method is not thread-safe and should only be used in single-threaded scenarios.
 * 
 * @param chr Chromosome name
 * @param pos Variant position (0-based)
 * @param deriveHP Derived haplotype
 * @param deriveHPsimilarity Similarity score
 */
void ReadHpDistriLog::recordChrDeriveHp(const std::string &chr, int &pos, int &deriveHP, float deriveHPsimilarity){
    //only use in single thread scenario
    chrVarReadHpResult[chr].recordDeriveHp(pos, deriveHP, deriveHPsimilarity);
}

/**
 * @brief Record alignment coverage region for a chromosome (single-thread use only)
 * 
 * Delegates to chrReadHpResult::recordAlignCoverRegion for the specified chromosome.
 * This method is not thread-safe and should only be used in single-threaded scenarios.
 * 
 * @param chr Chromosome name
 * @param pos Variant position (0-based)
 * @param startPos Start position of coverage region
 * @param endPos End position of coverage region
 */
void ReadHpDistriLog::recordChrAlignCoverRegion(const std::string &chr, int &pos, int &startPos, int &endPos){
    //only use in single thread scenario
    chrVarReadHpResult[chr].recordAlignCoverRegion(pos, startPos, endPos);
}

/**
 * @brief Write comprehensive read haplotype distribution report
 * 
 * Generates a detailed report containing haplotype distribution statistics
 * for each somatic SNP position across all specified chromosomes.
 * 
 * Report includes:
 * - Total somatic SNP count
 * - Per-position haplotype read counts and ratios
 * - Derived haplotype information and similarity scores
 * - Somatic variant read counts
 * 
 * @param logFileName Output file name
 * @param chrVec Vector of chromosome names to process
 */
void ReadHpDistriLog::writeReadHpDistriLog(const std::string logFileName, const std::vector<std::string> &chrVec){
    std::ofstream *readHpDistriLog=NULL;
    readHpDistriLog=new std::ofstream(logFileName);

    // Calculate total somatic SNP count across all chromosomes
    int somaticSnpCount = 0;
    for(auto chr: chrVec){
        if(!chrVarReadHpResult[chr].posReadHpResult.empty()){
            somaticSnpCount += chrVarReadHpResult[chr].posReadHpResult.size();
        }
    }

    if(!readHpDistriLog->is_open()){
        std::cerr<< "Fail to open write file: " << logFileName << "\n";
        exit(1);
    }else{
        (*readHpDistriLog) << "###################################################\n";
        (*readHpDistriLog) << "# Distribution of Read Haplotypes at Somatic SNPs #\n";
        (*readHpDistriLog) << "###################################################\n";
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

    // Process each chromosome and variant position
    for(auto chr: chrVec){
        std::map<int, ReadHpResult>::iterator curVarReadHpIter = chrVarReadHpResult[chr].posReadHpResult.begin();
        while(curVarReadHpIter != chrVarReadHpResult[chr].posReadHpResult.end()){
            int pos = (*curVarReadHpIter).first + 1;  // Convert to 1-based position for output
            
            // Extract read counts for different haplotype types
            int HP1readCount = (*curVarReadHpIter).second.readHpCounter[ReadHP::H1];
            int HP1_1readCount = (*curVarReadHpIter).second.readHpCounter[ReadHP::H1_1];

            int HP2readCount = (*curVarReadHpIter).second.readHpCounter[ReadHP::H2];
            int HP2_1readCount = (*curVarReadHpIter).second.readHpCounter[ReadHP::H2_1];

            int HP3readCount = (*curVarReadHpIter).second.readHpCounter[ReadHP::H3];

            int totaltagRead = HP1readCount + HP2readCount + HP3readCount + HP1_1readCount + HP2_1readCount;

            // Calculate haplotype ratios
            float HP1readRatio = (float)HP1readCount / (float)totaltagRead; 
            float HP1_1readRatio = (float)HP1_1readCount / (float)totaltagRead; 

            float HP2readRatio = (float)HP2readCount / (float)totaltagRead; 
            float HP2_1readRatio = (float)HP2_1readCount / (float)totaltagRead; 

            float HP3readRatio = (float)HP3readCount / (float)totaltagRead; 

            // Calculate mean derived haplotype similarity
            float meanDeriveHPsimilarity = 0.0;
            if((*curVarReadHpIter).second.deriveHPsimilarVec.size() != 0){
                float size = (*curVarReadHpIter).second.deriveHPsimilarVec.size();
                for(float deriveHPsimilarity: (*curVarReadHpIter).second.deriveHPsimilarVec){
                    meanDeriveHPsimilarity += deriveHPsimilarity;
                    // std::cout << deriveHPsimilarity << "\t";
                }
                meanDeriveHPsimilarity /= size;
            }
            
            // Write detailed statistics for this position
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

/**
 * @brief Write position coverage region report
 * 
 * Generates a report showing coverage regions for each somatic SNP position.
 * This report is useful for understanding the genomic regions covered by
 * reads that support each somatic variant.
 * 
 * @param logFileName Output file name
 * @param chrVec Vector of chromosome names to process
 */
void ReadHpDistriLog::writePosCoverRegionLog(const std::string logFileName, const std::vector<std::string> &chrVec){
    std::ofstream *posCoverRegionLog=NULL;
    posCoverRegionLog=new std::ofstream(logFileName);

    // Calculate total somatic SNP count
    int somaticSnpCount = 0;
    for(auto chr: chrVec){
        if(!chrVarReadHpResult[chr].posReadHpResult.empty()){
            somaticSnpCount += chrVarReadHpResult[chr].posReadHpResult.size();
        }
    }

    if(!posCoverRegionLog->is_open()){
        std::cerr<< "Fail to open write file: " << logFileName << "\n";
        exit(1);
    }else{
        (*posCoverRegionLog) << "############################\n";
        (*posCoverRegionLog) << "# Somatic SNP cover region #\n";
        (*posCoverRegionLog) << "############################\n";
        (*posCoverRegionLog) << "##SomaticSNP: " << somaticSnpCount << "\n";
        (*posCoverRegionLog) << "#Chr\t"
                              << "Pos\t"
                              << "Type\t"
                              << "StartPos\t"
                              << "EndPos\n";
    }

    // Write coverage region for each variant position
    for(auto chr: chrVec){
        std::map<int, ReadHpResult>::iterator curVarReadHpIter = chrVarReadHpResult[chr].posReadHpResult.begin();
        while(curVarReadHpIter != chrVarReadHpResult[chr].posReadHpResult.end()){
            int pos = (*curVarReadHpIter).first + 1;  // Convert to 1-based position

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

/**
 * @brief Write tagged read coverage region report with merged regions
 * 
 * Generates a comprehensive report showing merged coverage regions across
 * all somatic SNPs. This report includes:
 * - Coverage ratios for each chromosome
 * - Total genome coverage ratio
 * - Merged BED-format regions
 * 
 * The merging algorithm combines overlapping regions to create continuous
 * coverage blocks, which is useful for understanding the overall genomic
 * coverage of somatic variant analysis.
 * 
 * TODO: [PENDING] This function needs review and validation
 * 
 * @param logFileName Output file name
 * @param chrVec Vector of chromosome names to process
 * @param chrLength Map of chromosome names to their lengths
 */
void ReadHpDistriLog::writeTagReadCoverRegionLog(const std::string logFileName, const std::vector<std::string> &chrVec, std::map<std::string, int> &chrLength){
    std::ofstream *tagReadCoverRegionLog=NULL;
    tagReadCoverRegionLog=new std::ofstream(logFileName);

    std::map<std::string, std::vector<coverRegionInfo>> coverRegion;

    // Merge coverage regions from different SNPs
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

                // Check if regions overlap
                if (curEndPos < nextStartPos) {
                    // Regions do not overlap - save current region and start new one
                    coverRegionInfo regionInfo;
                    regionInfo.startPos = curStartPos;
                    regionInfo.endPos = curEndPos;
                    regionInfo.length = curEndPos - curStartPos + 1;
                    coverRegion[chr].emplace_back(regionInfo);
                    curStartPos = nextStartPos;
                    curEndPos = nextEndPos;
                }else{
                    // Regions overlap - merge them
                    curStartPos = std::min(curStartPos, nextStartPos);
                    curEndPos = std::max(curEndPos, nextEndPos);
                }
            // Handle last region
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

    // Calculate coverage ratios
    std::map<std::string, float> coverRegionRatio;
    long long totalChrLength = 0;
    long long totalChrCoverLength = 0;
    double totalChrCoverageRatio = 0.0;

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
        std::cerr<< "Fail to open write file: " << logFileName << "\n";
        exit(1);
    }else{
        // Write report header with coverage statistics
        (*tagReadCoverRegionLog) << "##################################\n";
        (*tagReadCoverRegionLog) << "# Somatic reads cover region bed #\n";
        (*tagReadCoverRegionLog) << "##################################\n";
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

    // Write merged coverage regions in BED format
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

/**
 * @brief Remove positions not derived by both H1 and H2 haplotypes
 * 
 * Filters out variant positions that are not derived from both H1 and H2
 * haplotypes. This is useful for quality control, ensuring that only
 * positions with sufficient haplotype support are retained for analysis.
 * 
 * @param chrVec Vector of chromosome names to process
 */
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
