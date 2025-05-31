#ifndef HAPLOTAG_LOGGING_H
#define HAPLOTAG_LOGGING_H

#include "HaplotagType.h"

struct chrReadHpResult{
    std::map<int, ReadHpResult> posReadHpResult;

    void recordReadHp(int &pos, int &hpResult, int &BaseHP);
    void recordDeriveHp(int &pos, int &deriveHP, float deriveHPsimilarity);
    void recordAlignCoverRegion(int& curVarPos, int &startPos, int &endPos);
};

class ReadHpDistriLog{
    private :
        struct coverRegionInfo{
            int startPos;
            int endPos;
            int length;

            coverRegionInfo(): startPos(0), endPos(0), length(0){}
        };
        // chr, variant position (0-base), reads HP 
        std::map<std::string, chrReadHpResult> chrVarReadHpResult;
    protected:

    public:
        ReadHpDistriLog();
        ~ReadHpDistriLog();

        // use in multi-thread scenario
        void loadChrKey(const std::string &chr);
        // Returns a pointer to ensure thread-safe access to chromosome results
        chrReadHpResult* getChrHpResultsPtr (const std::string &chr);

        // only use in single thread scenario
        void recordChrReadHp(const std::string &chr, int &pos, int &hpResult, int &BaseHP);
        void recordChrDeriveHp(const std::string &chr, int &pos, int &deriveHP, float deriveHPsimilarity);
        void recordChrAlignCoverRegion(const std::string &chr, int &pos, int &startPos, int &endPos);

        void writeReadHpDistriLog(const std::string logFileName, const std::vector<std::string> &chrVec);
        void writePosCoverRegionLog(const std::string logFileName, const std::vector<std::string> &chrVec);
        void writeTagReadCoverRegionLog(const std::string logFileName, const std::vector<std::string> &chrVec, std::map<std::string, int> &chrLength);
        void removeNotDeriveByH1andH2pos(const std::vector<std::string> &chrVec);
};

template<typename ParamType, typename LogType>
class HaplotagReadLog{
    protected:
        const ParamType& params;

        void checkStreamStatus() {
            if (!tagReadLog || !tagReadLog->is_open()) {
                throw std::runtime_error("Output stream is not valid");
            }
        }

        virtual void addParamsMessage() = 0;

        virtual void writeBasicColumns() = 0;

    public:
        std::ofstream* tagReadLog;

        HaplotagReadLog(const ParamType& params, std::string fileName) : params(params) {
            tagReadLog = new std::ofstream(fileName);
            if(!tagReadLog->is_open()){
                std::cerr<< "Fail to open write file: " << fileName << "\n";
                exit(1);
            }
        };

        virtual ~HaplotagReadLog(){
            tagReadLog->close();
            if(tagReadLog) delete tagReadLog;
        };

        virtual void writeHeader() {
            addParamsMessage();

            writeBasicColumns();
        }

        virtual void writeTagReadLog(LogType& data) = 0;
};

#endif
