#ifndef HAPLOTAG_LOGGING_H
#define HAPLOTAG_LOGGING_H

#include "HaplotagType.h"

struct chrReadHpResult{
    std::map<int, ReadHpResult> posReadHpResult;

    void recordReadHp(int &pos, int &hpResult, int &BaseHP);
    void recordDeriveHp(int &pos, int &deriveHP, float deriveHPsimilarity);
    void recordAlignCoverRegion(int& curVarPos, int &startPos, int &endPos);
};


class LogEntry {
    public:
        std::string key;
        std::string value;
        size_t order;

        LogEntry(const std::string& k, const std::string& v, size_t o) 
            : key(k), value(v), order(o) {}
};

class MessageManager{
    protected:
        std::vector<LogEntry> entries;
        size_t currentOrder;

        //Generally add entry (add at the end)
        template<typename T>
        void addEntry(const std::string& key, const T& value) {
            std::string valueStr = transformValueToString(value);

            entries.emplace_back(key, valueStr, currentOrder);
            // std::cerr << "Added entry: " << key << ":" << valueStr << " at index " << currentOrder << std::endl;
            currentOrder++;
        }

        // Insert entry after a specific order
        template<typename T>
        void insertByIndex(const std::string& key, const T& value, size_t insertIndex) {
            // move the order of the entries after the insertIndex
            for (auto& entry : entries) {
                if (entry.order >= insertIndex) {
                    entry.order ++;
                }
            }
            std::string valueStr = transformValueToString(value);
            entries.emplace_back(key, valueStr, insertIndex);
            // std::cerr << "Inserted entry: " << key << ":" << valueStr << " at index " << insertIndex << std::endl;
        }


        // Insert a new message after a specific key
        // If the target key is not found, append the message at the end
        template<typename T>
        void insertAfterKey(const std::string& newKey, const T& newValue, const std::string& targetKey) {
            // Find the position of the target key
            auto it = std::find_if(entries.begin(), entries.end(),
                [&targetKey](const LogEntry& entry) { return entry.key == targetKey; });
            
            if (it == entries.end()) {
                // If target key not found, append to the end
                std::cerr << "[ERROR](insertAfterKey) Target key not found: " << targetKey << std::endl;
                exit(1);
            }
            
            // Calculate insertion position (target key's order + 1)
            size_t insertIndex = it->order + 1;
            
            // Increment order for all entries after the insertion point
            for (auto& entry : entries) {
                if (entry.order >= insertIndex) {
                    entry.order++;
                }
            }
            
            // Insert the new entry
            std::string valueStr = transformValueToString(newValue);
            entries.emplace_back(newKey, valueStr, insertIndex);
        }

        void sortEntries(){
            std::sort(entries.begin(), entries.end(), [](const LogEntry& a, const LogEntry& b) {
                return a.order < b.order;
            });
        }

        void printMessage(){
            // Write all entries
            for (const auto& entry : entries) {
                std::cout << entry.key << ":" << entry.value << "\n";
            }
        }

        // basic type
        template<typename T>
        std::string transformValueToString(const T& value) {
            return std::to_string(value);
        }

        // Template Specialization for std::string
        std::string transformValueToString(const std::string& value) {
            return value;
        }

        // Template Specialization for bool
        std::string transformValueToString(const bool& value) {
            return value ? "true" : "false";
        }

        // Template Specialization for const char*
        std::string transformValueToString(const char* value) {
            return value;
        }

        // Template Specialization for char[]
        template<size_t N>
        std::string transformValueToString(const char (&value)[N]) {
            return std::string(value);
        }

        // Template Specialization for double
        std::string transformValueToString(const double& value) {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(2) << value;
            return oss.str();
        }

        // Template Specialization for float
        std::string transformValueToString(const float& value) {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(2) << value;
            return oss.str();
        }

    public:
        MessageManager() : currentOrder(0) {}

        virtual ~MessageManager() = default;

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

#endif
