#ifndef SOMATIC_HAPLOTAG_H
#define SOMATIC_HAPLOTAG_H

#include "Haplotag.h"
#include "SomaticHaplotagProcess.h"



class SomaticHaplotagHelpManager : public HaplotagHelpManager {
    public:
        SomaticHaplotagHelpManager(const std::string& program) : HaplotagHelpManager(program) {}

        virtual void buildMessage() override;
    
};

class SomaticHaplotagOptionManager : public HaplotagOptionManager {
    protected:
        virtual HelpMessageManager* createHelpManager(const std::string& program) override {
            return new SomaticHaplotagHelpManager(program);
        }

        virtual void initializeDefaultValues() override;

        virtual bool loadOptions(char& opt, std::istringstream& arg) override;

        virtual bool validateFiles() override;
        virtual bool validateNumericParameter() override;

    public:
        SomaticHaplotagOptionManager(const std::string& program) : HaplotagOptionManager(program) {}
        virtual void setOptions() override;
        virtual ~SomaticHaplotagOptionManager() = default;
};



// functions
int SomaticHaplotagMain(int argc, char** argv, std::string in_version);


#endif