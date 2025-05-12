#ifndef SOMATIC_HAPLOTAG_H
#define SOMATIC_HAPLOTAG_H

#include "Haplotag.h"
#include "SomaticHaplotagProcess.h"


class SomaticHaplotagHelpManager : public HaplotagHelpManager {
    public:
        SomaticHaplotagHelpManager(const std::string& program) : HaplotagHelpManager(program) {}

        virtual void modifyMessage() override;
    
};

class SomaticHaplotagOptionManager : public HaplotagOptionManager {
    protected:
        virtual HelpMessageManager* createHelpManager(const std::string& program) override {
            return new SomaticHaplotagHelpManager(program);
        }

        virtual void extendOptions() override;

        virtual bool loadExtendOptions(char& opt, std::istringstream& arg) override;
        virtual bool validateExtendFiles() override;

    public:
        SomaticHaplotagOptionManager(const std::string& program) : HaplotagOptionManager(program) {
            //extend somatic haplotag options
            extendOptions();
        }
        virtual ~SomaticHaplotagOptionManager() = default;
};



// functions
int SomaticHaplotagMain(int argc, char** argv, std::string in_version);


#endif