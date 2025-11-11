#ifndef PURITY_ESTIMATION_PROCESS_H
#define PURITY_ESTIMATION_PROCESS_H

#include "../haplotag/HaplotagProcess.h"
#include "SomaticVarCaller.h"
#include "TumorPurityEstimator.h"

struct PurityEstimParameters
{
    HaplotagParameters basic;
    std::string tumorBamFile;
    std::string tumorSnvFile;
};

class PurityEstimProcess : public HaplotagProcess {
    private:
        double tumorPurity;
        PurityEstimParameters& params;
        virtual void printParamsMessage() override;

    protected:
        virtual void parseVariantFiles(VcfParser& vcfParser) override;
        virtual void printExecutionReport() override;
    public:
        PurityEstimProcess(PurityEstimParameters& params);
        ~PurityEstimProcess();

        void estimatePurity();
};



#endif