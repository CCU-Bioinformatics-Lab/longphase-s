#ifndef PURITY_PREDICTION_PROCESS_H
#define PURITY_PREDICTION_PROCESS_H

#include "HaplotagProcess.h"
#include "SomaticVarCaller.h"
#include "TumorPurityPredictor.h"

struct PurityPredictionParameters
{
    HaplotagParameters basic;
    std::string tumorBamFile;
    std::string tumorSnpFile;
};

class PurityPredictionProcess : public HaplotagProcess {
    private:
        double tumorPurity;
        PurityPredictionParameters& params;
        virtual void printParamsMessage() override;

    protected:
        virtual void parseVariantFiles(VcfParser& vcfParser) override;
        virtual void printExecutionReport() override;
    public:
        PurityPredictionProcess(PurityPredictionParameters& params);
        ~PurityPredictionProcess();

        void predictPurity();
};



#endif