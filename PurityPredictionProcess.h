#ifndef PURITY_PREDICTION_PROCESS_H
#define PURITY_PREDICTION_PROCESS_H

#include "HaplotagProcess.h"
#include "PurityPrediction.h"
#include "SomaticVarCaller.h"
#include "TumorPurityPredictor.h"


class PurityPredictionProcess : public HaplotagProcess {
    private:
        PurityPredictionParameters& params;

    protected:
        virtual void parseVariantFiles(VcfParser& vcfParser) override;

    public:
        PurityPredictionProcess(PurityPredictionParameters& params);
        ~PurityPredictionProcess();

        void predictPurity();
};



#endif