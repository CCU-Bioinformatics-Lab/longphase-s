#ifndef SOMATIC_INDEL_MODEL_H
#define SOMATIC_INDEL_MODEL_H

#include <map>
#include <string>
#include <vector>

/**
 * @brief CNN post-filter for somatic INSERTION / DELETION calls.
 *
 * Hand-rolled forward pass of the FullAlignmentDualBranchCNN checkpoint
 * (output_multi_FA_event_thr_0.5_del_only/best_model.pt), so LongPhase-S gains
 * no ML runtime dependency.  Architecture, fixed by the checkpoint:
 *
 *   input   (6, 100, 99) = [tumor base/cigar/event, normal base/cigar/event]
 *           raw 0-255 values, NOT scaled to [0,1]
 *   BatchNorm2d(6, affine=False)
 *   tumor branch  <- channels 0..2  |  normal branch <- channels 3..5
 *       conv3x3(3->16)  BN ReLU  maxpool2
 *       conv3x3(16->32) BN ReLU  maxpool2
 *       conv3x3(32->64) BN ReLU  adaptiveAvgPool(4,4)   -> 1024
 *   scalars (6) -> BatchNorm1d(affine=False) -> 6->32 ReLU -> 32->32 ReLU
 *   concat(1024 + 1024 + 32 = 2080) -> 64 ReLU -> 32 ReLU -> 2 -> softmax
 *
 * Dropout / Dropout2d are identity at inference and are therefore absent.
 */
class SomaticIndelModel{
    public:
        SomaticIndelModel();

        /** Which built-in weight set to use. */
        enum VariantKind{ DELETION_MODEL, INSERTION_MODEL };

        /** Load one of the weight sets compiled into the binary. */
        bool loadEmbedded(VariantKind kind);

        /** Load an external weight blob ("LPSDELM1"), for testing / retraining. */
        bool load(const std::string &path);

        bool isLoaded() const { return loaded; }

        /**
         * @brief P(true positive) for one candidate.
         * @param image   6 * MAX_READS * WIDTH floats, channel-major, raw 0-255
         * @param scalars 6 floats, in SCALAR_FEATURE_NAMES order:
         *                purity, tumor_vaf, tumor_read_count, tumor_normal_ratio,
         *                event_length, normal_indel_frac
         */
        float predict(const float *image, const float *scalars) const;

        static const int IN_CHANNELS = 6;
        static const int PER_SAMPLE  = 3;
        static const int MAX_READS   = 100;
        static const int WIDTH       = 99;
        static const int N_SCALAR    = 6;

    private:
        struct Tensor{
            std::vector<int> dims;
            std::vector<float> data;
            int size() const { return (int)data.size(); }
        };

        std::map<std::string, Tensor> weights;
        bool loaded;

        const Tensor* get(const std::string &name) const;

        /** Parse a weight blob already in memory; shared by load / loadEmbedded. */
        bool parse(const unsigned char *blob, size_t length, const std::string &origin);

        /** One conv branch: 3 channels in -> 1024 features out. */
        void runBranch(const std::string &prefix, const std::vector<float> &in,
                       int inC, int inH, int inW, std::vector<float> &out) const;
};

#endif
