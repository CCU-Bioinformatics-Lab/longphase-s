#ifndef FULL_ALIGNMENT_ENCODER_H
#define FULL_ALIGNMENT_ENCODER_H

#include <map>
#include <string>
#include <vector>

#include "htslib/sam.h"

/**
 * @brief Per-read deletion geometry at one candidate site.
 *
 * Mirrors the (del_len, shift) pair the Python training pipeline read from
 * *_tumor_deletion_shift_info.log / *_normal_deletion_info.log:
 *   tumor  : delLen = TumorDelShiftInfo::CIGAR_Del_Length,
 *            shift  = TumorDelShiftInfo::shiftLength
 *   normal : delLen = NormalDelInfo::CIGAR_Del_Length,
 *            shift  = candidatePos(1-based) - NormalDelInfo::normalDelPos
 */
struct FaReadDel{
    int delLen;
    int shift;

    FaReadDel(): delLen(0), shift(0){}
    FaReadDel(int delLen, int shift): delLen(delLen), shift(shift){}
};

/** read id -> deletion geometry, for one candidate site. */
typedef std::map<std::string, FaReadDel> FaReadDelMap;

/**
 * @brief Builds the DeepSomatic-style full-alignment tensor consumed by
 *        SomaticDelModel.
 *
 * This is a C++ port of full_alignment.py (encode_one_read /
 * build_sample_channels / build_site_tensor) with extra=("cigar","event"),
 * verified byte-for-byte against the cached tensors the checkpoint was
 * trained on.
 *
 * Channel layout (6 x MAX_READS x WIDTH, channel major):
 *   0 tumor_base    A=250 C=30 G=180 T=100, indel anchor / no coverage = 0
 *   1 tumor_cigar   focused del-map: covered = 128 (M), deletion run = 16 (D);
 *                   falls back to the BAM CIGAR (M=128, I=255, D=16) when the
 *                   site has no per-read deletion entries
 *   2 tumor_event   deletion span = 16, inserted bases spread from the anchor
 *   3-5 the same three channels for the normal sample
 *
 * One instance owns its own BAM handles, so give each thread its own encoder.
 */
class FullAlignmentEncoder{
    public:
        FullAlignmentEncoder();
        ~FullAlignmentEncoder();

        /**
         * @brief Open the tumor / normal alignment pair (indexes required).
         * @param fastaFile reference FASTA, needed to decode CRAM input
         */
        bool open(const std::string &tumorBam, const std::string &normalBam,
                  const std::string &fastaFile);
        void close();
        bool isOpen() const { return opened; }

        /**
         * @brief Fill `image` (IN_CHANNELS * MAX_READS * WIDTH floats, 0-255).
         * @param chrom      chromosome name (chr-prefix mismatches are resolved)
         * @param pos1       candidate position, 1-based (VCF coordinate)
         * @param tumorMap   per-read del-map, or NULL to use the BAM CIGAR
         * @param normalMap  same for the normal sample
         */
        bool buildImage(const std::string &chrom, int pos1,
                        const FaReadDelMap *tumorMap,
                        const FaReadDelMap *normalMap,
                        std::vector<float> &image) const;

        /** Which candidate type the scalar rules should follow. */
        enum VariantKind{ DELETION_CANDIDATE, INSERTION_CANDIDATE };

        /**
         * @brief Fill the 6 scalar features, in SCALAR_FEATURE_NAMES order:
         *        purity, tumor_vaf, tumor_read_count, tumor_normal_ratio,
         *        event_length, normal_indel_frac.
         *
         * The two counts are read back off the cigar channels of `image`, which
         * is what the training pipeline did (data_loader.py
         * _build_channel_aligned_scalars), so they inherit the 100-read
         * truncation. INSERTION candidates count tumor support as insertions
         * anchored at the centre column and leave the normal deletion span
         * unbounded on the right; DELETION candidates use the interval rules.
         *
         * @param tumorDepth  LongPhase-S depth at the site, or 0 when the site
         *                    has no tumor per-read entries of this event type
         * @param normalDepth normal depth, or 0 when the site has no normal
         *                    per-read entries of this event type
         */
        static void buildScalars(const std::vector<float> &image, VariantKind kind,
                                 double purity, int tumorDepth, int normalDepth,
                                 int eventLen, float *scalars);

        static const int HALF        = 49;
        static const int WIDTH       = 99;
        static const int MAX_READS   = 100;
        static const int PER_SAMPLE  = 3;
        static const int IN_CHANNELS = 6;

    private:
        struct EncodedRead{
            std::vector<float> base;
            std::vector<float> cigar;
            std::vector<float> event;
            std::vector<float> strand;   // coverage mask; not an output channel
            int hp;
            std::string readid;
        };

        samFile *tumorFp;
        samFile *normalFp;
        hts_idx_t *tumorIdx;
        hts_idx_t *normalIdx;
        bam_hdr_t *tumorHdr;
        bam_hdr_t *normalHdr;
        bool opened;

        static bool encodeOneRead(const bam1_t *b, int winStart, EncodedRead &out);
        static int resolveTid(bam_hdr_t *hdr, const std::string &chrom);

        /** Fill PER_SAMPLE channels for one sample into `out` at `chOffset`. */
        bool buildSample(samFile *fp, hts_idx_t *idx, bam_hdr_t *hdr,
                         const std::string &chrom, int pos1,
                         const FaReadDelMap *delMap,
                         std::vector<float> &out, int chOffset) const;

        FullAlignmentEncoder(const FullAlignmentEncoder&);
        FullAlignmentEncoder& operator=(const FullAlignmentEncoder&);
};

#endif
