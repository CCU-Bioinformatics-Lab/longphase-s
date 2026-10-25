#include "FullAlignmentEncoder.h"

#include <algorithm>
#include <cmath>
#include <cerrno>
#include <cstdio>
#include <iostream>
#include <set>

namespace {

// Encoding constants; these mirror full_alignment.py and must not drift.
const float BASE_A = 250.0f;
const float BASE_C = 30.0f;
const float BASE_G = 180.0f;
const float BASE_T = 100.0f;
const float STRAND_FWD = 250.0f;
const float STRAND_REV = 100.0f;
const float FA_CIGAR_M = 128.0f;
const float FA_CIGAR_I = 255.0f;
const float FA_CIGAR_D = 16.0f;
const float FA_EVENT_DEL = 16.0f;
const float INS_ANCHOR_BASE = 0.0f;
// BASE_255.get(base, 60.0): only the inserted-base lookup has a fallback
const float INS_BASE_FALLBACK = 60.0f;

/** htslib 4-bit nt16 code -> base channel value ("=ACMGRSVTWYHKDBN"). */
float base255(uint8_t c){
    switch(c){
        case 1:  return BASE_A;
        case 2:  return BASE_C;
        case 4:  return BASE_G;
        case 8:  return BASE_T;
        default: return 0.0f;              // '=', N and ambiguity codes
    }
}

/** Inserted-base value for the event channel (N is a real 0, not the fallback). */
float insBase255(uint8_t c){
    switch(c){
        case 1:  return BASE_A;
        case 2:  return BASE_C;
        case 4:  return BASE_G;
        case 8:  return BASE_T;
        case 15: return 0.0f;              // 'N'
        default: return INS_BASE_FALLBACK;
    }
}

bool nearly(float a, float b){
    return std::fabs(a - b) <= 1e-5f;
}

/** Half-open runs of `value` within one row, like numpy's isclose grouping. */
void valueRuns(const float *row, int width, float value,
               std::vector<std::pair<int,int> > &runs){
    runs.clear();
    int start = -1;
    for(int i = 0; i < width; i++){
        const bool ok = nearly(row[i], value);
        if(ok && start < 0){
            start = i;
        }else if(!ok && start >= 0){
            runs.push_back(std::make_pair(start, i));
            start = -1;
        }
    }
    if(start >= 0) runs.push_back(std::make_pair(start, width));
}

/** Tumor alt support: rows whose deletion run hits the candidate interval. */
int countTumorDelRows(const float *cigarChannel, int maxReads, int width,
                      int center, int eventLen){
    if(eventLen <= 0) return 0;
    const int right = center + eventLen;
    int count = 0;
    std::vector<std::pair<int,int> > runs;
    for(int r = 0; r < maxReads; r++){
        valueRuns(cigarChannel + (size_t)r * width, width, FA_CIGAR_D, runs);
        for(size_t i = 0; i < runs.size(); i++){
            const int s = runs[i].first;
            const int e = runs[i].second;
            if(s < right && e > center && e <= right){
                count++;
                break;
            }
        }
    }
    return count;
}

/** Tumor INS support: rows whose insertion run STARTS at the centre column. */
int countInsAnchorRows(const float *cigarChannel, int maxReads, int width, int center){
    if(center < 0 || center >= width) return 0;
    int count = 0;
    for(int r = 0; r < maxReads; r++){
        const float *row = cigarChannel + (size_t)r * width;
        const bool atCenter = nearly(row[center], FA_CIGAR_I);
        const bool startsHere = (center > 0)
            ? (atCenter && !nearly(row[center - 1], FA_CIGAR_I))
            : atCenter;
        if(startsHere) count++;
    }
    return count;
}

/**
 * @brief Normal indel support: an insertion anchored at the centre column, or a
 *        deletion run crossing it.
 *
 * `boundDelRight` (DEL candidates only) additionally requires the deletion to
 * stop inside the candidate interval; INS candidates do not bound it.
 */
int countNormalIndelRows(const float *cigarChannel, int maxReads, int width,
                         int center, int eventLen, bool boundDelRight){
    std::set<int> rows;
    if(center >= 0 && center < width){
        for(int r = 0; r < maxReads; r++){
            const float *row = cigarChannel + (size_t)r * width;
            const bool atCenter = nearly(row[center], FA_CIGAR_I);
            const bool startsHere = (center > 0)
                ? (atCenter && !nearly(row[center - 1], FA_CIGAR_I))
                : atCenter;
            if(startsHere) rows.insert(r);
        }
    }
    const int right = center + eventLen;
    std::vector<std::pair<int,int> > runs;
    for(int r = 0; r < maxReads; r++){
        valueRuns(cigarChannel + (size_t)r * width, width, FA_CIGAR_D, runs);
        for(size_t i = 0; i < runs.size(); i++){
            const int s = runs[i].first;
            const int e = runs[i].second;
            const bool crossesPos = (s <= center && e > center);
            const bool withinRight = !boundDelRight || (eventLen > 0 && e <= right);
            if(crossesPos && withinRight){
                rows.insert(r);
                break;
            }
        }
    }
    return (int)rows.size();
}

}  // namespace

FullAlignmentEncoder::FullAlignmentEncoder()
    : tumorFp(NULL), normalFp(NULL), tumorIdx(NULL), normalIdx(NULL),
      tumorHdr(NULL), normalHdr(NULL), opened(false){}

FullAlignmentEncoder::~FullAlignmentEncoder(){
    close();
}

void FullAlignmentEncoder::close(){
    if(tumorHdr){ bam_hdr_destroy(tumorHdr); tumorHdr = NULL; }
    if(normalHdr){ bam_hdr_destroy(normalHdr); normalHdr = NULL; }
    if(tumorIdx){ hts_idx_destroy(tumorIdx); tumorIdx = NULL; }
    if(normalIdx){ hts_idx_destroy(normalIdx); normalIdx = NULL; }
    if(tumorFp){ sam_close(tumorFp); tumorFp = NULL; }
    if(normalFp){ sam_close(normalFp); normalFp = NULL; }
    opened = false;
}

bool FullAlignmentEncoder::open(const std::string &tumorBam, const std::string &normalBam,
                                const std::string &fastaFile){
    close();
    tumorFp = sam_open(tumorBam.c_str(), "r");
    normalFp = sam_open(normalBam.c_str(), "r");
    if(tumorFp == NULL || normalFp == NULL){
        std::cerr << "[FullAlignmentEncoder] cannot open BAM: "
                  << (tumorFp == NULL ? tumorBam : normalBam) << "\n";
        close();
        return false;
    }
    // CRAM input cannot decode read sequences without the reference
    if(!fastaFile.empty()){
        if(hts_set_fai_filename(tumorFp, fastaFile.c_str()) != 0 ||
           hts_set_fai_filename(normalFp, fastaFile.c_str()) != 0){
            std::cerr << "[FullAlignmentEncoder] cannot set reference: " << fastaFile << "\n";
            close();
            return false;
        }
    }
    // CRAM decodes every field except base qualities, which the encoder never reads (no effect on BAM input)
    const int requiredFields = SAM_QNAME | SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR |
                               SAM_RNEXT | SAM_PNEXT | SAM_TLEN | SAM_SEQ | SAM_AUX | SAM_RGAUX;
    if(hts_set_opt(tumorFp, CRAM_OPT_REQUIRED_FIELDS, requiredFields) != 0 ||
       hts_set_opt(normalFp, CRAM_OPT_REQUIRED_FIELDS, requiredFields) != 0){
        std::cerr << "[FullAlignmentEncoder] cannot set required fields\n";
        close();
        return false;
    }
    tumorIdx = sam_index_load(tumorFp, tumorBam.c_str());
    normalIdx = sam_index_load(normalFp, normalBam.c_str());
    if(tumorIdx == NULL || normalIdx == NULL){
        std::cerr << "[FullAlignmentEncoder] cannot load BAM index: "
                  << (tumorIdx == NULL ? tumorBam : normalBam) << "\n";
        close();
        return false;
    }
    tumorHdr = sam_hdr_read(tumorFp);
    normalHdr = sam_hdr_read(normalFp);
    if(tumorHdr == NULL || normalHdr == NULL){
        std::cerr << "[FullAlignmentEncoder] cannot read BAM header\n";
        close();
        return false;
    }
    opened = true;
    return true;
}

int FullAlignmentEncoder::resolveTid(bam_hdr_t *hdr, const std::string &chrom){
    int tid = bam_name2id(hdr, chrom.c_str());
    if(tid >= 0) return tid;
    if(chrom.size() > 3 && chrom.compare(0, 3, "chr") == 0){
        tid = bam_name2id(hdr, chrom.substr(3).c_str());
        if(tid >= 0) return tid;
    }
    return bam_name2id(hdr, ("chr" + chrom).c_str());
}

bool FullAlignmentEncoder::encodeOneRead(const bam1_t *b, int winStart, EncodedRead &out){
    const bam1_core_t *c = &b->core;
    if(c->n_cigar == 0 || c->l_qseq == 0) return false;

    out.base.assign(WIDTH, 0.0f);
    out.cigar.assign(WIDTH, 0.0f);
    out.event.assign(WIDTH, 0.0f);
    out.strand.assign(WIDTH, 0.0f);

    const float strandVal = (c->flag & BAM_FREVERSE) ? STRAND_REV : STRAND_FWD;
    const uint8_t *seq = bam_get_seq(b);
    const uint32_t *cig = bam_get_cigar(b);

    int refPos = c->pos;                       // 0-based
    int qPos = 0;
    bool touched = false;
    // Insertions are applied after the walk so a following M cannot overwrite
    // the spread; each entry is (anchor column, inserted nt16 bases).
    std::vector<std::pair<int, std::vector<uint8_t> > > insertions;

    for(uint32_t i = 0; i < c->n_cigar; i++){
        const int op = bam_cigar_op(cig[i]);
        const int len = bam_cigar_oplen(cig[i]);

        if(op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF){
            for(int k = 0; k < len; k++){
                const int col = (refPos + k) - winStart;
                if(col >= 0 && col < WIDTH){
                    out.base[col] = base255(bam_seqi(seq, qPos + k));
                    out.strand[col] = strandVal;
                    out.cigar[col] = FA_CIGAR_M;
                    touched = true;
                }
            }
            refPos += len;
            qPos += len;
        }else if(op == BAM_CINS){
            const int anchorCol = (refPos - 1) - winStart;
            if(anchorCol >= 0 && anchorCol < WIDTH){
                std::vector<uint8_t> insSeq(len);
                for(int k = 0; k < len; k++) insSeq[k] = bam_seqi(seq, qPos + k);
                insertions.push_back(std::make_pair(anchorCol, insSeq));
                touched = true;
            }
            qPos += len;
        }else if(op == BAM_CDEL || op == BAM_CREF_SKIP){
            // The left-anchor base is marked in the base channel (symmetric with
            // the insertion anchor); the cigar channel keeps it as a real match.
            const int anchorCol = (refPos - 1) - winStart;
            if(anchorCol >= 0 && anchorCol < WIDTH){
                out.base[anchorCol] = INS_ANCHOR_BASE;
                touched = true;
            }
            for(int k = 0; k < len; k++){
                const int col = (refPos + k) - winStart;
                if(col >= 0 && col < WIDTH){
                    out.strand[col] = strandVal;
                    out.cigar[col] = FA_CIGAR_D;
                    out.event[col] = FA_EVENT_DEL;
                    touched = true;
                }
            }
            refPos += len;
        }else if(op == BAM_CSOFT_CLIP){
            qPos += len;
        }
        // hard clip / pad consume neither reference nor query
    }

    for(size_t i = 0; i < insertions.size(); i++){
        const int anchorCol = insertions[i].first;
        const std::vector<uint8_t> &s = insertions[i].second;
        out.base[anchorCol] = INS_ANCHOR_BASE;
        const int n = std::min((int)s.size(), WIDTH - anchorCol);
        for(int k = 0; k < n; k++){
            out.event[anchorCol + k] = insBase255(s[k]);
            out.cigar[anchorCol + k] = FA_CIGAR_I;
        }
    }

    if(!touched) return false;

    out.hp = 0;
    uint8_t *hp = bam_aux_get(const_cast<bam1_t *>(b), "HP");
    if(hp != NULL){
        errno = 0;
        const int64_t v = bam_aux2i(hp);
        if(errno == 0) out.hp = (int)v;
    }
    out.readid = bam_get_qname(b);
    return true;
}

bool FullAlignmentEncoder::buildSample(samFile *fp, hts_idx_t *idx, bam_hdr_t *hdr,
                                       const std::string &chrom, int pos1,
                                       const FaReadDelMap *delMap,
                                       std::vector<float> &out, int chOffset) const{
    const int pos0 = pos1 - 1;
    const int winStart = pos0 - HALF;

    const int tid = resolveTid(hdr, chrom);
    if(tid < 0) return true;                       // leave the window all zeros

    hts_itr_t *iter = sam_itr_queryi(idx, tid, std::max(0, winStart), winStart + WIDTH);
    if(iter == NULL) return true;

    std::vector<EncodedRead> rows;
    bam1_t *b = bam_init1();
    while(sam_itr_next(fp, iter, b) >= 0){
        const uint16_t flag = b->core.flag;
        if((flag & BAM_FUNMAP) || (flag & BAM_FSECONDARY) ||
           (flag & BAM_FSUPPLEMENTARY) || (flag & BAM_FDUP)) continue;
        EncodedRead read;
        if(encodeOneRead(b, winStart, read)) rows.push_back(read);
    }
    bam_destroy1(b);
    hts_itr_destroy(iter);

    // Haplotype grouping: HP1, HP2, then everything else; stable so untagged
    // BAMs keep the coordinate order htslib returned.
    std::stable_sort(rows.begin(), rows.end(),
                     [](const EncodedRead &a, const EncodedRead &b2){
                         const int ka = (a.hp == 1 || a.hp == 2) ? a.hp : 9;
                         const int kb = (b2.hp == 1 || b2.hp == 2) ? b2.hp : 9;
                         return ka < kb;
                     });

    const int n = std::min((int)rows.size(), (int)MAX_READS);
    const size_t perCh = (size_t)MAX_READS * WIDTH;
    for(int i = 0; i < n; i++){
        float *pb = &out[(size_t)(chOffset + 0) * perCh + (size_t)i * WIDTH];
        float *pc = &out[(size_t)(chOffset + 1) * perCh + (size_t)i * WIDTH];
        float *pe = &out[(size_t)(chOffset + 2) * perCh + (size_t)i * WIDTH];
        std::copy(rows[i].base.begin(), rows[i].base.end(), pb);
        std::copy(rows[i].cigar.begin(), rows[i].cigar.end(), pc);
        std::copy(rows[i].event.begin(), rows[i].event.end(), pe);
    }

    // Focused del-map: replace the cigar channel with "covered = M, plus a
    // D-run of delLen starting at (centre - shift)" for every read the caller
    // recorded a deletion for.  Rows stay aligned to the BAM via the read id.
    if(delMap != NULL){
        for(int i = 0; i < n; i++){
            float *pc = &out[(size_t)(chOffset + 1) * perCh + (size_t)i * WIDTH];
            for(int c = 0; c < WIDTH; c++){
                pc[c] = (rows[i].strand[c] > 0.0f) ? FA_CIGAR_M : 0.0f;
            }
            FaReadDelMap::const_iterator it = delMap->find(rows[i].readid);
            if(it != delMap->end()){
                const int start = HALF - it->second.shift;
                for(int k = 0; k < it->second.delLen; k++){
                    const int col = start + k;
                    if(col >= 0 && col < WIDTH) pc[col] = FA_CIGAR_D;
                }
            }
        }
    }
    return true;
}

bool FullAlignmentEncoder::buildImage(const std::string &chrom, int pos1,
                                      const FaReadDelMap *tumorMap,
                                      const FaReadDelMap *normalMap,
                                      std::vector<float> &image) const{
    if(!opened) return false;
    image.assign((size_t)IN_CHANNELS * MAX_READS * WIDTH, 0.0f);
    bool ok = buildSample(tumorFp, tumorIdx, tumorHdr, chrom, pos1, tumorMap, image, 0);
    ok &= buildSample(normalFp, normalIdx, normalHdr, chrom, pos1, normalMap, image, PER_SAMPLE);
    return ok;
}

void FullAlignmentEncoder::buildScalars(const std::vector<float> &image, VariantKind kind,
                                        double purity, int tumorDepth, int normalDepth,
                                        int eventLen, float *scalars){
    const size_t perCh = (size_t)MAX_READS * WIDTH;
    const float *tumorCigar = &image[1 * perCh];
    const float *normalCigar = &image[(size_t)(PER_SAMPLE + 1) * perCh];

    const bool isIns = (kind == INSERTION_CANDIDATE);
    const int alt = isIns
        ? countInsAnchorRows(tumorCigar, MAX_READS, WIDTH, HALF)
        : countTumorDelRows(tumorCigar, MAX_READS, WIDTH, HALF, eventLen);
    const int normalAlt = countNormalIndelRows(normalCigar, MAX_READS, WIDTH, HALF,
                                               eventLen, !isIns);

    scalars[0] = (float)purity;
    scalars[1] = tumorDepth > 0 ? (float)((double)alt / tumorDepth) : 0.0f;
    scalars[2] = (float)(tumorDepth / 100.0);
    scalars[3] = (float)((double)tumorDepth / (normalDepth + 1) / 50.0);
    scalars[4] = (float)(eventLen / 50.0);
    scalars[5] = normalDepth > 0 ? (float)((double)normalAlt / normalDepth) : 0.0f;
}
