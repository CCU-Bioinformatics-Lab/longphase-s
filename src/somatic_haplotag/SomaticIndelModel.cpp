#include "SomaticIndelModel.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>

namespace {

const float BN_EPS = 1e-5f;

/** conv2d, stride 1, padding (kh/2, kw/2), kernel 3x3. */
void conv2d(const std::vector<float> &in, int inC, int inH, int inW,
            const std::vector<float> &w, const std::vector<float> &b,
            int outC, int kh, int kw, std::vector<float> &out) {
    const int padH = kh / 2, padW = kw / 2;
    const int outH = inH, outW = inW;                 // stride 1, same padding
    out.assign((size_t)outC * outH * outW, 0.0f);

    for (int oc = 0; oc < outC; ++oc) {
        float *op = &out[(size_t)oc * outH * outW];
        const float bias = b[oc];
        for (int i = 0; i < outH * outW; ++i) op[i] = bias;
        for (int ic = 0; ic < inC; ++ic) {
            const float *ip = &in[(size_t)ic * inH * inW];
            const float *wp = &w[((size_t)oc * inC + ic) * kh * kw];
            for (int kr = 0; kr < kh; ++kr) {
                for (int kc = 0; kc < kw; ++kc) {
                    const float wv = wp[kr * kw + kc];
                    if (wv == 0.0f) continue;
                    const int dr = kr - padH, dc = kc - padW;
                    const int r0 = std::max(0, -dr), r1 = std::min(outH, inH - dr);
                    const int c0 = std::max(0, -dc), c1 = std::min(outW, inW - dc);
                    for (int r = r0; r < r1; ++r) {
                        const float *irow = ip + (size_t)(r + dr) * inW + dc;
                        float *orow = op + (size_t)r * outW;
                        for (int c = c0; c < c1; ++c) orow[c] += wv * irow[c];
                    }
                }
            }
        }
    }
}

/** BatchNorm2d in eval mode, followed by ReLU when `relu` is set. */
void batchNormRelu(std::vector<float> &x, int C, int HW,
                   const std::vector<float> &gamma, const std::vector<float> &beta,
                   const std::vector<float> &mean, const std::vector<float> &var,
                   bool relu) {
    for (int c = 0; c < C; ++c) {
        const float inv = 1.0f / std::sqrt(var[c] + BN_EPS);
        const float g = gamma.empty() ? 1.0f : gamma[c];
        const float bt = beta.empty() ? 0.0f : beta[c];
        const float scale = g * inv;
        const float shift = bt - mean[c] * scale;
        float *p = &x[(size_t)c * HW];
        for (int i = 0; i < HW; ++i) {
            float v = p[i] * scale + shift;
            p[i] = (relu && v < 0.0f) ? 0.0f : v;
        }
    }
}

/** MaxPool2d(kernel=2, stride=2), floor mode (PyTorch default). */
void maxPool2(const std::vector<float> &in, int C, int inH, int inW,
              std::vector<float> &out, int &outH, int &outW) {
    outH = (inH - 2) / 2 + 1;
    outW = (inW - 2) / 2 + 1;
    out.assign((size_t)C * outH * outW, 0.0f);
    for (int c = 0; c < C; ++c) {
        const float *ip = &in[(size_t)c * inH * inW];
        float *op = &out[(size_t)c * outH * outW];
        for (int r = 0; r < outH; ++r) {
            for (int cc = 0; cc < outW; ++cc) {
                const int r0 = r * 2, c0 = cc * 2;
                float m = ip[(size_t)r0 * inW + c0];
                m = std::max(m, ip[(size_t)r0 * inW + c0 + 1]);
                m = std::max(m, ip[(size_t)(r0 + 1) * inW + c0]);
                m = std::max(m, ip[(size_t)(r0 + 1) * inW + c0 + 1]);
                op[(size_t)r * outW + cc] = m;
            }
        }
    }
}

/**
 * AdaptiveAvgPool2d((oh, ow)) with PyTorch's bin boundaries:
 *   start = floor(i * in / out), end = ceil((i+1) * in / out)
 * (bins overlap when `in` is not a multiple of `out`, e.g. 25 -> 4).
 */
void adaptiveAvgPool(const std::vector<float> &in, int C, int inH, int inW,
                     int oh, int ow, std::vector<float> &out) {
    out.assign((size_t)C * oh * ow, 0.0f);
    for (int c = 0; c < C; ++c) {
        const float *ip = &in[(size_t)c * inH * inW];
        float *op = &out[(size_t)c * oh * ow];
        for (int i = 0; i < oh; ++i) {
            const int r0 = (i * inH) / oh;
            const int r1 = ((i + 1) * inH + oh - 1) / oh;
            for (int j = 0; j < ow; ++j) {
                const int c0 = (j * inW) / ow;
                const int c1 = ((j + 1) * inW + ow - 1) / ow;
                double s = 0.0;
                for (int r = r0; r < r1; ++r)
                    for (int cc = c0; cc < c1; ++cc) s += ip[(size_t)r * inW + cc];
                op[(size_t)i * ow + j] = (float)(s / ((r1 - r0) * (c1 - c0)));
            }
        }
    }
}

/** y = W x + b, optionally followed by ReLU. */
void linear(const std::vector<float> &x, const std::vector<float> &w,
            const std::vector<float> &b, int inF, int outF,
            std::vector<float> &y, bool relu) {
    y.assign(outF, 0.0f);
    for (int o = 0; o < outF; ++o) {
        const float *wp = &w[(size_t)o * inF];
        double s = b[o];
        for (int i = 0; i < inF; ++i) s += (double)wp[i] * x[i];
        float v = (float)s;
        y[o] = (relu && v < 0.0f) ? 0.0f : v;
    }
}

}  // namespace

SomaticIndelModel::SomaticIndelModel() : loaded(false) {}

const SomaticIndelModel::Tensor* SomaticIndelModel::get(const std::string &name) const {
    std::map<std::string, Tensor>::const_iterator it = weights.find(name);
    return (it == weights.end()) ? NULL : &it->second;
}

// Weights compiled into the binary; SomaticIndelModelWeights.cpp .incbin's each
// blob between its pair of symbols.
extern "C" const unsigned char SOMATIC_DEL_MODEL_WEIGHTS[];
extern "C" const unsigned char SOMATIC_DEL_MODEL_WEIGHTS_END[];
extern "C" const unsigned char SOMATIC_INS_MODEL_WEIGHTS[];
extern "C" const unsigned char SOMATIC_INS_MODEL_WEIGHTS_END[];

bool SomaticIndelModel::loadEmbedded(VariantKind kind) {
    const unsigned char *begin = (kind == INSERTION_MODEL)
        ? SOMATIC_INS_MODEL_WEIGHTS : SOMATIC_DEL_MODEL_WEIGHTS;
    const unsigned char *end = (kind == INSERTION_MODEL)
        ? SOMATIC_INS_MODEL_WEIGHTS_END : SOMATIC_DEL_MODEL_WEIGHTS_END;
    const char *origin = (kind == INSERTION_MODEL) ? "<built-in ins>" : "<built-in del>";
    return parse(begin, (size_t)(end - begin), origin);
}

bool SomaticIndelModel::load(const std::string &path) {
    std::ifstream in(path.c_str(), std::ios::binary);
    if (!in.is_open()) {
        std::cerr << "[SomaticIndelModel] cannot open model file: " << path << "\n";
        loaded = false;
        return false;
    }
    std::vector<unsigned char> blob((std::istreambuf_iterator<char>(in)),
                                     std::istreambuf_iterator<char>());
    return parse(blob.empty() ? NULL : &blob[0], blob.size(), path);
}

bool SomaticIndelModel::parse(const unsigned char *blob, size_t length,
                            const std::string &origin) {
    loaded = false;
    weights.clear();

    size_t off = 0;
    // Every read is bounds-checked against `length`, so a truncated or corrupt
    // blob reports an error instead of running off the end.
    if (length < 12 || std::memcmp(blob, "LPSDELM1", 8) != 0) {
        std::cerr << "[SomaticIndelModel] bad magic in " << origin << " (expected LPSDELM1)\n";
        return false;
    }
    off = 8;

    int32_t n = 0;
    std::memcpy(&n, blob + off, 4);
    off += 4;

    for (int t = 0; t < n; ++t) {
        if (off + 4 > length) { std::cerr << "[SomaticIndelModel] truncated blob\n"; return false; }
        int32_t nameLen = 0;
        std::memcpy(&nameLen, blob + off, 4);
        off += 4;
        if (nameLen <= 0 || nameLen > 4096 || off + (size_t)nameLen > length) {
            std::cerr << "[SomaticIndelModel] corrupt tensor name in " << origin << "\n";
            return false;
        }
        std::string name(reinterpret_cast<const char *>(blob + off), (size_t)nameLen);
        off += nameLen;

        if (off + 4 > length) { std::cerr << "[SomaticIndelModel] truncated blob\n"; return false; }
        int32_t ndim = 0;
        std::memcpy(&ndim, blob + off, 4);
        off += 4;
        if (ndim < 0 || ndim > 8) {
            std::cerr << "[SomaticIndelModel] bad tensor rank for " << name << "\n";
            return false;
        }

        Tensor tensor;
        long long count = 1;
        for (int d = 0; d < ndim; ++d) {
            if (off + 4 > length) { std::cerr << "[SomaticIndelModel] truncated blob\n"; return false; }
            int32_t dim = 0;
            std::memcpy(&dim, blob + off, 4);
            off += 4;
            if (dim <= 0) {
                std::cerr << "[SomaticIndelModel] bad dimension for " << name << "\n";
                return false;
            }
            tensor.dims.push_back(dim);
            count *= dim;
        }
        if (ndim == 0) count = 1;                      // 0-d scalar buffers

        const size_t bytes = (size_t)count * 4;
        if (off + bytes > length) {
            std::cerr << "[SomaticIndelModel] truncated tensor: " << name << "\n";
            return false;
        }
        tensor.data.resize((size_t)count);
        std::memcpy(&tensor.data[0], blob + off, bytes);
        off += bytes;

        weights[name] = tensor;
    }

    // sanity: every tensor the forward pass reads must exist
    const char *required[] = {
        "input_norm.running_mean", "input_norm.running_var",
        "scalar_norm.running_mean", "scalar_norm.running_var",
        "scalar_mlp.0.weight", "scalar_mlp.0.bias",
        "scalar_mlp.2.weight", "scalar_mlp.2.bias",
        "fc1.weight", "fc1.bias", "fc2.weight", "fc2.bias",
        "fc3.weight", "fc3.bias", "branch_gate_fixed", NULL};
    for (int i = 0; required[i] != NULL; ++i) {
        if (get(required[i]) == NULL) {
            std::cerr << "[SomaticIndelModel] missing tensor: " << required[i] << "\n";
            return false;
        }
    }
    const char *branches[] = {"tumor_branch", "normal_branch", NULL};
    const int convIdx[] = {0, 5, 10};
    for (int b = 0; branches[b] != NULL; ++b) {
        for (int c = 0; c < 3; ++c) {
            char buf[128];
            snprintf(buf, sizeof(buf), "%s.net.%d.weight", branches[b], convIdx[c]);
            if (get(buf) == NULL) {
                std::cerr << "[SomaticIndelModel] missing tensor: " << buf << "\n";
                return false;
            }
        }
    }
    loaded = true;
    return true;
}

void SomaticIndelModel::runBranch(const std::string &prefix, const std::vector<float> &in,
                                int inC, int inH, int inW, std::vector<float> &out) const {
    // conv indices inside the Sequential: 0, 5, 10 (BN at +1, pool after 1st/2nd)
    const int convIdx[3] = {0, 5, 10};
    const int outC[3] = {16, 32, 64};

    std::vector<float> cur = in;
    int curC = inC, curH = inH, curW = inW;

    for (int stage = 0; stage < 3; ++stage) {
        char wn[128], bn[128], g[128], be[128], rm[128], rv[128];
        snprintf(wn, sizeof(wn), "%s.net.%d.weight", prefix.c_str(), convIdx[stage]);
        snprintf(bn, sizeof(bn), "%s.net.%d.bias",   prefix.c_str(), convIdx[stage]);
        snprintf(g,  sizeof(g),  "%s.net.%d.weight", prefix.c_str(), convIdx[stage] + 1);
        snprintf(be, sizeof(be), "%s.net.%d.bias",   prefix.c_str(), convIdx[stage] + 1);
        snprintf(rm, sizeof(rm), "%s.net.%d.running_mean", prefix.c_str(), convIdx[stage] + 1);
        snprintf(rv, sizeof(rv), "%s.net.%d.running_var",  prefix.c_str(), convIdx[stage] + 1);

        std::vector<float> conv;
        conv2d(cur, curC, curH, curW, get(wn)->data, get(bn)->data,
               outC[stage], 3, 3, conv);
        curC = outC[stage];
        batchNormRelu(conv, curC, curH * curW, get(g)->data, get(be)->data,
                      get(rm)->data, get(rv)->data, true);
        if (stage < 2) {
            std::vector<float> pooled;
            int nh = 0, nw = 0;
            maxPool2(conv, curC, curH, curW, pooled, nh, nw);
            cur.swap(pooled);
            curH = nh; curW = nw;
        } else {
            cur.swap(conv);
        }
    }
    adaptiveAvgPool(cur, curC, curH, curW, 4, 4, out);   // 64*4*4 = 1024
}

float SomaticIndelModel::predict(const float *image, const float *scalars) const {
    if (!loaded) return 1.0f;                 // fail open: keep the call

    const int HW = MAX_READS * WIDTH;

    // ---- input_norm: BatchNorm2d(6, affine=False) --------------------------
    const std::vector<float> &imean = get("input_norm.running_mean")->data;
    const std::vector<float> &ivar  = get("input_norm.running_var")->data;
    std::vector<float> x((size_t)IN_CHANNELS * HW);
    for (int c = 0; c < IN_CHANNELS; ++c) {
        const float inv = 1.0f / std::sqrt(ivar[c] + BN_EPS);
        const float m = imean[c];
        const float *src = image + (size_t)c * HW;
        float *dst = &x[(size_t)c * HW];
        for (int i = 0; i < HW; ++i) dst[i] = (src[i] - m) * inv;
    }

    // ---- conv branches ----------------------------------------------------
    std::vector<float> tumorIn(x.begin(), x.begin() + (size_t)PER_SAMPLE * HW);
    std::vector<float> normalIn(x.begin() + (size_t)PER_SAMPLE * HW, x.end());
    std::vector<float> tFeat, nFeat;
    runBranch("tumor_branch",  tumorIn,  PER_SAMPLE, MAX_READS, WIDTH, tFeat);
    runBranch("normal_branch", normalIn, PER_SAMPLE, MAX_READS, WIDTH, nFeat);

    // fixed branch gate (checkpoint: [1.0, 1.0], i.e. identity)
    const std::vector<float> &gate = get("branch_gate_fixed")->data;
    for (size_t i = 0; i < tFeat.size(); ++i) tFeat[i] *= gate[0];
    for (size_t i = 0; i < nFeat.size(); ++i) nFeat[i] *= gate[1];

    // ---- scalar branch ----------------------------------------------------
    const std::vector<float> &smean = get("scalar_norm.running_mean")->data;
    const std::vector<float> &svar  = get("scalar_norm.running_var")->data;
    std::vector<float> s(N_SCALAR);
    for (int i = 0; i < N_SCALAR; ++i)
        s[i] = (scalars[i] - smean[i]) / std::sqrt(svar[i] + BN_EPS);
    std::vector<float> h1, sFeat;
    linear(s, get("scalar_mlp.0.weight")->data, get("scalar_mlp.0.bias")->data,
           N_SCALAR, 32, h1, true);
    linear(h1, get("scalar_mlp.2.weight")->data, get("scalar_mlp.2.bias")->data,
           32, 32, sFeat, true);

    // ---- head -------------------------------------------------------------
    std::vector<float> z;
    z.reserve(tFeat.size() + nFeat.size() + sFeat.size());
    z.insert(z.end(), tFeat.begin(), tFeat.end());
    z.insert(z.end(), nFeat.begin(), nFeat.end());
    z.insert(z.end(), sFeat.begin(), sFeat.end());

    std::vector<float> a, b, logits;
    linear(z, get("fc1.weight")->data, get("fc1.bias")->data, (int)z.size(), 64, a, true);
    linear(a, get("fc2.weight")->data, get("fc2.bias")->data, 64, 32, b, true);
    linear(b, get("fc3.weight")->data, get("fc3.bias")->data, 32, 2, logits, false);

    const float m = std::max(logits[0], logits[1]);
    const float e0 = std::exp(logits[0] - m), e1 = std::exp(logits[1] - m);
    return e1 / (e0 + e1);                       // P(class 1 = true positive)
}
