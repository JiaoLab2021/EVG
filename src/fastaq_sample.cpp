#include "../include/fastaq_sample.hpp"

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>

#include <htslib/khash.h>
#include <htslib/kseq.h>

#include "../include/hts_file.hpp"
#include "../include/logger.hpp"
#include "../include/save.hpp"

KSEQ_INIT(InputFile*, read_input)

namespace {

KHASH_SET_INIT_INT64(64)

typedef uint64_t krint64_t;

// Append one record as FASTA/FASTQ text (seqtk stk_printseq without line wrapping) and write the text through
// output once it holds 1 MB.
void write_record(SAVE& output, std::string& text, const kseq_t *s) {
    text += s->qual.l ? '@' : '>';
    text += s->name.s;
    if (s->comment.l) {
        text += ' ';
        text += s->comment.s;
    }
    text += '\n';
    text += s->seq.s;
    text += '\n';
    if (s->qual.l) {
        text += "+\n";
        text += s->qual.s;
        text += '\n';
    }
    if (text.size() >= (1 << 20)) {
        output.write(text);
        text.clear();
    }
}

/*
64-bit Mersenne Twister pseudorandom number generator. Adapted from:
    http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/mt19937-64.c
which was written by Takuji Nishimura and Makoto Matsumoto and released
under the 3-clause BSD license.
*/

struct _krand_t;
typedef struct _krand_t krand_t;

#define KR_NN 312
#define KR_MM 156
#define KR_UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define KR_LM 0x7FFFFFFFULL /* Least significant 31 bits */

struct _krand_t {
    int mti;
    krint64_t mt[KR_NN];
};

// Initialise the generator state from a seed.
static void kr_srand0(krint64_t seed, krand_t *kr) {
    kr->mt[0] = seed;
    for (kr->mti = 1; kr->mti < KR_NN; ++kr->mti)
        kr->mt[kr->mti] = 6364136223846793005ULL * (kr->mt[kr->mti - 1] ^ (kr->mt[kr->mti - 1] >> 62)) + kr->mti;
}

// Allocate and seed a generator.
krand_t *kr_srand(krint64_t seed) {
    krand_t *kr;
    kr = (krand_t *)malloc(sizeof(krand_t));
    kr_srand0(seed, kr);
    return kr;
}

// Return the next 64-bit random number.
krint64_t kr_rand(krand_t *kr) {
    krint64_t x;
    static const krint64_t mag01[2] = { 0, 0xB5026F5AA96619E9ULL };
    if (kr->mti >= KR_NN) {
        int i;
        if (kr->mti == KR_NN + 1) kr_srand0(5489ULL, kr);
        for (i = 0; i < KR_NN - KR_MM; ++i) {
            x = (kr->mt[i] & KR_UM) | (kr->mt[i+1] & KR_LM);
            kr->mt[i] = kr->mt[i + KR_MM] ^ (x>>1) ^ mag01[(int)(x&1)];
        }
        for (; i < KR_NN - 1; ++i) {
            x = (kr->mt[i] & KR_UM) | (kr->mt[i+1] & KR_LM);
            kr->mt[i] = kr->mt[i + (KR_MM - KR_NN)] ^ (x>>1) ^ mag01[(int)(x&1)];
        }
        x = (kr->mt[KR_NN - 1] & KR_UM) | (kr->mt[0] & KR_LM);
        kr->mt[KR_NN - 1] = kr->mt[KR_MM - 1] ^ (x>>1) ^ mag01[(int)(x&1)];
        kr->mti = 0;
    }
    x = kr->mt[kr->mti++];
    x ^= (x >> 29) & 0x5555555555555555ULL;
    x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
    x ^= (x << 37) & 0xFFF7EEE000000000ULL;
    x ^= (x >> 43);
    return x;
}

#define kr_drand(_kr) ((kr_rand(_kr) >> 11) * (1.0/9007199254740992.0))

// Copy a non-empty string, growing the destination buffer if needed.
static void cpy_kstr(kstring_t *dst, const kstring_t *src) {
    if (src->l == 0) return;
    if (src->l + 1 > dst->m) {
        dst->m = src->l + 1;
        kroundup32(dst->m);
        dst->s = (char *)realloc(dst->s, dst->m);
    }
    dst->l = src->l;
    memcpy(dst->s, src->s, src->l + 1);
}

// Copy name, sequence, quality and comment of a record.
static void cpy_kseq(kseq_t *dst, const kseq_t *src) {
    cpy_kstr(&dst->name, &src->name);
    cpy_kstr(&dst->seq,  &src->seq);
    cpy_kstr(&dst->qual, &src->qual);
    cpy_kstr(&dst->comment, &src->comment);
}

}

// Write a random fraction of the records, or `frac` records when frac >= 1 (seqtk sample), to outputFileName or stdout;
// throw on any error.
void stk_sample(const std::string& inputFileName, const std::string& outputFileName, double frac, const long int& seedNum,
                const std::string& mode, uint32_t threads) {
    int twopass = 0;
    uint64_t i, num = 0, n_seqs = 0;
    kseq_t *seq;
    int status;
    krand_t *kr = 0;

    kr = kr_srand(seedNum);

    if (mode == "2") {
        twopass = 1;
    }

    if (frac >= 1.0) num = (uint64_t)(frac + .499), frac = 0.;
    else if (twopass) {
        warning_stream() << "when sampling a fraction, -m 2 is ignored.\n";
        twopass = 0;
    }

    if (!twopass) {
        // ------ One pass: print by fraction, or keep a reservoir of num records ------
        kseq_t *buf = 0;
        if (num > 0) buf = (kseq_t *)calloc(num, sizeof(kseq_t));
        if (num > 0 && buf == NULL) {
            free(kr);
            throw std::runtime_error("could not allocate enough memory for " + std::to_string(num) + " sequences.");
        }

        InputFile input(inputFileName, threads);
        SAVE output(outputFileName, threads);
        std::string text;
        seq = kseq_init(&input);
        n_seqs = 0;
        while ((status = kseq_read(seq)) >= 0) {
            double r = kr_drand(kr);
            ++n_seqs;
            if (num) {
                uint64_t y = n_seqs - 1 < num? n_seqs - 1 : (uint64_t)(r * n_seqs);
                if (y < num) cpy_kseq(&buf[y], seq);
            } else if (r < frac) write_record(output, text, seq);
        }
        output.write(text);
        text.clear();
        if (status < -1) throw std::runtime_error("'" + inputFileName + "' -> truncated or damaged FASTA/FASTQ file.");
        for (i = 0; i < num; ++i) {
            kseq_t *p = &buf[i];
            if (p->seq.l) write_record(output, text, p);
            free(p->seq.s); free(p->qual.s); free(p->name.s);
        }
        output.write(text);
        if (buf != NULL) free(buf);
        kseq_destroy(seq);
    } else {
        // ------ Two passes: pick record indices, then print those records ------
        uint64_t *buf;
        khash_t(64) *hash;
        int absent;

        buf = (uint64_t *)calloc(num, 8);
        if (buf == NULL) {
            free(kr);
            throw std::runtime_error("could not allocate enough memory for " + std::to_string(num) + " sequences.");
        }
        for (i = 0; i < num; ++i) buf[i] = UINT64_MAX;
        {
            InputFile input1(inputFileName, threads);
            seq = kseq_init(&input1);
            n_seqs = 0;
            while ((status = kseq_read(seq)) >= 0) {
                double r = kr_drand(kr);
                uint64_t y;
                ++n_seqs;
                y = n_seqs - 1 < num? n_seqs - 1 : (uint64_t)(r * n_seqs);
                if (y < num) buf[y] = n_seqs;
            }
            if (status < -1) throw std::runtime_error("'" + inputFileName + "' -> truncated or damaged FASTA/FASTQ file.");
            kseq_destroy(seq);
        }
        hash = kh_init(64);
        for (i = 0; i < num; ++i) kh_put(64, hash, buf[i], &absent);
        free(buf);

        InputFile input2(inputFileName, threads);
        SAVE output(outputFileName, threads);
        std::string text;
        seq = kseq_init(&input2);
        n_seqs = 0;
        while ((status = kseq_read(seq)) >= 0)
            if (kh_get(64, hash, ++n_seqs) != kh_end(hash))
                write_record(output, text, seq);
        output.write(text);
        kh_destroy(64, hash);
        if (status < -1) throw std::runtime_error("'" + inputFileName + "' -> truncated or damaged FASTA/FASTQ file.");
        kseq_destroy(seq);
    }

    free(kr);
}
