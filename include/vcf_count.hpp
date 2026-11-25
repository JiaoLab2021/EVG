#pragma once

#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <unordered_map>
#include <cmath>

#include <cstdlib>

#include "strip_split_join.hpp"
#include "vcf_open.hpp"

using namespace std;

// True if a length lies in [d1, d2].
class f_mod_c
{
private:
    int64_t dv1;
    int64_t dv2;

public:
    f_mod_c(int64_t d1 = 1, int64_t d2 = 2) : dv1(d1), dv2(d2) {}

    bool operator() (int64_t x) {return  dv1 <= x && x <= dv2;}
};

// Count alleles by variant type and length bin.
class VCFCount
{
private:
    string vcfFileName_;
    uint32_t threads_;

    vector<string> lengthVec_;

    uint64_t snpNum_ = 0;
    uint64_t indelNum_ = 0;
    uint64_t insNum_ = 0;
    uint64_t delNum_ = 0;
    uint64_t invNum_ = 0;
    uint64_t dupNum_ = 0;
    uint64_t otherNum_ = 0;

    uint64_t snpBiAllelicNum_ = 0;
    uint64_t snpMuAllelicNum_ = 0;
    uint64_t indelBiAllelicNum_ = 0;
    uint64_t indelMuAllelicNum_ = 0;
    uint64_t insBiAllelicNum_ = 0;
    uint64_t insMuAllelicNum_ = 0;
    uint64_t delBiAllelicNum_ = 0;
    uint64_t delMuAllelicNum_ = 0;
    uint64_t invBiAllelicNum_ = 0;
    uint64_t invMuAllelicNum_ = 0;
    uint64_t dupBiAllelicNum_ = 0;
    uint64_t dupMuAllelicNum_ = 0;
    uint64_t otherBiAllelicNum_ = 0;
    uint64_t otherMuAllelicNum_ = 0;

    vector<int32_t> svLenVec_;
    vector<int32_t> InvLenVec_;
    vector<int32_t> DupLenVec_;

    vector<uint64_t> svNumVec_;
    vector<uint64_t> svTotalLenVec_;
    vector<uint64_t> InvNumVec_;
    vector<uint64_t> InvTotalLenVec_;
    vector<uint64_t> DupNumVec_;
    vector<uint64_t> DupTotalLenVec_;

public:
    VCFCount(const string & vcfFileName, uint32_t threads);

    void count_variants();

    void count_length_bins();

    void print_result();
};
