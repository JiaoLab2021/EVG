#pragma once

#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <vector>
#include <map>
#include <tuple>
#include <unordered_map>
#include <algorithm>
#include <iomanip>
#include <cstdlib>
#include <regex>
#include <sstream>

#include "vcf_open.hpp"
#include "save.hpp"
#include "strip_split_join.hpp"
#include "sort_find.hpp"

using namespace std;

// Index of the truth VCF.
struct vcfStructure
{
    map<string,vector<uint32_t> > refStartVecMap;  // map<chromosome, vector<start> >

    map<string, map<uint32_t, tuple<uint32_t, vector<uint32_t>, string, int32_t, vector<int> > > > chrStartLenInfoGtTupMap;  // map<chromosome, map<start, tuple<refLen, qryLenVec, vcfInfo, svLen, gtVec> > >

    vector<int64_t> allLengthList;
};

// True if a length lies in [d1, d2].
class f_mod_r
{
private:
    int64_t dv1;
    int64_t dv2;

public:
    f_mod_r(int64_t d1 = 1, int64_t d2 = 2) : dv1(d1), dv2(d2) {}

    bool operator() (int64_t x) {return  dv1 <= x && x <= dv2;}
};

// Compare the calls of a genotyped VCF with a truth VCF and write evaluation tables.
class VCFRecall
{
private:
    string trueVCFFileName_;
    string evaluateFileName_;
    string rocKey_;
    string rocType_;
    uint32_t threads_;

    vcfStructure trueVCFStrcuture_;

    int rocColNum_;
    vector<string> rocKeyVec_;
    regex rocInfoReg_;  // <key>=(\d+) for an INFO ROC field

    vector<string> lengthVec_;

    static const int32_t bufferSize_ = 10 * 1024 * 1024;
    map<float, vector<int32_t> > rocCallMap_;  // map<score, vector<length> > of all calls
    map<float, vector<int32_t> > rocRecallMap_;  // map<score, vector<length> > of true calls

public:
    VCFRecall(
        const string& trueVCFFileName,
        const string& evaluateFileName,
        const string& rocKey,
        const string& rocType,
        uint32_t threads
    );

    void build_true_index();

    void index_true_record(const VCFINFOSTRUCT & INFOSTRUCTTMP);

    void convert_gramtools();

    void evaluate_genotype();

    void evaluate_recall();

    vector<int> get_gt(
        const vector<string> & informationsVec,
        int sampleIdx = 0
    ) const;

    float get_roc(
        const VCFINFOSTRUCT& INFOSTRUCTTMP
    );

    int32_t sv_length_select(
        const uint32_t & refLen,
        const vector<uint32_t> & qryLenVec,
        const vector<int> & gtVec
    ) const;

    vector<int64_t> count_num(
        vector<int64_t> length_list
    );

    int saveFailCall(
        map<string, map<uint32_t, tuple<uint32_t, vector<uint32_t>, string, int32_t, vector<int> > > > & chrStartLenInfoGtTupMap,
        const string & outFileName
    );

    void roc_calculate(
        const vector<int64_t> & allLengthList
    );
};
