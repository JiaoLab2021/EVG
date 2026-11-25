#pragma once

#include <fstream>
#include <string>
#include <iostream>
#include <numeric>
#include <cstring>
#include <map>
#include <cmath>

#include "vcf_open.hpp"
#include "save.hpp"
#include "strip_split_join.hpp"
#include "sort_find.hpp"

using namespace std;

// Count the variants of every sample in sliding windows, per group.
class VCFWin
{
private:
    uint32_t windowSize_, stepSize_;

    string chrLenFilename_;
    string vcfFilename_;
    string waysFilename_;

    string outputFileName_;

    string mode_;

    uint32_t threads_;

    map<string, uint32_t> chrLenMap_;  // map<chromosome, length>

    map<string, vector<uint32_t> > winStartMap_;  // map<chromosome, vector<window start> >
    map<string, vector<uint32_t> > winEndMap_;  // map<chromosome, vector<window end> >

    map<string,string> waysGroupMap_;  // map<sample, group>

    map<string,map<uint32_t, vector<uint32_t> > > chrWinStartLenMap_;  // map<chromosome, map<window start, vector<length> > >
    map<string,map<uint32_t, vector<string> > > chrWinStartGroupMap_;  // map<chromosome, map<window start, vector<group> > >

    vector<int> get_gt(
        const vector<string>& lineVec,
        int sampleIdx = 0
    );
public:
    VCFWin(
        const uint32_t& windowSize,
        const uint32_t& stepSize,
        const string& chrLenFilename,
        const string& vcfFilename,
        const string& waysFilename,
        const string& outputFileName,
        const string& mode,
        uint32_t threads
    );

    void load_chromosome_lengths();

    void build_windows();

    void load_groups();

    void count_window_variants();

    void save_result();
};
