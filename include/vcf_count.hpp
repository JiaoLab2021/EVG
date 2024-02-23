#ifndef VCF_COUNT_HPP
#define VCF_COUNT_HPP

#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include "zlib.h"
#include <unordered_map>
#include <cmath>

#include "strip_split_join.hpp"
#include "get_time.hpp"
#include "sort_find.hpp"
#include "vcf_open.hpp"
#include <getopt.h>
#include <cstdlib>


using namespace std;


void help_count(char** argv);
int main_count(int argc, char** argv);

// Bool function that performs statistics on the variation length
class f_mod_c
{
private:
    int64_t dv1;
    int64_t dv2;

public:
    f_mod_c(int64_t d1 = 1, int64_t d2 = 2) : dv1(d1), dv2(d2) {}

    bool operator() (int64_t x) {return  dv1 <= x && x <= dv2;}
};

class VCFCount
{
private:
    string vcfFileName_;

    // Variable length range
    vector<string> lengthVec_;

    // Record the number of variations
    uint64_t snpNum_ = 0;
    uint64_t indelNum_ = 0;
    uint64_t insNum_ = 0;
    uint64_t delNum_ = 0;
    uint64_t invNum_ = 0;
    uint64_t dupNum_ = 0;
    uint64_t otherNum_ = 0;

    // Allelic Number
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

    // Record the length of variations
    vector<int32_t> svLenVec_;
    vector<int32_t> InvLenVec_;
    vector<int32_t> DupLenVec_;

    // Length statistics
    vector<uint64_t> svNumVec_;
    vector<uint64_t> svTotalLenVec_;
    vector<uint64_t> InvNumVec_;
    vector<uint64_t> InvTotalLenVec_;
    vector<uint64_t> DupNumVec_;
    vector<uint64_t> DupTotalLenVec_;

public:
    /**
	 * init
	 *
	 * @param vcfFileName       input VCF  file name
     * 
	**/
    VCFCount(const string & vcfFileName);

    void count();

    void count_num_len();

    void get_result();
};

#endif