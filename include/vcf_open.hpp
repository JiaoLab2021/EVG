#ifndef VCF_OPEN_HPP
#define VCF_OPEN_HPP

#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <map>
#include <tuple>
#include <iomanip>
#include "zlib.h"
#include <sstream>
#include <vector>
#include <iterator>


#include <unordered_map>
#include "strip_split_join.hpp"
#include "get_time.hpp"

using namespace std;

struct VCFINFOSTRUCT
{
    public:
    // general information
    string line;
    vector<string> lineVec;
    vector<string> ALTVec;  // [4]

    string CHROM;  // [0]
    uint32_t POS;  // [1]
    string ID;  // [2]
    string REF;  // [3]
    string ALT;  // [4]
    double QUAL;  // [5]
    string FILTER;  // [6]
    string INFO;  // [7]
    string FORMAT;  // [8]

    uint32_t LEN;  // REF length
    uint32_t END;  // ALT length

    VCFINFOSTRUCT() : line(""), CHROM(""), POS(0), ID(""), REF(""), ALT(""), QUAL(0.0), FILTER(""), INFO(""), FORMAT("") {}

    // clear
    void clear();
};


/**
 * @brief Open the vcf file
 * 
 * @param vcfFileName   the output of vcf file
 * 
 * @return
**/
class VCFOPEN
{
private:
    // vcf file
    string vcfFileName_;

    // Input file stream
    gzFile gzfpI;

    // buffer size
    uint32_t bufferSize_ ;
    char *line_;
public:
    VCFOPEN(
        const string & vcfFileName
    );

    ~VCFOPEN();

    bool read(
        VCFINFOSTRUCT & INFOSTRUCTTMP
    );

    /**
     * Get the type of variant
     * 
     * @param refLen  the length of REF variant
     * @param ALTVec  vector<qrySeq>
     * 
     * @return string   TYPE: SNP, InDel, Deletion, Insertion, Inversion, Duplication, Other
    **/
    string get_TYPE(
        const uint32_t & refLen,
        const vector<string> & ALTVec
    );

    /**
     * Get a list of locus genotypes.
     *
     * @param lineVec  lineVec
     * 
     * 
     * @return GTVecMap   map<int, vector<string> >,  map<idx, vector<GTString> >
    **/
    map<int, vector<string> > get_gt(
        const vector<string> & lineVec
    );

    /**
     * split gt
     *
     * @param gtTxt
     * 
     * 
     * @return vector<gt>
    **/
    vector<string> gt_split(const string & gtTxt);

    
    /**
     * Calculate MAF and MISSRATE.
     *
     * @param GTVecMap    map<int, vector<string> >,  map<idx, vector<GTString> >
     * @param sampleNum   total number of samples
     * 
     * 
     * @return tuple<double, double>   tuple<MAF, MISSRATE>
    **/
    tuple<double, double> calculate(
        const map<int, vector<string> > & GTVecMap, 
        uint32_t sampleNum
    );
};

#endif