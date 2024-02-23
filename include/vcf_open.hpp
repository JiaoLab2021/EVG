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

    VCFINFOSTRUCT() : POS(0), QUAL(0.0), LEN(0), END(0) {}

    // clear
    void clear() {
        VCFINFOSTRUCT temp;
        swap(temp);
    }

    // swap
    void swap(VCFINFOSTRUCT& other) {
        using std::swap;
        swap(line, other.line);
        swap(lineVec, other.lineVec);
        swap(ALTVec, other.ALTVec);
        swap(CHROM, other.CHROM);
        swap(POS, other.POS);
        swap(ID, other.ID);
        swap(REF, other.REF);
        swap(ALT, other.ALT);
        swap(QUAL, other.QUAL);
        swap(FILTER, other.FILTER);
        swap(INFO, other.INFO);
        swap(FORMAT, other.FORMAT);
        swap(LEN, other.LEN);
        swap(END, other.END);
    }

    // assignment operator overloading
    VCFINFOSTRUCT& operator=(const VCFINFOSTRUCT& other) {
        if (this != &other) { // Prevent self-assignment
            // Copy member variables one by one
            line = other.line;
            lineVec = other.lineVec;
            ALTVec = other.ALTVec;
            CHROM = other.CHROM;
            POS = other.POS;
            ID = other.ID;
            REF = other.REF;
            ALT = other.ALT;
            QUAL = other.QUAL;
            FILTER = other.FILTER;
            INFO = other.INFO;
            FORMAT = other.FORMAT;
            LEN = other.LEN;
            END = other.END;
        }
        return *this;
    }

    // Move constructor
    VCFINFOSTRUCT(VCFINFOSTRUCT&& other) noexcept
        : line(std::move(other.line)),
        lineVec(std::move(other.lineVec)),
        ALTVec(std::move(other.ALTVec)),
        CHROM(std::move(other.CHROM)),
        POS(other.POS),
        ID(std::move(other.ID)),
        REF(std::move(other.REF)),
        ALT(std::move(other.ALT)),
        QUAL(other.QUAL),
        FILTER(std::move(other.FILTER)),
        INFO(std::move(other.INFO)),
        FORMAT(std::move(other.FORMAT)),
        LEN(other.LEN),
        END(other.END) {
        // Optionally clear or reset the state of 'other'
    }

    // Move assignment operator
    VCFINFOSTRUCT& operator=(VCFINFOSTRUCT&& other) noexcept {
        if (this != &other) {
            line = std::move(other.line);
            lineVec = std::move(other.lineVec);
            ALTVec = std::move(other.ALTVec);
            CHROM = std::move(other.CHROM);
            POS = other.POS;
            ID = std::move(other.ID);
            REF = std::move(other.REF);
            ALT = std::move(other.ALT);
            QUAL = other.QUAL;
            FILTER = std::move(other.FILTER);
            INFO = std::move(other.INFO);
            FORMAT = std::move(other.FORMAT);
            LEN = other.LEN;
            END = other.END;
            // Optionally clear or reset the state of 'other'
        }
        return *this;
    }
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
    size_t bufferSize_ ;
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
     * @return GTVecMap, misSampleNum   map<idx, vector<GTString> >
    **/
    tuple<unordered_map<int, vector<string> >, uint32_t> get_gt(
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
    vector<string> gt_split(string & gtTxt);

    
    /**
     * Calculate MAF and MISSRATE.
     *
     * @param GTVecMap      map<int, vector<string> >,  map<idx, vector<GTString> >
     * @param misSampleNum  number of missing samples
     * 
     * 
     * @return tuple<double, double>   tuple<MAF, MISSRATE>
    **/
    tuple<double, double> calculate(
        const unordered_map<int, vector<string> > & GTVecMap, 
        uint32_t misSampleNum
    );
};

#endif