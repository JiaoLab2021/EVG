#ifndef VCF_RECALL_HPP
#define VCF_RECALL_HPP

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
#include "zlib.h"
#include <getopt.h>
#include <sstream>

#include "vcf_open.hpp"
#include "save.hpp"
#include "strip_split_join.hpp"
#include "sort_find.hpp"
#include "get_time.hpp"
#include "check_vcf_sort.hpp"


using namespace std;



void help_recall(char** argv);
int main_recall(int argc, char** argv);


struct vcfStructure
{
    map<string,vector<uint32_t> > refStartVecMap;  // map<chr,vector<refStart> >
    
    map<string, map<uint32_t, tuple<uint32_t, vector<uint32_t>, string, int32_t, vector<int> > > > chrStartLenInfoGtTupMap;  // unordered_map<chr, map<refStart, tuple<refLen, qryLenVec, vcfInfo, svLen, gtVec> > >
    
    vector<int64_t> allLengthList;  // All the variations in length
};

// Bool function that performs statistics on the variation length
class f_mod_r
{
private:
    int64_t dv1;
    int64_t dv2;

public:
    f_mod_r(int64_t d1 = 1, int64_t d2 = 2) : dv1(d1), dv2(d2) {}

    bool operator() (int64_t x) {return  dv1 <= x && x <= dv2;}
};


class VCFRecall
{
private:
    // input
    string trueVCFFileName_;
    string evaluateFileName_;
    string rocKey_;
    string rocType_;

    vcfStructure trueVCFStrcuture_;

    // ROC -> save the number of DP or GQ
    int rocColNum_;
    vector<string> rocKeyVec_;

    // Variable length range
    vector<string> lengthVec_;

    // output
    static const int32_t bufferSize_ = 10 * 1024 * 1024;
    map<float, vector<int32_t> > rocCallMap_;  // map<DP/GQ, vector<length>> The software identifies all variants
    map<float, vector<int32_t> > rocRecallMap_;  // map<DP/GQ, vector<length>> The software identifies the correct variants
    
public:
    /**
	 * init
	 *
	 * @param trueVCFFileName    true set
     * @param evaluateFileName   evaluation set
     * @param rocKey             Keywords used to extract roc score
     * @param rocType            roc calculation rules (recall/genotype)
     * 
	**/
    VCFRecall(
        const string& trueVCFFileName, 
        const string& evaluateFileName, 
        const string& rocKey, 
        const string& rocType
    );


    /**
     * build the index of true VCF file
     * 
     * @return void
    **/
    void build_true_index();


    /**
	 * Convert the results of gramtools
     * 
     * @return void
	**/
    void gramtools_convert();


    /**
	 * genotype evaluation function
     * 
     * @return void
	**/
    void evulate_gt();


    /**
	 * genotype evaluation function
     * 
     * @return void
	**/
    void evulate_recall();


    /**
	 * Get a list of loci genotypes.
	 *
	 * @param informationsVec  vcfInfoList
     * @param sampleIdx        Index of the sample genotype, with the default value 0 representing the last column
     * 
     * 
     * @return gtVec           vector <int>
	**/
    vector<int> get_gt(
        const vector<string> & informationsVec, 
        int sampleIdx = 0
    );


    /**
	 * Get a list of loci genotypes.
	 *
     * @param INFOSTRUCTTMP    line information
     * 
     * 
     * @return rocNum
	**/
    float get_roc(
        const VCFINFOSTRUCT& INFOSTRUCTTMP
    );


    /**
	 * Get the length information of the variation
	 *
     * @param refLen            ref length
	 * @param qryLenVec         qry Length list
     * @param gtVec             Genotype list
     * 
     * 
     * @return int32_t              svLength
	**/
    int32_t sv_length_select(
        const uint32_t & refLen, 
        const vector<uint32_t> & qryLenVec, 
        const vector<int> & gtVec
    );


    /**
	 * Statistical variation length information
	 *
	 * @param length_list         Length list
     * 
     * 
     * @return vector<int64_t>   The length of each interval
	**/
    vector<int64_t> count_num(
        vector<int64_t> length_list
    );


    /**
	 * Get the length of the haplotype.
	 *
     * @param svType                         Variation type
	 * @param refSeq                         ref column information
     * @param qrySeqs                        qry Column information
     * @param gtVec                          Typing information of loci
     * @param lenType                        Take only the length corresponding to the haplotype or all the lengths (hap/all)
     * 
     * 
     * @return tuple<int, vector<int> >      tuple<refLen, vector<qryLen> >
	**/
    tuple<uint32_t, vector<uint32_t> > get_hap_len(
        const string & svType, 
        const string & refSeq, 
        const string & qrySeqs, 
        const vector<int> & gtVec, 
        const string & lenType
    );
    

    /**
	 * Save the result of failCall in recall
	 *
	 * @param chrStartLenInfoGtTupMap     The rest of the vcf information is used in the true set
     * @param outFileName                 Output file name
     * 
     * 
     * @return int             0
	**/
    int saveFailCall(
        map<string, map<uint32_t, tuple<uint32_t, vector<uint32_t>, string, int32_t, vector<int> > > > & chrStartLenInfoGtTupMap, 
        const string & outFileName
    );


    /**
	 * Save the result of failCall in recall
	 *
     * @param allLengthList     All variation lengths
     * 
     * 
     * @return int             0
	**/
    void roc_calculate(
        const vector<int64_t> & allLengthList
    );
};

#endif