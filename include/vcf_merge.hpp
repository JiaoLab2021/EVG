#ifndef VCF_MERGE_HPP
#define VCF_MERGE_HPP

#include <getopt.h>
#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <vector>
#include <numeric>
#include <map>
#include <malloc.h>
#include "zlib.h"
#include <cmath>
#include <sstream>

#include "vcf_open.hpp"
#include "save.hpp"
#include "strip_split_join.hpp"
#include "sort_find.hpp"
#include "get_time.hpp"
#include "check_vcf_sort.hpp"


using namespace std;

void help_merge(char** argv);

int main_merge(int argc, char** argv);


 struct baseVcfStruct
{
    string headInfo; // Comment lines of the vcf

    map<string, vector<uint32_t> > startMap;  // map<chromosome, vector<start>>
    map<string, vector<uint32_t> > refLenMap;  // map<chromosome, vector<refLen>>
    map<string, vector<vector<uint32_t> > > qryLenVecMap;  // map<chromosome, vector<vector<qryLen> > >, When a locus has multiple alleles, this is used to preserve the length of the qry

    // Store the vcfInfo of baseVcf
    map<string, map<uint32_t, string> > BaseInfoMap;  // map<chromosome, map<start, vcfInfo>>

    // Store the results after binary search method, save the software name and the corresponding typing results
    map<string, map<uint32_t, map<string, tuple<vector<float>, vector<int> > > > > recallSoftwareGtDepVecMap;  // map<chr, map<start, map<software, tuple<vector<depth>, vector<gt> > > > >
    
    map<string, tuple<float, float, float>> depthMap;  // map<software, pair<aveDepth, variance, sd>> The last two are the variance and the standard deviation
    map<string, int> softwateSampleIdx;  // map<software, index>, Index of the sample corresponding to the software file

    // Records the number of columns used to find elements after that column when Paragraph is extracted
    int colNum;

    baseVcfStruct() : colNum(0) {}
};

struct softwareVcfStruct
{
    string software;  // Software name

    map<string, vector<uint32_t> > startMap;  // map<chromosome, vector<start>>
    map<string, vector<uint32_t> > refLenMap;  // map<chromosome, vector<refLen>>
    map<string, vector<vector<uint32_t> > > qryLenVecMap;  // map<chromosome, vector<vector<qryLen> > >, probably 1/2, etc. Therefore, the length of each allele should be stored. If the length was 1, it represented homozygous variation
    map<string, vector<vector<int> > > gtVecMap;  // map<chromosome, vector<vector<genotype>>>
    map<string, vector<vector<float> > > depthVecMap;  // map<chromosome, vector<vector<depth>>>

    vector<float> depthVec;  // vector<depth>
};



class VCFMerge
{
private:
    // EVG running mode
    string mode_;

    string trueVcf_;

    string ParagraphVcf_;
    string GraphTyper2Vcf_;
    string BayesTyperVcf_;
    string MAPVcf_;
    string GiraffeVcf_;
    string GraphAlignerVcf_;
    string PanGenieVcf_;

    string sampleName_;

    uint16_t softwareNum_;  // Records the amount of software used to filter the results

    map<string, map<int, string > > outChrStartInfoMap_;  // outMap<chromosome, map<refStart, vcfInfo> >

    baseVcfStruct mergeVcfStruct_;  // recore the index of all software

    string outputFileName_;


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
     * @return tuple<uint32_t, vector<uint32_t> >      tuple<refLen, vector<qryLen> >
    **/
    tuple<uint32_t, vector<uint32_t> > get_hap_len(
        const string & svType, 
        const string & refSeq, 
        const string & qrySeqs, 
        const vector<int> & gtVec, 
        const string & lenType
    );


    /**
     * The results of the software are filtered and indexed, and the files are the output files of the software.
     *
     * @date 2023/07/09
     * 
     * @param vcfFileName        Software output vcf file
     * @param software           Software name
     * 
     * 
     * @return void
    **/
    void build_softwarefile_index(
        const string & vcfFileName, 
        const string & software
    );


    // paragraph DP -> 1
    // graphaligner, vg-map- vg-giraffe DP -> 6
    // bayestyper MAC -> -1,1.87332
    // graphtyper DP -> 6
    // pangenie KC -> 4
    /**
     * Get site DP.
     *
     * @param informationsVec     vcfInfoVec
     * @param sampleIdx           sample corresponding column
     * 
     * 
     * @return depthVec        vector<depth>
    **/
    vector<float> get_depth(
        const vector<string> & informationsVec, 
        const uint32_t& sampleIdx
    );


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
        uint32_t sampleIdx
    );



    /**
     * vcf_merge Function added to the total map after recall
     *
     * @date 2023/07/09
     * 
     * @param chromosome         Chromosome number
     * @param trueRefStart       The actual ref starting position
     * @param software           software
     * @param qryLenVec          The haplotype of software typing corresponds to the length
     * @param gtVec              Classification gt list of software
     * @param gt                 The current cycle of gt
     * @param depth              The depth corresponding to the site
     * @param j                  Index of multiple haplotype loops of software
     * @param k                  Binary search method index of two coordinate loops
     * @param l                  Cyclic index of true site haplotype
     * @param indexLeft          Left index of binary search method
     * @param indexRight         Right index of binary search method
     * 
     * 
     * @return tuple<int, int>   tuple<state, j>  0 -> Add correctly; -1 -> The next loop is required
    **/
    tuple<int, int> recall_push(
        const string chromosome, 
        const uint32_t trueRefStart, 
        const string software, 
        const vector<uint32_t> qryLenVec, 
        const vector<int> gtVec, 
        const int gt, 
        const float depth, 
        uint32_t j, 
        uint32_t k,
        uint32_t l,
        int64_t & indexLeft, 
        int64_t & indexRight
    );


    /**
     * The results of the software are filtered and indexed, and the files are the output files of the software.
     *
     * @date 2023/07/09
     * 
     * @param softvcfStructure   build_softwarefile_index Indicates the index of the built software vcf
     * 
     * @return void
    **/
    void vcf_merge(
        softwareVcfStruct & softvcfStructure
    );


    /**
     * Get the length information of the variation
     *
     * @param vcfInfo           vcfInfo
     * @param gtVec             Genotype list
     * 
     * 
     * @return int              svLength
    **/
    int64_t sv_length_select(
        const string & vcfInfo, 
        const vector<int> & gtVec
    );


    /**
     * Filter result according to the depth and number of variants.
     *
     * @param gtAllVec        vector of all gt sites
     * @param softwareVec     Locus vector of all software
     * @param depthNorVec     The vector of all normalized depths at the site
     * @param vcfInfo         Variation information at the site
     * 
     * @return software
    **/
    string filter_rule(
        vector<vector<int>> & gtAllVec, 
        vector<string> & softwareVec, 
        vector<float> & depthNorVec, 
        const string & vcfInfo
    );


    /**
     * Calculate the variance and standard deviation
     *
     * @param data     List with site depth
     * 
     * @return tuple   make_tuple(mean, variance, std_deviation)
    **/
    tuple<float, float, float> cal_var_sd(const vector<float>& data);
public:
    /**
	 * init
	 *
     * @param mode
	 * @param trueVcf
     * @param ParagraphVcf
     * @param GraphTyper2Vcf
     * @param BayesTyperVcf
     * @param MAPVcf
     * @param GiraffeVcf
     * @param GraphAlignerVcf
     * @param PanGenieVcf
     * @param sampleName       Sample names to be merged
     * @param outputFileName
     * 
	**/
    VCFMerge(
        const string& mode,
        const string& trueVcf, 
        const string& ParagraphVcf, 
        const string& GraphTyper2Vcf, 
        const string& BayesTyperVcf, 
        const string& MAPVcf, 
        const string& GiraffeVcf, 
        const string& GraphAlignerVcf, 
        const string& PanGenieVcf, 
        const string& sampleName, 
        const string& outputFileName
    );


    /**
     * The base file index is built without filtering the vcf and is used to create the base index. The file is a real vcf file
     * 
     * 
     * @return void
    **/
    void build_basefile_index();


    /**
	 * Building index for different software.
     * 
     * 
     * @return softwareNum_
	**/
    void run_index_merge();


    /**
     * Select software with coverage close to 0
     * mergeVcfStruct_ -> struct after software merger
    **/

    /**
     * Filter and merge result according to the depth and number of variants.
     *
     * @param mode_               EVG running mode
     * 
     * @return void
    **/
   void vcf_merge_filter();


    /**
     * save result
     * mergeVcfStruct_ -> struct after software merger
     * outChrStartInfoMap_ -> vcf_merge_filter Output result   outChrStartInfoMap_[chromosome][refStart][refLen][qryLen]
     * prefix -> Output file name prefix
    **/

    /**
	 * save result
	 *
     * 
     * @return void
	**/
    void result_save();

};

#endif