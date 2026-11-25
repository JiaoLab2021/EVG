#pragma once

#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <vector>
#include <numeric>
#include <map>
#include <cmath>
#include <sstream>
#include <tuple>
#include <unordered_map>

#include "vcf_open.hpp"
#include "save.hpp"
#include "strip_split_join.hpp"
#include "sort_find.hpp"

using namespace std;

// Index of the population VCF and the matched tool calls.
struct baseVcfStruct {
    string headInfo;

    map<string, vector<uint32_t> > startMap;  // map<chromosome, vector<start> >
    map<string, vector<uint32_t> > refLenMap;  // map<chromosome, vector<refLen> >
    map<string, vector<vector<uint32_t> > > qryLenVecVecMap;  // map<chromosome, vector<vector<qryLen> > >

    vector<string> sampleNameVec;

    unordered_map<string, unordered_map<uint32_t, tuple<string, vector<string> > > > BaseInfoTupMap;  // map<chromosome, map<start, tuple<vcfInfo, vector<gt> > > >

    map<string, map<uint32_t, map<string, tuple<vector<float>, vector<int> > > > > recallSoftwareGtDepVecMap;  // map<chromosome, map<start, map<software, tuple<vector<depth>, vector<gt> > > > >

    map<string, tuple<float, float, float>> depthMap;  // map<software, tuple<mean, variance, sd> >
    map<string, int> softwateSampleIdx;  // map<software, sample column>

    int colNum;

    baseVcfStruct() : colNum(0) {}
};

// Index of the calls of one tool, per chromosome.
struct softwareVcfStruct {
    string software;

    map<string, vector<uint32_t> > startMap;  // map<chromosome, vector<start> >
    map<string, vector<uint32_t> > refLenMap;  // map<chromosome, vector<refLen> >
    map<string, vector<vector<uint32_t> > > qryLenVecVecMap;  // map<chromosome, vector<vector<qryLen of each haplotype> > >
    map<string, vector<vector<int> > > gtVecMap;  // map<chromosome, vector<vector<genotype> > >
    map<string, vector<vector<float> > > depthVecMap;  // map<chromosome, vector<vector<depth> > >

    vector<float> depthVec;

    int sampleIdx = 0;  // column of the sample in the tool VCF
    tuple<float, float, float> depthStats;  // tuple<mean, variance, sd> of depthVec
};

// Merge the calls of several genome graph tools for one sample on the sites of the population VCF.
class VCFMerge
{
private:
    string trueVcf_;

    string ParagraphVcf_;
    string GraphTyper2Vcf_;
    string BayesTyperVcf_;
    string MAPVcf_;
    string GiraffeVcf_;
    string GraphAlignerVcf_;
    string PanGenieVcf_;

    string sampleName_;

    uint32_t threads_;

    uint16_t softwareNum_ = 0;
    bool paragraphInputBool_ = false;

    map<string, map<int, string > > outChrStartInfoMap_;  // map<chromosome, map<start, vcfInfo> >

    baseVcfStruct mergeVcfStruct_;

    uint32_t mostSimilarSampleIdx_ = 0;

    Outcome<baseVcfStruct> build_basefile_index() const;

    void index_base_header(const VCFINFOSTRUCT & INFOSTRUCTTMP, baseVcfStruct & index) const;

    void index_base_record(const VCFINFOSTRUCT & INFOSTRUCTTMP, baseVcfStruct & index) const;

    vector<string> get_gt_all(
        const vector<string> & lineVec
    ) const;

    Outcome<softwareVcfStruct> build_softwarefile_index(
        const string & vcfFileName,
        const string & software
    ) const;

    void index_software_header(const VCFINFOSTRUCT & INFOSTRUCTTMP, softwareVcfStruct & softvcfStructure) const;

    void index_software_record(const VCFINFOSTRUCT & INFOSTRUCTTMP, softwareVcfStruct & softvcfStructure) const;

    vector<float> get_depth(
        const vector<string> & informationsVec,
        const uint32_t& sampleIdx
    ) const;

    vector<int> get_gt_sample(
        const vector<string> & informationsVec,
        uint32_t sampleIdx
    ) const;

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

    void vcf_merge(
        softwareVcfStruct & softvcfStructure
    );

    int64_t sv_length_select(
        const string & vcfInfo,
        const vector<int> & gtVec
    );

    string filter_rule(
        vector<vector<int>> & gtAllVec,
        vector<string> & softwareVec,
        vector<float> & depthNorVec,
        const string & vcfInfo
    );

    tuple<float, float, float> cal_var_sd(const vector<float>& data) const;
public:
    string mMode;
    string mOutputFileName;

    VCFMerge(
        const string& trueVcf,
        const string& ParagraphVcf,
        const string& GraphTyper2Vcf,
        const string& BayesTyperVcf,
        const string& MAPVcf,
        const string& GiraffeVcf,
        const string& GraphAlignerVcf,
        const string& PanGenieVcf,
        const string& sampleName,
        const string& mode,
        const string& outputFileName,
        uint32_t threads
    );

    void run_index_merge();

    bool vcf_merge_filter(bool secondBool = false);

    void result_save();

    void extract_directory_FileName(const std::string& path, std::string& directory, std::string& fileName);
};
