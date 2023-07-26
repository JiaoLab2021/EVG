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
    
    vector<uint32_t> allLengthList;  // 所有变异的长度
};

// 对变异长度进行统计的Bool函数
class f_mod
{
private:
    uint32_t dv1;
    uint32_t dv2;

public:
    f_mod(uint32_t d1 = 1, uint32_t d2 = 2) : dv1(d1), dv2(d2) {}

    bool operator() (uint32_t x) {return  dv1 <= x && x <= dv2;}
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
    int bufferSize_;
    map<float, vector<uint32_t> > rocCallMap_;  // map<DP/GQ, vector<length>> The software identifies all variants
    map<float, vector<uint32_t> > rocRecallMap_;  // map<DP/GQ, vector<length>> The software identifies the correct variants
    
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
	 * 获取位点基因型列表.
	 *
	 * @param informationsVec  vcfInfoList
     * @param sampleIdx        sample基因型的索引,默认值0代表最后一列
     * 
     * 
     * @return gtVec           vector <int>
	**/
    vector<int> get_gt(
        const vector<string> & informationsVec, 
        int sampleIdx = 0
    );


    /**
	 * 获取位点基因型列表.
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
	 * 获取变异的长度信息
	 *
     * @param refLen            ref长度
	 * @param qryLenVec         qry长度列表
     * @param gtVec             基因型列表
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
	 * 统计变异长度信息
	 *
	 * @param length_list         长度列表
     * 
     * 
     * @return vector<uint64_t>   每个区间的长度
	**/
    vector<uint64_t> count_num(
        vector<uint32_t> length_list
    );


    /**
	 * 获取单倍型对应的长度信息.
	 *
     * @param svType                         变异类型
	 * @param refSeq                         ref列信息
     * @param qrySeqs                        qry列信息
     * @param gtVec                          位点的分型信息
     * @param lenType                        只取单倍型对应的长度还是所有的长度(hap/all)
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
	 * 保存recall中failCall的结果
	 *
	 * @param chrStartLenInfoGtTupMap     真集中用剩下的vcf信息
     * @param outFileName                 输出文件名
     * 
     * 
     * @return int             0
	**/
    int saveFailCall(
        map<string, map<uint32_t, tuple<uint32_t, vector<uint32_t>, string, int32_t, vector<int> > > > & chrStartLenInfoGtTupMap, 
        const string & outFileName
    );


    /**
	 * 保存recall中failCall的结果
	 *
     * @param allLengthList     所有变异长度
     * 
     * 
     * @return int             0
	**/
    void roc_calculate(
        const vector<uint32_t> & allLengthList
    );
};

#endif