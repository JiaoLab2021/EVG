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
#include "strip_split_join.hpp"
#include "sort_find.hpp"
#include "get_time.hpp"
#include "check_vcf_sort.hpp"


using namespace std;


// 全局变量
// rocKey  
extern string rocKey;  // (FORMAT.<name>/INFO.<name>) [FORMAT.DP]
// rocType判断是用recall来计算roc，还是只用genotype计算roc  recall/genotype [genotype]
extern string rocType;


void help_recall(char** argv);
int main_recall(int argc, char** argv);


namespace RECALL
{
    struct vcfStructure
    {
        map<string,vector<int> > refStartVecMap;  // map<chr,vector<refStart> >
        
        map<string, map<int, tuple<int, vector<int>, string, int, vector<int> > > > chrStartLenInfoGtVecMap;  // unordered_map<chr, map<refStart, tuple<refLen, qryLenVec, vcfInfo, svLen, gtVec> > >
        
        vector<int> allLengthList;  // 所有变异的长度
    };


    // 对变异长度进行统计的Bool函数
    class f_mod
    {
    private:
        int dv1;
        int dv2;

    public:
        f_mod(int d1 = 1, int d2 = 2) : dv1(d1), dv2(d2) {}
    
        bool operator() (int x) {return  dv1 <= x && x <= dv2;}
    };


    /**
	 * 对vcf文件进行排序并构建索引
	 *
	 * @param genotype_filename    输入文件
     * 
     * 
     * @return outStructure        vcfStructure
	**/
    vcfStructure build_index(
        string genotype_filename
    );


    /**
	 * 对gramtools的结果进行转换
	 *
	 * @param evaluateFilename    输入文件
     * 
     * 
     * @return outFileName        输出文件
	**/
    string gramtools_convert(
        string evaluateFilename
    );


    /**
	 * genotype评估函数
	 *
	 * @param evaluateFilename            输入文件（待评估）
     * @param chrStartLenInfoGtVecMap     真集所有的vcf信息
     * @param all_length_list             真集所有vcf长度
     * 
     * 
     * @return int             0
	**/
    int evulate_gt(
        string evaluateFilename, 
        map<string, map<int, tuple<int, vector<int>, string, int, vector<int> > > > & chrStartLenInfoGtVecMap, 
        vector<int> & all_length_list
    );


    /**
	 * recall评估函数
	 *
	 * @param evaluateFilename            输入文件（待评估）
     * @param trueVcfStructure            build_index()输出结果
     * 
     * 
     * @return int             0
	**/
    int evulate_recall(
        const string & evaluateFilename, 
        vcfStructure & trueVcfStructure
    );


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
	 * 获取变异的长度信息
	 *
     * @param refLen            ref长度
	 * @param qryLenVec         qry长度列表
     * @param gtVec             基因型列表
     * 
     * 
     * @return int              svLength
	**/
    int sv_length_select(
        const int & refLen, 
        const vector<int> & qryLenVec, 
        const vector<int> & gtVec
    );


    /**
	 * 统计变异长度信息
	 *
     * @param sv_length         划分的区间
	 * @param length_list       长度列表
     * 
     * 
     * @return vector<int>      每个区间的长度
	**/
    vector<int> count_num(
        vector<string> sv_length, 
        vector<int> length_list
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
    tuple<int, vector<int> > get_hap_len(
        const string & svType, 
        const string & refSeq, 
        const string & qrySeqs, 
        const vector<int> & gtVec, 
        const string & lenType
    );
    

    /**
	 * 保存recall中failCall的结果
	 *
	 * @param chrStartLenInfoGtVecMap     真集中用剩下的vcf信息
     * @param outFileName                 输出文件名
     * 
     * 
     * @return int             0
	**/
    int saveFailCall(
        map<string, map<int, tuple<int, vector<int>, string, int, vector<int> > > > & chrStartLenInfoGtVecMap, 
        const string & outFileName
    );


    /**
	 * 保存recall中failCall的结果
	 *
	 * @param rocCallMap
     * @param rocRecallMap
     * @param all_length_list     所有变异长度
     * 
     * 
     * @return int             0
	**/
    void roc_calculate(
        const map<float, vector<int> > & rocCallMap, 
        const map<float, vector<int> > & rocRecallMap, 
        const vector<int> & all_length_list
    );
    
} // namespace RECALL

#endif