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
    string headInfo; // vcf的注释行

    map<string, vector<uint32_t> > startMap;  // map<chromosome, vector<start>>
    map<string, vector<uint32_t> > refLenMap;  // map<chromosome, vector<refLen>>
    map<string, vector<vector<uint32_t> > > qryLenVecMap;  // map<chromosome, vector<vector<qryLen> > >, 一个位点有多个等位基因的时候，用这个来保存qry的长度

    // 存储baseVcf的vcfInfo
    map<string, map<uint32_t, string> > BaseInfoMap;  // map<chromosome, map<start, vcfInfo>>

    // 存储二分查找法后的结果，保存软件名字和对应的分型结果
    map<string, map<uint32_t, map<string, tuple<vector<float>, vector<int> > > > > recallSoftwareGtDepVecMap;  // map<chr, map<start, map<software, tuple<vector<depth>, vector<gt> > > > >
    
    map<string, tuple<float, float, float>> depthMap;  // map<software, pair<aveDepth, variance, sd>>  最后两个为方差和标准差
    map<string, int> softwateSampleIdx;  // map<software, index>, 软件文件对应的sample的索引

    // 记录列数，用于Paragraph提取时，从该列之后找元素
    int colNum;

    baseVcfStruct() : colNum(0) {}
};

struct softwareVcfStruct
{
    string software;  // 软件名字

    map<string, vector<uint32_t> > startMap;  // map<chromosome, vector<start>>
    map<string, vector<uint32_t> > refLenMap;  // map<chromosome, vector<refLen>>
    map<string, vector<vector<uint32_t> > > qryLenVecMap;  // map<chromosome, vector<vector<qryLen> > >，可能是1/2等情况，因此要存储每个等位基因的长度，如果长度为1代表为纯合变异
    map<string, vector<vector<int> > > gtVecMap;  // map<chromosome, vector<vector<genotype>>>
    map<string, vector<vector<float> > > depthVecMap;  // map<chromosome, vector<vector<depth>>>

    vector<float> depthVec;  // vector<depth>
};



class VCFMerge
{
private:
    // EVG的运行模式
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

    uint16_t softwareNum_;  // 记录软件的数量，用于对结果过滤

    map<string, map<int, string > > outChrStartInfoMap_;  // outMap<chromosome, map<refStart, vcfInfo> >

    baseVcfStruct mergeVcfStruct_;  // recore the index of all software

    string outputFileName_;


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
	 * 对软件的结果进行过滤并构建索引，文件为软件的输出文件。
	 *
     * @date 2023/07/09
     * 
	 * @param vcfFileName        软件输出的vcf文件
     * @param software           软件名
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
	 * 获取位点DP.
	 *
	 * @param informationsVec     vcfInfoVec
     * @param sampleIdx           sample对应的列
     * 
     * 
     * @return depthVec        vector<depth>
	**/
    vector<float> get_depth(
        const vector<string> & informationsVec, 
        const uint32_t& sampleIdx
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
        uint32_t sampleIdx
    );



    /**
	 * vcf_merge中recall后添加到总图中函数
	 *
     * @date 2023/07/09
     * 
	 * @param chromosome         染色体号
     * @param trueRefStart       真实的ref起始位置
     * @param software           软件
     * @param qryLenVec          软件分型的单倍型对应长度
     * @param gtVec              软件的分型gt列表
     * @param gt                 当前循环的gt
     * @param depth              位点对应的深度
     * @param j                  软件多个单倍型循环的索引
     * @param k                  二分查找法两个坐标循环的索引
     * @param l                  真实位点单倍型的循环索引
     * @param indexLeft          二分查找法的左索引
     * @param indexRight         二分查找法的右索引
     * 
     * 
     * @return tuple<int, int>   tuple<state, j>  0 -> 正确添加；-1 -> 需要进行下一个循环
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
	 * 对软件的结果进行过滤并构建索引，文件为软件的输出文件。
	 *
     * @date 2023/07/09
     * 
	 * @param softvcfStructure   build_softwarefile_index构建的软件vcf索引
     * 
     * 
     * @return void
	**/
    void vcf_merge(
        softwareVcfStruct & softvcfStructure
    );


    /**
	 * 获取变异的长度信息
	 *
     * @param vcfInfo           vcfInfo
     * @param gtVec             基因型列表
     * 
     * 
     * @return int64_t          svLength
	**/
    int64_t sv_length_select(
        const string & vcfInfo, 
        const vector<int> & gtVec
    );


	/**
	 * Filter result according to the depth and number of variants.
	 *
	 * @param gtAllVec        位点所有的gt的vector
	 * @param softwareVec     位点所有的software的vector
	 * @param depthNorVec     位点所有标准化后的depth的vector
     * @param vcfInfo         位点的变异信息
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
	 * 计算方差和标准差
	 *
	 * @param data     含有位点深度的列表
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
	 * 构件base文件索引，不对vcf进行过滤，用于创建基准index。文件为真实的vcf文件
	 *
     * 
     * @return outStruct
	**/
    void build_basefile_index();


    /**
	 * Building index for different software.
     * 
     * 
     * @return softwareNum_
	**/
    void run_index_merge();


    /*
        vcf过滤
        mergeVcfStruct -> 软件合并后的struct
    */

    /**
	 * Filter and merge result according to the depth and number of variants.
     * 
     * @return void
	**/
   void vcf_merge_filter();


    /*
        save result
        mergeVcfStruct -> 软件合并后的struct
        outMap -> vcf_merge_filter输出结果   outMap[chromosome][refStart][refLen][qryLen]
        prefix -> 输出文件名前缀
    */

    /**
	 * save result
	 *
     * 
     * @return void
	**/
    void result_save();

};

#endif