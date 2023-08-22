#ifndef VCF_WIN_HPP
#define VCF_WIN_HPP

#include <fstream>
#include <string>
#include <iostream>
#include <numeric>
#include <getopt.h>
#include "zlib.h"
#include <thread>
#include <cstring>
#include <map>
#include <cmath>

#include "vcf_open.hpp"
#include "save.hpp"
#include "strip_split_join.hpp"
#include "sort_find.hpp"
#include "get_time.hpp"

using namespace std;


void help_win(char** argv);
int main_win(int argc, char** argv);


class VCFWin
{
private:
    uint32_t windowSize_, stepSize_; 

    string chrLenFilename_;
    string vcfFilename_;
    string waysFilename_;

    string outputFileName_;
    
    string mode_;

    map<string, uint32_t> chrLenMap_;  // map<chr, length>
    vector<string> vcfFileVec_;

     // 步长字典
    map<string, vector<uint32_t> > winStartMap_;  // map<string, vector<start> >
    map<string, vector<uint32_t> > winEndMap_;  // map<string, vector<end> >

    // ways分类字典
    map<string,string> waysGroupMap_;

    map<string,map<uint32_t, vector<uint32_t> > > chrWinStartLenMap_;  // <chromosome, map<winStart, vector<length> > >
    map<string,map<uint32_t, vector<string> > > chrWinStartGroupMap_;  // <chromosome, map<winStart, vector<group> > >

    /**
     * 获取位点基因型列表.
     *
     * @param lineVec          lineVec
     * @param sampleIdx        sample基因型的索引,默认值0代表最后一列
     * 
     * 
     * @return gtVec           vector <int>
    **/
    vector<int> get_gt(
        const vector<string>& lineVec, 
        int sampleIdx = 0
    );
public:
    /**
     * init
     *
     * @param windowSize
     * @param stepSize
     * @param chrLenFilename
     * @param vcfFilename
     * @param waysFilename
     * @param outputFileName
     * @param mode
     * 
    **/
    VCFWin(
        const uint32_t& windowSize, 
        const uint32_t& stepSize, 
        const string& chrLenFilename, 
        const string& vcfFilename, 
        const string& waysFilename, 
        const string& outputFileName, 
        const string& mode
    );

    // 统计染色体条数和名字
    void count_chrName();
    

    // 步长字典
    void step_count();
    

    // ways的分类
    void ways_group();


    // 计算窗口内变异的长度
    void window_len_count();


    // save
    void save_result();
};

#endif