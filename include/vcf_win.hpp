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

#include "strip_split_join.hpp"
#include "sort_find.hpp"
#include "get_time.hpp"

using namespace std;


void help_win(char** argv);
int main_win(int argc, char** argv);


namespace VCFWIN
{
    // 统计染色体条数和名字
    vector<string> count_chrName(string chrLenFilename);
    

    // 按染色体拆分vcf文件
    vector<string> split_file(string vcfFilename, vector<string> chrNameVector);
    

    // 步长字典
    pair<map<string,vector<long int>>, map<string,vector<long int>>> step_count(string chrLenFilename, int windowSize, int stepSize);
    

    // ways的分类
    map<string,string> ways_group(string waysFilename);


    // 计算窗口内变异的长度
    int window_len_count(string vcfFilename, map<string,vector<long int>> & winStartMap, map<string,vector<long int>> & winEndMap, int windowSize, map<string,string> & waysGroupMap, string mode);
    

} // namespace VCFWIN

#endif