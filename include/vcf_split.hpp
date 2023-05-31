#ifndef VCF_SPLIT_HPP
#define VCF_SPLIT_HPP

#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include "zlib.h"
#include <map>
#include <unordered_map>

#include "strip_split_join.hpp"
#include "get_time.hpp"
#include "sort_find.hpp"


using namespace std;

void help_split(char** argv);

int main_split(int argc, char** argv);

namespace VCFSPLIT
{
    // 存储vcf信息
    struct svInfo
    {
        string information;
        uint32_t refEnd;
    };


    struct snpIndelInfo
    {
        // 保存snp和indel的起始和终止位置
        vector<uint32_t> refStartVec;
        vector<uint32_t> refEndVec;
    };
    

    // 对GT进行拆分
    vector<string> gt_split(const string & gtTxt);
    

    // 按类型分类
    int vcf_split_type(const string & vcfFileName, const string & prefix, const string & filterBool);
    

    // 按附近snp和indel数量分类
    int vcf_split_number(
        const string & vcfFileName, 
        string & baseVcfFileName, 
        const int & length, 
        const string & prefix, 
        const string & filterBool
    );
} // namespace VCFSPLIT

#endif