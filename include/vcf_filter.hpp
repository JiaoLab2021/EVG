#ifndef vcf_filter_hpp
#define vcf_filter_hpp

#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include "zlib.h"
#include <unordered_map>
#include <cmath>
#include "strip_split_join.hpp"
#include "get_time.hpp"
#include "sort_find.hpp"
#include "vcf_open.hpp"
#include "save.hpp"
#include <getopt.h>
#include <cstdlib>


using namespace std;


void help_filter(char** argv);
int main_filter(int argc, char** argv);


namespace VCFFILTER
{
    /**
     * @brief 根据次等位基因频率和缺失率过滤SNPs
     * 
     * @param vcfFileName    输入vcf文件
     * @param outputFileName 输出文件名
     * @param MAF            次等位基因频率
     * @param MISSRATE       缺失率
     * 
     * @return 0
    **/
    int vcf_filter(
        const string & vcfFileName, 
        const string & outputFileName, 
        const double & MAF, 
        const double & MISSRATE
    );

} // namespace VCFFILTER

#endif