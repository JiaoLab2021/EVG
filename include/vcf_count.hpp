#ifndef VCF_COUNT_HPP
#define VCF_COUNT_HPP

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
#include <getopt.h>
#include <cstdlib>


using namespace std;


void help_count(char** argv);
int main_count(int argc, char** argv);


namespace VCFCOUNT
{
    /**
     * @brief 统计vcf文件
     * 
     * @param vcfFileName  输入vcf文件
     * 
     * @return 0
    **/
    int count(
        const string & vcfFileName
    );

} // namespace VCFCOUNT

#endif