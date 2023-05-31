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
    // �洢vcf��Ϣ
    struct svInfo
    {
        string information;
        uint32_t refEnd;
    };


    struct snpIndelInfo
    {
        // ����snp��indel����ʼ����ֹλ��
        vector<uint32_t> refStartVec;
        vector<uint32_t> refEndVec;
    };
    

    // ��GT���в��
    vector<string> gt_split(const string & gtTxt);
    

    // �����ͷ���
    int vcf_split_type(const string & vcfFileName, const string & prefix, const string & filterBool);
    

    // ������snp��indel��������
    int vcf_split_number(
        const string & vcfFileName, 
        string & baseVcfFileName, 
        const int & length, 
        const string & prefix, 
        const string & filterBool
    );
} // namespace VCFSPLIT

#endif