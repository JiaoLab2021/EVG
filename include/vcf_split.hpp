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
#include "vcf_open.hpp"
#include "save.hpp"
#include "get_time.hpp"
#include "sort_find.hpp"


using namespace std;

void help_split(char** argv);

int main_split(int argc, char** argv);


class VCFSplit
{
private:
    // 存储vcf信息
    struct svInfo
    {
        string information;
        uint32_t refEnd;

        svInfo() : information(""), refEnd(0) {}
    };

    struct snpIndelInfo
    {
        // 保存snp和indel的起始和终止位置
        vector<uint32_t> refStartVec;
        vector<uint32_t> refEndVec;
    };

    string vcfFileName_;  // input VCF file name
    string baseVcfFileName_;  // VCF file for extracting snp and indel locations
    int length_;  // Count the number of SNPs and indels within the 'length' near the breakpoint
    bool filterBool_;  // Whether to filter according to the 'FILTER' column
    string prefix_;  // output file prefix

public:
    /**
     * @brief Convert vcf to the format required by graph genome tools
     * 
     * @param vcfFileName      input VCF file name
     * @param baseVcfFileName  VCF file for extracting snp and indel locations
     * @param length           Count the number of SNPs and indels within the 'length' near the breakpoint
     * @param filterBool       Whether to filter according to the 'FILTER' column
     * @param prefix           output file prefix
     * 
    **/
    VCFSplit(
        string vcfFileName,
        string baseVcfFileName,
        int length, 
        bool filterBool,
        string prefix
    );


    /**
     * @brief Get the length information of variant
     * 
     * @param INFOSTRUCTTMP   variant information
     * @param VCFOPENCLASS    VCFOPEN
     * @param refLen          store the length of REF
     * @param qryLen          store the length of ALT
     * 
    **/
    bool get_len(
        const VCFINFOSTRUCT& INFOSTRUCTTMP, 
        VCFOPEN& VCFOPENCLASS, 
        uint32_t& refLen, 
        uint32_t& qryLen
    );

    // by type
    void vcf_split_type();
    

    // Classified by the number of nearby SNPs and indels
    void vcf_split_number();
};

#endif