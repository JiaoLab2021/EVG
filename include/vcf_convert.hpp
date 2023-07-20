#ifndef VCF_CONVERT_HPP
#define VCF_CONVERT_HPP

#include <fstream>
#include <string>
#include <iostream>
#include <regex>
#include <map>
#include <algorithm>
#include "zlib.h"
#include <getopt.h>

#include "strip_split_join.hpp"
#include "kseq.h"
#include "get_time.hpp"
#include "check_vcf_sort.hpp"
#include "save.hpp"
#include "vcf_open.hpp"

using namespace std;

void help_convert(char** argv);

int main_convert(int argc, char** argv);


struct refIndexStruct
{
    map<string, string> sequenceMap;  // map<chromosome, sequence>
    string chrLenTxt;  // Store the contig length and save it to the vcf header

    refIndexStruct() : chrLenTxt("") {}
};


class Convert
{
private:
    string refFileName_;  // reference genome
    string vcfFileName_;  // input VCF file name

    int readLen_;  // read length

    double MAF_;  // Minimal allele frequency
    double MISSRATE_;  // Missing rate

    refIndexStruct refIndexS_;  // the index of reference genome

    string outFileName_;  // output file name

    // open file
    KSEQ_INIT(gzFile, gzread)

public:
    /**
     * @brief Convert vcf to the format required by graph genome tools
     * 
     * @param refFileName     reference genome
     * @param vcfFileName     input VCF file name
     * @param readLen         read length
     * @param MAF             Minimal allele frequency
     * @param MISSRATE        Missing rate
     * @param outFileName     output file name
     * 
    **/
    Convert(string refFileName, string vcfFileName, int readLen, double MAF, double MISSRATE, string outFileName);

    /**
     * @brief build the reference genome index
     * 
     * @return void
    **/
    void build_reference_index();

    /*
     a. 表头加染色体长度
     b. 第一个变异要大于read length
     c. 序列必须与reference对应
     d. 对vcf的第八列 'END=' 进行替换
     e. refSeq的第一个碱基要和qrySeq的第一个碱基一样
     f. 检查ref和qry序列是否一样，一样的位点跳过
     g. 检查refSeq和qrySeq中是否含有atgcnATGCN外的字符，含有的话跳过该位点
     h. 检查替换后的qry有没有相同的，有的话跳过该位点 --> e.
     i. 检查是否有位置重复的变异
     j. 将基因型中的.转为.|.
     k. 将基因型中的/转为|
     l. 只保留基因型中的二倍体变异
     m. 检查GT是不是比qry的序列还多
    */
    /**
     * @brief Convert vcf to the format required by graph genome tools
     * 
     * @return void
    **/
    void vcf_convert();
};

#endif