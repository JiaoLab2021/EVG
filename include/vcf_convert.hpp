#ifndef VCF_CONVERT_HPP
#define VCF_CONVERT_HPP

#include <fstream>
#include <string>
#include <iostream>
#include <regex>
#include <map>
#include <algorithm>
#include "zlib.h"
#include "strip_split_join.hpp"
#include "kseq.h"
#include "get_time.hpp"
#include "check_vcf_sort.hpp"
#include "save.hpp"
#include "vcf_open.hpp"

using namespace std;

void help_convert(char** argv);

int main_convert(int argc, char** argv);

namespace CONVERT{
    // kseq.h 打开文件
    KSEQ_INIT(gzFile, gzread)

    struct refIndexStruct
    {
        map<string, string> sequenceMap;  // map<chromosome, sequence>
        string chrLenTxt;  // 存储contig长度，保存到vcf头
    };


    /*
        打开fasta文件
        inputFileName -> refgenome
        refIndexS -> 保存contig长度和序列信息
    */
    void build_reference_index(string inputFileName, refIndexStruct & refIndexS);


    /*
        对vcf文件进行替换
        vcfFilename -> 需要转换的vcf文件
        readLen -> 测序文件的read长度
        refIndexS -> contig长度和序列信息
        outFilename -> 输出文件名
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
     * @brief 将vcf转换为 graph genome tools 需要的格式
     * 
     * @param vcfFileName     输入vcf文件
     * @param readLen         读长
     * @param MAF             次等位基因频率
     * @param MISSRATE        缺失率
     * @param refIndexS       参考基因组长度和序列信息
     * @param outFilename     输出文件名
     * 
     * @return 0
    **/
    void vcf_convert(
        const string & vcfFilename, 
        const int & readLen, 
        const double & MAF, 
        const double & MISSRATE, 
        const refIndexStruct & refIndexS, 
        const string & outFilename
    );

}  // namespace CONVERT

#endif