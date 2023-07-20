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
     a. ��ͷ��Ⱦɫ�峤��
     b. ��һ������Ҫ����read length
     c. ���б�����reference��Ӧ
     d. ��vcf�ĵڰ��� 'END=' �����滻
     e. refSeq�ĵ�һ�����Ҫ��qrySeq�ĵ�һ�����һ��
     f. ���ref��qry�����Ƿ�һ����һ����λ������
     g. ���refSeq��qrySeq���Ƿ���atgcnATGCN����ַ������еĻ�������λ��
     h. ����滻���qry��û����ͬ�ģ��еĻ�������λ�� --> e.
     i. ����Ƿ���λ���ظ��ı���
     j. ���������е�.תΪ.|.
     k. ���������е�/תΪ|
     l. ֻ�����������еĶ��������
     m. ���GT�ǲ��Ǳ�qry�����л���
    */
    /**
     * @brief Convert vcf to the format required by graph genome tools
     * 
     * @return void
    **/
    void vcf_convert();
};

#endif