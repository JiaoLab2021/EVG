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
    // kseq.h ���ļ�
    KSEQ_INIT(gzFile, gzread)

    struct refIndexStruct
    {
        map<string, string> sequenceMap;  // map<chromosome, sequence>
        string chrLenTxt;  // �洢contig���ȣ����浽vcfͷ
    };


    /*
        ��fasta�ļ�
        inputFileName -> refgenome
        refIndexS -> ����contig���Ⱥ�������Ϣ
    */
    void build_reference_index(string inputFileName, refIndexStruct & refIndexS);


    /*
        ��vcf�ļ������滻
        vcfFilename -> ��Ҫת����vcf�ļ�
        readLen -> �����ļ���read����
        refIndexS -> contig���Ⱥ�������Ϣ
        outFilename -> ����ļ���
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
     * @brief ��vcfת��Ϊ graph genome tools ��Ҫ�ĸ�ʽ
     * 
     * @param vcfFileName     ����vcf�ļ�
     * @param readLen         ����
     * @param MAF             �ε�λ����Ƶ��
     * @param MISSRATE        ȱʧ��
     * @param refIndexS       �ο������鳤�Ⱥ�������Ϣ
     * @param outFilename     ����ļ���
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