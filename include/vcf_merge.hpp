#ifndef VCF_MERGE_HPP
#define VCF_MERGE_HPP

#include <getopt.h>
#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <vector>
#include <numeric>
#include <map>
#include <malloc.h>
#include "zlib.h"
#include <cmath>

#include "vcf_open.hpp"
#include "save.hpp"
#include "strip_split_join.hpp"
#include "sort_find.hpp"
#include "get_time.hpp"
#include "check_vcf_sort.hpp"


using namespace std;

void help_merge(char** argv);

int main_merge(int argc, char** argv);


 struct baseVcfStruct
{
    string headInfo; // vcf��ע����

    map<string, vector<uint32_t> > startMap;  // map<chromosome, vector<start>>
    map<string, vector<uint32_t> > refLenMap;  // map<chromosome, vector<refLen>>
    map<string, vector<vector<uint32_t> > > qryLenVecMap;  // map<chromosome, vector<vector<qryLen> > >, һ��λ���ж����λ�����ʱ�������������qry�ĳ���

    // �洢baseVcf��vcfInfo
    map<string, map<uint32_t, string> > BaseInfoMap;  // map<chromosome, map<start, vcfInfo>>

    // �洢���ֲ��ҷ���Ľ��������������ֺͶ�Ӧ�ķ��ͽ��
    map<string, map<uint32_t, map<string, tuple<vector<float>, vector<int> > > > > recallSoftwareGtDepVecMap;  // map<chr, map<start, map<software, tuple<vector<depth>, vector<gt> > > > >
    
    map<string, tuple<float, float, float>> depthMap;  // map<software, pair<aveDepth, variance, sd>>  �������Ϊ����ͱ�׼��
    map<string, int> softwateSampleIdx;  // map<software, index>, ����ļ���Ӧ��sample������

    // ��¼����������Paragraph��ȡʱ���Ӹ���֮����Ԫ��
    int colNum;

    baseVcfStruct() : colNum(0) {}
};

struct softwareVcfStruct
{
    string software;  // �������

    map<string, vector<uint32_t> > startMap;  // map<chromosome, vector<start>>
    map<string, vector<uint32_t> > refLenMap;  // map<chromosome, vector<refLen>>
    map<string, vector<vector<uint32_t> > > qryLenVecMap;  // map<chromosome, vector<vector<qryLen> > >��������1/2����������Ҫ�洢ÿ����λ����ĳ��ȣ��������Ϊ1����Ϊ���ϱ���
    map<string, vector<vector<int> > > gtVecMap;  // map<chromosome, vector<vector<genotype>>>
    map<string, vector<vector<float> > > depthVecMap;  // map<chromosome, vector<vector<depth>>>

    vector<float> depthVec;  // vector<depth>
};



class VCFMerge
{
private:
    // EVG������ģʽ
    string mode_;

    string trueVcf_;

    string ParagraphVcf_;
    string GraphTyper2Vcf_;
    string BayesTyperVcf_;
    string MAPVcf_;
    string GiraffeVcf_;
    string GraphAlignerVcf_;
    string PanGenieVcf_;

    string sampleName_;

    uint16_t softwareNum_;  // ��¼��������������ڶԽ������

    map<string, map<int, string > > outChrStartInfoMap_;  // outMap<chromosome, map<refStart, vcfInfo> >

    baseVcfStruct mergeVcfStruct_;  // recore the index of all software

    string outputFileName_;


    /**
     * ��ȡ�����Ͷ�Ӧ�ĳ�����Ϣ.
     *
     * @param svType                         ��������
     * @param refSeq                         ref����Ϣ
     * @param qrySeqs                        qry����Ϣ
     * @param gtVec                          λ��ķ�����Ϣ
     * @param lenType                        ֻȡ�����Ͷ�Ӧ�ĳ��Ȼ������еĳ���(hap/all)
     * 
     * 
     * @return tuple<uint32_t, vector<uint32_t> >      tuple<refLen, vector<qryLen> >
    **/
    tuple<uint32_t, vector<uint32_t> > get_hap_len(
        const string & svType, 
        const string & refSeq, 
        const string & qrySeqs, 
        const vector<int> & gtVec, 
        const string & lenType
    );


    /**
	 * ������Ľ�����й��˲������������ļ�Ϊ���������ļ���
	 *
     * @date 2023/07/09
     * 
	 * @param vcfFileName        ��������vcf�ļ�
     * @param software           �����
     * 
     * 
     * @return void
	**/
    void build_softwarefile_index(
        const string & vcfFileName, 
        const string & software
    );


    // paragraph DP -> 1
    // graphaligner, vg-map- vg-giraffe DP -> 6
    // bayestyper MAC -> -1,1.87332
    // graphtyper DP -> 6
    // pangenie KC -> 4
    /**
	 * ��ȡλ��DP.
	 *
	 * @param informationsVec     vcfInfoVec
     * @param sampleIdx           sample��Ӧ����
     * 
     * 
     * @return depthVec        vector<depth>
	**/
    vector<float> get_depth(
        const vector<string> & informationsVec, 
        const uint32_t& sampleIdx
    );


    /**
	 * ��ȡλ��������б�.
	 *
	 * @param informationsVec  vcfInfoList
     * @param sampleIdx        sample�����͵�����,Ĭ��ֵ0�������һ��
     * 
     * 
     * @return gtVec           vector <int>
	**/
    vector<int> get_gt(
        const vector<string> & informationsVec, 
        uint32_t sampleIdx
    );



    /**
	 * vcf_merge��recall����ӵ���ͼ�к���
	 *
     * @date 2023/07/09
     * 
	 * @param chromosome         Ⱦɫ���
     * @param trueRefStart       ��ʵ��ref��ʼλ��
     * @param software           ���
     * @param qryLenVec          ������͵ĵ����Ͷ�Ӧ����
     * @param gtVec              ����ķ���gt�б�
     * @param gt                 ��ǰѭ����gt
     * @param depth              λ���Ӧ�����
     * @param j                  ������������ѭ��������
     * @param k                  ���ֲ��ҷ���������ѭ��������
     * @param l                  ��ʵλ�㵥���͵�ѭ������
     * @param indexLeft          ���ֲ��ҷ���������
     * @param indexRight         ���ֲ��ҷ���������
     * 
     * 
     * @return tuple<int, int>   tuple<state, j>  0 -> ��ȷ��ӣ�-1 -> ��Ҫ������һ��ѭ��
	**/
    tuple<int, int> recall_push(
        const string chromosome, 
        const uint32_t trueRefStart, 
        const string software, 
        const vector<uint32_t> qryLenVec, 
        const vector<int> gtVec, 
        const int gt, 
        const float depth, 
        uint32_t j, 
        uint32_t k,
        uint32_t l,
        int64_t & indexLeft, 
        int64_t & indexRight
    );


    /**
	 * ������Ľ�����й��˲������������ļ�Ϊ���������ļ���
	 *
     * @date 2023/07/09
     * 
	 * @param softvcfStructure   build_softwarefile_index���������vcf����
     * 
     * 
     * @return void
	**/
    void vcf_merge(
        softwareVcfStruct & softvcfStructure
    );


    /**
	 * ��ȡ����ĳ�����Ϣ
	 *
     * @param vcfInfo           vcfInfo
     * @param gtVec             �������б�
     * 
     * 
     * @return int64_t          svLength
	**/
    int64_t sv_length_select(
        const string & vcfInfo, 
        const vector<int> & gtVec
    );


	/**
	 * Filter result according to the depth and number of variants.
	 *
	 * @param gtAllVec        λ�����е�gt��vector
	 * @param softwareVec     λ�����е�software��vector
	 * @param depthNorVec     λ�����б�׼�����depth��vector
     * @param vcfInfo         λ��ı�����Ϣ
     * 
     * @return software
	**/
    string filter_rule(
        vector<vector<int>> & gtAllVec, 
        vector<string> & softwareVec, 
        vector<float> & depthNorVec, 
        const string & vcfInfo
    );


    /**
	 * ���㷽��ͱ�׼��
	 *
	 * @param data     ����λ����ȵ��б�
     * 
     * @return tuple   make_tuple(mean, variance, std_deviation)
	**/
    tuple<float, float, float> cal_var_sd(const vector<float>& data);
public:
    /**
	 * init
	 *
     * @param mode
	 * @param trueVcf
     * @param ParagraphVcf
     * @param GraphTyper2Vcf
     * @param BayesTyperVcf
     * @param MAPVcf
     * @param GiraffeVcf
     * @param GraphAlignerVcf
     * @param PanGenieVcf
     * @param sampleName       Sample names to be merged
     * @param outputFileName
     * 
	**/
    VCFMerge(
        const string& mode,
        const string& trueVcf, 
        const string& ParagraphVcf, 
        const string& GraphTyper2Vcf, 
        const string& BayesTyperVcf, 
        const string& MAPVcf, 
        const string& GiraffeVcf, 
        const string& GraphAlignerVcf, 
        const string& PanGenieVcf, 
        const string& sampleName, 
        const string& outputFileName
    );


    /**
	 * ����base�ļ�����������vcf���й��ˣ����ڴ�����׼index���ļ�Ϊ��ʵ��vcf�ļ�
	 *
     * 
     * @return outStruct
	**/
    void build_basefile_index();


    /**
	 * Building index for different software.
     * 
     * 
     * @return softwareNum_
	**/
    void run_index_merge();


    /*
        vcf����
        mergeVcfStruct -> ����ϲ����struct
    */

    /**
	 * Filter and merge result according to the depth and number of variants.
     * 
     * @return void
	**/
   void vcf_merge_filter();


    /*
        save result
        mergeVcfStruct -> ����ϲ����struct
        outMap -> vcf_merge_filter������   outMap[chromosome][refStart][refLen][qryLen]
        prefix -> ����ļ���ǰ׺
    */

    /**
	 * save result
	 *
     * 
     * @return void
	**/
    void result_save();

};

#endif