#ifndef VCF_MERGE_HPP
#define VCF_MERGE_HPP

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
#include "strip_split_join.hpp"
#include "sort_find.hpp"
#include "get_time.hpp"
#include "check_vcf_sort.hpp"


using namespace std;

void help_merge(char** argv);

int main_merge(int argc, char** argv);

namespace MERGE
{ 
    struct baseVcfStruct
    {
        string headInfo; // vcf��ע����

        map<string, vector<int>> startMap;  // map<chromosome, vector<start>>
        map<string, vector<int>> refLenMap;  // map<chromosome, vector<refLen>>
        map<string, vector<vector<int>>> qryLenVecMap;  // map<chromosome, vector<vector<qryLen> > >, һ��λ���ж����λ�����ʱ�������������qry�ĳ���

        // �洢baseVcf��vcfInfo
        map<string, map<int, string>> BaseInfoMap;  // map<chromosome, map<start, vcfInfo>>

        // �洢���ֲ��ҷ���Ľ��������������ֺͶ�Ӧ�ķ��ͽ��
        map<string, map<int, map<string, tuple<vector<float>, vector<int> > > > > recallSoftwareGtDepVecMap;  // map<chr, map<start, map<software, tuple<vector<depth>, vector<gt> > > > >
        
        map<string, tuple<float, float, float>> depthMap;  // map<software, pair<aveDepth, variance, sd>>  �������Ϊ����ͱ�׼��
        map<string, int> softwateSampleIdx;  // map<software, index>, ����ļ���Ӧ��sample������

        // ��¼����������Paragraph��ȡʱ���Ӹ���֮����Ԫ��
        int colNum;
    };

    struct softwareVcfStruct
    {
        string software;  // �������

        map<string, vector<int>> startMap;  // map<chromosome, vector<start>>
        map<string, vector<int>> refLenMap;  // map<chromosome, vector<refLen>>
        map<string, vector<vector<int> > > qryLenVecMap;  // map<chromosome, vector<vector<qryLen> > >��������1/2����������Ҫ�洢ÿ����λ����ĳ��ȣ��������Ϊ1����Ϊ���ϱ���
        map<string, vector<vector<int>>> gtVecMap;  // map<chromosome, vector<vector<genotype>>>
        map<string, vector<vector<float>>> depthVecMap;  // map<chromosome, vector<vector<depth>>>

        
        vector<float> depthVec;  // vector<depth>
    };

    /********************************************************************************************/

   	/**
	 * Building index for different software.
	 *
	 * @param trueVcf              λ�����е�gt��vector
	 * @param ParagraphVcf         Paragraph������
	 * @param GraphTyper2Vcf       GraphTyper2������
	 * @param BayesTyperVcf        BayesTyper������
     * @param MAPVcf               VG-MAP������
     * @param GiraffeVcf           VG-Giraffe������
     * @param GraphAlignerVcf      GraphAligner������
     * @param PanGenieVcf          PanGenie������
     * @param mergeVcfStruct       �ϲ����struct
     * 
     * 
     * @return softwareNum
	**/

    int run_index_merge(
        string trueVcf, 
        string ParagraphVcf, 
        string GraphTyper2Vcf, 
        string BayesTyperVcf, 
        string MAPVcf, 
        string GiraffeVcf, 
        string GraphAlignerVcf, 
        string & PanGenieVcf, 
        baseVcfStruct & mergeVcfStruct, 
        string sample_name
    );


    /*
        ����ļ��Ƿ����
        fileName -> ��Ҫ�����ļ�
    */
    bool check_file(string fileName);



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
        const int & sampleIdx
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
        int sampleIdx
    );


    /**
	 * ����base�ļ�����������vcf���й��ˣ����ڴ�����׼index���ļ�Ϊ��ʵ��vcf�ļ�
	 *
	 * @param vcfFileName     base��vcf�ļ�
	 * @param sample_name     ������
     * 
     * 
     * @return outStruct
	**/
    baseVcfStruct build_basefile_index(
        const string & vcfFileName, 
        const string & sample_name
    );
    

    /**
	 * ������Ľ�����й��˲������������ļ�Ϊ���������ļ���
	 *
     * @date 2023/3/6
	 * @param mergeVcfStruct     build_basefile_index������base����
	 * @param vcfFileName        ��������vcf�ļ�
     * @param software           �����
     * @param sample_name        ��������������ȡ��������Ϣ
     * 
     * 
     * @return 0
	**/
    int build_softwarefile_index(
        baseVcfStruct & mergeVcfStruct, 
        const string & vcfFileName, 
        const string & software, 
        const string & sample_name
    );


    /**
	 * vcf_merge��recall����ӵ���ͼ�к���
	 *
     * @date 2023/3/7
	 * @param mergeVcfStruct     build_basefile_index������base����
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
        baseVcfStruct & mergeVcfStruct, 
        const string chromosome, 
        const int trueRefStart, 
        const string software, 
        const vector<int> qryLenVec, 
        const vector<int> gtVec, 
        const int gt, 
        const float depth, 
        int j, 
        int k,
        int l,
        int & indexLeft, 
        int & indexRight
    );


    /**
	 * ������Ľ�����й��˲������������ļ�Ϊ���������ļ���
	 *
     * @date 2023/3/6
	 * @param mergeVcfStruct     build_basefile_index������base����
	 * @param softvcfStructure   build_softwarefile_index���������vcf����
     * 
     * 
     * @return 0
	**/
    int vcf_merge(
        baseVcfStruct & mergeVcfStruct, 
        softwareVcfStruct & softvcfStructure
    );


    /**
	 * ��ȡ����ĳ�����Ϣ
	 *
     * @param vcfInfo           vcfInfo
     * @param gtVec             �������б�
     * 
     * 
     * @return int              svLength
	**/
    int sv_length_select(
        const string & vcfInfo, 
        const vector<int> & gtVec
    );


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
     * @return tuple<int, vector<int> >      tuple<refLen, vector<qryLen> >
	**/
    tuple<int, vector<int> > get_hap_len(
        const string & svType, 
        const string & refSeq, 
        const string & qrySeqs, 
        const vector<int> & gtVec, 
        const string & lenType
    );


	/**
	 * Filter result according to the depth and number of variants.
	 *
	 * @param gtAllVec        λ�����е�gt��vector
	 * @param softwareVec     λ�����е�software��vector
	 * @param depthNorVec     λ�����б�׼�����depth��vector
     * @param softwareNum     �û�������������
     * @param vcfInfo         λ��ı�����Ϣ
     * @param mode            EVG���е�ģʽ
     * @param mergeVcfStruct  ��ͼ
     * 
     * @return software
	**/
    string filter_rule(
        vector<vector<int>> & gtAllVec, 
        vector<string> & softwareVec, 
        vector<float> & depthNorVec, 
        const string & vcfInfo, 
        const int & softwareNum, 
        const string & mode, 
        baseVcfStruct & mergeVcfStruct
    );


    /*
        vcf����
        mergeVcfStruct -> ����ϲ����struct
    */

   /**
	 * Filter and merge result according to the depth and number of variants.
	 *
	 * @param mergeVcfStruct     ����ϲ����struct
	 * @param softwareNum        �û��ύ���������
	 * @param mode               EVG���е�ģʽ
     * 
     * @return outMap            outMap<chromosome, map<refStart, vcfInfo> >
	**/
   map<string, map<int, string > > vcf_merge_filter(
        baseVcfStruct & mergeVcfStruct, 
        const int & softwareNum, 
        const string & mode
    );
    

    /**
	 * ���㷽��ͱ�׼��
	 *
	 * @param data     ����λ����ȵ��б�
     * 
     * @return tuple   make_tuple(mean, variance, std_deviation)
	**/
    tuple<float, float, float> cal_var_sd(const vector<float>& data);


    /*
        ������
        mergeVcfStruct -> ����ϲ����struct
        outMap -> vcf_merge_filter������   outMap[chromosome][refStart][refLen][qryLen]
        prefix -> ����ļ���ǰ׺
    */

    /**
	 * ������
	 *
     * @param headInfo           vcfע����
	 * @param outMap             vcf_merge_filter������   outMap[chromosome][refStart]
	 * @param outputFileName     ����ļ���
     * 
     * @return tuple   make_tuple(mean, variance, std_deviation)
	**/
    int result_save(
        string & headInfo, 
        map<string, map<int, string > > & outMap, 
        const string & outputFileName
    );

}  // namespace MERGE

#endif