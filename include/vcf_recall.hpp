#ifndef VCF_RECALL_HPP
#define VCF_RECALL_HPP

#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <vector>
#include <map>
#include <tuple>
#include <unordered_map>
#include <algorithm>
#include <iomanip>
#include <cstdlib>
#include "zlib.h"
#include "strip_split_join.hpp"
#include "sort_find.hpp"
#include "get_time.hpp"
#include "check_vcf_sort.hpp"


using namespace std;


// ȫ�ֱ���
// rocKey  
extern string rocKey;  // (FORMAT.<name>/INFO.<name>) [FORMAT.DP]
// rocType�ж�����recall������roc������ֻ��genotype����roc  recall/genotype [genotype]
extern string rocType;


void help_recall(char** argv);
int main_recall(int argc, char** argv);


namespace RECALL
{
    struct vcfStructure
    {
        map<string,vector<int> > refStartVecMap;  // map<chr,vector<refStart> >
        
        map<string, map<int, tuple<int, vector<int>, string, int, vector<int> > > > chrStartLenInfoGtVecMap;  // unordered_map<chr, map<refStart, tuple<refLen, qryLenVec, vcfInfo, svLen, gtVec> > >
        
        vector<int> allLengthList;  // ���б���ĳ���
    };


    // �Ա��쳤�Ƚ���ͳ�Ƶ�Bool����
    class f_mod
    {
    private:
        int dv1;
        int dv2;

    public:
        f_mod(int d1 = 1, int d2 = 2) : dv1(d1), dv2(d2) {}
    
        bool operator() (int x) {return  dv1 <= x && x <= dv2;}
    };


    /**
	 * ��vcf�ļ��������򲢹�������
	 *
	 * @param genotype_filename    �����ļ�
     * 
     * 
     * @return outStructure        vcfStructure
	**/
    vcfStructure build_index(
        string genotype_filename
    );


    /**
	 * ��gramtools�Ľ������ת��
	 *
	 * @param evaluateFilename    �����ļ�
     * 
     * 
     * @return outFileName        ����ļ�
	**/
    string gramtools_convert(
        string evaluateFilename
    );


    /**
	 * genotype��������
	 *
	 * @param evaluateFilename            �����ļ�����������
     * @param chrStartLenInfoGtVecMap     �漯���е�vcf��Ϣ
     * @param all_length_list             �漯����vcf����
     * 
     * 
     * @return int             0
	**/
    int evulate_gt(
        string evaluateFilename, 
        map<string, map<int, tuple<int, vector<int>, string, int, vector<int> > > > & chrStartLenInfoGtVecMap, 
        vector<int> & all_length_list
    );


    /**
	 * recall��������
	 *
	 * @param evaluateFilename            �����ļ�����������
     * @param trueVcfStructure            build_index()������
     * 
     * 
     * @return int             0
	**/
    int evulate_recall(
        const string & evaluateFilename, 
        vcfStructure & trueVcfStructure
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
        int sampleIdx = 0
    );


    /**
	 * ��ȡ����ĳ�����Ϣ
	 *
     * @param refLen            ref����
	 * @param qryLenVec         qry�����б�
     * @param gtVec             �������б�
     * 
     * 
     * @return int              svLength
	**/
    int sv_length_select(
        const int & refLen, 
        const vector<int> & qryLenVec, 
        const vector<int> & gtVec
    );


    /**
	 * ͳ�Ʊ��쳤����Ϣ
	 *
     * @param sv_length         ���ֵ�����
	 * @param length_list       �����б�
     * 
     * 
     * @return vector<int>      ÿ������ĳ���
	**/
    vector<int> count_num(
        vector<string> sv_length, 
        vector<int> length_list
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
	 * ����recall��failCall�Ľ��
	 *
	 * @param chrStartLenInfoGtVecMap     �漯����ʣ�µ�vcf��Ϣ
     * @param outFileName                 ����ļ���
     * 
     * 
     * @return int             0
	**/
    int saveFailCall(
        map<string, map<int, tuple<int, vector<int>, string, int, vector<int> > > > & chrStartLenInfoGtVecMap, 
        const string & outFileName
    );


    /**
	 * ����recall��failCall�Ľ��
	 *
	 * @param rocCallMap
     * @param rocRecallMap
     * @param all_length_list     ���б��쳤��
     * 
     * 
     * @return int             0
	**/
    void roc_calculate(
        const map<float, vector<int> > & rocCallMap, 
        const map<float, vector<int> > & rocRecallMap, 
        const vector<int> & all_length_list
    );
    
} // namespace RECALL

#endif