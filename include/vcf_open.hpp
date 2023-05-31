#ifndef VCF_OPEN_HPP
#define VCF_OPEN_HPP

#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <map>
#include <tuple>
#include <iomanip>
#include "zlib.h"
#include <unordered_map>
#include "strip_split_join.hpp"
#include "get_time.hpp"

using namespace std;

struct VCFINFOSTRUCT
{
    public:
    // �ܵ���Ϣ
    string INFO = "";
    vector<string> INFOVec;

    string CHROM = "";  // CHROM  [0]
    uint32_t POS = 0;  // POS  [1]
    string REF = "";  // REF  [3]
    vector<string> ALTVec;  // ALT �б�  [4]
    double QUAL = 0;  // QUAL  [5]
    string FILTER = ""; // FILTER  [6]

    string TYPE = "";  // ��������
    uint32_t LEN = 0;  // ref����
    uint32_t END = 0;  // ref����

    vector<string> GTVec; // �������б�
    double MAF = 0.0;  // �ε�λ����Ƶ��
    double MISSRATE = 0.0;  // ȱʧ����

    // ��սṹ��
    void clear();
};


/**
 * @brief ��vcf�ļ�
 * 
 * @param vcfFileName_   the output of vcf file
 * 
 * @return
**/
class VCFOPEN
{
private:
    // vcf�ļ�
    string vcfFileName;

    // �����ļ���
    gzFile gzfpI;
public:
    void init(
        const string & vcfFileName_
    );

    int open();

    bool read(
        VCFINFOSTRUCT & INFOSTRUCTTMP_
    );

    string get_TYPE(
        const uint32_t & LEN,
        const vector<string> & ALTVec
    );

    // ������line�ķ��ͽ��
    map<int, vector<string> > get_gt(
        const vector<string> & INFOVec
    );

    // ����MAF��MISSRATE
    tuple<double, double> calculate(
        const map<int, vector<string> > & GTVecMap, 
        uint32_t sampleNum
    );

    int close();
};

#endif