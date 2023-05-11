#ifndef vcf_open_hpp
#define vcf_open_hpp

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

namespace VCFOPEN
{


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

        void clear();
    };

    // ��սṹ��
    void VCFINFOSTRUCT::clear(

    )
    {
        // �ܵ���Ϣ
         INFO = "";
         INFOVec.clear();
         vector<string>().swap(INFOVec);

         CHROM = "";  // CHROM  [0]
         POS = 0;  // POS  [1]
         REF = "";  // REF  [3]
         ALTVec.clear();  // ALT �б�  [4]
         vector<string>().swap(ALTVec);
         QUAL = 0;  // QUAL  [5]
         FILTER = ""; // FILTER  [6]

         TYPE = "";  // ��������
         LEN = 0;  // ref����
         END = 0;  // ref����

         GTVec.clear(); // �������б�
         vector<string>().swap(GTVec);
         MAF = 0.0;  // �ε�λ����Ƶ��
         MISSRATE = 0.0;  // ȱʧ����
    }



    
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

    /**
     * @brief ��ʼ��
     * 
     * @param vcfFileName_   the output of vcf file
     * 
     * @return
    **/
    void VCFOPEN::init(
        const string & vcfFileName_
    )
    {
        vcfFileName = vcfFileName_;
    }

    /**
     * @brief ���ļ�
     * 
     * @return 0
    **/
   int VCFOPEN::open(

   )
   {
        // �����ļ���
        gzfpI = gzopen(vcfFileName.c_str(), "rb");
        if(!gzfpI)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'" << vcfFileName << "': No such file or directory." << endl;
            exit(1);
        }

        return 0;
    }

    /**
     * @brief �����ļ�
     * 
     * @param INFOSTRUCTTMP_  �洢���е�����
     * 
     * @return bool
    **/
    bool VCFOPEN::read(
        VCFINFOSTRUCT & INFOSTRUCTTMP_
    )
    {
        INFOSTRUCTTMP_.clear();

        string info = ""; // ��ʱ�ַ���
        char line[1024]; // һ��ֻ��1024�ֽڵ�����

        if(gzgets(gzfpI, line, 1024))
        {
            info += line;

            // ��� line
            memset(line, '\0', sizeof(line));

            // ����������з���˵��һ��û������������
            while (info.find("\n") == string::npos && gzgets(gzfpI, line, 1024))
            {
                info += line;

                // ��� line
                memset(line, '\0', sizeof(line));
            }
            
            // ȥ�����з�
            if (info.size() > 0)
            {
                info = strip(info, '\n');
            }
            else  // ���� false
            {
                return false;
            }
            
            // ��ֵ
            INFOSTRUCTTMP_.INFO = info;

            // ע����
            if (info.find("#") != string::npos)
            {
                return true;
            }

            INFOSTRUCTTMP_.INFOVec = split(info, "\t");  // �ָ�

            INFOSTRUCTTMP_.CHROM = INFOSTRUCTTMP_.INFOVec[0];  // CHROM
            INFOSTRUCTTMP_.POS = stoi(INFOSTRUCTTMP_.INFOVec[1]);  // POS
            INFOSTRUCTTMP_.REF = INFOSTRUCTTMP_.INFOVec[3];  // REF
            INFOSTRUCTTMP_.ALTVec = split(INFOSTRUCTTMP_.INFOVec[4], ",");  // 'ALT' �� 'vec'
            if (isdigit(INFOSTRUCTTMP_.INFOVec[5][0]))  // ��������ֵĻ���ֵ QUAL ,������ 0.0
            {
                INFOSTRUCTTMP_.QUAL = stod(INFOSTRUCTTMP_.INFOVec[5]);  // QUAL
            }
            INFOSTRUCTTMP_.FILTER = INFOSTRUCTTMP_.INFOVec[6];  // FILTER

            INFOSTRUCTTMP_.LEN = INFOSTRUCTTMP_.REF.size();  // ref���г���
            INFOSTRUCTTMP_.END = INFOSTRUCTTMP_.POS + INFOSTRUCTTMP_.LEN - 1;  // ��ֹλ��
        }
        else  // ������ֱ�ӷ��� false
        {
            return false;
        }

        return true;
    }

    /**
     * @brief �ر��ļ�
     * 
     * @return 0
    **/
    int VCFOPEN::close(

    )
    {
        // �ر��ļ�
        gzclose(gzfpI);

        return 0;
    }


    /**
	 * ����MAF��MISSRATE.
     * 
     * @param LEN     ref�ĳ���
     * @param ALTVec  vector<qrySeq>
     * 
     * @return string   TYPE: SNP, InDel, Deletion, Insertion, Other
	**/
    string VCFOPEN::get_TYPE(
        const uint32_t & LEN,
        const vector<string> & ALTVec
    )
    {
        string TYPE = "";

        uint32_t qryLEN = 0;
        for (auto iter1 : ALTVec)
        {
            uint32_t qryLENTmp = iter1.size();
            qryLEN = max(qryLEN, qryLENTmp);  // ����� qrySeq
        }

        uint32_t size = qryLEN - LEN;
        if (size == 0)
        {
            TYPE = "SNP";
        }
        else if (size >= -49 && size <= 49)
        {
            TYPE = "InDel";
        }
        else if (size < -49)
        {
            TYPE = "Deletion";
        }
        else if (size > 49)
        {
            TYPE = "Insertion";
        }
        else
        {
            TYPE = "Other";
        }
        
        return TYPE;
    }
    

    /**
	 * ��ȡλ��������б�.
	 *
	 * @param INFOVec  INFOVec
     * 
     * 
     * @return GTVecMap   map<int, vector<string> >,  map<idx, vector<GTString> >
	**/
    map<int, vector<string> > VCFOPEN::get_gt(
        const vector<string> & INFOVec
    )
    {
        map<int, vector<string> > GTVecMap;  // ����Line���͵�map,  map<idx, vector<GTString> >

        vector<string> formatVec;  // FORMAT�ַ����
        int formatIndex = 8; // FORMAT������
        formatVec = split(INFOVec[formatIndex], ":");

        int gtIndex = distance(formatVec.begin(), find(formatVec.begin(), formatVec.end(), "GT"));  // ��ȡGT������λ��

        if (gtIndex == formatVec.size())  // �ж�index�Ƿ���ڣ������ڵĻ����ػ����Ͷ�Ϊ0��
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: no [GT] information in FORMAT -> " << INFOVec[0] << ":" << INFOVec[1] << endl;
            return GTVecMap;
        }
        else  // ������ڣ�����б���
        {
            string GTString;  // �洢�������ֶ�

            for (size_t i = 9; i < INFOVec.size(); i++)  // �ӵھ��п�ʼѭ��
            {
                GTString = split(INFOVec[i], ":")[gtIndex];  // gt�ֶ�

                // ���� '.' ��λ��ֱ������
                if (GTString.find(".") != string::npos)
                {
                    continue;
                }
                
                string splitStr;  // gt�еķָ���
                if (GTString.find("/") != string::npos)  // �жϡ�/���ָ���
                {
                    splitStr = "/";
                }
                else if (GTString.find("|") != string::npos)  // �жϡ�|��Ϊ�ָ���
                {
                    splitStr = "|";
                }
                
                // ���֪���ָ���ٸ�ֵ
                if (splitStr.size() > 0)
                {
                    GTVecMap[i - 9] = split(GTString, splitStr);  // ��ֵ
                }
            }
        }

        return GTVecMap;
    }


    /**
	 * ����MAF��MISSRATE.
	 *
	 * @param GTVecMap    map<int, vector<string> >,  map<idx, vector<GTString> >
     * @param sampleNum   �ܵ�sample����
     * 
     * 
     * @return tuple<double, double>   tuple<MAF, MISSRATE>
	**/
    tuple<double, double> VCFOPEN::calculate(
        const map<int, vector<string> > & GTVecMap, 
        uint32_t sampleNum
    )
    {
        double MAFTMP;  // ��С��λ����Ƶ��
        double MISSRATETMP;  // ȱʧ��

        /* ******************************** ���� MAF ******************************** */ 
        map<string, uint32_t> GTFreMapTmp;  // ��¼������Ƶ��  map<gt, fre>
        uint16_t ploidy = 1;  // ��¼����

        for (auto iter1 : GTVecMap)  // map<idx, vector<GTString> >
        {
            for (auto iter2 : iter1.second)  // vector<GTString>
            {
                GTFreMapTmp[iter2]++;  // �����Ͷ�Ӧ��Ƶ�ʼ�1
            }
            // �����1�ټ���
            if (ploidy == 1)
            {
                ploidy = iter1.second.size();  //��¼����
            }
        }

        map<uint32_t, string> GTFreMapTmpTmp;  // mapת��, map<fre, gt>
        for (auto iter1 : GTFreMapTmp)  // map<gt, fre>
        {
            GTFreMapTmpTmp[iter1.second] = iter1.first;
        }

        if (GTFreMapTmpTmp.size() > 1)
        {
            // �����ڶ��Ǵε�λ����
            auto iter1 = GTFreMapTmpTmp.end();
            iter1--; iter1--;
            // frequence/(ploidy*N)
            MAFTMP = iter1->first/(double)(ploidy*sampleNum);
        }
        
        /* ******************************** ���� MISSRATE ******************************** */ 
        // 1- number/allNumber
        MISSRATETMP = 1 - (GTVecMap.size()/(double)sampleNum);

        return make_tuple(MAFTMP, MISSRATETMP);
    }

} // namespace VCFOPEN

#endif