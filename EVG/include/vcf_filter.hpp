#ifndef vcf_filter_hpp
#define vcf_filter_hpp

#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include "zlib.h"
#include <unordered_map>
#include <cmath>
#include "strip_split_join.hpp"
#include "get_time.hpp"
#include "sort_find.hpp"
#include "vcf_open.hpp"
#include "save.hpp"

using namespace std;

namespace VCFFILTER
{
    /**
     * @brief ���ݴε�λ����Ƶ�ʺ�ȱʧ�ʹ���SNPs
     * 
     * @param vcfFileName    ����vcf�ļ�
     * @param outputFileName ����ļ���
     * @param MAF            �ε�λ����Ƶ��
     * @param MISSRATE       ȱʧ��
     * 
     * @return 0
    **/
    int vcf_filter(
        const string & vcfFileName, 
        const string & outputFileName, 
        const double & MAF, 
        const double & MISSRATE
    )
    {
        // �����ļ���
        // �洢vcf��Ϣ
        VCFOPEN::VCFINFOSTRUCT INFOSTRUCTTMP;
        VCFOPEN::VCFOPEN VCFOPENCLASS;
        VCFOPENCLASS.init(
            vcfFileName
        );
        // ��vcf�ļ�
        VCFOPENCLASS.open();


        // ����ļ���
        SAVE::SAVE SAVECLASS;
        SAVECLASS.init(
            outputFileName
        );
        SAVECLASS.open();


        // ��ʱ�洢����ַ���
        string outTxt = "";

        // ���û�б����꣬����
        while (VCFOPENCLASS.read(INFOSTRUCTTMP))
        {
            // ���ע���У�ֱ�ӱ���
            if (INFOSTRUCTTMP.INFO.find("#") != string::npos)
            {
                outTxt += INFOSTRUCTTMP.INFO + "\n";
                continue;
            }

            // ��ȡ��������
            INFOSTRUCTTMP.TYPE = VCFOPENCLASS.get_TYPE(
                INFOSTRUCTTMP.LEN, 
                INFOSTRUCTTMP.ALTVec
            );

            if (INFOSTRUCTTMP.TYPE == "SNP")  // SNPʱ���ж��Ƿ����
            {
                double MAFTMP;  // ��С��λ����Ƶ��
                double MISSRATETMP;  // ȱʧ��

                // ��ȡ���еĻ�����   map<idx, vector<gtString>>
                map<int, vector<string> > GTVecMapTmp = VCFOPENCLASS.get_gt(
                    INFOSTRUCTTMP.INFOVec
                );

                // ���ֻ��һ�������ͣ�������
                if (GTVecMapTmp.size() <= 1)
                {
                    continue;
                }
                
                // ������С��λ����Ƶ�ʺ�ȱʧ��
                tie(MAFTMP, MISSRATETMP) = VCFOPENCLASS.calculate(
                    GTVecMapTmp, 
                    INFOSTRUCTTMP.INFOVec.size() - 9
                );

                // ͨ����ֵ��ֱ�ӱ���
                if (MAFTMP >= MAF && MISSRATETMP <= MISSRATE)
                {
                    outTxt += INFOSTRUCTTMP.INFO + "\n";
                }
            }
            else  // �������͵ı���ֱ�ӱ���
            {
                outTxt += INFOSTRUCTTMP.INFO + "\n";
            }

            if (outTxt.size() >= 10000000)  // ÿ10mдһ��
            {
                // ����ļ���
                SAVECLASS.save(
                    outTxt
                );

                // ���
                outTxt.clear();
                string().swap(outTxt);
            }
        }

        if (outTxt.size() >= 0)  // ���дһ��
        {
            // �������͵ı���ֱ�ӱ���
            SAVECLASS.save(
                outTxt
            );

            // ���
            outTxt.clear();
            string().swap(outTxt);
        }

        // �ر��ļ�
        VCFOPENCLASS.close();
        SAVECLASS.close();
        
        return 0;
    }

} // namespace VCFFILTER

#endif