#ifndef vcf_count_hpp
#define vcf_count_hpp

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

using namespace std;

namespace VCFCOUNT
{
    /**
     * @brief ͳ��vcf�ļ�
     * 
     * @param vcfFileName  ����vcf�ļ�
     * 
     * @return 0
    **/
    int vcf_count(
        const string & vcfFileName
    )
    {
        // ��¼���ֱ��������
        uint32_t snpNum = 0;
        uint32_t indelNum = 0;
        uint32_t insNum = 0;
        uint32_t delNum = 0;
        // uint32_t invNum = 0;
        // uint32_t dupNum = 0;
        uint32_t otherNum = 0;

        // �洢vcf��Ϣ
        VCFOPEN::VCFINFOSTRUCT INFOSTRUCTTMP;

        // ��ʼ��
        VCFOPEN::VCFOPEN VCFOPENCLASS;
        VCFOPENCLASS.init(
            vcfFileName
        );

        // ��vcf�ļ�
        VCFOPENCLASS.open();

        // ���û�б����꣬����
        while (VCFOPENCLASS.read(INFOSTRUCTTMP))
        {
            // �г�����Ϣ�ټ���
            if (INFOSTRUCTTMP.POS > 0)
            {
                // ����qry�����б�
                for (auto qrySeq : INFOSTRUCTTMP.ALTVec)
                {
                    uint32_t qryLen = qrySeq.size();

                    int32_t svLen = qryLen -INFOSTRUCTTMP.LEN;

                    if (svLen == 0)
                    {
                        snpNum++;
                    }
                    else if (svLen <= 49 && svLen >= -49)
                    {
                        indelNum++;
                    }
                    else if (svLen < -49)
                    {
                        delNum++;
                    }
                    else if (svLen > 49)
                    {
                        insNum++;
                    }
                    else
                    {
                        otherNum++;
                    }
                }
            }
        }

        // �ر�vcf�ļ�
        VCFOPENCLASS.close();

        // ��ӡ���
        cout << "SNP\tInDels\tDeletion\tInsertion\tOther\n";
        cout << snpNum << "\t" << indelNum << "\t" << delNum << "\t" << insNum << "\t" << otherNum << endl;
        
        return 0;
    }

} // namespace VCFCOUNT

#endif