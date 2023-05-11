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
     * @brief 统计vcf文件
     * 
     * @param vcfFileName  输入vcf文件
     * 
     * @return 0
    **/
    int vcf_count(
        const string & vcfFileName
    )
    {
        // 记录各种变异的数量
        uint32_t snpNum = 0;
        uint32_t indelNum = 0;
        uint32_t insNum = 0;
        uint32_t delNum = 0;
        // uint32_t invNum = 0;
        // uint32_t dupNum = 0;
        uint32_t otherNum = 0;

        // 存储vcf信息
        VCFOPEN::VCFINFOSTRUCT INFOSTRUCTTMP;

        // 初始化
        VCFOPEN::VCFOPEN VCFOPENCLASS;
        VCFOPENCLASS.init(
            vcfFileName
        );

        // 打开vcf文件
        VCFOPENCLASS.open();

        // 如果没有遍历完，继续
        while (VCFOPENCLASS.read(INFOSTRUCTTMP))
        {
            // 有长度信息再计算
            if (INFOSTRUCTTMP.POS > 0)
            {
                // 遍历qry序列列表
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

        // 关闭vcf文件
        VCFOPENCLASS.close();

        // 打印结果
        cout << "SNP\tInDels\tDeletion\tInsertion\tOther\n";
        cout << snpNum << "\t" << indelNum << "\t" << delNum << "\t" << insNum << "\t" << otherNum << endl;
        
        return 0;
    }

} // namespace VCFCOUNT

#endif