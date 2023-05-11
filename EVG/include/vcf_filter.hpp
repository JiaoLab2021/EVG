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
     * @brief 根据次等位基因频率和缺失率过滤SNPs
     * 
     * @param vcfFileName    输入vcf文件
     * @param outputFileName 输出文件名
     * @param MAF            次等位基因频率
     * @param MISSRATE       缺失率
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
        // 输入文件流
        // 存储vcf信息
        VCFOPEN::VCFINFOSTRUCT INFOSTRUCTTMP;
        VCFOPEN::VCFOPEN VCFOPENCLASS;
        VCFOPENCLASS.init(
            vcfFileName
        );
        // 打开vcf文件
        VCFOPENCLASS.open();


        // 输出文件流
        SAVE::SAVE SAVECLASS;
        SAVECLASS.init(
            outputFileName
        );
        SAVECLASS.open();


        // 临时存储输出字符串
        string outTxt = "";

        // 如果没有遍历完，继续
        while (VCFOPENCLASS.read(INFOSTRUCTTMP))
        {
            // 如果注释行，直接保存
            if (INFOSTRUCTTMP.INFO.find("#") != string::npos)
            {
                outTxt += INFOSTRUCTTMP.INFO + "\n";
                continue;
            }

            // 获取变异类型
            INFOSTRUCTTMP.TYPE = VCFOPENCLASS.get_TYPE(
                INFOSTRUCTTMP.LEN, 
                INFOSTRUCTTMP.ALTVec
            );

            if (INFOSTRUCTTMP.TYPE == "SNP")  // SNP时再判断是否过滤
            {
                double MAFTMP;  // 最小等位基因频率
                double MISSRATETMP;  // 缺失率

                // 获取所有的基因型   map<idx, vector<gtString>>
                map<int, vector<string> > GTVecMapTmp = VCFOPENCLASS.get_gt(
                    INFOSTRUCTTMP.INFOVec
                );

                // 如果只有一个基因型，跳过。
                if (GTVecMapTmp.size() <= 1)
                {
                    continue;
                }
                
                // 计算最小等位基因频率和缺失率
                tie(MAFTMP, MISSRATETMP) = VCFOPENCLASS.calculate(
                    GTVecMapTmp, 
                    INFOSTRUCTTMP.INFOVec.size() - 9
                );

                // 通过阈值了直接保存
                if (MAFTMP >= MAF && MISSRATETMP <= MISSRATE)
                {
                    outTxt += INFOSTRUCTTMP.INFO + "\n";
                }
            }
            else  // 其它类型的变异直接保存
            {
                outTxt += INFOSTRUCTTMP.INFO + "\n";
            }

            if (outTxt.size() >= 10000000)  // 每10m写一次
            {
                // 输出文件流
                SAVECLASS.save(
                    outTxt
                );

                // 清空
                outTxt.clear();
                string().swap(outTxt);
            }
        }

        if (outTxt.size() >= 0)  // 最后写一次
        {
            // 其它类型的变异直接保存
            SAVECLASS.save(
                outTxt
            );

            // 清空
            outTxt.clear();
            string().swap(outTxt);
        }

        // 关闭文件
        VCFOPENCLASS.close();
        SAVECLASS.close();
        
        return 0;
    }

} // namespace VCFFILTER

#endif