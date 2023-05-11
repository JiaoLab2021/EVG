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
        // 总的信息
        string INFO = "";
        vector<string> INFOVec;

        string CHROM = "";  // CHROM  [0]
        uint32_t POS = 0;  // POS  [1]
        string REF = "";  // REF  [3]
        vector<string> ALTVec;  // ALT 列表  [4]
        double QUAL = 0;  // QUAL  [5]
        string FILTER = ""; // FILTER  [6]

        string TYPE = "";  // 变异类型
        uint32_t LEN = 0;  // ref长度
        uint32_t END = 0;  // ref结束

        vector<string> GTVec; // 基因型列表
        double MAF = 0.0;  // 次等位基因频率
        double MISSRATE = 0.0;  // 缺失比例

        void clear();
    };

    // 清空结构体
    void VCFINFOSTRUCT::clear(

    )
    {
        // 总的信息
         INFO = "";
         INFOVec.clear();
         vector<string>().swap(INFOVec);

         CHROM = "";  // CHROM  [0]
         POS = 0;  // POS  [1]
         REF = "";  // REF  [3]
         ALTVec.clear();  // ALT 列表  [4]
         vector<string>().swap(ALTVec);
         QUAL = 0;  // QUAL  [5]
         FILTER = ""; // FILTER  [6]

         TYPE = "";  // 变异类型
         LEN = 0;  // ref长度
         END = 0;  // ref结束

         GTVec.clear(); // 基因型列表
         vector<string>().swap(GTVec);
         MAF = 0.0;  // 次等位基因频率
         MISSRATE = 0.0;  // 缺失比例
    }



    
    /**
     * @brief 打开vcf文件
     * 
     * @param vcfFileName_   the output of vcf file
     * 
     * @return
    **/
    class VCFOPEN
    {
    private:
        // vcf文件
        string vcfFileName;

        // 输入文件流
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
        // 找所有line的分型结果
        map<int, vector<string> > get_gt(
            const vector<string> & INFOVec
        );
        // 计算MAF和MISSRATE
        tuple<double, double> calculate(
            const map<int, vector<string> > & GTVecMap, 
            uint32_t sampleNum
        );
        int close();
    };

    /**
     * @brief 初始化
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
     * @brief 打开文件
     * 
     * @return 0
    **/
   int VCFOPEN::open(

   )
   {
        // 输入文件流
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
     * @brief 遍历文件
     * 
     * @param INFOSTRUCTTMP_  存储该行的内容
     * 
     * @return bool
    **/
    bool VCFOPEN::read(
        VCFINFOSTRUCT & INFOSTRUCTTMP_
    )
    {
        INFOSTRUCTTMP_.clear();

        string info = ""; // 临时字符串
        char line[1024]; // 一次只读1024字节的数据

        if(gzgets(gzfpI, line, 1024))
        {
            info += line;

            // 清空 line
            memset(line, '\0', sizeof(line));

            // 如果不含换行符，说明一行没结束，继续读
            while (info.find("\n") == string::npos && gzgets(gzfpI, line, 1024))
            {
                info += line;

                // 清空 line
                memset(line, '\0', sizeof(line));
            }
            
            // 去除换行符
            if (info.size() > 0)
            {
                info = strip(info, '\n');
            }
            else  // 返回 false
            {
                return false;
            }
            
            // 赋值
            INFOSTRUCTTMP_.INFO = info;

            // 注释行
            if (info.find("#") != string::npos)
            {
                return true;
            }

            INFOSTRUCTTMP_.INFOVec = split(info, "\t");  // 分割

            INFOSTRUCTTMP_.CHROM = INFOSTRUCTTMP_.INFOVec[0];  // CHROM
            INFOSTRUCTTMP_.POS = stoi(INFOSTRUCTTMP_.INFOVec[1]);  // POS
            INFOSTRUCTTMP_.REF = INFOSTRUCTTMP_.INFOVec[3];  // REF
            INFOSTRUCTTMP_.ALTVec = split(INFOSTRUCTTMP_.INFOVec[4], ",");  // 'ALT' 的 'vec'
            if (isdigit(INFOSTRUCTTMP_.INFOVec[5][0]))  // 如果是数字的话赋值 QUAL ,否则是 0.0
            {
                INFOSTRUCTTMP_.QUAL = stod(INFOSTRUCTTMP_.INFOVec[5]);  // QUAL
            }
            INFOSTRUCTTMP_.FILTER = INFOSTRUCTTMP_.INFOVec[6];  // FILTER

            INFOSTRUCTTMP_.LEN = INFOSTRUCTTMP_.REF.size();  // ref序列长度
            INFOSTRUCTTMP_.END = INFOSTRUCTTMP_.POS + INFOSTRUCTTMP_.LEN - 1;  // 终止位置
        }
        else  // 遍历完直接返回 false
        {
            return false;
        }

        return true;
    }

    /**
     * @brief 关闭文件
     * 
     * @return 0
    **/
    int VCFOPEN::close(

    )
    {
        // 关闭文件
        gzclose(gzfpI);

        return 0;
    }


    /**
	 * 计算MAF和MISSRATE.
     * 
     * @param LEN     ref的长度
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
            qryLEN = max(qryLEN, qryLENTmp);  // 找最长的 qrySeq
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
	 * 获取位点基因型列表.
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
        map<int, vector<string> > GTVecMap;  // 所有Line分型的map,  map<idx, vector<GTString> >

        vector<string> formatVec;  // FORMAT字符拆分
        int formatIndex = 8; // FORMAT所在列
        formatVec = split(INFOVec[formatIndex], ":");

        int gtIndex = distance(formatVec.begin(), find(formatVec.begin(), formatVec.end(), "GT"));  // 获取GT的索引位置

        if (gtIndex == formatVec.size())  // 判断index是否存在，不存在的话返回基因型都为0。
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: no [GT] information in FORMAT -> " << INFOVec[0] << ":" << INFOVec[1] << endl;
            return GTVecMap;
        }
        else  // 如果存在，则进行保存
        {
            string GTString;  // 存储基因型字段

            for (size_t i = 9; i < INFOVec.size(); i++)  // 从第九列开始循环
            {
                GTString = split(INFOVec[i], ":")[gtIndex];  // gt字段

                // 含有 '.' 的位点直接跳过
                if (GTString.find(".") != string::npos)
                {
                    continue;
                }
                
                string splitStr;  // gt中的分隔符
                if (GTString.find("/") != string::npos)  // 判断‘/’分隔符
                {
                    splitStr = "/";
                }
                else if (GTString.find("|") != string::npos)  // 判断‘|’为分隔符
                {
                    splitStr = "|";
                }
                
                // 如果知道分割符再赋值
                if (splitStr.size() > 0)
                {
                    GTVecMap[i - 9] = split(GTString, splitStr);  // 赋值
                }
            }
        }

        return GTVecMap;
    }


    /**
	 * 计算MAF和MISSRATE.
	 *
	 * @param GTVecMap    map<int, vector<string> >,  map<idx, vector<GTString> >
     * @param sampleNum   总的sample数量
     * 
     * 
     * @return tuple<double, double>   tuple<MAF, MISSRATE>
	**/
    tuple<double, double> VCFOPEN::calculate(
        const map<int, vector<string> > & GTVecMap, 
        uint32_t sampleNum
    )
    {
        double MAFTMP;  // 最小等位基因频率
        double MISSRATETMP;  // 缺失率

        /* ******************************** 计算 MAF ******************************** */ 
        map<string, uint32_t> GTFreMapTmp;  // 记录基因型频率  map<gt, fre>
        uint16_t ploidy = 1;  // 记录倍型

        for (auto iter1 : GTVecMap)  // map<idx, vector<GTString> >
        {
            for (auto iter2 : iter1.second)  // vector<GTString>
            {
                GTFreMapTmp[iter2]++;  // 基因型对应的频率加1
            }
            // 如果是1再计算
            if (ploidy == 1)
            {
                ploidy = iter1.second.size();  //记录倍型
            }
        }

        map<uint32_t, string> GTFreMapTmpTmp;  // map转换, map<fre, gt>
        for (auto iter1 : GTFreMapTmp)  // map<gt, fre>
        {
            GTFreMapTmpTmp[iter1.second] = iter1.first;
        }

        if (GTFreMapTmpTmp.size() > 1)
        {
            // 倒数第二是次等位基因
            auto iter1 = GTFreMapTmpTmp.end();
            iter1--; iter1--;
            // frequence/(ploidy*N)
            MAFTMP = iter1->first/(double)(ploidy*sampleNum);
        }
        
        /* ******************************** 计算 MISSRATE ******************************** */ 
        // 1- number/allNumber
        MISSRATETMP = 1 - (GTVecMap.size()/(double)sampleNum);

        return make_tuple(MAFTMP, MISSRATETMP);
    }

} // namespace VCFOPEN

#endif