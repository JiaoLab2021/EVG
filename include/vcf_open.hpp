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

    // 清空结构体
    void clear();
};


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

#endif