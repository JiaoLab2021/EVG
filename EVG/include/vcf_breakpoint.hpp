#ifndef vcf_breakpoint_hpp
#define vcf_breakpoint_hpp

#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include "zlib.h"
#include "strip_split_join.hpp"
#include "get_time.hpp"

using namespace std;

namespace BREAKPOINT
{
    int vcf_breakpoint(const string & vcfFileName, const string & prefix, const int & breakpointErrorSize)
    {
        // 输入文件流
        gzFile gzfpI = gzopen(vcfFileName.c_str(), "rb");
        if(!gzfpI)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "'" << vcfFileName << "': No such file or directory." << endl;
            exit(1);
        }

        // 记录转换后的结果
        string outTxt;

        // 输出文件名
        string outFileName = prefix + ".vcf.gz";

        // 输出文件流
        gzFile outFile = gzopen(outFileName.c_str(), "wb");

        if(!outFile)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                 << "'"
                 << outFileName
                 << "': No such file or directory." 
                 << endl;
            exit(1);
        }

        // 打开vcf文件进行转换
        string information = ""; // 临时字符串
        char line[1024]; // 一次只读1024字节的数据

        while(gzgets(gzfpI, line, 1024))
        {
            information += line;
            if (information.find("\n") != string::npos) // 一行结束
            {
                if(line != nullptr && line[0] == '\0')
                {
                    continue;
                }

                information = strip(information, '\n'); // 去掉换行符

                if (information.find("#") != string::npos) // 检查是否是注释行
                {
                    string headLine = information + "\n";
                    gzwrite(outFile, headLine.c_str(), headLine.length());
                }
                else
                {
                    vector<string> informationVec = split(information, "\t");

                    // vcf的各种信息
                    string refChr = informationVec[0];
                    long long int refStart = stol(informationVec[1]);

                    string refSeq = informationVec[3];
                    string qrySeqs = informationVec[4];
                    vector<string> qrySeqsVec = split(qrySeqs, ",");
                    long long int refLen = refSeq.length();
                    long long int refEnd = refStart + refLen - 1;
                    long long int svLen = qrySeqs.length() - refLen;


                    // 搜索INFO里边的end信息
                    smatch endResult;
                    regex endReg("END=(\\d+)");

                    string::const_iterator iterStart = informationVec[7].begin();
                    string::const_iterator iterEnd = informationVec[7].end();
                    regex_search(iterStart, iterEnd, endResult, endReg);

                    string endResultTmp = endResult[1];
                    if (endResultTmp.empty())
                    {
                        if (refEnd == stol(endResult[1]))
                        {
                            continue;
                        }
                        else
                        {
                            cerr << "[" << __func__ << "::" << getTime() << "] " <<  information << endl;
                        }
                    }

                    // 添加误差
                    refStart -= breakpointErrorSize;
                    informationVec[1] = to_string(refStart);
                    refEnd -= breakpointErrorSize;
                    informationVec[7] = regex_replace(informationVec[7], endReg, "END=" + to_string(refEnd));

                    // 添加到输出字符串中
                    for (size_t i = 0; i < informationVec.size(); i++)
                    {
                        if (i > 0)
                        {
                            outTxt += "\t";
                        }
                        outTxt += informationVec[i];
                    }
                    outTxt += "\n";

                    // 保存结果
                    if (outTxt.size() > 10000000) // 每10Mb写入一次
                    {
                        gzwrite(outFile, outTxt.c_str(), outTxt.length());

                        // 清空字符串
                        outTxt.clear();
                        string().swap(outTxt);
                    }
                }

                // 清空字符串
                information.clear();
                string().swap(information);
            }
        }
        // 关闭文件
        gzclose(gzfpI);

        // 最后写入一次
        gzwrite(outFile, outTxt.c_str(), outTxt.length());

        // 清空字符串
        outTxt.clear();
        string().swap(outTxt);

        gzclose(outFile);

        return 0;
    }
} // namespace BREAKPOINT

#endif