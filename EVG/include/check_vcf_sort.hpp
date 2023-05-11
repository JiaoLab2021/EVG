#ifndef check_vcf_sort_hpp
#define check_vcf_sort_hpp
#include <string>
#include <iostream>
#include "zlib.h"
#include "get_time.hpp"
#include "strip_split_join.hpp"

using namespace std;

int check_vcf_sort(string inputVcf)
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Checking sort: '"<< inputVcf << "'\n";

    // 输入文件流
    gzFile gzfp = gzopen(inputVcf.c_str(), "rb");

    // 检查是否是注释行
    string::size_type idx;

    // 记录上一个变异的start和染色体号
    string preChromosome;
    long long int preRefStart = 0;

    if(!gzfp)
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "'" << inputVcf << "': No such file or directory." << endl;
        exit(1);
    }
    else
    {
        string informations;
        char line[1024]; // 一次只读1024字节的数据

        while(gzgets(gzfp, line, 1024))
        {
            informations += line;
            if (informations.find("\n") != string::npos) // 一行结束
            {
                if(informations.empty())
                {
                    continue;
                }

                informations = strip(informations, '\n'); // 去掉换行符

                idx = informations.find("#");
                // 跳过注释行
                if (idx == string::npos)
                {
                    vector<string> informationsVec = split(informations, "\t");
                    string chromosome = informationsVec[0];
                    long long int refStart = stol(informationsVec[1]);

                    // 如果是新的染色体，则将起始位置归零
                    if (chromosome != preChromosome)
                    {
                        preChromosome = chromosome;
                        preRefStart = 0;
                    }

                    // 比上一个变异大或者等于
                    if (refStart >= preRefStart)
                    {
                        preRefStart = refStart;
                    }
                    else // 比上一个变异小，代表vcf没排序，提醒用户进行排序
                    {
                        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: unsorted file -> "<< inputVcf << " " 
                            << chromosome << " " << preRefStart << ">" << refStart << endl;
                        exit(1);
                    }
                }
                // 清空字符串
                informations.clear();
                string().swap(informations);
            }
        }
    }
    // 释放内存，关闭文件
    gzclose(gzfp);

    return 0;
}

#endif