#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <getopt.h>
#include "get_time.hpp"
#include "strip_split_join.hpp"

using namespace std;

namespace sample
{
    void save(string outputFileName, string outTxt);

    // 打开fastq.gz文件
    void sample(string inputFileName1, string inputFileName2, const vector<long long int> & randVecSel, string prefix)
    {
        if (inputFileName1.size() > 0 && inputFileName2.size() > 0) // 双端测序
        {
            // 输出文件名
            string outputFileName1 = prefix + "." + split(inputFileName1, "/")[split(inputFileName1, "/").size() - 1];
            string outputFileName2 = prefix + "." + split(inputFileName2, "/")[split(inputFileName2, "/").size() - 1];
            if (outputFileName1.rfind(".gz") == string::npos && outputFileName1.rfind(".GZ") == string::npos)
            {
                outputFileName1 += ".gz";
            }
            if (outputFileName2.rfind(".gz") == string::npos && outputFileName2.rfind(".GZ") == string::npos)
            {
                outputFileName2 += ".gz";
            }

            // 看文件是否存在，存在则清空文件
            // 输出文件流
            gzFile gzfpO1 = gzopen(outputFileName1.c_str(), "wb");
            gzFile gzfpO2 = gzopen(outputFileName2.c_str(), "wb");
            gzclose(gzfpO1);
            gzclose(gzfpO2);
            
            // 输入文件流
            gzFile gzfpI1 = gzopen(inputFileName1.c_str(), "rb");
            gzFile gzfpI2 = gzopen(inputFileName2.c_str(), "rb");

            // 记录read数
            long long int readNum = 0;

            // 输出字符串
            string outTxt1;
            string outTxt2;

            // 记录在vec中的index
            long long int index = 0;

            // 打开文件
            if(!gzfpI1 || !gzfpI2)
            {
                cerr << "[" << __func__ << "::" << getTime() << "] '" 
                     << inputFileName1 << "' or '" << inputFileName2
                     << "': No such file or directory.\n";
                exit(1);
            }
            else
            {
                kseq_t *ks1;
                kseq_t *ks2;
                ks1 = kseq_init(gzfpI1);
                ks2 = kseq_init(gzfpI2);
            
                while(kseq_read(ks1) >= 0 && kseq_read(ks2) >= 0)
                {
                    if (randVecSel[index] == readNum) // 如果是选择的read，则添加到outTxt里
                    {
                        string readName1 = ks1->name.s;
                        readName1 = "@" + readName1;
                        string readSeq1 = ks1->seq.s;
                        string readQual1 = ks1->qual.s;

                        string readName2 = ks2->name.s;
                        readName2 = "@" + readName2;
                        string readSeq2 = ks2->seq.s;
                        string readQual2 = ks2->qual.s;

                        // 判断fastq有没有接头序列信息
                        string readNameLast1;
                        string readNameLast2;
                        if (ks1->comment.s)
                        {
                            readNameLast1 = ' ' + string(ks1->comment.s);
                            readName1 += readNameLast1;
                        }
                        else
                        {
                            readNameLast1 = "";
                            readName1 += " 1";
                        }

                        if (ks2->comment.s)
                        {
                            readNameLast2 = ' ' + string(ks2->comment.s);
                            readName2 += readNameLast2;
                        }
                        else
                        {
                            readNameLast2 = "";
                            readName2 += " 1";
                        }

                        outTxt1 += readName1 + 
                                readNameLast1 + 
                                "\n" + 
                                readSeq1 + 
                                "\n+\n" + 
                                readQual1 +
                                "\n";

                        outTxt2 += readName2 + 
                                readNameLast2 + 
                                "\n" + 
                                readSeq2 + 
                                "\n+\n" + 
                                readQual2 +
                                "\n";

                        if (outTxt1.length() > 100000000) // 每100Mb写一次
                        {
                            // 写入磁盘
                            save(outputFileName1, outTxt1);
                            save(outputFileName2, outTxt2);

                            // 清空字符串
                            outTxt1.clear();
                            outTxt2.clear();
                            string().swap(outTxt1);
                            string().swap(outTxt2);
                        }

                        // index指向下一个元素
                        index++;
                    }
                    // read数+1
                    readNum++;
                }

                // 写入到磁盘中
                save(outputFileName1, outTxt1);
                save(outputFileName2, outTxt2);

                // 清空字符串
                outTxt1.clear();
                outTxt2.clear();
                string().swap(outTxt1);
                string().swap(outTxt2);

                // 释放内存，关闭文件
                kseq_destroy(ks1);
                kseq_destroy(ks2);
                gzclose(gzfpI1);
                gzclose(gzfpI2);
            }
        }
        else if (inputFileName1.size() > 0 && inputFileName2.size() == 0) // 单端测序
        {
            // 输出文件名
            string outputFileName1 = prefix + "." + split(inputFileName1, "/")[split(inputFileName1, "/").size() - 1];
            if (outputFileName1.rfind(".gz") == string::npos && outputFileName1.rfind(".GZ") == string::npos)
            {
                outputFileName1 += ".gz";
            }

            // 看文件是否存在，存在则清空文件
            // 输出文件流
            gzFile gzfpO1 = gzopen(outputFileName1.c_str(), "wb");
            gzclose(gzfpO1);
            
            // 输入文件流
            gzFile gzfpI1 = gzopen(inputFileName1.c_str(), "rb");

            // 记录read数
            long long int readNum = 0;

            // 输出字符串
            string outTxt1;

            // 记录在vec中的index
            long long int index = 0;

            // 打开文件
            if(!gzfpI1)
            {
                cerr << "[" << __func__ << "::" << getTime() << "] '" 
                     << inputFileName1
                     << "': No such file or directory.\n";
                    << endl;
                exit(1);
            }
            else
            {
                kseq_t *ks1;
                ks1 = kseq_init(gzfpI1);
            
                while(kseq_read(ks1) >= 0)
                {
                    if (randVecSel[index] == readNum) // 如果是选择的read，则添加到outTxt里
                    {
                        string readName1 = ks1->name.s;
                        readName1 = "@" + readName1;
                        string readSeq1 = ks1->seq.s;
                        string readQual1 = ks1->qual.s;

                        // 判断fastq有没有接头序列信息
                        string readNameLast1;
                        if (ks1->comment.s)
                        {
                            readNameLast1 = ' ' + string(ks1->comment.s);
                            readName1 += readNameLast1;
                        }
                        else
                        {
                            readNameLast1 = "";
                            readName1 += " 1";
                        }

                        outTxt1 += readName1 + 
                                readNameLast1 + 
                                "\n" + 
                                readSeq1 + 
                                "\n+\n" + 
                                readQual1 +
                                "\n";

                        if (outTxt1.length() > 100000000) // 每100Mb写一次
                        {
                            // 写入磁盘
                            save(outputFileName1, outTxt1);

                            // 清空字符串
                            outTxt1.clear();
                            string().swap(outTxt1);
                        }

                        // index指向下一个元素
                        index++;
                    }
                    // read数+1
                    readNum++;
                }

                // 写入到磁盘中
                save(outputFileName1, outTxt1);

                // 清空字符串
                outTxt1.clear();
                string().swap(outTxt1);

                // 释放内存，关闭文件
                kseq_destroy(ks1);
                gzclose(gzfpI1);
            }
        }
        else
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                    << " please input sequencing file." 
                    << endl;
            exit(1);
        }
    }

    // 保存fastq.gz文件
    void save(string outputFileName, string outTxt)
    {
        // 输出文件流
        gzFile gzfp = gzopen(outputFileName.c_str(), "ab");

        // 打开文件
        if(!gzfp)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] '" 
                 << outputFileName
                 << "': No such file or directory.\n";
            exit(1);
        }
        else
        {
            gzwrite(gzfp, outTxt.c_str(), outTxt.length());
        }

        // 释放内存，关闭文件
        gzclose(gzfp);
    }

    // 随机一组数据
    vector<long long int> randVector(long long int num)
    {
        vector<long long int> result;
        result.clear();
        result.reserve(num);
        srand((long long int)time(0));
        for (size_t i = 0; i < num; i++)
        {
            result.push_back(i);
        }
        long long int p1;
        long long int p2;
        long long int temp;
        long long int count = num;

        while (--num)
        {
            p1 = num;
            p2 = rand() % num;
            temp = result[p1];
            result[p1] = result[p2];
            result[p2] = temp;
        }
        return result;
    }
}