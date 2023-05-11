#ifndef count_hpp
#define count_hpp
#include <fstream>
#include <string>
#include <iostream>
#include "zlib.h"
#include "kseq.h"
#include "get_time.hpp"
#include <getopt.h>
#include <mutex>
#include <malloc.h>
#include "ThreadPool.hpp"

std::mutex mtx;

using namespace std;

// kseq.h 打开文件
KSEQ_INIT(gzFile, gzread)

namespace count
{
    struct countStruct
    {
        long long int readNum = 0;
        long long int readBase = 0;
    };

    int fastq_a_count_run(vector<long long int> reads1Vec, 
                          vector<long long int> reads2Vec, 
                          countStruct & countOut);

    // 打开fastq/a.gz文件
    void fastq_a_count(string inputFileName1, 
                       string inputFileName2, 
                       countStruct & countOut, 
                       const string & outputFileName, 
                       const int & threadsNum, 
                       const int & readSplitNum)
    {
        // 记录测序序列信息，用于多线程提交
        vector<long long int> reads1Vec;
        vector<long long int> reads2Vec;

        // 进程池
        ThreadPool pool(threadsNum);

        // 初始化线程池
        pool.init();

        if (inputFileName1.length() > 0 && inputFileName2.length() > 0) // 双端测序
        {
            // 输入文件流
            gzFile gzfp1 = gzopen(inputFileName1.c_str(), "rb");
            gzFile gzfp2 = gzopen(inputFileName2.c_str(), "rb");

            // 打开文件
            if(!gzfp1 || !gzfp2)
            {
                cerr << "[" << __func__ << "::" << getTime() << "] " 
                     << inputFileName1 << "' or '" << inputFileName2
                     << "': No such file or directory.\n";
                exit(1);
            }
            else
            {
                kseq_t *ks1;
                kseq_t *ks2;
                ks1 = kseq_init(gzfp1);
                ks2 = kseq_init(gzfp2);
            
                while(kseq_read(ks1) >= 0 && kseq_read(ks2) >= 0)
                {
                    long long int length1 = ks1->seq.l;
                    long long int length2 = ks2->seq.l;
                    reads1Vec.push_back(length1);
                    reads2Vec.push_back(length2);
                    countOut.readNum += 2;

                    if (reads1Vec.size() >= readSplitNum)
                    {
                        pool.submit(fastq_a_count_run, 
                                    reads1Vec, 
                                    reads2Vec, 
                                    ref(countOut));
                        // 清空vector
                        reads1Vec.clear();
                        reads2Vec.clear();
                        vector<long long int>().swap(reads1Vec);
                        vector<long long int>().swap(reads2Vec);
                    }

                    // 检查任务队列是否超过阈值，超过了等待，以防把数据一次性加载到内存中
                    while (pool.get_queue() >= threadsNum*100)
                    {
                        // 每隔0.5秒检查一次
                        sleep(0.5);
                    }
                }

                // 最后提交一次任务
                if (reads1Vec.size() > 0)
                {
                    pool.submit(fastq_a_count_run, 
                                reads1Vec, 
                                reads2Vec, 
                                ref(countOut));
                    // 清空vector
                    reads1Vec.clear();
                    reads2Vec.clear();
                    vector<long long int>().swap(reads1Vec);
                    vector<long long int>().swap(reads2Vec);
                }

                // 检查任务队列是否执行完，执行完则关闭线程池，否则每隔0.5s检查一次
                while (pool.get_queue() > 0)
                {
                    // 每隔0.5秒检查一次
                    sleep(0.5);
                }

                // 关闭线程池
                pool.shutdown();

                // 释放内存，关闭文件
                kseq_destroy(ks1);
                kseq_destroy(ks2);
                gzclose(gzfp1);
                gzclose(gzfp2);
            }
        }
        else if (inputFileName1.length() > 0 && inputFileName2.length() == 0) // 单端测序
        {
            // 输入文件流
            gzFile gzfp = gzopen(inputFileName1.c_str(), "rb");

            // 打开文件
            if(!gzfp)
            {
                cerr << "[" << __func__ << "::" << getTime() << "] " 
                     << inputFileName1
                     << "': No such file or directory.\n";
                exit(1);
            }
            else
            {
                kseq_t *ks;
                ks = kseq_init(gzfp);
            
                while(kseq_read(ks) >= 0)
                {
                    long long int length = ks->seq.l;
                    reads1Vec.push_back(length);
                    countOut.readNum++;

                    if (reads1Vec.size() >= readSplitNum)
                    {
                        pool.submit(fastq_a_count_run, 
                                    reads1Vec, 
                                    reads2Vec, 
                                    ref(countOut));
                        // 清空vector
                        reads1Vec.clear();
                        reads2Vec.clear();
                        vector<long long int>().swap(reads1Vec);
                        vector<long long int>().swap(reads2Vec);
                    }

                    // 检查任务队列是否超过阈值，超过了等待，以防把数据一次性加载到内存中
                    while (pool.get_queue() >= threadsNum*100)
                    {
                        // 每隔0.5秒检查一次
                        sleep(0.5);
                    }
                }

                // 最后提交一次任务
                if (reads1Vec.size() > 0)
                {
                    pool.submit(fastq_a_count_run, 
                                reads1Vec, 
                                reads2Vec, 
                                ref(countOut));
                    // 清空vector
                    reads1Vec.clear();
                    reads2Vec.clear();
                    vector<long long int>().swap(reads1Vec);
                    vector<long long int>().swap(reads2Vec);
                }

                // 检查任务队列是否执行完，执行完则关闭线程池，否则每隔0.5s检查一次
                while (pool.get_queue() > 0)
                {
                    // 每隔0.5秒检查一次
                    sleep(0.5);
                }

                // 关闭线程池
                pool.shutdown();

                // 释放内存，关闭文件
                kseq_destroy(ks);
                gzclose(gzfp);
            }
        }

        // 输出结果
        if (outputFileName.empty()) // 如果没有指定输出文件名，则打印到标准输出
        {
            cout << "readBase:" << countOut.readBase << "\n" 
                 << "readNum:" << countOut.readNum << "\n" 
                 << "readLen:" << countOut.readBase/countOut.readNum 
                 << endl;
        }
        else // 保存到文件
        {
            // 输出文件流
            ofstream outputFile;
            outputFile.open(outputFileName, ios::out);

            

            // 打开文件
            if(!outputFile)
            {
                cerr << "[" << __func__ << "::" << getTime() << "] " 
                    << "'"
                    << outputFileName 
                    << "': No such file or directory." 
                    << endl;
                outputFile.close();
                exit(1);
            }
            else
            {
                outputFile << "readBase:" << countOut.readBase << "\n" 
                        << "readNum:" << countOut.readNum << "\n" 
                        << "readLen:" << countOut.readBase/countOut.readNum 
                        << endl;
            }

            // 关闭文件
            outputFile.close();
        }
    }

    // fastq/a.gz多线程函数
    int fastq_a_count_run(vector<long long int> reads1Vec, 
                          vector<long long int> reads2Vec, 
                          countStruct & countOut)
    {
        // 临时参数，减少线程锁时间
        long long int readBaseTmp = 0;

        for (size_t i = 0; i < reads1Vec.size(); i++)
        {
            if( reads1Vec[i] >= 0 )
            {
                readBaseTmp += reads1Vec[i];
            }

            if (reads2Vec.size() > 0)
            {
                if ( reads2Vec[i] >= 0 )
                {
                    readBaseTmp += reads2Vec[i];
                }
            }
        }

        // 清空内存
        vector<long long int>().swap(reads1Vec);
        vector<long long int>().swap(reads1Vec);

        // 多线程数据锁
        std::lock_guard<std::mutex> mtx_locker(mtx);
        countOut.readBase += readBaseTmp;

        // 释放内存
        malloc_trim(0);	// 0 is for heap memory

        return 0;
    }
}

#endif