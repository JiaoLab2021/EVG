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

// kseq.h ���ļ�
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

    // ��fastq/a.gz�ļ�
    void fastq_a_count(string inputFileName1, 
                       string inputFileName2, 
                       countStruct & countOut, 
                       const string & outputFileName, 
                       const int & threadsNum, 
                       const int & readSplitNum)
    {
        // ��¼����������Ϣ�����ڶ��߳��ύ
        vector<long long int> reads1Vec;
        vector<long long int> reads2Vec;

        // ���̳�
        ThreadPool pool(threadsNum);

        // ��ʼ���̳߳�
        pool.init();

        if (inputFileName1.length() > 0 && inputFileName2.length() > 0) // ˫�˲���
        {
            // �����ļ���
            gzFile gzfp1 = gzopen(inputFileName1.c_str(), "rb");
            gzFile gzfp2 = gzopen(inputFileName2.c_str(), "rb");

            // ���ļ�
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
                        // ���vector
                        reads1Vec.clear();
                        reads2Vec.clear();
                        vector<long long int>().swap(reads1Vec);
                        vector<long long int>().swap(reads2Vec);
                    }

                    // �����������Ƿ񳬹���ֵ�������˵ȴ����Է�������һ���Լ��ص��ڴ���
                    while (pool.get_queue() >= threadsNum*100)
                    {
                        // ÿ��0.5����һ��
                        sleep(0.5);
                    }
                }

                // ����ύһ������
                if (reads1Vec.size() > 0)
                {
                    pool.submit(fastq_a_count_run, 
                                reads1Vec, 
                                reads2Vec, 
                                ref(countOut));
                    // ���vector
                    reads1Vec.clear();
                    reads2Vec.clear();
                    vector<long long int>().swap(reads1Vec);
                    vector<long long int>().swap(reads2Vec);
                }

                // �����������Ƿ�ִ���ִ꣬������ر��̳߳أ�����ÿ��0.5s���һ��
                while (pool.get_queue() > 0)
                {
                    // ÿ��0.5����һ��
                    sleep(0.5);
                }

                // �ر��̳߳�
                pool.shutdown();

                // �ͷ��ڴ棬�ر��ļ�
                kseq_destroy(ks1);
                kseq_destroy(ks2);
                gzclose(gzfp1);
                gzclose(gzfp2);
            }
        }
        else if (inputFileName1.length() > 0 && inputFileName2.length() == 0) // ���˲���
        {
            // �����ļ���
            gzFile gzfp = gzopen(inputFileName1.c_str(), "rb");

            // ���ļ�
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
                        // ���vector
                        reads1Vec.clear();
                        reads2Vec.clear();
                        vector<long long int>().swap(reads1Vec);
                        vector<long long int>().swap(reads2Vec);
                    }

                    // �����������Ƿ񳬹���ֵ�������˵ȴ����Է�������һ���Լ��ص��ڴ���
                    while (pool.get_queue() >= threadsNum*100)
                    {
                        // ÿ��0.5����һ��
                        sleep(0.5);
                    }
                }

                // ����ύһ������
                if (reads1Vec.size() > 0)
                {
                    pool.submit(fastq_a_count_run, 
                                reads1Vec, 
                                reads2Vec, 
                                ref(countOut));
                    // ���vector
                    reads1Vec.clear();
                    reads2Vec.clear();
                    vector<long long int>().swap(reads1Vec);
                    vector<long long int>().swap(reads2Vec);
                }

                // �����������Ƿ�ִ���ִ꣬������ر��̳߳أ�����ÿ��0.5s���һ��
                while (pool.get_queue() > 0)
                {
                    // ÿ��0.5����һ��
                    sleep(0.5);
                }

                // �ر��̳߳�
                pool.shutdown();

                // �ͷ��ڴ棬�ر��ļ�
                kseq_destroy(ks);
                gzclose(gzfp);
            }
        }

        // ������
        if (outputFileName.empty()) // ���û��ָ������ļ��������ӡ����׼���
        {
            cout << "readBase:" << countOut.readBase << "\n" 
                 << "readNum:" << countOut.readNum << "\n" 
                 << "readLen:" << countOut.readBase/countOut.readNum 
                 << endl;
        }
        else // ���浽�ļ�
        {
            // ����ļ���
            ofstream outputFile;
            outputFile.open(outputFileName, ios::out);

            

            // ���ļ�
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

            // �ر��ļ�
            outputFile.close();
        }
    }

    // fastq/a.gz���̺߳���
    int fastq_a_count_run(vector<long long int> reads1Vec, 
                          vector<long long int> reads2Vec, 
                          countStruct & countOut)
    {
        // ��ʱ�����������߳���ʱ��
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

        // ����ڴ�
        vector<long long int>().swap(reads1Vec);
        vector<long long int>().swap(reads1Vec);

        // ���߳�������
        std::lock_guard<std::mutex> mtx_locker(mtx);
        countOut.readBase += readBaseTmp;

        // �ͷ��ڴ�
        malloc_trim(0);	// 0 is for heap memory

        return 0;
    }
}

#endif