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

    // �����ļ���
    gzFile gzfp = gzopen(inputVcf.c_str(), "rb");

    // ����Ƿ���ע����
    string::size_type idx;

    // ��¼��һ�������start��Ⱦɫ���
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
        char line[1024]; // һ��ֻ��1024�ֽڵ�����

        while(gzgets(gzfp, line, 1024))
        {
            informations += line;
            if (informations.find("\n") != string::npos) // һ�н���
            {
                if(informations.empty())
                {
                    continue;
                }

                informations = strip(informations, '\n'); // ȥ�����з�

                idx = informations.find("#");
                // ����ע����
                if (idx == string::npos)
                {
                    vector<string> informationsVec = split(informations, "\t");
                    string chromosome = informationsVec[0];
                    long long int refStart = stol(informationsVec[1]);

                    // ������µ�Ⱦɫ�壬����ʼλ�ù���
                    if (chromosome != preChromosome)
                    {
                        preChromosome = chromosome;
                        preRefStart = 0;
                    }

                    // ����һ���������ߵ���
                    if (refStart >= preRefStart)
                    {
                        preRefStart = refStart;
                    }
                    else // ����һ������С������vcfû���������û���������
                    {
                        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: unsorted file -> "<< inputVcf << " " 
                            << chromosome << " " << preRefStart << ">" << refStart << endl;
                        exit(1);
                    }
                }
                // ����ַ���
                informations.clear();
                string().swap(informations);
            }
        }
    }
    // �ͷ��ڴ棬�ر��ļ�
    gzclose(gzfp);

    return 0;
}

#endif