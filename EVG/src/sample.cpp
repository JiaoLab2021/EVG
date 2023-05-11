// g++ sample.cpp -o sample -lz
#include <fstream>
#include <string>
#include <iostream>
#include "zlib.h"
#include <vector>
#include <algorithm>
#include <getopt.h>
#include "../include/kseq.h"
#include "../include/get_time.hpp"
#include "../include/sample.hpp"

using namespace std;

void help_sample(char** argv);

int main_sample(int argc, char** argv)
{
    // �����ļ�
	string inputFastqFile1;
    string inputFastqFile2;

    // ����ļ�ǰ׺
    string prefix = "sample";

    // �²����ı���
    float proportion = 0;

    // read����
    unsigned long long int raedNumber = 0;

    // �������
    int c;
    
    while (true)
    {
        static const struct option long_options[] = 
		{
			{"fastq", required_argument, 0, 'i'},
            {"proportion", required_argument, 0, 'p'},

            {"prefix", required_argument, 0, '1'},
            {"number", required_argument, 0, 'n'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "i:p:1:n:h", long_options, &option_index);

        if (c == -1)
            break;
        
        switch (c)
        {
        case 'i':
            if (inputFastqFile1.empty())
            {
                inputFastqFile1 = optarg;
            }
            else
            {
                inputFastqFile2 = optarg;
            }
            break;
        case 'p':
            proportion = stof(optarg);
            break;
        case '1':
            prefix = optarg;
            break;
        case 'n':
            raedNumber = stol(optarg);
            break;
        case 'h':
        case '?':
            help_sample(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (inputFastqFile1.empty() && inputFastqFile2.empty())
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
             << "parameter error: -i.\n";
        help_sample(argv);
        return 1;
    }

    // log
    cerr << "[" << __func__ << "::" << getTime() << "] " 
         << "Running.\n";

    // �ж�proportion�Ƿ�Ϊ0
    if (proportion == 0 || proportion >= 1)
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
             << "Skip: proportion=" << proportion << ".\n";
        exit(0);
    }

    // ����raed���������жϲ�������û�����˵Ļ�ֱ������
    if (raedNumber == 0)
    {
        count::countStruct countOut;
        count::fastq_a_count(inputFastqFile1, 
                             "", 
                             countOut, 
                             "", 
                             10, 
                             500);
        raedNumber = countOut.readNum;
    }
    
    // ���������
    vector<long long int> randVec = sample::randVector(raedNumber);

    // ����ȡ���������
    vector<long long int> randVecSel;
    long long int selNum = raedNumber*proportion;
    for (long long int i = 0; i < selNum; i++)
    {
        randVecSel.push_back(randVec[i]);
    }

    // vector����
    sort(begin(randVecSel), end(randVecSel));
    
    sample::sample(inputFastqFile1, inputFastqFile2, randVecSel, prefix);

    // log
    cerr << "[" << __func__ << "::" << getTime() << "] " 
         << "Done.\n";
    
    return 0;
}

// �����ĵ�
void help_sample(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -p -i [options]" << endl
       << "sample sequences by proportion" << endl
       << endl
       << "required arguments:" << endl
	   << "    -i, --fastq        FILE     input fastq, possibly compressed, two are allowed, one for each mate" << endl
       << "    -p, --proportion   FLOAT    sample by proportion (0-1)" << endl
       << endl
	   << "optional arguments:" << endl
       << "    --prefix           STRING   output filename prefix [sample]" << endl
       << "    -n, --number       INT      number of reads [auto]" << endl
       << "    -h, --help                  print this help document" << endl;
}