// g++ vcf_split.cpp -o vcf_split -lz
#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <getopt.h>
#include <iomanip>
#include <cstdlib>
#include "../include/get_time.hpp"
#include "../include/vcf_count.hpp"

using namespace std;

void help_count(char** argv);

int main_count(int argc, char** argv)
{
    // �����ļ�
    string vcfFileName;

    // �������
    int c;
    while (true)
    {
        static const struct option long_options[] = {
            {"vcf", required_argument, 0, 'v'},

            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "v:h", long_options, &option_index);

        if (c == -1)
            break;

        switch (c)
        {
        case 'v':
            vcfFileName = optarg;
            break;
        case 'h':
        case '?':
            help_count(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    // �������Ƿ���ȷ
    if (vcfFileName.empty())
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error.\n";
        help_count(argv);
        return 1;
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running.\n";
    
    // ͳ�ƽ��
    VCFCOUNT::vcf_count(vcfFileName);
    
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done.\n";

    return 0;
}

// �����ĵ�
void help_count(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -v [options]" << endl
       << "count the number of alleles" << endl
       << endl
       << "required arguments:" << endl
       << "    -v, --vcf           FILE       vcf file to be converted" << endl
       << endl
       << "    -h, --help                     print this help document" << endl;
}