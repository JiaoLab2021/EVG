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
#include "../include/vcf_split.hpp"

using namespace std;

void help_split(char** argv);

int main_split(int argc, char** argv)
{
    // �����ļ�
    string vcfFileName;

    // ������ȡsnp��indelλ�õ�vcf�ļ�
    string baseVcfFileName;
    
    // ��ֵ�ģʽ
    string mode = "type";

    // snp��indel��sv�ľ���
    int length = 100;

    // ����ļ�ǰ׺
    string prefix = "split";

    // �Ƿ����FILTER�н��й���
    string filterBool = "true";

    // �������
    int c;
    while (true)
    {
        static const struct option long_options[] = {
            {"vcf", required_argument, 0, 'v'},

            {"mode", required_argument, 0, 'm'},
            {"length", required_argument, 0, 'l'},
            {"basevcf", required_argument, 0, 'V'},
            {"prefix", required_argument, 0, 'p'},
            {"filter", required_argument, 0, 'f'},

            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "v:m:l:V:p:f:h", long_options, &option_index);

        if (c == -1)
            break;

        switch (c)
        {
        case 'v':
            vcfFileName = optarg;
            break;
        case 'm':
            mode = optarg;
            break;
        case 'l':
            length = stoi(optarg);
            break;
        case 'V':
            baseVcfFileName = optarg;
            break;
        case 'p':
            prefix = optarg;
            break;;
        case 'f':
            filterBool = optarg;
            break;
        case 'h':
        case '?':
            help_split(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    // �������Ƿ���ȷ
    if (vcfFileName.empty() || (mode != "type" && mode != "number"))
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error.\n";
        help_split(argv);
        return 1;
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running.\n";
    
    if (mode == "type")
    {
        VCFSPLIT::vcf_split_type(vcfFileName, prefix, filterBool);
    }
    else if (mode == "number")
    {
        VCFSPLIT::vcf_split_number(vcfFileName, baseVcfFileName, length, prefix, filterBool);
    }
    else
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error.\n";
        help_split(argv);
        return 1;
    }
    
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done.\n";

    return 0;
}

// �����ĵ�
void help_split(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -v [options]" << endl
       << "split variation by type or number of nearby SNP+Indels" << endl
       << endl
       << "required arguments (Note: vcf files must be sorted):" << endl
       << "    -v, --vcf           FILE       vcf file to be converted" << endl
       << endl
	   << "optional arguments:" << endl
       << "    -m, --mode          STRING     split standard: (type -> SNP+Indel+Ins+Del / number -> number of nearby SNP+Indel) [type]" << endl
       << "    -l, --length        INT        distance of snp and indel from sv [100]" << endl      
       << "    -V, --basevcf       INT        vcf file used to extract snp and indel locations [same as -v]" << endl      
       << "    -p, --prefix        STRING     output prefix [split]" << endl
       << "    -f, --filter        BOOL       filter using the FILTER column (true/false) [true]" << endl
       << endl
       << "    -h, --help                     print this help document" << endl;
}