// g++ vcf_convert.cpp -o vcf_convert -lz
#include <fstream>
#include <string>
#include <iostream>
#include <regex>
#include <map>
#include <getopt.h>
#include <algorithm>
#include "../include/get_time.hpp"
#include "../include/vcf_convert.hpp"

using namespace std;

void help_convert(char** argv);

int main_convert(int argc, char** argv)
{
    // �����ļ��Ͳ���
    string referenceFilename;
    string vcfFilename;
    string outFilename = "vcfConvert.out.vcf.gz";
    int readLen = 350; // ����Ķ���

    // ������ֵ
    double MAF = 0.;  // ��С��λ����Ƶ��
    double MISSRATE = 1.0;  // ȱʧ��

    // �������
    int c;
    while (true)
    {
        static const struct option long_options[] = {
            {"reference", required_argument, 0, 'r'},
            {"vcf", required_argument, 0, 'v'},

            {"length", required_argument, 0, 'l'},
            {"maf", required_argument, 0, '1'},
            {"geno", required_argument, 0, '2'},
            {"out", required_argument, 0, 'o'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "r:v:l:1:2:o:h", long_options, &option_index);

        if (c == -1)
            break;

        switch (c)
        {
        case 'r':
            referenceFilename = optarg;
            break;
        case 'v':
            vcfFilename = optarg;
            break;
        case 'l':
            readLen = stoi(optarg);
            break;
        case '1':
            MAF = stod(optarg);
            break;
        case '2':
            MISSRATE = stod(optarg);
            break;
        case 'o':
            outFilename = optarg;
            break;
        case 'h':
        case '?':
            help_convert(argv);
            exit(1);
            break;
        default:
            abort();
        }
    }

    if (argc <= 2) {
        help_convert(argv);
        return 1;
    }

    // ��ӡlog
    cerr << "[" << __func__ << "::" << getTime() << "] " 
         << "Running.\n";

    // ����reference������
    CONVERT::refIndexStruct refIndexS;
    CONVERT::build_reference_index(
        referenceFilename, 
        refIndexS
    );

    // ���raed�ĳ��ȴ���1000������read length���ˣ���Ϊ������paragraph
    if (readLen > 1000)
    {
        readLen = 0;
    }

    // vcf�ļ�ת��
    CONVERT::vcf_convert(
        vcfFilename, 
        readLen, 
        MAF, 
        MISSRATE, 
        refIndexS, 
        outFilename
    );

    // ��ӡlog
    cerr << "[" << __func__ << "::" << getTime() << "] " 
         << "Done.\n";

    return 0;
}

// �����ĵ�
void help_convert(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -r -v [options]" << endl
       << "convert vcf files merged by vcftools to the format required by genome graph." << endl
       << endl
       << "required arguments:" << endl
       << "    -r, --reference   FILE     input FASTA reference" << endl
       << "    -v, --vcf         FILE     input VCF" << endl
       << endl
	   << "optional arguments:" << endl
       << "    -l, --length      INT      read length [350]" << endl
       << "    --maf             FLOAT    exclude variants with minor allele frequency lower than threshold [0.0]" << endl
       << "    --geno            FLOAT    exclude variants with missing call frequencies greater than threshold [1.0]" << endl
       << "    -o, --out         FILE     output file name [vcfConvert.out.vcf.gz]" << endl
       << endl
       << "    -h, --help                 print this help document" << endl;
}