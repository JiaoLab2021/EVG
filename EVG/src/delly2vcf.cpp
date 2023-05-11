#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <regex>
#include <getopt.h>
#include "../include/delly2vcf.hpp"

using namespace std;

void help_delly2vcf(char** argv);

int main_delly2vcf(int argc, char** argv)
{
    string dellyFilename;
    string refgenomeFilename;
    string outputFilename;

    // �������
    int c;
    while (true)
    {
        static const struct option long_options[] = 
        {
            {"vcf", required_argument, 0, 'v'},
            {"refgenome", required_argument, 0, 'r'},
            {"out", no_argument, 0, 'o'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "v:r:o:h", long_options, &option_index);

        if (c == -1)
            break;

        switch (c)
        {
        case 'v':
            dellyFilename = optarg;
            break;
        case 'r':
            refgenomeFilename = optarg;
            break;
        case 'o':
            outputFilename = optarg;
            break;
        case 'h':
        case '?':
            help_delly2vcf(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc <= 2)
    {
        help_delly2vcf(argv);
        return 1;
    }

    // ��������ļ�������
    if (outputFilename.empty())
    {
        outputFilename = split(dellyFilename, ".")[0] + ".convert.vcf";
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " 
         << "Running.\n";

    // ���������ֵ�
    map<string,string> seq_dict = DELLY2VCF::build_dict(& refgenomeFilename);

    // vcfת��
    DELLY2VCF::vcf_convert(& dellyFilename, seq_dict, & outputFilename);

    cerr << "[" << __func__ << "::" << getTime() << "] " 
         << "Done.\n";
    return 0;
}

// �����ĵ�
void help_delly2vcf(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -v -r [options]" << endl
       << "convert the format of delly to vcf" << endl
       << endl
       << "required arguments:" << endl
       << "    -v, --vcf       FILE    vcf file output by delly" << endl
       << "    -r, --refgenome FILE    refgenome genome (fasta)" << endl
       << endl
	   << "optional arguments:" << endl
       << "    -o, --out       FILE    output file name [xxx.convert.vcf]" << endl
       << endl
       << "    -h, --help              print this help document" << endl;
}