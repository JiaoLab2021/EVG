#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <regex>
#include <getopt.h>
#include "../include/count.hpp"
#include "../include/get_time.hpp"

using namespace std;

void help_count(char** argv);

int main_count(int argc, char** argv)
{
    // 输入文件
    string fastqFileName1 = "";
    string fastqFileName2 = "";
    string fastaFileName = "";

    // 输出文件
    string outputFileName = "";

    // read分割数目
    int readSplitNum = 500;

    // 线程数
	int threadsNum = 10;

    // 输入参数
    int c;
    while (true)
    {
        static const struct option long_options[] = 
        {
            {"fastq", required_argument, 0, 'i'},
            {"fasta", required_argument, 0, 'I'},
            {"out", required_argument, 0, 'o'},
            {"number", required_argument, 0, 'n'},
            {"threads", required_argument, 0, 't'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "i:I:o:n:t:h", long_options, &option_index);

        if (c == -1)
            break;

        switch (c)
        {
        case 'i':
            if (fastqFileName1.empty())
            {
                fastqFileName1 = optarg;
            }
            else
            {
                fastqFileName2 = optarg;
            }
            break;
        case 'I':
            fastaFileName = optarg;
            break;
        case 'o':
            outputFileName = optarg;
            break;
        case 'n':
            readSplitNum = stoi(optarg);
            break;
        case 't':
            threadsNum = stoi(optarg);
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

    if (argc <= 2)
    {
        help_count(argv);
        return 1;
    }

    // print log
    cerr << "[" << __func__ << "::" << getTime() << "] " 
         << "Running.\n";
    
    count::countStruct countOut;
    if (fastqFileName1.length() > 0)
    {
        count::fastq_a_count(fastqFileName1, 
                             fastqFileName2, 
                             countOut, 
                             outputFileName, 
                             threadsNum, 
                             readSplitNum);
    }
    else if (fastaFileName.length() > 0)
    {
        count::fastq_a_count(fastaFileName, 
                             "", 
                             countOut, 
                             outputFileName, 
                             threadsNum, 
                             readSplitNum);
    }
    else
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
             << "parameter error: -i / -I.\n";
        exit(1);
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " 
         << "Done.\n";

    return 0;
}

// 帮助文档
void help_count(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -i / -I [options]" << endl
       << "calculate the number of bases and reads of fasta/q" << endl
       << endl
       << "required arguments:" << endl
       << "    -i, --fastq     FILE    input fastq, possibly compressed, two are allowed, one for each mate" << endl
       << "    -I, --fasta     FILE    input fasta, possibly compressed" << endl
       << endl
	   << "optional arguments:" << endl
       << "    -o, --out       FILE    output file name [stdout]" << endl
       << "    -n, --number    INT     number of reads used by each thread [500]" << endl
       << "    -t, --threads   INT     number of compute threads to use [10]" << endl
       << "    -h, --help              print this help document" << endl;
}
