// g++ sample_fast.cpp -o sample -lz
#include <fstream>
#include <string>
#include <iostream>
#include "zlib.h"
#include <vector>
#include <algorithm>
#include <getopt.h>
#include "../include/kseq.h"
#include "../include/get_time.hpp"
#include "../include/sample_fast.hpp"
#include "../include/get_time.hpp"

using namespace std;

void help_sample(char** argv);

int main_sample(int argc, char** argv)
{
    // 输入文件
	string imputFileName;

    // 下采样的比例
    double frac = 0;

    // seedNum
	long int seedNum = 11;

    // mode
    string mode = "2";

    // 输入参数
    int c;
    
    while (true)
    {
        static const struct option long_options[] = 
		{
			{"input", required_argument, 0, 'i'},
            {"frac", required_argument, 0, 'f'},

            {"seed", required_argument, 0, 's'},
            {"mode", required_argument, 0, 'm'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "i:f:s:m:h", long_options, &option_index);

        if (c == -1)
            break;
        
        switch (c)
        {
        case 'i':
            imputFileName = optarg;
            break;
        case 'f':
            frac = stod(optarg);
            break;
        case 's':
            seedNum = stol(optarg);
            break;
        case 'm':
            mode = optarg;
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

    if (imputFileName.empty())
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
             << "parameter error: -i.\n";
        help_sample(argv);
        exit(1);
    }

    // print log
    cerr << "[" << __func__ << "::" << getTime() << "] " 
         << "Running.\n";

    // 判断frc是否为0
    if (frac == 0 || frac == 1)
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
             << "Skip: frac=" << frac << ".\n";
        exit(1);
    }

    sample_fast::stk_sample(
        imputFileName, 
        frac, 
        seedNum, 
        mode
    );

    // log
    cerr << "[" << __func__ << "::" << getTime() << "] " 
         << "Done.\n";
    
    return 0;
}

// 帮助文档
void help_sample(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -i -f [options]" << endl
       << "sample sequences by proportion" << endl
       << endl
       << "required arguments:" << endl
	   << "    -i, --input     FILE     input file, possibly compressed" << endl
       << "    -f, --frac      FLOAT    frac(0-1) / number(>1)" << endl
       << endl
	   << "optional arguments:" << endl
       << "    -s, --seed      INT      RNG seed [11]" << endl
       << "    -m, --mode      INT      1-fast, 2-reduced memory [2]" << endl
       << "    -h, --help               print this help document" << endl;
}