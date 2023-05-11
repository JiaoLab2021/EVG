// g++ vcf_breakpoint.cpp -o vcf_breakpoint -lz
#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <getopt.h>
#include <iomanip>
#include <cstdlib>
#include "../include/get_time.hpp"
#include "../include/vcf_breakpoint.hpp"

using namespace std;

void help_breakpoint(char** argv);

int main_breakpoint(int argc, char** argv)
{
    // 输入文件
    string vcfFileName;
    
    // 输出文件前缀
    string prefix = "breakpoint";

    // 断点误差大小
    int breakpointErrorSize = 1;

    // 输入参数
    int c;
    while (true)
    {
        static const struct option long_options[] = {
            {"vcf", required_argument, 0, 'v'},

            {"prefix", required_argument, 0, 'p'},
            {"size", required_argument, 0, 's'},

            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "v:p:s:h", long_options, &option_index);

        if (c == -1)
            break;

        switch (c)
        {
        case 'v':
            vcfFileName = optarg;
            break;
        case 'p':
            prefix = optarg;
            break;;
        case 's':
            breakpointErrorSize = stoi(optarg);
            break;
        case 'h':
        case '?':
            help_breakpoint(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    // 检查参数是否正确
    if (vcfFileName.empty())
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
             << "Parameter error.\n";
        help_breakpoint(argv);
        return 1;
    }
    
    cerr << "[" << __func__ << "::" << getTime() << "] " 
         << "Running.\n";
    
    BREAKPOINT::vcf_breakpoint(vcfFileName, prefix, breakpointErrorSize);
    
    cerr << "[" << __func__ << "::" << getTime() << "] " 
         << "Done.\n";

    return 0;
}

// 帮助文档
void help_breakpoint(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -v [options]" << endl
       << "set error for breakpoint of variations" << endl
       << endl
       << "required arguments (Note: vcf files must be sorted):" << endl
       << "    -v, --vcf           FILE       vcf file to set breakpoint error" << endl
       << endl
	   << "optional arguments:" << endl
       << "    -p, --prefix        STRING     output prefix [breakpoint]" << endl
       << "    -s, --size          INT        breakpoint error size [1]" << endl
       << endl
       << "    -h, --help                     print this help document" << endl;
}