// g++ vcf_merge_depth.cpp -o vcf_merge_depth -lz
#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <vector>
#include <numeric>
#include <map>
#include <malloc.h>
#include <getopt.h>
#include <cmath>
#include "../include/get_time.hpp"
#include "../include/vcf_merge.hpp"


using namespace std;

void help_merge(char** argv);

int main_merge(int argc, char** argv)
{
    // EVG的运行模式
    string mode = "specific";

    string trueVcf;
    string ParagraphVcf;
    string GraphTyper2Vcf;
    string BayesTyperVcf;
    string mapVcf;
    string giraffeVcf;
    string GraphAlignerVcf;
    string PanGenieVcf;

    // 样本名称
    string sample_name = "";

    // 输出文件名字
    string outputFileName = "";

    // 输入参数
    int c;
    
    while (true)
    {
        static const struct option long_options[] = {
            {"vcf", required_argument, 0, 'v'},
            {"Paragraph", required_argument, 0, '1'},
            {"GraphTyper2", required_argument, 0, '2'},
            {"BayesTyper", required_argument, 0, '3'},
            {"VG-MAP", required_argument, 0, '4'},
            {"VG-Giraffe", required_argument, 0, '5'},
            {"GraphAligner", required_argument, 0, '6'},
            {"PanGenie", required_argument, 0, '7'},
            {"name", required_argument, 0, 'n'},

            {"mode", required_argument, 0, 'm'},
            {"out", required_argument, 0, 'o'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "v:1:2:3:4:5:6:7:n:m:o:h", long_options, &option_index);

        if (c == -1)
            break;
        
        switch (c)
        {
        case 'v':
            trueVcf = optarg;
            break;
        case '1':
            ParagraphVcf = optarg;
            break;
        case '2':
            GraphTyper2Vcf = optarg;
            break;
        case '3':
            BayesTyperVcf = optarg;
            break;
        case '4':
            mapVcf = optarg;
            break;
        case '5':
            giraffeVcf = optarg;
            break;
        case '6':
            GraphAlignerVcf = optarg;
            break;
        case '7':
            PanGenieVcf = optarg;
            break;
        case 'n':
            sample_name = optarg;
            break;
        case 'm':
            mode = optarg;
            break;
        case 'o':
            outputFileName = optarg;
            break;
        case 'h':
        case '?':
            help_merge(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc <= 2) {
        help_merge(argv);
        return 1;
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running.\n";

    // 检查参数
    if (mode != "specific" && mode != "all")
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Mode error: specific/all." << endl;
        help_merge(argv);
        exit(1);
    }
    if (sample_name.size() == 0)
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Empty sample name: -n." << endl;
        help_merge(argv);
        exit(1);
    }

    MERGE::baseVcfStruct mergeVcfStruct;

    // 运行程序构建索引和合并
    int softwareNum = MERGE::run_index_merge(
        trueVcf, 
        ParagraphVcf, 
        GraphTyper2Vcf, 
        BayesTyperVcf, 
        mapVcf, 
        giraffeVcf, 
        GraphAlignerVcf, 
        PanGenieVcf, 
        mergeVcfStruct, 
        sample_name
    );

    // 对合并的vcf进行过滤
    map<string, map<int, string > > outMap = MERGE::vcf_merge_filter(
        mergeVcfStruct, 
        softwareNum, 
        mode
    );

    // 保存结果
    MERGE::result_save(
        mergeVcfStruct.headInfo, 
        outMap, 
        outputFileName
    );

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done.\n";

    return 0;
}

// 帮助文档
void help_merge(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -v [options]" << endl
       << "merge the output of genomes graph software" << endl
       << endl
       << "required arguments (Note: vcf files must be sorted):" << endl
       << "    -v, --vcf            FILE      base vcf file for genotyping" << endl
       << "    --Paragraph          FILE      output file of Paragraph" << endl
       << "    --GraphTyper2        FILE      output file of GraphTyper2" << endl
       << "    --BayesTyper         FILE      output file of BayesTyper" << endl
       << "    --VG-MAP             FILE      output file of VG-MAP" << endl
       << "    --VG-Giraffe         FILE      output file of VG-Giraffe" << endl
       << "    --GraphAligner       FILE      output file of GraphAligner" << endl
       << "    --PanGenie           FILE      output file of PanGenie" << endl
       << "    -n, --name           STRING    used to extract genotyping result, after the FORMAT of vcf" << endl
       << endl
       << "optional arguments:" << endl
       << "    -m, --mode           STRING    software mode (specific/all) [specific]" << endl
       << "    -o, --out            FILE      output filename [stdout]" << endl
       << endl
       << "    -h, --help                     print this help document" << endl;
}