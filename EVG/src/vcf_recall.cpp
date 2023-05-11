// g++ vcf_recall.cpp -o vcf_recall -lz
#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <getopt.h>
#include <iomanip>
#include <cstdlib>
#include "../include/get_time.hpp"
#include "../include/vcf_recall.hpp"

using namespace std;

void help_recall(char** argv);

int main_recall(int argc, char** argv)
{
    string genotype_filename;
    string evaluateFilename;
    string model;

    // roc计算规则
    string rocKey = "FORMAT.DP";
    string rocType = "genotype";

    // 输入参数
    int c;
    while (true)
    {
        static const struct option long_options[] = {
            {"vcf1", required_argument, 0, 't'},
            {"vcf2", required_argument, 0, 'e'},         
            {"model", required_argument, 0, 'm'},

            {"score-field", required_argument, 0, 's'},
            {"roc", required_argument, 0, 'r'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "t:e:m:s:r:h", long_options, &option_index);

        if (c == -1)
            break;

        switch (c)
        {
        case 't':
            genotype_filename = optarg;
            break;
        case 'e':
            evaluateFilename = optarg;
            break;
        case 'm':
            model = optarg;
            break;
        case 's':
            rocKey = optarg;
            break;
        case 'r':
            rocType = optarg;
            break;
        case 'h':
        case '?':
            help_recall(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    // 检查参数是否正确
    if (argc <= 6 || genotype_filename.empty() || evaluateFilename.empty() || model.empty())
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
             << "Parameter error.\n";
        help_recall(argv);
        return 1;
    }

    // 检查rocType是否正确
    if (rocType != "recall" && rocType != "genotype")
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
             << "Parameter error: -r\n";
        help_recall(argv);
        return 1;
    }
    
    cerr << "[" << __func__ << "::" << getTime() << "] " 
         << "Running.\n";

    // 命名空间
    
    // 修改rocKey and rocType全局变量
    RECALL::change_roc(
        rocKey, 
        rocType
    );

    // 构建索引
    RECALL::vcfStructure trueVcfStructure = RECALL::build_index(
        genotype_filename
    );

    // 评估
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Start evaluation.\n";

    if (model == "vg" || model == "recall" || model == "GraphAligner")
    {
        RECALL::evulate_recall(
            evaluateFilename, 
            trueVcfStructure
        );
    }
    else if (model == "bayestyper" || model == "graphtyper" || model == "paragraph" || model == "genotype")
    {
        RECALL::evulate_gt(
            evaluateFilename, 
            trueVcfStructure.chrStartLenInfoGtVecMap, 
            trueVcfStructure.allLengthList
        );
    }
    else if (model == "gramtools")
    {
        string outFileName = RECALL::gramtools_convert(evaluateFilename);
        RECALL::evulate_gt(
            outFileName, 
            trueVcfStructure.chrStartLenInfoGtVecMap, 
            trueVcfStructure.allLengthList
        );
    }
    else
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "The model is incorrect. [vg/paragraph/graphtyper/bayestyper/GraphAligner/gramtools/pangenie]" << endl;
        exit(1);
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " 
         << "Done.\n";

    return 0;
}

// 帮助文档
void help_recall(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -t -e -m [options]" << endl
       << "assessing the results of genomes graph software" << endl
       << endl
       << "required arguments (Note: vcf files must be sorted):" << endl
       << "    -t, --vcf1          FILE       vcf file containing truth values" << endl
       << "    -e, --vcf2          FILE       output vcf file by genomes graph software" << endl
       << "    -m, --model         STRING     the mode in which the software runs (vg/paragraph/graphtyper/bayestyper/GraphAligner/gramtools/pangenie)" << endl
       << endl
	   << "optional arguments:" << endl
       << "    -s, --score-field   STRING     the name of the VCF FORMAT/INFO field to use as the ROC score (FORMAT.<name>/INFO.<name>) [FORMAT.DP]" << endl
       << "    -r, --roc           STRING     roc calculation rules (recall/genotype) [genotype]" << endl
       << endl
       << "    -h, --help                     print this help document" << endl;
}