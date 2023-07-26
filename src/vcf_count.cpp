// g++ vcf_count.cpp -o vcf_count -lz -O2
#include "../include/vcf_count.hpp"


using namespace std;


int main_count(int argc, char** argv)
{
    // 输入文件
    string vcfFileName;

    // 输入参数
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

    // 检查参数是否正确
    if (vcfFileName.empty())
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error.\n";
        help_count(argv);
        return 1;
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running ...\n";
    
    // init
    VCFCount VCFCountClass(vcfFileName);
    // count
    VCFCountClass.count();
    
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done ...\n";

    return 0;
}

// 帮助文档
void help_count(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -v [options]" << endl
       << "calculate the number of alleles." << endl
       << endl
       << "required arguments:" << endl
       << "    -v, --vcf           FILE       vcf file to be converted" << endl
       << endl
       << "    -h, --help                     print this help document" << endl;
}


/**
 * init
 *
 * @param vcfFileName       input VCF  file name
 * 
**/
VCFCount::VCFCount(
    const string & vcfFileName
) : vcfFileName_(vcfFileName) {}


void VCFCount::count()
{
    // 记录各种变异的数量
    uint64_t snpNum = 0;
    uint64_t indelNum = 0;
    uint64_t insNum = 0;
    uint64_t delNum = 0;
    uint64_t invNum = 0;
    uint64_t dupNum = 0;
    uint64_t otherNum = 0;

    // 存储vcf信息
    VCFINFOSTRUCT INFOSTRUCTTMP;

    // 初始化
    VCFOPEN VCFOPENCLASS(vcfFileName_);

    // 如果没有遍历完，继续
    while (VCFOPENCLASS.read(INFOSTRUCTTMP)) {
        // empty or comment line, skip
        if (INFOSTRUCTTMP.line.empty() || INFOSTRUCTTMP.line.find("#") != string::npos) {
            continue;
        }
        
        // 遍历qry序列列表
        for (const auto& qrySeq : INFOSTRUCTTMP.ALTVec) {
            // Record the amount of variation
            uint32_t qryLen = qrySeq.size();
            int32_t svLen = qryLen - INFOSTRUCTTMP.LEN;
            double lengthRatio = qryLen / float(INFOSTRUCTTMP.LEN);

            if (svLen == 0 && INFOSTRUCTTMP.LEN == 1 && qryLen == 1) {
                snpNum++;
            } else if (svLen <= 49 && svLen >= -49 && INFOSTRUCTTMP.LEN <= 49 && qryLen <= 49) {
                indelNum++;
            } else if (svLen >= -2 && svLen <= 2 && INFOSTRUCTTMP.LEN > 49 && qryLen > 49 ) {
                invNum++;
            } else if (lengthRatio >= 1.8 && lengthRatio <= 2.2 && INFOSTRUCTTMP.LEN > 49 && qryLen > 49) {
                dupNum++;
            } else if (svLen < 0) {
                delNum++;
            } else if (svLen > 0) {
                insNum++;
            } else {
                otherNum++;
            }
        }
    }

    // 打印结果
    cout << "           - SNP: " << snpNum << endl
        << "           - InDels: " << indelNum << endl
        << "           - Insertion: " << insNum << endl
        << "           - Deletion: " << delNum << endl
        << "           - Inversion: " << invNum << endl
        << "           - Duplication: " << dupNum << endl
        << "           - Other: " << otherNum << endl;
}