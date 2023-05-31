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

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running.\n";
    
    // 统计结果
    VCFCOUNT::count(vcfFileName);
    
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done.\n";

    return 0;
}

// 帮助文档
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


/**
 * @brief 统计vcf文件
 * 
 * @param vcfFileName  输入vcf文件
 * 
 * @return 0
**/
int VCFCOUNT::count(
    const string & vcfFileName
)
{
    // 记录各种变异的数量
    uint32_t snpNum = 0;
    uint32_t indelNum = 0;
    uint32_t insNum = 0;
    uint32_t delNum = 0;
    // uint32_t invNum = 0;
    // uint32_t dupNum = 0;
    uint32_t otherNum = 0;

    // 存储vcf信息
    VCFINFOSTRUCT INFOSTRUCTTMP;

    // 初始化
    VCFOPEN VCFOPENCLASS;
    VCFOPENCLASS.init(
        vcfFileName
    );

    // 打开vcf文件
    VCFOPENCLASS.open();

    // 如果没有遍历完，继续
    while (VCFOPENCLASS.read(INFOSTRUCTTMP))
    {
        // 有长度信息再计算
        if (INFOSTRUCTTMP.POS > 0)
        {
            // 遍历qry序列列表
            for (auto qrySeq : INFOSTRUCTTMP.ALTVec)
            {
                uint32_t qryLen = qrySeq.size();

                int32_t svLen = qryLen -INFOSTRUCTTMP.LEN;

                if (svLen == 0)
                {
                    snpNum++;
                }
                else if (svLen <= 49 && svLen >= -49)
                {
                    indelNum++;
                }
                else if (svLen < -49)
                {
                    delNum++;
                }
                else if (svLen > 49)
                {
                    insNum++;
                }
                else
                {
                    otherNum++;
                }
            }
        }
    }

    // 关闭vcf文件
    VCFOPENCLASS.close();

    // 打印结果
    cout << "SNP\tInDels\tDeletion\tInsertion\tOther\n";
    cout << snpNum << "\t" << indelNum << "\t" << delNum << "\t" << insNum << "\t" << otherNum << endl;
    
    return 0;
}