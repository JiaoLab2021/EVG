// g++ vcf_split.cpp -o vcf_split -lz -O2
#include "../include/vcf_filter.hpp"

using namespace std;

int main_filter(int argc, char** argv)
{
    // 输入文件
    string vcfFileName;

    // 过滤阈值
    double MAF = 0.01;  // 最小等位基因频率
    double MISSRATE = 0.1;  // 缺失率

    // 输出文件名
    string outputFileName = "";

    // 输入参数
    int c;
    while (true)
    {
        static const struct option long_options[] = {
            {"vcf", required_argument, 0, 'v'},

            {"maf", required_argument, 0, '1'},
            {"geno", required_argument, 0, '2'},
            {"out", required_argument, 0, 'o'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "v:1:2:o:h", long_options, &option_index);

        if (c == -1)
            break;

        switch (c)
        {
        case 'v':
            vcfFileName = optarg;
            break;
        case '1':
            MAF = stod(optarg);
            break;
        case '2':
            MISSRATE = stod(optarg);
            break;
        case 'o':
            outputFileName = optarg;
            break;
        case 'h':
        case '?':
            help_filter(argv);
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
        help_filter(argv);
        return 1;
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running.\n";
    
    // 统计结果
    VCFFILTER::vcf_filter(
        vcfFileName, 
        outputFileName, 
        MAF, 
        MISSRATE
    );
    
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done.\n";

    return 0;
}

// 帮助文档
void help_filter(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -v [options]" << endl
       << "filter SNPs by maf and missing rate" << endl
       << endl
       << "required arguments:" << endl
       << "    -v, --vcf        FILE      vcf file to be converted" << endl
       << endl
       << "optional arguments:" << endl
       << "    --maf            FLOAT     exclude variants with minor allele frequency lower than threshold [0.01]" << endl
       << "    --geno           FLOAT     exclude variants with missing call frequencies greater than threshold [0.1]" << endl
       << "    -o, --out        FILE      output filename [stdout]" << endl
       << endl
       << "    -h, --help                 print this help document" << endl;
}


/**
 * @brief 根据次等位基因频率和缺失率过滤SNPs
 * 
 * @param vcfFileName    输入vcf文件
 * @param outputFileName 输出文件名
 * @param MAF            次等位基因频率
 * @param MISSRATE       缺失率
 * 
 * @return 0
**/
int VCFFILTER::vcf_filter(
    const string & vcfFileName, 
    const string & outputFileName, 
    const double & MAF, 
    const double & MISSRATE
)
{
    // 输入文件流
    // 存储vcf信息
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS;
    VCFOPENCLASS.init(
        vcfFileName
    );
    // 打开vcf文件
    VCFOPENCLASS.open();


    // 输出文件流
    SAVE::SAVE SAVECLASS;
    SAVECLASS.init(
        outputFileName
    );
    SAVECLASS.open();


    // 临时存储输出字符串
    string outTxt = "";

    // 如果没有遍历完，继续
    while (VCFOPENCLASS.read(INFOSTRUCTTMP))
    {
        // 如果注释行，直接保存
        if (INFOSTRUCTTMP.INFO.find("#") != string::npos)
        {
            outTxt += INFOSTRUCTTMP.INFO + "\n";
            continue;
        }

        // 获取变异类型
        INFOSTRUCTTMP.TYPE = VCFOPENCLASS.get_TYPE(
            INFOSTRUCTTMP.LEN, 
            INFOSTRUCTTMP.ALTVec
        );

        if (INFOSTRUCTTMP.TYPE == "SNP")  // SNP时再判断是否过滤
        {
            double MAFTMP;  // 最小等位基因频率
            double MISSRATETMP;  // 缺失率

            // 获取所有的基因型   map<idx, vector<gtString>>
            map<int, vector<string> > GTVecMapTmp = VCFOPENCLASS.get_gt(
                INFOSTRUCTTMP.INFOVec
            );

            // 如果只有一个基因型，跳过。
            if (GTVecMapTmp.size() <= 1)
            {
                continue;
            }
            
            // 计算最小等位基因频率和缺失率
            tie(MAFTMP, MISSRATETMP) = VCFOPENCLASS.calculate(
                GTVecMapTmp, 
                INFOSTRUCTTMP.INFOVec.size() - 9
            );

            // 通过阈值了直接保存
            if (MAFTMP >= MAF && MISSRATETMP <= MISSRATE)
            {
                outTxt += INFOSTRUCTTMP.INFO + "\n";
            }
        }
        else  // 其它类型的变异直接保存
        {
            outTxt += INFOSTRUCTTMP.INFO + "\n";
        }

        if (outTxt.size() >= 10000000)  // 每10m写一次
        {
            // 输出文件流
            SAVECLASS.save(
                outTxt
            );

            // 清空
            outTxt.clear();
            string().swap(outTxt);
        }
    }

    if (outTxt.size() >= 0)  // 最后写一次
    {
        // 其它类型的变异直接保存
        SAVECLASS.save(
            outTxt
        );

        // 清空
        outTxt.clear();
        string().swap(outTxt);
    }

    // 关闭文件
    VCFOPENCLASS.close();
    SAVECLASS.close();
    
    return 0;
}