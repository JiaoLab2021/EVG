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

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running ...\n";
    
    // init
    VCFFilter VCFFilterClass(
        vcfFileName, 
        outputFileName, 
        MAF, 
        MISSRATE
    );
    VCFFilterClass.vcf_filter();
    
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done ...\n";

    return 0;
}

// 帮助文档
void help_filter(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -v FILE [options]" << endl
       << "filter SNPs by maf and missing rate." << endl
       << endl
       << "required arguments:" << endl
       << "    -v, --vcf        FILE      vcf file to be converted" << endl
       << endl
       << "optional arguments:" << endl
       << "    --maf            FLOAT     exclude SNPs with minor allele frequency lower than threshold [0.01]" << endl
       << "    --geno           FLOAT     exclude SNPs with missing call frequencies greater than threshold [0.1]" << endl
       << "    -o, --out        FILE      output filename [stdout]" << endl
       << endl
       << "    -h, --help                 print this help document" << endl;
}


/**
 * init
 *
 * @param vcfFileName         input VCF  file name
 * @param outputFileName      output file name
 * @param MAF                 MAF
 * @param MISSRATE            MISSRATE
 * 
**/
VCFFilter::VCFFilter(
    const string & vcfFileName, 
    const string & outputFileName, 
    const double & MAF, 
    const double & MISSRATE
) : vcfFileName_(vcfFileName), outputFileName_(outputFileName), MAF_(MAF), MISSRATE_(MISSRATE) {}



void VCFFilter::vcf_filter()
{
    // input file stream
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS(vcfFileName_);


    // output file stream
    SAVE SAVECLASS(outputFileName_);


    // 临时存储输出字符串
    string outTxt = "";

    // 如果没有遍历完，继续
    while (VCFOPENCLASS.read(INFOSTRUCTTMP)) {
        // empty line, skip
        if (INFOSTRUCTTMP.line.empty()) {
            continue;
        }
        
        // 如果注释行，直接保存
        if (INFOSTRUCTTMP.line.find("#") != string::npos) {
            outTxt += INFOSTRUCTTMP.line + "\n";
            continue;
        }

        // 获取变异类型
        INFOSTRUCTTMP.ID = VCFOPENCLASS.get_TYPE(
            INFOSTRUCTTMP.LEN, 
            INFOSTRUCTTMP.ALTVec
        );

        if (INFOSTRUCTTMP.ID == "SNP") {  // SNP时再判断是否过滤
            double MAFTMP;  // 最小等位基因频率
            double MISSRATETMP;  // 缺失率

            // 获取所有的基因型   map<idx, vector<gtString>>
            map<int, vector<string> > GTVecMapTmp = VCFOPENCLASS.get_gt(
                INFOSTRUCTTMP.lineVec
            );

            // 如果只有一个基因型，跳过。
            if (GTVecMapTmp.size() <= 1) {
                continue;
            }
            
            // 计算最小等位基因频率和缺失率
            tie(MAFTMP, MISSRATETMP) = VCFOPENCLASS.calculate(
                GTVecMapTmp, 
                INFOSTRUCTTMP.lineVec.size() - 9
            );

            // 通过阈值了直接保存
            if (MAFTMP >= MAF_ && MISSRATETMP <= MISSRATE_) {
                outTxt += INFOSTRUCTTMP.line + "\n";
            }
        } else {  // 其它类型的变异直接保存
            outTxt += INFOSTRUCTTMP.line + "\n";
        }

        if (outTxt.size() >= 10 * 1024 * 1024) {  // 每10m写一次
            // 输出文件流
            SAVECLASS.save(
                outTxt
            );

            // 清空
            outTxt.clear();
        }
    }

    if (outTxt.size() > 0) {  // 最后写一次
        // 其它类型的变异直接保存
        SAVECLASS.save(
            outTxt
        );

        // 清空
        string().swap(outTxt);
    }
}