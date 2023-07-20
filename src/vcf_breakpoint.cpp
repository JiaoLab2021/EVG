// g++ vcf_breakpoint.cpp -o vcf_breakpoint -lz -O2
#include "../include/vcf_breakpoint.hpp"

using namespace std;

int main_breakpoint(int argc, char** argv)
{
    // 输入文件
    string vcfFileName;
    
    // 输出文件前缀
    string prefix = "breakpoint";

    // 断点误差大小
    int breakpointErrorSize_ = 1;

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
            breakpointErrorSize_ = stoi(optarg);
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
    
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running ...\n";

    // init
    VCFBreakpoint VCFBreakpointClass(vcfFileName, prefix, breakpointErrorSize_);
    // set breakpoint
    VCFBreakpointClass.vcf_breakpoint();
    
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done ...\n";

    return 0;
}

// 帮助文档
void help_breakpoint(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -v [options]" << endl
       << "set error for variants breakpoints." << endl
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


/**
 * init
 *
 * @param vcfFileName                     tinput VCF  file name
 * @param prefix                          output file prefix
 * @param breakpointErrorSize             breakpoint error
 * 
**/
VCFBreakpoint::VCFBreakpoint(
    const string & vcfFileName, 
    const string & prefix, 
    const int & breakpointErrorSize
) : vcfFileName_(vcfFileName), prefix_(prefix), breakpointErrorSize_(breakpointErrorSize) {}


void VCFBreakpoint::vcf_breakpoint()
{
    // input file stream
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS(vcfFileName_);

    // 记录转换后的结果
    string outTxt;

    // ouyput file stream
    SAVE outFile(prefix_ + ".vcf.gz");

    while(VCFOPENCLASS.read(INFOSTRUCTTMP)) {
        if (INFOSTRUCTTMP.line.find("#") != string::npos) {  // 检查是否是注释行
            string headLine = INFOSTRUCTTMP.line + "\n";
            outFile.save(headLine);
        } else {
            // 使用正则表达式提取长度信息
            std::regex reg("END=(\\d+)");
            
            // 添加误差
            INFOSTRUCTTMP.lineVec[1] = (INFOSTRUCTTMP.POS > breakpointErrorSize_) ? to_string(INFOSTRUCTTMP.POS - breakpointErrorSize_) : "1";

            // end position
            uint32_t refEnd = stoul(INFOSTRUCTTMP.lineVec[1]) + INFOSTRUCTTMP.LEN - 1;
            INFOSTRUCTTMP.lineVec[7] = regex_replace(string(INFOSTRUCTTMP.lineVec[7]), regex(reg), string("END=" + to_string(refEnd)));

            // 添加到输出字符串中
            outTxt += join(INFOSTRUCTTMP.lineVec, "\t") + "\n";

            // save result
            if (outTxt.size() > 10 * 1024 * 1024) {  // 每10Mb写入一次
                outFile.save(outTxt);

                // 清空字符串
                outTxt.clear();
            }
        }
    }

    // 最后写入一次
    outFile.save(outTxt);

    // 清空字符串
    string().swap(outTxt);
}