// g++ vcf_split.cpp -o vcf_split -lz -O2
#include <getopt.h>
#include <cstdlib>

#include "../include/vcf_split.hpp"


using namespace std;

int main_split(int argc, char** argv)
{
    // 输入文件
    string vcfFileName;

    // 输入提取snp和indel位置的vcf文件
    string baseVcfFileName;
    
    // 拆分的模式
    string mode = "type";

    // snp和indel离sv的距离
    int length = 100;

    // 输出文件前缀
    string prefix = "split";

    // 是否根据FILTER列进行过滤
    bool filterBool = false;

    // 输入参数
    int c;
    while (true)
    {
        static const struct option long_options[] = {
            {"vcf", required_argument, 0, 'v'},

            {"mode", required_argument, 0, 'm'},
            {"length", required_argument, 0, 'l'},
            {"basevcf", required_argument, 0, 'V'},
            {"prefix", required_argument, 0, 'p'},
            {"filter", no_argument, 0, 'f'},

            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "v:m:l:V:p:fh", long_options, &option_index);

        if (c == -1)
            break;

        switch (c)
        {
        case 'v':
            vcfFileName = optarg;
            break;
        case 'm':
            mode = optarg;
            break;
        case 'l':
            length = stoi(optarg);
            break;
        case 'V':
            baseVcfFileName = optarg;
            break;
        case 'p':
            prefix = optarg;
            break;;
        case 'f':
            filterBool = true;
            break;
        case 'h':
        case '?':
            help_split(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    // 检查参数是否正确
    if (vcfFileName.empty() || (mode != "type" && mode != "number"))
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error.\n";
        help_split(argv);
        return 1;
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running ...\n";

    // init
    VCFSplit VCFSplitClass(vcfFileName, baseVcfFileName, length, filterBool, prefix);
    
    if (mode == "type")
    {
        VCFSplitClass.vcf_split_type();
    }
    else if (mode == "number")
    {
        VCFSplitClass.vcf_split_number();
    }
    else
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error.\n";
        help_split(argv);
        return 1;
    }
    
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done ...\n";

    return 0;
}

// 帮助文档
void help_split(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -v [options]" << endl
       << "split VCF files by type or number of nearby SNP+Indels." << endl
       << endl
       << "required arguments (Note: vcf files must be sorted):" << endl
       << "    -v, --vcf           FILE       vcf file to be converted" << endl
       << endl
	   << "optional arguments:" << endl
       << "    -m, --mode          STRING     split standard: (type -> SNP+Indel+Ins+Del / number -> number of nearby SNP+Indel) [type]" << endl
       << "    -l, --length        INT        distance of snp and indel from sv [100]" << endl      
       << "    -V, --basevcf       INT        vcf file used to extract snp and indel locations [same as -v]" << endl      
       << "    -p, --prefix        STRING     output prefix [split]" << endl
       << "    -f, --filter        BOOL       filter using the FILTER column" << endl
       << endl
       << "    -h, --help                     print this help document" << endl;
}



/**
 * @brief Convert vcf to the format required by graph genome tools
 * 
 * @param vcfFileName      input VCF file name
 * @param baseVcfFileName  VCF file for extracting snp and indel locations
 * @param length           Count the number of SNPs and indels within the 'length' near the breakpoint
 * @param filterBool       Whether to filter according to the 'FILTER' column
 * @param prefix           output file prefix
 * 
**/
VCFSplit::VCFSplit(
    string vcfFileName,
    string baseVcfFileName,
    int length, 
    bool filterBool,
    string prefix
) : vcfFileName_(vcfFileName), baseVcfFileName_(baseVcfFileName), length_(length), filterBool_(filterBool), prefix_(prefix) {}


/**
 * @brief Get the length information of variant
 * 
 * @param INFOSTRUCTTMP   variant information
 * @param VCFOPENCLASS    VCFOPEN
 * @param refLen          store the length of REF
 * @param qryLen          store the length of ALT
 * 
**/
bool VCFSplit::get_len(
    const VCFINFOSTRUCT& INFOSTRUCTTMP, 
    VCFOPEN& VCFOPENCLASS, 
    uint32_t& refLen, 
    uint32_t& qryLen
)
{
    // 最后一列信息
    vector<string> lastColVec = split(INFOSTRUCTTMP.lineVec[INFOSTRUCTTMP.lineVec.size()-1], ":");

    // 找FORMAT字段中gt的位置
    vector<string> formatVec = split(strip(INFOSTRUCTTMP.lineVec[8], '\n'), ":");
    vector<string>::iterator gtItera = find(formatVec.begin(), formatVec.end(), "GT");

    string gt;
    int gtIndex = 0;
    if (gtItera != formatVec.end()) {  // FORMAT中有GT
        // GT的index
        gtIndex = distance(formatVec.begin(), gtItera);
        gt = strip(lastColVec[gtIndex], '\n');
    } else {  // 没有的话退出代码
        cerr << "[" << __func__ << "::" << getTime() << "] "
            << "Error: GT not in FORMAT column -> " 
            << INFOSTRUCTTMP.line << endl;
        exit(1);
    }

    if (gt == "0/0" || gt == "." || gt == "./." || gt == ".|." || gt == "0|0") {  // 如果是(0/0, .)格式的则直接跳过，不计数。
        return false;
    }

    // 根据FILTER字段进行过滤
    if (INFOSTRUCTTMP.FILTER != "PASS" && filterBool_) {  // 如果不是(PASS)则直接跳过，如果不过滤，跳过该判断
        return false;
    }

    // 软件的分型结果
    vector<string> evaluate_gt;
    if (filterBool_) {
        evaluate_gt = VCFOPENCLASS.gt_split(gt);
    }
    else {  // 如果不过滤，则这里也不判断，因为gramtools的GT格式是1，没有分隔符
        evaluate_gt = {"1", "1"};
    }
    
    if (INFOSTRUCTTMP.ALT.find(">") != string::npos) {  // graphtyper软件的结果
        // 使用正则表达式提取长度信息
        std::regex reg("SVSIZE=(\\d+)");
        std::smatch match;
        // 检查qry序列中有没有长度信息，没有的话跳过该位点
        if (!std::regex_search(INFOSTRUCTTMP.ALT, match, reg)) {  // <INS:SVSIZE=97:BREAKPOINT1>
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: No length information in the fourth column, skip this site -> " << INFOSTRUCTTMP.line << endl; 
            return false;
        }

        if (INFOSTRUCTTMP.ALT.find("<INS") != string::npos) {  // 字符段中包含插入的字段（graphtyper）<INS:SVSIZE=97:BREAKPOINT1>
            refLen = INFOSTRUCTTMP.LEN;
            qryLen = std::stoul(match[1].str());
        } else if (INFOSTRUCTTMP.ALT.find("<DEL") != string::npos) {  // 字符段中包含缺失的字段（graphtyper） <DEL:SVSIZE=233:COVERAGE>
            refLen = std::stoul(match[1].str());
            qryLen = INFOSTRUCTTMP.LEN;
        } else if (INFOSTRUCTTMP.ALT.find("<DUP") != string::npos) {  // 字符段中包含重复的字段（graphtyper）<DUP:SVSIZE=2806:BREAKPOINT1>
            refLen = std::stoul(match[1].str());
            qryLen = refLen * 2;
        } else {  // 不认识的字段
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: unknown string -> ref_seq:" << INFOSTRUCTTMP.REF << " qey_seq:" << INFOSTRUCTTMP.ALT << endl;
            refLen = INFOSTRUCTTMP.LEN;

            qryLen = std::stoul(match[1].str());
        }
    } else {  // 正常的变异
        refLen = INFOSTRUCTTMP.LEN;

        string qrySeq;
        
        uint32_t maxGtNum = stoi(*max_element(evaluate_gt.begin(), evaluate_gt.end()));
        if (maxGtNum > INFOSTRUCTTMP.ALTVec.size()) {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: genotyping and ALT sequence numbers do not match -> " << INFOSTRUCTTMP.line << endl;
            exit(1);
        } else {
            qrySeq = INFOSTRUCTTMP.ALTVec[maxGtNum-1]; // GT号减去1是序列的索引
        }

        qryLen = qrySeq.size();
    }

    // 判断bayestyper的结果为duplication，且ref_len为1-2，qry_seq却很长时
    if ((INFOSTRUCTTMP.ID.find("Duplication") != string::npos) && (refLen <= 2)) {  // bayestyper会把duplication变成插入，所以重新计算长度
        refLen = qryLen;
        qryLen *= 2;
    }

    return true;
}


void VCFSplit::vcf_split_type()
{
    // 输入文件流
    // 存储vcf信息
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS(vcfFileName_);

    // 记录每种变异和注释行
    string simpleSVs, invSVs, dupSVs, otherSVs;

    // 输出文件流
    SAVE SimpleSave(prefix_ + ".snp.indel.ins.del.vcf.gz");
    SAVE InvSave(prefix_ + ".inv.vcf.gz");
    SAVE DupSave(prefix_ + ".dup.vcf.gz");
    SAVE OtherSave(prefix_ + ".other.vcf.gz");

    // If not traversed, continue
    while (VCFOPENCLASS.read(INFOSTRUCTTMP)) {
        // empty line, skip
        if (INFOSTRUCTTMP.line.empty()) {
            continue;
        }

        if (INFOSTRUCTTMP.line.find("#") != string::npos) {  // 检查是否是注释行
            string headLine = INFOSTRUCTTMP.line + "\n";
            SimpleSave.save(headLine);
            InvSave.save(headLine);
            DupSave.save(headLine);
            OtherSave.save(headLine);
        } else {
            // Get the length information of variant. If encountering an unrecognized string, continue to the next loop.
            uint32_t refLen, qryLen;
            if (!get_len(INFOSTRUCTTMP, VCFOPENCLASS, refLen, qryLen)) {
                continue;
            }

            // 判断变异类型，并保存到字符串中
            if ((refLen < 50 && qryLen < 50) || (refLen >= 50 && qryLen < 50) || (refLen < 50 && qryLen >= 50)) {  // snp+indel+del+ins
                simpleSVs += INFOSTRUCTTMP.line + "\n";
            } else if (refLen >= 50 && qryLen >= 50) {
                if (refLen*2 <= qryLen+10 && refLen*2 >= qryLen-10) {  // dup
                    dupSVs += INFOSTRUCTTMP.line + "\n";
                } else if (refLen <= qryLen+10 && refLen >= qryLen-10) {  // inv
                    invSVs += INFOSTRUCTTMP.line + "\n";
                } else {
                    otherSVs += INFOSTRUCTTMP.line + "\n";
                }
            } else {
                otherSVs += INFOSTRUCTTMP.line + "\n";
            }

            // save result
            if (simpleSVs.size() > 10 * 1024 * 1024) {  // 每10Mb写入一次
                SimpleSave.save(simpleSVs);
                InvSave.save(invSVs);
                DupSave.save(dupSVs);
                OtherSave.save(otherSVs);

                // 清空字符串
                simpleSVs.clear();
                invSVs.clear();
                dupSVs.clear();
                otherSVs.clear();
            }
        }
    }

    // 最后写入一次
    SimpleSave.save(simpleSVs);
    InvSave.save(invSVs);
    DupSave.save(dupSVs);
    OtherSave.save(otherSVs);

    // 清空字符串
    simpleSVs.clear();
    invSVs.clear();
    dupSVs.clear();
    otherSVs.clear();
}


// Classified by the number of nearby SNPs and indels
void VCFSplit::vcf_split_number()
{
    if (baseVcfFileName_.empty()) {
        baseVcfFileName_ = vcfFileName_;
    }

    // 记录sv拆分的结果
    string SVs0, SVs1, SVs2, SVs4, SVs6, SVs8, SVsMore;

    // 输出文件流
    SAVE SVs0Save(prefix_ + ".0.vcf.gz");
    SAVE SVs1Save(prefix_ + ".1.vcf.gz");
    SAVE SVs2Save(prefix_ + ".2.vcf.gz");
    SAVE SVs4Save(prefix_ + ".4.vcf.gz");
    SAVE SVs6Save(prefix_ + ".6.vcf.gz");
    SAVE SVs8Save(prefix_ + ".8.vcf.gz");
    SAVE SVsMoreSave(prefix_ + ".more.vcf.gz");

    // snp+indel的索引
    VCFINFOSTRUCT INFOSTRUCTBase;
    VCFOPEN VCFOPENCLASSBase(baseVcfFileName_);

    // Record the location information of snp and indel
    unordered_map<string, snpIndelInfo> snpIndelInfoMap;  // unordered_map<chromosome, snpIndelInfo>

    // If not traversed, continue
    while (VCFOPENCLASSBase.read(INFOSTRUCTBase)) {
        // empty line, skip
        if (INFOSTRUCTBase.line.empty()) {
            continue;
        }

        if (INFOSTRUCTBase.line.find("#") == string::npos) {
            // Get the length information of variant. If encountering an unrecognized string, continue to the next loop.
            uint32_t refLen, qryLen;
            if (!get_len(INFOSTRUCTBase, VCFOPENCLASSBase, refLen, qryLen)) {
                continue;
            }

            uint32_t refEnd;
            refEnd = INFOSTRUCTBase.POS + refLen - 1;

            // 判断变异类型，并保存到哈希表中
            if ((refLen < 50 && qryLen < 50)) {  // 保存snp+indel的起始和终止位置
                snpIndelInfoMap[INFOSTRUCTBase.CHROM].refStartVec.push_back(INFOSTRUCTBase.POS);
                snpIndelInfoMap[INFOSTRUCTBase.CHROM].refEndVec.push_back(refEnd);
            }
        }
    }


    // 输入文件流
    // 存储vcf信息
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS(vcfFileName_);

    // 构建哈希表 （结构变异）
    map<string, map<uint32_t, svInfo> > svInfoMap; // map<chromosome, map<refStart, svInfo>>

    // If not traversed, continue
    while (VCFOPENCLASS.read(INFOSTRUCTTMP)) {
        // empty line, skip
        if (INFOSTRUCTTMP.line.empty()) {
            continue;
        }

        if (INFOSTRUCTTMP.line.find("#") != string::npos) {  // 检查是否是注释行
            string headLine = INFOSTRUCTTMP.line + "\n";
            SVs0Save.save(headLine);
            SVs1Save.save(headLine);
            SVs2Save.save(headLine);
            SVs4Save.save(headLine);
            SVs6Save.save(headLine);
            SVs8Save.save(headLine);
            SVsMoreSave.save(headLine);
        } else {
            // Get the length information of variant. If encountering an unrecognized string, continue to the next loop.
            uint32_t refLen, qryLen;
            if (!get_len(INFOSTRUCTTMP, VCFOPENCLASS, refLen, qryLen)) {
                continue;
            }

            uint32_t refEnd;
            refEnd = INFOSTRUCTTMP.POS + refLen - 1;

            // 判断变异类型，并保存到哈希表中
            if (refLen >= 50 || qryLen >= 50) {  // 保存sv的终止和变异信息
                svInfoMap[INFOSTRUCTTMP.CHROM][INFOSTRUCTTMP.POS].information = INFOSTRUCTTMP.line;
                svInfoMap[INFOSTRUCTTMP.CHROM][INFOSTRUCTTMP.POS].refEnd = refEnd;
            }
        }
    }


    // 找sv附近的snp和indel数量并进行拆分
    // 对sv哈希表进行循环
    for (const auto& [chromosome, refStartInfoMap] : svInfoMap) {
        // snp+indel起始和终止的列表
        const vector<uint32_t>& snpIndelRefStartVec = snpIndelInfoMap[chromosome].refStartVec;
        const vector<uint32_t>& snpIndelRefEndVec = snpIndelInfoMap[chromosome].refEndVec;

        for (const auto& [refStart, Info] : refStartInfoMap) {
            // 二分查找法找起始和终止200bp内的snp和indel数量
            // 寻找左右索引
            int64_t startLeftIdxTmp = search_Binary_right(snpIndelRefEndVec, refStart - length_);
            int64_t startRightIdxTmp = search_Binary_left(snpIndelRefEndVec, refStart);
            int64_t endLeftIdxTmp = search_Binary_right(snpIndelRefStartVec, Info.refEnd);
            int64_t endRightIdxTmp = search_Binary_left(snpIndelRefStartVec, Info.refEnd + length_);

            // sv两端snp+indel的数量
            int64_t leftSnpIndelNum = startRightIdxTmp - startLeftIdxTmp;
            if (leftSnpIndelNum < 0) {  // 没有
                leftSnpIndelNum = 0;
            } else {  // 有一个或多个
                leftSnpIndelNum++;
            }
            
            int64_t rightSnpIndelNum = endRightIdxTmp - endLeftIdxTmp;
            if (rightSnpIndelNum < 0) {  // 没有
                rightSnpIndelNum = 0;
            } else {  // 有一个或多个
                rightSnpIndelNum++;
            }

            int64_t snpIndelNum = leftSnpIndelNum + rightSnpIndelNum;
            if (snpIndelNum == 0) {
                SVs0 += Info.information + "\n";
            } else if (snpIndelNum == 1) {
                SVs1 += Info.information + "\n";
            } else if (snpIndelNum == 2) {
                SVs2 += Info.information + "\n";
            } else if (snpIndelNum == 4) {
                SVs4 += Info.information + "\n";
            } else if (snpIndelNum == 6) {
                SVs6 += Info.information + "\n";
            } else if (snpIndelNum == 8) {
                SVs8 += Info.information + "\n";
            } else {
                SVsMore += Info.information + "\n";
            }
            
            // save result
            if (SVs0.size() > 10 * 1024 * 1024) // 每10Mb写入一次
            {
                SVs0Save.save(SVs0);
                SVs1Save.save(SVs1);
                SVs2Save.save(SVs2);
                SVs4Save.save(SVs4);
                SVs6Save.save(SVs6);
                SVs8Save.save(SVs8);
                SVsMoreSave.save(SVsMore);

                // 清空字符串
                SVs0.clear();
                SVs1.clear();
                SVs2.clear();
                SVs4.clear();
                SVs6.clear();
                SVs8.clear();
                SVsMore.clear();
            }
        }
    }
    
    SVs0Save.save(SVs0);
    SVs1Save.save(SVs1);
    SVs2Save.save(SVs2);
    SVs4Save.save(SVs4);
    SVs6Save.save(SVs6);
    SVs8Save.save(SVs8);
    SVsMoreSave.save(SVsMore);

    // 清空字符串
    SVs0.clear();
    SVs1.clear();
    SVs2.clear();
    SVs4.clear();
    SVs6.clear();
    SVs8.clear();
    SVsMore.clear();
}