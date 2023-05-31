// g++ vcf_convert.cpp -o vcf_convert -lz
#include <getopt.h>

#include "../include/vcf_convert.hpp"

using namespace std;


int main_convert(int argc, char** argv)
{
    // 输入文件和参数
    string referenceFilename;
    string vcfFilename;
    string outFilename = "vcfConvert.out.vcf.gz";
    int readLen = 350; // 测序的读长

    // 过滤阈值
    double MAF = 0.;  // 最小等位基因频率
    double MISSRATE = 1.0;  // 缺失率

    // 输入参数
    int c;
    while (true)
    {
        static const struct option long_options[] = {
            {"reference", required_argument, 0, 'r'},
            {"vcf", required_argument, 0, 'v'},

            {"length", required_argument, 0, 'l'},
            {"maf", required_argument, 0, '1'},
            {"geno", required_argument, 0, '2'},
            {"out", required_argument, 0, 'o'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "r:v:l:1:2:o:h", long_options, &option_index);

        if (c == -1)
            break;

        switch (c)
        {
        case 'r':
            referenceFilename = optarg;
            break;
        case 'v':
            vcfFilename = optarg;
            break;
        case 'l':
            readLen = stoi(optarg);
            break;
        case '1':
            MAF = stod(optarg);
            break;
        case '2':
            MISSRATE = stod(optarg);
            break;
        case 'o':
            outFilename = optarg;
            break;
        case 'h':
        case '?':
            help_convert(argv);
            exit(1);
            break;
        default:
            abort();
        }
    }

    if (argc <= 2) {
        help_convert(argv);
        return 1;
    }

    // 打印log
    cerr << "[" << __func__ << "::" << getTime() << "] " 
         << "Running.\n";

    // 构建reference的索引
    CONVERT::refIndexStruct refIndexS;
    CONVERT::build_reference_index(
        referenceFilename, 
        refIndexS
    );

    // 如果raed的长度大于1000，则不用read length过滤，因为不用跑paragraph
    if (readLen > 1000)
    {
        readLen = 0;
    }

    // vcf文件转换
    CONVERT::vcf_convert(
        vcfFilename, 
        readLen, 
        MAF, 
        MISSRATE, 
        refIndexS, 
        outFilename
    );

    // 打印log
    cerr << "[" << __func__ << "::" << getTime() << "] " 
         << "Done.\n";

    return 0;
}

// 帮助文档
void help_convert(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -r -v [options]" << endl
       << "convert vcf files merged by vcftools to the format required by genome graph." << endl
       << endl
       << "required arguments:" << endl
       << "    -r, --reference   FILE     input FASTA reference" << endl
       << "    -v, --vcf         FILE     input VCF" << endl
       << endl
	   << "optional arguments:" << endl
       << "    -l, --length      INT      read length [350]" << endl
       << "    --maf             FLOAT    exclude variants with minor allele frequency lower than threshold [0.0]" << endl
       << "    --geno            FLOAT    exclude variants with missing call frequencies greater than threshold [1.0]" << endl
       << "    -o, --out         FILE     output file name [vcfConvert.out.vcf.gz]" << endl
       << endl
       << "    -h, --help                 print this help document" << endl;
}


/*
    打开fasta文件
    inputFileName -> refgenome
    refIndexS -> 保存contig长度和序列信息
*/
void CONVERT::build_reference_index(string inputFileName, refIndexStruct & refIndexS)
{
    // 染色体名称输出，graphtyper需要
    ofstream outFile;
    outFile.open("CHROMOSOME.NAME", ios::out);
    if(!outFile.is_open())
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'CHROMOSOME.NAME': No such file or directory." 
                << endl;
        exit(1);
    }

    // 输入文件流
    gzFile gzfp = gzopen(inputFileName.c_str(), "rb");

    // 打开文件
    if(!gzfp)
    {
        cerr << "[" << __func__ << "::" << getTime() << "] "
                << "'"
                << inputFileName 
                << "': No such file or directory." 
                << endl;
        exit(1);
    }
    else
    {
        kseq_t *ks;
        ks = kseq_init(gzfp);
    
        while( kseq_read(ks) >= 0 )
        {
            string chromosome = ks->name.s;
            long long int chrLen = ks->seq.l;
            string sequence = ks->seq.s;

            refIndexS.chrLenTxt += "##contig=<ID=" + chromosome + ",length=" + to_string(chrLen) + ">\n";
            refIndexS.sequenceMap[chromosome] = sequence;

            // 输出染色体名称
            outFile << chromosome + "\n";
        }

        // 释放内存，关闭文件
        kseq_destroy(ks);
        gzclose(gzfp);
    }

    // 释放内存
    outFile.close();
}


/*
    对vcf文件进行替换
    vcfFilename -> 需要转换的vcf文件
    readLen -> 测序文件的read长度
    refIndexS -> contig长度和序列信息
    outFilename -> 输出文件名
    a. 表头加染色体长度
    b. 第一个变异要大于read length
    c. 序列必须与reference对应
    d. 对vcf的第八列 'END=' 进行替换
    e. refSeq的第一个碱基要和qrySeq的第一个碱基一样
    f. 检查ref和qry序列是否一样，一样的位点跳过
    g. 检查refSeq和qrySeq中是否含有atgcnATGCN外的字符，含有的话跳过该位点
    h. 检查替换后的qry有没有相同的，有的话跳过该位点 --> e.
    i. 检查是否有位置重复的变异
    j. 将基因型中的.转为.|.
    k. 将基因型中的/转为|
    l. 只保留基因型中的二倍体变异
    m. 检查GT是不是比qry的序列还多
*/
/**
 * @brief 将vcf转换为 graph genome tools 需要的格式
 * 
 * @param vcfFileName     输入vcf文件
 * @param readLen         读长
 * @param MAF            次等位基因频率
 * @param MISSRATE       缺失率
 * @param refIndexS       参考基因组长度和序列信息
 * @param outFilename     输出文件名
 * 
 * @return 0
**/
void CONVERT::vcf_convert(
    const string & vcfFilename, 
    const int & readLen, 
    const double & MAF, 
    const double & MISSRATE, 
    const refIndexStruct & refIndexS, 
    const string & outFilename
)
{
    // 检查文件是否排序
    check_vcf_sort(vcfFilename);

    // 基因组上的序列信息
    map<string, string> seqMap = refIndexS.sequenceMap;

    // 保存转换后的结果
    string outTxt;


    // 输出文件流
    SAVE::SAVE SAVECLASS;
    SAVECLASS.init(
        outFilename
    );
    SAVECLASS.open();


    // 记录上一个变异的start和染色体号
    string preChromosome;
    long long int preRefStart = 0;


    // 输入文件流
    // 存储vcf信息
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS;
    VCFOPENCLASS.init(
        vcfFilename
    );
    // 打开vcf文件
    VCFOPENCLASS.open();


    // 如果没有遍历完，继续
    while (VCFOPENCLASS.read(INFOSTRUCTTMP))
    {
        // 注释行
        if (INFOSTRUCTTMP.INFO.find("#") != string::npos)
        {
            // 如果vcf中含有染色体长度信息则跳过
            if (INFOSTRUCTTMP.INFO.find("#CHROM") != string::npos)
            {
                // 保存染色体长度信息
                outTxt += refIndexS.chrLenTxt;

                // 保存表头
                outTxt += INFOSTRUCTTMP.INFO + "\n";
            }
            else if (INFOSTRUCTTMP.INFO.find(",length") == string::npos)
            {
                outTxt += INFOSTRUCTTMP.INFO + "\n";
            }

            continue;
        }

        // 非注释行
        if (INFOSTRUCTTMP.INFOVec.size() < 9) // 先检查文件对不对，不对则跳出代码
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Error '"
                    << vcfFilename 
                    << "' -> The number of vcf columns is less than 9." 
                    << " (" << strip(INFOSTRUCTTMP.INFO, '\n') << ")"
                    << endl;
            exit(1);
        }



        /* ************************ filter SNPs by maf and missing rate ************************ */
        // 判断是否需要过滤
        if (MAF > 0 && MISSRATE < 1)
        {
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

                // 没通过阈值，直接下一个循环
                if (MAFTMP < MAF || MISSRATETMP > MISSRATE)
                {
                    continue;
                }
            }
        }



        // 基因组上的ref序列信息
        // 先看提交的基因组中有没有对应的染色体信息，没有的话退出程序
        if (seqMap.find(INFOSTRUCTTMP.CHROM) == seqMap.end())
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Error: "
                    << INFOSTRUCTTMP.CHROM 
                    << " -> Chromosome does not exist in refgenome." 
                    << endl;
            exit(1);
        }

        // 检查qrySeq有没有<之类的信息，有的话该位点直接写出 (<INS>;<DUP>)
        if (INFOSTRUCTTMP.INFOVec[4].find("<") != string::npos)
        {
            outTxt += INFOSTRUCTTMP.INFO + "\n";

            continue;
        }
        

        // 1. 检查第一个变异是否大于read length，如果大于的话再进行转换，否则跳过该变异
        if (readLen > INFOSTRUCTTMP.POS) // paragraph的要求
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Warning: start of variation is less than read length, skip this site -> "
                    << INFOSTRUCTTMP.CHROM << "\t" 
                    << INFOSTRUCTTMP.POS << endl;

            continue;
        }
        

        // 2. vcf对应到reference上的序列
        if (seqMap[INFOSTRUCTTMP.CHROM].size() < (INFOSTRUCTTMP.END)) // 检查染色体长度是否正确
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Chromosome length error: " 
                    << INFOSTRUCTTMP.CHROM << ".size()=" << seqMap[INFOSTRUCTTMP.CHROM].size() << endl;
            exit(1);
        }
        
        string trueRefSeq = seqMap[INFOSTRUCTTMP.CHROM].substr(INFOSTRUCTTMP.POS - 1, INFOSTRUCTTMP.LEN);
        if (trueRefSeq != INFOSTRUCTTMP.REF) // 如果和参考基因组上的序列不一样，则替换成参考基因组的序列
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Warning: sequence difference between refgenome and vcf, replace by refgenome sequence -> "
                    << INFOSTRUCTTMP.CHROM << " " 
                    << INFOSTRUCTTMP.POS << " " 
                    << INFOSTRUCTTMP.REF << "(" << trueRefSeq << ")\n";

            INFOSTRUCTTMP.INFOVec[3] = trueRefSeq;
            INFOSTRUCTTMP.REF = trueRefSeq;
        }


        // 3. 对vcf的第八列 'END=' 进行替换
        // paragrpah (raise Exception("{}:{} error in adding ref support.".format(start, end)))
        // 正则表达式对END进行重新替换
        std::regex endReg("END=\\d+");
        INFOSTRUCTTMP.INFOVec[7] = regex_replace(INFOSTRUCTTMP.INFOVec[7], endReg, "END=" + to_string(INFOSTRUCTTMP.END));


        // 4. Different padding base for REF and ALT. (paragraph要求SV的ALT第一个字母要和REF一样)
        for (int i = 0; i < INFOSTRUCTTMP.ALTVec.size(); i++)
        {
            string qrySeq = INFOSTRUCTTMP.ALTVec[i];

            // 第一个碱基不一样的时候，对qrySeq的序列加上refSrq的前一个碱基
            if (qrySeq[0] != INFOSTRUCTTMP.REF[0] && (qrySeq.length() > 1 || INFOSTRUCTTMP.REF.length() > 1))
            {
                // 坐标往前挪1
                INFOSTRUCTTMP.POS = INFOSTRUCTTMP.POS - 1;
                INFOSTRUCTTMP.INFOVec[1] = to_string(INFOSTRUCTTMP.POS);
                INFOSTRUCTTMP.LEN = INFOSTRUCTTMP.LEN + 1;
                INFOSTRUCTTMP.REF = seqMap[INFOSTRUCTTMP.CHROM].substr(INFOSTRUCTTMP.POS - 1, INFOSTRUCTTMP.LEN); // 重新提取序列信息
                INFOSTRUCTTMP.INFOVec[3] = INFOSTRUCTTMP.REF; // 对Vector进行赋值

                // 将 'INFOSTRUCTTMP.ALTVec' 中所有序列都加上 'refSeq[0]'
                for (size_t j = 0; j < INFOSTRUCTTMP.ALTVec.size(); j++)
                {
                    INFOSTRUCTTMP.ALTVec[j] = INFOSTRUCTTMP.REF[0] + INFOSTRUCTTMP.ALTVec[j];
                }

                qrySeq = INFOSTRUCTTMP.ALTVec[i]; // qrySeq 重新赋值
            }

            // 5. 检查refSeq和qrySeq是否一样，一样则跳过
            if (qrySeq == INFOSTRUCTTMP.REF)
            {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "Warning: sequence same in REF and ALT, skip this site -> "
                        << " (" << strip(INFOSTRUCTTMP.INFO, '\n') << ")"
                        << endl;

                continue;
            }

            // 6. 检查refSeq和qrySeq中是否含有atgcnATGCN外的字符，含有的话跳过该位点
            smatch results;
            std::regex atgcReg("[^ATGCNatgcn]");
            if (regex_search(qrySeq, results, atgcReg) != 0)
            {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "Warning: sequence contains non-ATGCNatgcn characters, skip this site -> "
                        << " (" << strip(INFOSTRUCTTMP.INFO, '\n') << ")"
                        << endl;

                continue;
            }
        }
        // 对qry序列进行替换
        INFOSTRUCTTMP.INFOVec[4] = join(INFOSTRUCTTMP.ALTVec, ",");


        // 7. 检查替换后的qry有没有相同的，有的话跳过该位点
        // bayestyper要求 (A    ATG,TGT,ATG)
        // Assertion `count(alt_alleles.begin() + i + 1, alt_alleles.end(), alt_alleles.at(i)) == 0' failed.
        int qrySeqNum = 0;
        for (auto it : INFOSTRUCTTMP.ALTVec)
        {
            int qrySeqNumTmp = count(INFOSTRUCTTMP.ALTVec.begin(), INFOSTRUCTTMP.ALTVec.end(), it);
            qrySeqNum = max(qrySeqNumTmp, qrySeqNum);
        }
        if (qrySeqNum > 1)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Warning: allelic repeat, skip -> "
                    << " (" << strip(INFOSTRUCTTMP.INFO, '\n') << ")"
                    << endl;

            continue;
        }


        // 8. 检查变异是否重复，重复的话跳过
        // 如果是新的染色体，则将起始位置归零
        if (INFOSTRUCTTMP.CHROM != preChromosome)
        {
            preChromosome = INFOSTRUCTTMP.CHROM;
            preRefStart = 0;
        }
        if (INFOSTRUCTTMP.POS == preRefStart)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Warning: multiple variants observed on position " 
                    << INFOSTRUCTTMP.POS 
                    << " on " 
                    << INFOSTRUCTTMP.CHROM
                    << ", skip this site"
                    << endl;

            continue;
        }


        // 9. 对基因型列进行检查。只能有二倍体变异、'/'替换为'|'、'.'替换为'.|.'。PanGenie软件
        // 找FORMAT字段中gt的位置
        vector<string> formatVec = split(strip(INFOSTRUCTTMP.INFOVec[8], '\n'), ":");
        vector<string>::iterator gtItera = find(formatVec.begin(), formatVec.end(), "GT");
        int maxGT = 0; // 记录最大的GT，看是否比qry的序列还多，多的话跳过该位点
        // 只要GT字段
        INFOSTRUCTTMP.INFOVec[8] = "GT";

        int gtIndex = 0;
        if (gtItera != formatVec.end()) // FORMAT中有GT
        {
            // GT的index
            gtIndex = distance(formatVec.begin(), gtItera);

            // 对GT列进行循环
            for (int i = 9; i < INFOSTRUCTTMP.INFOVec.size(); i++)
            {
                // 找到基因型
                string gt = strip(split(INFOSTRUCTTMP.INFOVec[i], ":")[gtIndex], '\n');

                // 看是不是'.'，是的话替换为'.|.'
                if (gt == ".")
                {
                    gt = ".|.";
                }

                // 将'/'替换为'|'
                if (gt.find("/") != string::npos)
                {
                    std::regex reg("/");
                    gt = regex_replace(string(gt), regex(reg), string("|"));
                }
                INFOSTRUCTTMP.INFOVec[i] = gt;

                // 循环找最大的GT
                if (gt.find("/") != string::npos)
                {
                    for (auto it1 : split(gt, "/"))
                    {
                        try
                        {
                            if (stoi(it1) > maxGT)
                            {
                                maxGT = stoi(it1);
                            }
                        }
                        catch(const std::invalid_argument& e)
                        {
                            continue;
                        }
                    }
                }
                else if (gt.find("|") != string::npos)
                {
                    for (auto it1 : split(gt, "|"))
                    {
                        try
                        {
                            if (stoi(it1) > maxGT)
                            {
                                maxGT = stoi(it1);
                            }
                        }
                        catch(const std::invalid_argument& e)
                        {
                            continue;
                        }
                    }
                }
                else
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] "
                            << "Error: invalid delimiter: " 
                            << " (" << strip(INFOSTRUCTTMP.INFO, '\n') << ")"
                            << endl;
                    exit(1);
                }
            }
        }
        else // 没有的话退出代码
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Error: GT not in FORMAT column -> " 
                    << " (" << strip(INFOSTRUCTTMP.INFO, '\n') << ")"
                    << endl;
            exit(1);
        }


        // 10. 检查GT是不是比qry的序列还多
        // pangenie要求：VariantReader::VariantReader: invalid genotype in VCF.
        if (INFOSTRUCTTMP.ALTVec.size() < maxGT)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Warning: there are more GTs than qry sequences: " 
                    << " (" << strip(INFOSTRUCTTMP.INFO, '\n') << ")"
                    << endl;

            continue;
        }
        

        // 给preRefStart重新赋值
        preRefStart = INFOSTRUCTTMP.POS;

        // 检查information是否被清空，清空的话跳过
        if (INFOSTRUCTTMP.INFO.empty())
        {
            continue;
        }

        // 将替换后的字符串加到outTxt中
        outTxt += join(INFOSTRUCTTMP.INFOVec, "\t") + "\n";

        if (outTxt.size() > 10000000) // 每10Mb写入一次，减少磁盘IO
        {
            SAVECLASS.save(
                outTxt
            );

            // 清空字符串
            outTxt.clear();
            string().swap(outTxt);
        }
    }


    if (outTxt.size() >= 0)  // 最后写一次
    {
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
}