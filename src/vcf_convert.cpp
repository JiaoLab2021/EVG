// g++ vcf_convert.cpp -o vcf_convert -lz

#include "../include/vcf_convert.hpp"

using namespace std;


int main_convert(int argc, char** argv)
{
    // 输入文件和参数
    string refFileName;
    string vcfFileName;
    string outFileName = "vcfConvert.out.vcf.gz";
    int readLen = 350; // 测序的读长

    // 过滤阈值
    double MAF = 0;  // 最小等位基因频率
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
            refFileName = optarg;
            break;
        case 'v':
            vcfFileName = optarg;
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
            outFileName = optarg;
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
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running ...\n";

    // If the length of "raed" is greater than 1000, there's no need to filter by read length because we don't have to process paragraph
    readLen = (readLen > 1000) ? 0 : readLen;
    

    // init
    Convert ConvertClass(refFileName, vcfFileName, readLen, MAF, MISSRATE, outFileName);
    // build reference index
    ConvertClass.build_reference_index();
    // Convert VCF file
    ConvertClass.vcf_convert();

    // 打印log
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done ...\n";

    return 0;
}

// 帮助文档
void help_convert(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -r FILE -v FILE [options]" << endl
       << "convert VCF files to the required format for genome graph software." << endl
       << endl
       << "required arguments:" << endl
       << "    -r, --reference   FILE     input FASTA reference" << endl
       << "    -v, --vcf         FILE     input VCF" << endl
       << endl
	   << "optional arguments:" << endl
       << "    -l, --length      INT      read length [350]" << endl
       << "    --maf             FLOAT    exclude SNPs with minor allele frequency lower than threshold [0.0]" << endl
       << "    --geno            FLOAT    exclude SNPs with missing call frequencies greater than threshold [1.0]" << endl
       << "    -o, --out         FILE     output file name [vcfConvert.out.vcf.gz]" << endl
       << endl
       << "    -h, --help                 print this help document" << endl;
}



/**
 * @brief Convert vcf to the format required by graph genome tools
 * 
 * @param refFileName     reference genome
 * @param vcfFileName     input VCF file name
 * @param readLen         read length
 * @param MAF             Minimal allele frequency
 * @param MISSRATE        Missing rate
 * @param outFileName     output file name
 * 
**/
Convert::Convert(
    string refFileName, 
    string vcfFileName, 
    int readLen, 
    double MAF, 
    double MISSRATE, 
    string outFileName
) : refFileName_(refFileName), vcfFileName_(vcfFileName), readLen_(readLen), MAF_(MAF), MISSRATE_(MISSRATE), outFileName_(outFileName) {}


/**
 * @brief build the reference genome index
 * 
 * @return void
**/
void Convert::build_reference_index()
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
    gzFile gzfp = gzopen(refFileName_.c_str(), "rb");

    // 打开文件
    if(!gzfp)
    {
        cerr << "[" << __func__ << "::" << getTime() << "] "
            << "'" << refFileName_ << "': No such file or directory." << endl;
        exit(1);
    }
    else
    {
        kseq_t *ks;
        ks = kseq_init(gzfp);
    
        while( kseq_read(ks) >= 0 )
        {
            string chromosome = ks->name.s;
            uint32_t chrLen = ks->seq.l;
            string sequence = ks->seq.s;

            refIndexS_.chrLenTxt += "##contig=<ID=" + chromosome + ",length=" + to_string(chrLen) + ">\n";
            refIndexS_.sequenceMap[chromosome] = sequence;

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
 1. Add Chromosome Length to Header
 2. Check if qrySeq contains any '<'. If it does, skip it. (<INS>;<DUP>)
 3. First mutation should be greater than read length
 4. Sequences must correspond to the reference
 5. Replace the 'END=' in the eighth column of the VCF
 6. Different padding base for REF and ALT. (paragraph要求SV的ALT第一个字母要和REF一样)
 7. 检查ref和qry序列是否一样，一样的位点跳过
 8. 检查refSeq和qrySeq中是否含有atgcnATGCN外的字符，含有的话跳过该位点
 9. 检查替换后的qry有没有相同的，有的话跳过该位点 --> e.
 10. 检查是否有位置重复的变异
 11. 将基因型中的.转为.|. -> PanGenie
 12. 将基因型中的/转为|
 13. 只保留基因型中的二倍体变异
 14. 检查GT是不是比qry的序列还多
*/
/**
 * @brief Convert vcf to the format required by graph genome tools
 * 
 * @return void
**/
void Convert::vcf_convert()
{
    // Whether check file is sorted
    check_vcf_sort(vcfFileName_);

    // 基因组上的序列信息
    const map<string, string>& seqMap = refIndexS_.sequenceMap;

    // 保存转换后的结果
    string outTxt;

    // 输出文件流
    SAVE SAVECLASS(outFileName_);

    // 记录上一个变异的start和染色体号
    string preChromosome;
    uint32_t preRefStart = 0;


    // 输入文件流
    // 存储vcf信息
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS(vcfFileName_);

    // If not traversed, continue
    while (VCFOPENCLASS.read(INFOSTRUCTTMP))
    {
        // empty line, skip
        if (INFOSTRUCTTMP.line.empty())
        {
            continue;
        }
        
        // comment line
        if (INFOSTRUCTTMP.line.find("#") != string::npos)
        {
            if (INFOSTRUCTTMP.line.find("#CHROM") != string::npos)
            {
                // 1. Add Chromosome Length to Header
                outTxt += refIndexS_.chrLenTxt;

                // save header
                outTxt += INFOSTRUCTTMP.line + "\n";
            }
            // Skip if the vcf contains chromosome length information
            else if (INFOSTRUCTTMP.line.find(",length") == string::npos)
            {
                outTxt += INFOSTRUCTTMP.line + "\n";
            }

            continue;
        }


        /* ************************ filter SNPs by maf and missing rate ************************ */
        // 判断是否需要过滤
        if (MAF_ > 0 && MISSRATE_ < 1)
        {
            // 获取变异类型
            INFOSTRUCTTMP.ID = VCFOPENCLASS.get_TYPE(
                INFOSTRUCTTMP.LEN, 
                INFOSTRUCTTMP.ALTVec
            );

            if (INFOSTRUCTTMP.ID == "SNP")  // SNP时再判断是否过滤
            {
                double MAFTMP;  // 最小等位基因频率
                double MISSRATETMP;  // 缺失率

                // 获取所有的基因型   map<idx, vector<gtString>>
                map<int, vector<string> > GTVecMapTmp = VCFOPENCLASS.get_gt(
                    INFOSTRUCTTMP.lineVec
                );

                // 如果只有一个基因型，跳过。
                if (GTVecMapTmp.empty())
                {
                    continue;
                }
                
                // 计算最小等位基因频率和缺失率
                tie(MAFTMP, MISSRATETMP) = VCFOPENCLASS.calculate(
                    GTVecMapTmp, 
                    INFOSTRUCTTMP.lineVec.size() - 9
                );

                // 没通过阈值，直接下一个循环
                if (MAFTMP < MAF_ || MISSRATETMP > MISSRATE_)
                {
                    continue;
                }
            }
        }


        // Check if there is corresponding chromosome information in the submitted genome
        if (seqMap.find(INFOSTRUCTTMP.CHROM) == seqMap.end())
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Error: '"
                << INFOSTRUCTTMP.CHROM 
                << "' is not present in '" << refFileName_ << "'" 
                << endl;
            exit(1);
        }

        // 2. Check if qrySeq contains any '<'. If it does, skip it. (<INS>;<DUP>)
        if (INFOSTRUCTTMP.lineVec[4].find("<") != string::npos)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: The query sequence contains the '>' symbol, skip this site -> " 
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << "\t" 
                << INFOSTRUCTTMP.lineVec[4] << endl;
            continue;
        }
        

        // 3. First mutation should be greater than read length
        if (static_cast<uint32_t>(readLen_) > INFOSTRUCTTMP.POS) // paragraph的要求
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: Start of variation is less than read length, skip this site -> "
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }
        

        // 4.Sequences must correspond to the reference
        if (seqMap.at(INFOSTRUCTTMP.CHROM).size() < (INFOSTRUCTTMP.END)) // 检查染色体长度是否正确
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Error: The variant end position is greater than the chromosome length -> " 
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << endl;
            exit(1);
        }
        
        string trueRefSeq = seqMap.at(INFOSTRUCTTMP.CHROM).substr(INFOSTRUCTTMP.POS - 1, INFOSTRUCTTMP.LEN);
        if (trueRefSeq != INFOSTRUCTTMP.REF) // 如果和reference genome上的序列不一样，则替换成reference genome的序列
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: sequence difference between refgenome and vcf, replace by refgenome sequence -> "
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << endl;

            INFOSTRUCTTMP.lineVec[3] = trueRefSeq;
            INFOSTRUCTTMP.REF = trueRefSeq;
        }


        // 5. Replace the 'END=' in the eighth column of the VCF
        // paragrpah (raise Exception("{}:{} error in adding ref support.".format(start, end)))
        // 正则表达式对END进行重新替换
        std::regex endReg("END=\\d+");
        INFOSTRUCTTMP.lineVec[7] = regex_replace(INFOSTRUCTTMP.lineVec[7], endReg, "END=" + to_string(INFOSTRUCTTMP.END));


        // 6. Different padding base for REF and ALT. (paragraph要求SV的ALT第一个字母要和REF一样)
        for (size_t i = 0; i < INFOSTRUCTTMP.ALTVec.size(); i++)
        {
            string qrySeq = INFOSTRUCTTMP.ALTVec[i];

            // When the first base is different, add the preceding base of refSeq to the sequence of qrySeq
            if (qrySeq[0] != INFOSTRUCTTMP.REF[0] && (qrySeq.length() > 1 || INFOSTRUCTTMP.REF.length() > 1))
            {
                // Move the coordinates one step backward
                INFOSTRUCTTMP.POS = INFOSTRUCTTMP.POS - 1;
                INFOSTRUCTTMP.lineVec[1] = to_string(INFOSTRUCTTMP.POS);
                INFOSTRUCTTMP.LEN = INFOSTRUCTTMP.LEN + 1;
                INFOSTRUCTTMP.REF = seqMap.at(INFOSTRUCTTMP.CHROM).substr(INFOSTRUCTTMP.POS - 1, INFOSTRUCTTMP.LEN); // 重新提取序列信息
                INFOSTRUCTTMP.lineVec[3] = INFOSTRUCTTMP.REF; // 对Vector进行赋值

                // Add 'refSeq[0]' to all sequences in 'INFOSTRUCTTMP.ALTVec'
                for (size_t j = 0; j < INFOSTRUCTTMP.ALTVec.size(); j++)
                {
                    INFOSTRUCTTMP.ALTVec[j] = INFOSTRUCTTMP.REF[0] + INFOSTRUCTTMP.ALTVec[j];
                }

                qrySeq = INFOSTRUCTTMP.ALTVec[i]; // qrySeq 重新赋值
            }

            // 7. 检查ref和qry序列是否一样，一样的位点跳过
            if (qrySeq == INFOSTRUCTTMP.REF)
            {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Warning: Sequence same in REF and ALT, skip this site -> "
                    << INFOSTRUCTTMP.CHROM << "\t" 
                    << INFOSTRUCTTMP.POS << endl;
                continue;
            }

            // 8. 检查refSeq和qrySeq中是否含有atgcnATGCN外的字符，含有的话跳过该位点
            smatch results;
            std::regex atgcReg("[^ATGCNatgcn]");
            if (!regex_search(qrySeq, results, atgcReg))
            {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Warning: Sequence contains non-ATGCNatgcn characters, skip this site -> "
                    << INFOSTRUCTTMP.CHROM << "\t" 
                    << INFOSTRUCTTMP.POS << endl;
                continue;
            }
        }
        // 对qry序列进行替换
        INFOSTRUCTTMP.lineVec[4] = join(INFOSTRUCTTMP.ALTVec, ",");


        // 9. 检查替换后的qry有没有相同的，有的话跳过该位点 --> e.
        // bayestyper要求 (A    ATG,TGT,ATG)
        // Assertion `count(alt_alleles.begin() + i + 1, alt_alleles.end(), alt_alleles.at(i)) == 0' failed.
        for (const auto& it : INFOSTRUCTTMP.ALTVec)
        {
            if (count(INFOSTRUCTTMP.ALTVec.begin(), INFOSTRUCTTMP.ALTVec.end(), it) > 1)
            {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Warning: Allelic repeat, skip this site -> "
                    << INFOSTRUCTTMP.CHROM << "\t" 
                    << INFOSTRUCTTMP.POS << endl;
                continue;
            }
        }


        // 10. 检查是否有位置重复的变异
        // 如果是新的染色体，则将起始位置归零
        if (INFOSTRUCTTMP.CHROM != preChromosome)
        {
            preChromosome = INFOSTRUCTTMP.CHROM;
            preRefStart = 0;
        }
        if (INFOSTRUCTTMP.POS == preRefStart)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: multiple variants observed, skip this site -> " 
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }


        // 找FORMAT字段中gt的位置
        vector<string> formatVec = split(INFOSTRUCTTMP.FORMAT, ":");
        vector<string>::iterator gtItera = find(formatVec.begin(), formatVec.end(), "GT");
        uint32_t maxGT = 0; // 记录最大的GT，看是否比qry的序列还多，多的话跳过该位点
        // 只要GT字段
        INFOSTRUCTTMP.lineVec[8] = "GT";

        if (gtItera != formatVec.end()) // FORMAT中有GT
        {
            // GT的index
            uint32_t gtIndex = distance(formatVec.begin(), gtItera);
            
            // 对GT列进行循环
            for (size_t i = 9; i < INFOSTRUCTTMP.lineVec.size(); i++) {
                // 找到基因型
                string gt = split(INFOSTRUCTTMP.lineVec[i], ":")[gtIndex];

                // 11. 将基因型中的.转为.|. -> PanGenie
                if (gt == ".")
                {
                    gt = ".|.";
                }
                else if (gt == "0")
                {
                    gt = "0|0";
                }
                

                // 12. 将基因型中的/转为| -> PanGenie
                if (gt.find("/") != string::npos)
                {
                    std::regex reg("/");
                    gt = regex_replace(string(gt), regex(reg), string("|"));
                }

                // 循环找最大的GT
                vector<string> gtVec = VCFOPENCLASS.gt_split(gt);
                for (auto it1 : gtVec) {
                    try {
                        if (stoul(it1) > maxGT) {
                            maxGT = stoi(it1);
                        }
                    } catch(const std::invalid_argument& e) {
                        continue;
                    }
                }

                // 13. 只保留基因型中的二倍体变异
                if (gtVec.size() == 1)
                {
                    gt = gtVec[0] + "|0";
                }
                else if (gtVec.size() > 2)
                {
                    gt = gtVec[0] + "|" + gtVec[1];
                }
                
                INFOSTRUCTTMP.lineVec[i] = gt;
            }
        }
        else // 跳过该位点
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: GT not in FORMAT column, skip this site -> " 
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }


        // 14. 检查GT是不是比qry的序列还多
        // pangenie要求：VariantReader::VariantReader: invalid genotype in VCF.
        if (INFOSTRUCTTMP.ALTVec.size() < maxGT)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: there are more GTs than qry sequences: " 
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }
        

        // 给preRefStart重新赋值
        preRefStart = INFOSTRUCTTMP.POS;

        // 检查information是否被清空，清空的话跳过
        if (INFOSTRUCTTMP.line.empty())
        {
            continue;
        }

        // 将替换后的字符串加到outTxt中
        outTxt += join(INFOSTRUCTTMP.lineVec, "\t") + "\n";

        if (outTxt.size() > 10 * 1024 * 1024) // 每10Mb写入一次，减少磁盘IO
        {
            SAVECLASS.save(outTxt);

            // 清空字符串
            outTxt.clear();
            string().swap(outTxt);
        }
    }


    if (outTxt.size() > 0)  // 最后写一次
    {
        SAVECLASS.save(outTxt);

        // 清空
        outTxt.clear();
        string().swap(outTxt);
    }
}