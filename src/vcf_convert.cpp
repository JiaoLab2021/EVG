// g++ vcf_convert.cpp -o vcf_convert -lz

#include "../include/vcf_convert.hpp"

using namespace std;


int main_convert(int argc, char** argv)
{
    // �����ļ��Ͳ���
    string refFileName;
    string vcfFileName;
    string outFileName = "vcfConvert.out.vcf.gz";
    int readLen = 350; // ����Ķ���

    // ������ֵ
    double MAF = 0;  // ��С��λ����Ƶ��
    double MISSRATE = 1.0;  // ȱʧ��

    // �������
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

    // ��ӡlog
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running ...\n";

    // If the length of "raed" is greater than 1000, there's no need to filter by read length because we don't have to process paragraph
    readLen = (readLen > 1000) ? 0 : readLen;
    

    // init
    Convert ConvertClass(refFileName, vcfFileName, readLen, MAF, MISSRATE, outFileName);
    // build reference index
    ConvertClass.build_reference_index();
    // Convert VCF file
    ConvertClass.vcf_convert();

    // ��ӡlog
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done ...\n";

    return 0;
}

// �����ĵ�
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
    // Ⱦɫ�����������graphtyper��Ҫ
    ofstream outFile;
    outFile.open("CHROMOSOME.NAME", ios::out);
    if(!outFile.is_open())
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
            << "'CHROMOSOME.NAME': No such file or directory." 
            << endl;
        exit(1);
    }

    // �����ļ���
    gzFile gzfp = gzopen(refFileName_.c_str(), "rb");

    // ���ļ�
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

            // ���Ⱦɫ������
            outFile << chromosome + "\n";
        }

        // �ͷ��ڴ棬�ر��ļ�
        kseq_destroy(ks);
        gzclose(gzfp);
    }

    // �ͷ��ڴ�
    outFile.close();
}


/*
 1. Add Chromosome Length to Header
 2. Check if qrySeq contains any '<'. If it does, skip it. (<INS>;<DUP>)
 3. First mutation should be greater than read length
 4. Sequences must correspond to the reference
 5. Replace the 'END=' in the eighth column of the VCF
 6. Different padding base for REF and ALT. (paragraphҪ��SV��ALT��һ����ĸҪ��REFһ��)
 7. ���ref��qry�����Ƿ�һ����һ����λ������
 8. ���refSeq��qrySeq���Ƿ���atgcnATGCN����ַ������еĻ�������λ��
 9. ����滻���qry��û����ͬ�ģ��еĻ�������λ�� --> e.
 10. ����Ƿ���λ���ظ��ı���
 11. ���������е�.תΪ.|. -> PanGenie
 12. ���������е�/תΪ|
 13. ֻ�����������еĶ��������
 14. ���GT�ǲ��Ǳ�qry�����л���
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

    // �������ϵ�������Ϣ
    const map<string, string>& seqMap = refIndexS_.sequenceMap;

    // ����ת����Ľ��
    string outTxt;

    // ����ļ���
    SAVE SAVECLASS(outFileName_);

    // ��¼��һ�������start��Ⱦɫ���
    string preChromosome;
    uint32_t preRefStart = 0;


    // �����ļ���
    // �洢vcf��Ϣ
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
        // �ж��Ƿ���Ҫ����
        if (MAF_ > 0 && MISSRATE_ < 1)
        {
            // ��ȡ��������
            INFOSTRUCTTMP.ID = VCFOPENCLASS.get_TYPE(
                INFOSTRUCTTMP.LEN, 
                INFOSTRUCTTMP.ALTVec
            );

            if (INFOSTRUCTTMP.ID == "SNP")  // SNPʱ���ж��Ƿ����
            {
                double MAFTMP;  // ��С��λ����Ƶ��
                double MISSRATETMP;  // ȱʧ��

                // ��ȡ���еĻ�����   map<idx, vector<gtString>>
                map<int, vector<string> > GTVecMapTmp = VCFOPENCLASS.get_gt(
                    INFOSTRUCTTMP.lineVec
                );

                // ���ֻ��һ�������ͣ�������
                if (GTVecMapTmp.empty())
                {
                    continue;
                }
                
                // ������С��λ����Ƶ�ʺ�ȱʧ��
                tie(MAFTMP, MISSRATETMP) = VCFOPENCLASS.calculate(
                    GTVecMapTmp, 
                    INFOSTRUCTTMP.lineVec.size() - 9
                );

                // ûͨ����ֵ��ֱ����һ��ѭ��
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
        if (static_cast<uint32_t>(readLen_) > INFOSTRUCTTMP.POS) // paragraph��Ҫ��
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: Start of variation is less than read length, skip this site -> "
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }
        

        // 4.Sequences must correspond to the reference
        if (seqMap.at(INFOSTRUCTTMP.CHROM).size() < (INFOSTRUCTTMP.END)) // ���Ⱦɫ�峤���Ƿ���ȷ
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Error: The variant end position is greater than the chromosome length -> " 
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << endl;
            exit(1);
        }
        
        string trueRefSeq = seqMap.at(INFOSTRUCTTMP.CHROM).substr(INFOSTRUCTTMP.POS - 1, INFOSTRUCTTMP.LEN);
        if (trueRefSeq != INFOSTRUCTTMP.REF) // �����reference genome�ϵ����в�һ�������滻��reference genome������
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
        // ������ʽ��END���������滻
        std::regex endReg("END=\\d+");
        INFOSTRUCTTMP.lineVec[7] = regex_replace(INFOSTRUCTTMP.lineVec[7], endReg, "END=" + to_string(INFOSTRUCTTMP.END));


        // 6. Different padding base for REF and ALT. (paragraphҪ��SV��ALT��һ����ĸҪ��REFһ��)
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
                INFOSTRUCTTMP.REF = seqMap.at(INFOSTRUCTTMP.CHROM).substr(INFOSTRUCTTMP.POS - 1, INFOSTRUCTTMP.LEN); // ������ȡ������Ϣ
                INFOSTRUCTTMP.lineVec[3] = INFOSTRUCTTMP.REF; // ��Vector���и�ֵ

                // Add 'refSeq[0]' to all sequences in 'INFOSTRUCTTMP.ALTVec'
                for (size_t j = 0; j < INFOSTRUCTTMP.ALTVec.size(); j++)
                {
                    INFOSTRUCTTMP.ALTVec[j] = INFOSTRUCTTMP.REF[0] + INFOSTRUCTTMP.ALTVec[j];
                }

                qrySeq = INFOSTRUCTTMP.ALTVec[i]; // qrySeq ���¸�ֵ
            }

            // 7. ���ref��qry�����Ƿ�һ����һ����λ������
            if (qrySeq == INFOSTRUCTTMP.REF)
            {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Warning: Sequence same in REF and ALT, skip this site -> "
                    << INFOSTRUCTTMP.CHROM << "\t" 
                    << INFOSTRUCTTMP.POS << endl;
                continue;
            }

            // 8. ���refSeq��qrySeq���Ƿ���atgcnATGCN����ַ������еĻ�������λ��
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
        // ��qry���н����滻
        INFOSTRUCTTMP.lineVec[4] = join(INFOSTRUCTTMP.ALTVec, ",");


        // 9. ����滻���qry��û����ͬ�ģ��еĻ�������λ�� --> e.
        // bayestyperҪ�� (A    ATG,TGT,ATG)
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


        // 10. ����Ƿ���λ���ظ��ı���
        // ������µ�Ⱦɫ�壬����ʼλ�ù���
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


        // ��FORMAT�ֶ���gt��λ��
        vector<string> formatVec = split(INFOSTRUCTTMP.FORMAT, ":");
        vector<string>::iterator gtItera = find(formatVec.begin(), formatVec.end(), "GT");
        uint32_t maxGT = 0; // ��¼����GT�����Ƿ��qry�����л��࣬��Ļ�������λ��
        // ֻҪGT�ֶ�
        INFOSTRUCTTMP.lineVec[8] = "GT";

        if (gtItera != formatVec.end()) // FORMAT����GT
        {
            // GT��index
            uint32_t gtIndex = distance(formatVec.begin(), gtItera);
            
            // ��GT�н���ѭ��
            for (size_t i = 9; i < INFOSTRUCTTMP.lineVec.size(); i++) {
                // �ҵ�������
                string gt = split(INFOSTRUCTTMP.lineVec[i], ":")[gtIndex];

                // 11. ���������е�.תΪ.|. -> PanGenie
                if (gt == ".")
                {
                    gt = ".|.";
                }
                else if (gt == "0")
                {
                    gt = "0|0";
                }
                

                // 12. ���������е�/תΪ| -> PanGenie
                if (gt.find("/") != string::npos)
                {
                    std::regex reg("/");
                    gt = regex_replace(string(gt), regex(reg), string("|"));
                }

                // ѭ��������GT
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

                // 13. ֻ�����������еĶ��������
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
        else // ������λ��
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: GT not in FORMAT column, skip this site -> " 
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }


        // 14. ���GT�ǲ��Ǳ�qry�����л���
        // pangenieҪ��VariantReader::VariantReader: invalid genotype in VCF.
        if (INFOSTRUCTTMP.ALTVec.size() < maxGT)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: there are more GTs than qry sequences: " 
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }
        

        // ��preRefStart���¸�ֵ
        preRefStart = INFOSTRUCTTMP.POS;

        // ���information�Ƿ���գ���յĻ�����
        if (INFOSTRUCTTMP.line.empty())
        {
            continue;
        }

        // ���滻����ַ����ӵ�outTxt��
        outTxt += join(INFOSTRUCTTMP.lineVec, "\t") + "\n";

        if (outTxt.size() > 10 * 1024 * 1024) // ÿ10Mbд��һ�Σ����ٴ���IO
        {
            SAVECLASS.save(outTxt);

            // ����ַ���
            outTxt.clear();
            string().swap(outTxt);
        }
    }


    if (outTxt.size() > 0)  // ���дһ��
    {
        SAVECLASS.save(outTxt);

        // ���
        outTxt.clear();
        string().swap(outTxt);
    }
}