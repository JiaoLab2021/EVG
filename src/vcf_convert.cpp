// g++ vcf_convert.cpp -o vcf_convert -lz
#include <getopt.h>

#include "../include/vcf_convert.hpp"

using namespace std;


int main_convert(int argc, char** argv)
{
    // �����ļ��Ͳ���
    string referenceFilename;
    string vcfFilename;
    string outFilename = "vcfConvert.out.vcf.gz";
    int readLen = 350; // ����Ķ���

    // ������ֵ
    double MAF = 0.;  // ��С��λ����Ƶ��
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

    // ��ӡlog
    cerr << "[" << __func__ << "::" << getTime() << "] " 
         << "Running.\n";

    // ����reference������
    CONVERT::refIndexStruct refIndexS;
    CONVERT::build_reference_index(
        referenceFilename, 
        refIndexS
    );

    // ���raed�ĳ��ȴ���1000������read length���ˣ���Ϊ������paragraph
    if (readLen > 1000)
    {
        readLen = 0;
    }

    // vcf�ļ�ת��
    CONVERT::vcf_convert(
        vcfFilename, 
        readLen, 
        MAF, 
        MISSRATE, 
        refIndexS, 
        outFilename
    );

    // ��ӡlog
    cerr << "[" << __func__ << "::" << getTime() << "] " 
         << "Done.\n";

    return 0;
}

// �����ĵ�
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
    ��fasta�ļ�
    inputFileName -> refgenome
    refIndexS -> ����contig���Ⱥ�������Ϣ
*/
void CONVERT::build_reference_index(string inputFileName, refIndexStruct & refIndexS)
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
    gzFile gzfp = gzopen(inputFileName.c_str(), "rb");

    // ���ļ�
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
    ��vcf�ļ������滻
    vcfFilename -> ��Ҫת����vcf�ļ�
    readLen -> �����ļ���read����
    refIndexS -> contig���Ⱥ�������Ϣ
    outFilename -> ����ļ���
    a. ��ͷ��Ⱦɫ�峤��
    b. ��һ������Ҫ����read length
    c. ���б�����reference��Ӧ
    d. ��vcf�ĵڰ��� 'END=' �����滻
    e. refSeq�ĵ�һ�����Ҫ��qrySeq�ĵ�һ�����һ��
    f. ���ref��qry�����Ƿ�һ����һ����λ������
    g. ���refSeq��qrySeq���Ƿ���atgcnATGCN����ַ������еĻ�������λ��
    h. ����滻���qry��û����ͬ�ģ��еĻ�������λ�� --> e.
    i. ����Ƿ���λ���ظ��ı���
    j. ���������е�.תΪ.|.
    k. ���������е�/תΪ|
    l. ֻ�����������еĶ��������
    m. ���GT�ǲ��Ǳ�qry�����л���
*/
/**
 * @brief ��vcfת��Ϊ graph genome tools ��Ҫ�ĸ�ʽ
 * 
 * @param vcfFileName     ����vcf�ļ�
 * @param readLen         ����
 * @param MAF            �ε�λ����Ƶ��
 * @param MISSRATE       ȱʧ��
 * @param refIndexS       �ο������鳤�Ⱥ�������Ϣ
 * @param outFilename     ����ļ���
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
    // ����ļ��Ƿ�����
    check_vcf_sort(vcfFilename);

    // �������ϵ�������Ϣ
    map<string, string> seqMap = refIndexS.sequenceMap;

    // ����ת����Ľ��
    string outTxt;


    // ����ļ���
    SAVE::SAVE SAVECLASS;
    SAVECLASS.init(
        outFilename
    );
    SAVECLASS.open();


    // ��¼��һ�������start��Ⱦɫ���
    string preChromosome;
    long long int preRefStart = 0;


    // �����ļ���
    // �洢vcf��Ϣ
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS;
    VCFOPENCLASS.init(
        vcfFilename
    );
    // ��vcf�ļ�
    VCFOPENCLASS.open();


    // ���û�б����꣬����
    while (VCFOPENCLASS.read(INFOSTRUCTTMP))
    {
        // ע����
        if (INFOSTRUCTTMP.INFO.find("#") != string::npos)
        {
            // ���vcf�к���Ⱦɫ�峤����Ϣ������
            if (INFOSTRUCTTMP.INFO.find("#CHROM") != string::npos)
            {
                // ����Ⱦɫ�峤����Ϣ
                outTxt += refIndexS.chrLenTxt;

                // �����ͷ
                outTxt += INFOSTRUCTTMP.INFO + "\n";
            }
            else if (INFOSTRUCTTMP.INFO.find(",length") == string::npos)
            {
                outTxt += INFOSTRUCTTMP.INFO + "\n";
            }

            continue;
        }

        // ��ע����
        if (INFOSTRUCTTMP.INFOVec.size() < 9) // �ȼ���ļ��Բ��ԣ���������������
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
        // �ж��Ƿ���Ҫ����
        if (MAF > 0 && MISSRATE < 1)
        {
            // ��ȡ��������
            INFOSTRUCTTMP.TYPE = VCFOPENCLASS.get_TYPE(
                INFOSTRUCTTMP.LEN, 
                INFOSTRUCTTMP.ALTVec
            );

            if (INFOSTRUCTTMP.TYPE == "SNP")  // SNPʱ���ж��Ƿ����
            {
                double MAFTMP;  // ��С��λ����Ƶ��
                double MISSRATETMP;  // ȱʧ��

                // ��ȡ���еĻ�����   map<idx, vector<gtString>>
                map<int, vector<string> > GTVecMapTmp = VCFOPENCLASS.get_gt(
                    INFOSTRUCTTMP.INFOVec
                );

                // ���ֻ��һ�������ͣ�������
                if (GTVecMapTmp.size() <= 1)
                {
                    continue;
                }
                
                // ������С��λ����Ƶ�ʺ�ȱʧ��
                tie(MAFTMP, MISSRATETMP) = VCFOPENCLASS.calculate(
                    GTVecMapTmp, 
                    INFOSTRUCTTMP.INFOVec.size() - 9
                );

                // ûͨ����ֵ��ֱ����һ��ѭ��
                if (MAFTMP < MAF || MISSRATETMP > MISSRATE)
                {
                    continue;
                }
            }
        }



        // �������ϵ�ref������Ϣ
        // �ȿ��ύ�Ļ���������û�ж�Ӧ��Ⱦɫ����Ϣ��û�еĻ��˳�����
        if (seqMap.find(INFOSTRUCTTMP.CHROM) == seqMap.end())
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Error: "
                    << INFOSTRUCTTMP.CHROM 
                    << " -> Chromosome does not exist in refgenome." 
                    << endl;
            exit(1);
        }

        // ���qrySeq��û��<֮�����Ϣ���еĻ���λ��ֱ��д�� (<INS>;<DUP>)
        if (INFOSTRUCTTMP.INFOVec[4].find("<") != string::npos)
        {
            outTxt += INFOSTRUCTTMP.INFO + "\n";

            continue;
        }
        

        // 1. ����һ�������Ƿ����read length��������ڵĻ��ٽ���ת�������������ñ���
        if (readLen > INFOSTRUCTTMP.POS) // paragraph��Ҫ��
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Warning: start of variation is less than read length, skip this site -> "
                    << INFOSTRUCTTMP.CHROM << "\t" 
                    << INFOSTRUCTTMP.POS << endl;

            continue;
        }
        

        // 2. vcf��Ӧ��reference�ϵ�����
        if (seqMap[INFOSTRUCTTMP.CHROM].size() < (INFOSTRUCTTMP.END)) // ���Ⱦɫ�峤���Ƿ���ȷ
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Chromosome length error: " 
                    << INFOSTRUCTTMP.CHROM << ".size()=" << seqMap[INFOSTRUCTTMP.CHROM].size() << endl;
            exit(1);
        }
        
        string trueRefSeq = seqMap[INFOSTRUCTTMP.CHROM].substr(INFOSTRUCTTMP.POS - 1, INFOSTRUCTTMP.LEN);
        if (trueRefSeq != INFOSTRUCTTMP.REF) // ����Ͳο��������ϵ����в�һ�������滻�ɲο������������
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Warning: sequence difference between refgenome and vcf, replace by refgenome sequence -> "
                    << INFOSTRUCTTMP.CHROM << " " 
                    << INFOSTRUCTTMP.POS << " " 
                    << INFOSTRUCTTMP.REF << "(" << trueRefSeq << ")\n";

            INFOSTRUCTTMP.INFOVec[3] = trueRefSeq;
            INFOSTRUCTTMP.REF = trueRefSeq;
        }


        // 3. ��vcf�ĵڰ��� 'END=' �����滻
        // paragrpah (raise Exception("{}:{} error in adding ref support.".format(start, end)))
        // ������ʽ��END���������滻
        std::regex endReg("END=\\d+");
        INFOSTRUCTTMP.INFOVec[7] = regex_replace(INFOSTRUCTTMP.INFOVec[7], endReg, "END=" + to_string(INFOSTRUCTTMP.END));


        // 4. Different padding base for REF and ALT. (paragraphҪ��SV��ALT��һ����ĸҪ��REFһ��)
        for (int i = 0; i < INFOSTRUCTTMP.ALTVec.size(); i++)
        {
            string qrySeq = INFOSTRUCTTMP.ALTVec[i];

            // ��һ�������һ����ʱ�򣬶�qrySeq�����м���refSrq��ǰһ�����
            if (qrySeq[0] != INFOSTRUCTTMP.REF[0] && (qrySeq.length() > 1 || INFOSTRUCTTMP.REF.length() > 1))
            {
                // ������ǰŲ1
                INFOSTRUCTTMP.POS = INFOSTRUCTTMP.POS - 1;
                INFOSTRUCTTMP.INFOVec[1] = to_string(INFOSTRUCTTMP.POS);
                INFOSTRUCTTMP.LEN = INFOSTRUCTTMP.LEN + 1;
                INFOSTRUCTTMP.REF = seqMap[INFOSTRUCTTMP.CHROM].substr(INFOSTRUCTTMP.POS - 1, INFOSTRUCTTMP.LEN); // ������ȡ������Ϣ
                INFOSTRUCTTMP.INFOVec[3] = INFOSTRUCTTMP.REF; // ��Vector���и�ֵ

                // �� 'INFOSTRUCTTMP.ALTVec' ���������ж����� 'refSeq[0]'
                for (size_t j = 0; j < INFOSTRUCTTMP.ALTVec.size(); j++)
                {
                    INFOSTRUCTTMP.ALTVec[j] = INFOSTRUCTTMP.REF[0] + INFOSTRUCTTMP.ALTVec[j];
                }

                qrySeq = INFOSTRUCTTMP.ALTVec[i]; // qrySeq ���¸�ֵ
            }

            // 5. ���refSeq��qrySeq�Ƿ�һ����һ��������
            if (qrySeq == INFOSTRUCTTMP.REF)
            {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "Warning: sequence same in REF and ALT, skip this site -> "
                        << " (" << strip(INFOSTRUCTTMP.INFO, '\n') << ")"
                        << endl;

                continue;
            }

            // 6. ���refSeq��qrySeq���Ƿ���atgcnATGCN����ַ������еĻ�������λ��
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
        // ��qry���н����滻
        INFOSTRUCTTMP.INFOVec[4] = join(INFOSTRUCTTMP.ALTVec, ",");


        // 7. ����滻���qry��û����ͬ�ģ��еĻ�������λ��
        // bayestyperҪ�� (A    ATG,TGT,ATG)
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


        // 8. �������Ƿ��ظ����ظ��Ļ�����
        // ������µ�Ⱦɫ�壬����ʼλ�ù���
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


        // 9. �Ի������н��м�顣ֻ���ж�������졢'/'�滻Ϊ'|'��'.'�滻Ϊ'.|.'��PanGenie���
        // ��FORMAT�ֶ���gt��λ��
        vector<string> formatVec = split(strip(INFOSTRUCTTMP.INFOVec[8], '\n'), ":");
        vector<string>::iterator gtItera = find(formatVec.begin(), formatVec.end(), "GT");
        int maxGT = 0; // ��¼����GT�����Ƿ��qry�����л��࣬��Ļ�������λ��
        // ֻҪGT�ֶ�
        INFOSTRUCTTMP.INFOVec[8] = "GT";

        int gtIndex = 0;
        if (gtItera != formatVec.end()) // FORMAT����GT
        {
            // GT��index
            gtIndex = distance(formatVec.begin(), gtItera);

            // ��GT�н���ѭ��
            for (int i = 9; i < INFOSTRUCTTMP.INFOVec.size(); i++)
            {
                // �ҵ�������
                string gt = strip(split(INFOSTRUCTTMP.INFOVec[i], ":")[gtIndex], '\n');

                // ���ǲ���'.'���ǵĻ��滻Ϊ'.|.'
                if (gt == ".")
                {
                    gt = ".|.";
                }

                // ��'/'�滻Ϊ'|'
                if (gt.find("/") != string::npos)
                {
                    std::regex reg("/");
                    gt = regex_replace(string(gt), regex(reg), string("|"));
                }
                INFOSTRUCTTMP.INFOVec[i] = gt;

                // ѭ��������GT
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
        else // û�еĻ��˳�����
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Error: GT not in FORMAT column -> " 
                    << " (" << strip(INFOSTRUCTTMP.INFO, '\n') << ")"
                    << endl;
            exit(1);
        }


        // 10. ���GT�ǲ��Ǳ�qry�����л���
        // pangenieҪ��VariantReader::VariantReader: invalid genotype in VCF.
        if (INFOSTRUCTTMP.ALTVec.size() < maxGT)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Warning: there are more GTs than qry sequences: " 
                    << " (" << strip(INFOSTRUCTTMP.INFO, '\n') << ")"
                    << endl;

            continue;
        }
        

        // ��preRefStart���¸�ֵ
        preRefStart = INFOSTRUCTTMP.POS;

        // ���information�Ƿ���գ���յĻ�����
        if (INFOSTRUCTTMP.INFO.empty())
        {
            continue;
        }

        // ���滻����ַ����ӵ�outTxt��
        outTxt += join(INFOSTRUCTTMP.INFOVec, "\t") + "\n";

        if (outTxt.size() > 10000000) // ÿ10Mbд��һ�Σ����ٴ���IO
        {
            SAVECLASS.save(
                outTxt
            );

            // ����ַ���
            outTxt.clear();
            string().swap(outTxt);
        }
    }


    if (outTxt.size() >= 0)  // ���дһ��
    {
        SAVECLASS.save(
            outTxt
        );

        // ���
        outTxt.clear();
        string().swap(outTxt);
    }

    // �ر��ļ�
    VCFOPENCLASS.close();
    SAVECLASS.close();
}