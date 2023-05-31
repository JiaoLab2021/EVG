// g++ vcf_split.cpp -o vcf_split -lz -O2
#include "../include/vcf_filter.hpp"

using namespace std;

int main_filter(int argc, char** argv)
{
    // �����ļ�
    string vcfFileName;

    // ������ֵ
    double MAF = 0.01;  // ��С��λ����Ƶ��
    double MISSRATE = 0.1;  // ȱʧ��

    // ����ļ���
    string outputFileName = "";

    // �������
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

    // �������Ƿ���ȷ
    if (vcfFileName.empty())
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error.\n";
        help_filter(argv);
        return 1;
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running.\n";
    
    // ͳ�ƽ��
    VCFFILTER::vcf_filter(
        vcfFileName, 
        outputFileName, 
        MAF, 
        MISSRATE
    );
    
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done.\n";

    return 0;
}

// �����ĵ�
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
 * @brief ���ݴε�λ����Ƶ�ʺ�ȱʧ�ʹ���SNPs
 * 
 * @param vcfFileName    ����vcf�ļ�
 * @param outputFileName ����ļ���
 * @param MAF            �ε�λ����Ƶ��
 * @param MISSRATE       ȱʧ��
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
    // �����ļ���
    // �洢vcf��Ϣ
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS;
    VCFOPENCLASS.init(
        vcfFileName
    );
    // ��vcf�ļ�
    VCFOPENCLASS.open();


    // ����ļ���
    SAVE::SAVE SAVECLASS;
    SAVECLASS.init(
        outputFileName
    );
    SAVECLASS.open();


    // ��ʱ�洢����ַ���
    string outTxt = "";

    // ���û�б����꣬����
    while (VCFOPENCLASS.read(INFOSTRUCTTMP))
    {
        // ���ע���У�ֱ�ӱ���
        if (INFOSTRUCTTMP.INFO.find("#") != string::npos)
        {
            outTxt += INFOSTRUCTTMP.INFO + "\n";
            continue;
        }

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

            // ͨ����ֵ��ֱ�ӱ���
            if (MAFTMP >= MAF && MISSRATETMP <= MISSRATE)
            {
                outTxt += INFOSTRUCTTMP.INFO + "\n";
            }
        }
        else  // �������͵ı���ֱ�ӱ���
        {
            outTxt += INFOSTRUCTTMP.INFO + "\n";
        }

        if (outTxt.size() >= 10000000)  // ÿ10mдһ��
        {
            // ����ļ���
            SAVECLASS.save(
                outTxt
            );

            // ���
            outTxt.clear();
            string().swap(outTxt);
        }
    }

    if (outTxt.size() >= 0)  // ���дһ��
    {
        // �������͵ı���ֱ�ӱ���
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
    
    return 0;
}