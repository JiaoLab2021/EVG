// g++ vcf_count.cpp -o vcf_count -lz -O2
#include "../include/vcf_count.hpp"


using namespace std;


int main_count(int argc, char** argv)
{
    // �����ļ�
    string vcfFileName;

    // �������
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

    // �������Ƿ���ȷ
    if (vcfFileName.empty())
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error.\n";
        help_count(argv);
        return 1;
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running.\n";
    
    // ͳ�ƽ��
    VCFCOUNT::count(vcfFileName);
    
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done.\n";

    return 0;
}

// �����ĵ�
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
 * @brief ͳ��vcf�ļ�
 * 
 * @param vcfFileName  ����vcf�ļ�
 * 
 * @return 0
**/
int VCFCOUNT::count(
    const string & vcfFileName
)
{
    // ��¼���ֱ��������
    uint32_t snpNum = 0;
    uint32_t indelNum = 0;
    uint32_t insNum = 0;
    uint32_t delNum = 0;
    // uint32_t invNum = 0;
    // uint32_t dupNum = 0;
    uint32_t otherNum = 0;

    // �洢vcf��Ϣ
    VCFINFOSTRUCT INFOSTRUCTTMP;

    // ��ʼ��
    VCFOPEN VCFOPENCLASS;
    VCFOPENCLASS.init(
        vcfFileName
    );

    // ��vcf�ļ�
    VCFOPENCLASS.open();

    // ���û�б����꣬����
    while (VCFOPENCLASS.read(INFOSTRUCTTMP))
    {
        // �г�����Ϣ�ټ���
        if (INFOSTRUCTTMP.POS > 0)
        {
            // ����qry�����б�
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

    // �ر�vcf�ļ�
    VCFOPENCLASS.close();

    // ��ӡ���
    cout << "SNP\tInDels\tDeletion\tInsertion\tOther\n";
    cout << snpNum << "\t" << indelNum << "\t" << delNum << "\t" << insNum << "\t" << otherNum << endl;
    
    return 0;
}