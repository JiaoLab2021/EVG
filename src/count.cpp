#include "../include/count.hpp"

using namespace std;

int main_count(int argc, char** argv)
{
    // �����ļ�
    string FileName1 = "";
    string FileName2 = "";
    string fastaFileName = "";

    // ����ļ�
    string outputFileName = "";

    // �������
    int c;
    while (true)
    {
        static const struct option long_options[] = 
        {
            {"input", required_argument, 0, 'i'},
            {"out", required_argument, 0, 'o'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "i:o:h", long_options, &option_index);

        if (c == -1)
            break;

        switch (c)
        {
        case 'i':
            if (FileName1.empty())
            {
                FileName1 = optarg;
            }
            else
            {
                FileName2 = optarg;
            }
            break;
        case 'o':
            outputFileName = optarg;
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

    if (argc <= 2)
    {
        help_count(argv);
        return 1;
    }

    // print log
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running.\n";
    
    count::countStruct countOut;
    if (FileName1.length() > 0)
    {
        count::fastq_a_count(
            FileName1, 
            FileName2, 
            countOut, 
            outputFileName
        );
    }
    else
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
             << "parameter error: -i.\n";
        exit(1);
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done.\n";

    return 0;
}

// help document
void help_count(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -i FILE [options]" << endl
       << "calculate the number of bases and reads of fasta/q" << endl
       << endl
       << "required arguments:" << endl
       << "    -i, --input     FILE    input fastq/a, possibly compressed, two are allowed, one for each mate" << endl
       << endl
	   << "optional arguments:" << endl
       << "    -o, --out       FILE    output file name [stdout]" << endl
       << "    -h, --help              print this help document" << endl;
}


// ��fastq/a.gz�ļ�
void count::fastq_a_count(
    string inputFileName1, 
    string inputFileName2, 
    countStruct & countOut, 
    const string & outputFileName
)
{
    if (inputFileName1.length() > 0 && inputFileName2.length() > 0) // ˫�˲���
    {
        // �����ļ���
        gzFile gzfp1 = gzopen(inputFileName1.c_str(), "rb");
        gzFile gzfp2 = gzopen(inputFileName2.c_str(), "rb");

        // ���ļ�
        if(!gzfp1 || !gzfp2)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                    << "'" << inputFileName1 << "' or '" << inputFileName2
                    << "': No such file or directory.\n";
            exit(1);
        }
        else
        {
            kseq_t *ks1;
            kseq_t *ks2;
            ks1 = kseq_init(gzfp1);
            ks2 = kseq_init(gzfp2);
        
            while(kseq_read(ks1) >= 0 && kseq_read(ks2) >= 0)
            {
                countOut.readBase +=  ks1->seq.l;
                countOut.readBase += ks2->seq.l;
                countOut.readNum += 2;
            }

            // �ͷ��ڴ棬�ر��ļ�
            kseq_destroy(ks1);
            kseq_destroy(ks2);
            gzclose(gzfp1);
            gzclose(gzfp2);
        }
    }
    else if (inputFileName1.length() > 0 && inputFileName2.length() == 0) // ���˲���
    {
        // �����ļ���
        gzFile gzfp = gzopen(inputFileName1.c_str(), "rb");

        // ���ļ�
        if(!gzfp)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                    << inputFileName1
                    << "': No such file or directory.\n";
            exit(1);
        }
        else
        {
            kseq_t *ks;
            ks = kseq_init(gzfp);
        
            while(kseq_read(ks) >= 0)
            {
                countOut.readBase +=  ks1->seq.l;
                countOut.readNum++;
            }

            // �ͷ��ڴ棬�ر��ļ�
            kseq_destroy(ks);
            gzclose(gzfp);
        }
    }

    // ������
    if (outputFileName.empty()) // ���û��ָ������ļ��������ӡ����׼���
    {
        cout << "readBase:" << countOut.readBase << "\n" 
            << "readNum:" << countOut.readNum << "\n" 
            << "readLen:" << countOut.readBase/countOut.readNum 
            << endl;
    }
    else // ���浽�ļ�
    {
        // ����ļ���
        ofstream outputFile;
        outputFile.open(outputFileName, ios::out);

        // ���ļ�
        if(!outputFile)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'"
                << outputFileName 
                << "': No such file or directory." 
                << endl;
            outputFile.close();
            exit(1);
        }
        else
        {
            outputFile << "readBase:" << countOut.readBase << "\n" 
                    << "readNum:" << countOut.readNum << "\n" 
                    << "readLen:" << countOut.readBase/countOut.readNum 
                    << endl;
        }

        // �ر��ļ�
        outputFile.close();
    }
}