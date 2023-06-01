#include "../include/count.hpp"

using namespace std;

int main_count(int argc, char** argv)
{
    // 输入文件
    string FileName1 = "";
    string FileName2 = "";
    string fastaFileName = "";

    // 输出文件
    string outputFileName = "";

    // 输入参数
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


// 打开fastq/a.gz文件
void count::fastq_a_count(
    string inputFileName1, 
    string inputFileName2, 
    countStruct & countOut, 
    const string & outputFileName
)
{
    if (inputFileName1.length() > 0 && inputFileName2.length() > 0) // 双端测序
    {
        // 输入文件流
        gzFile gzfp1 = gzopen(inputFileName1.c_str(), "rb");
        gzFile gzfp2 = gzopen(inputFileName2.c_str(), "rb");

        // 打开文件
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

            // 释放内存，关闭文件
            kseq_destroy(ks1);
            kseq_destroy(ks2);
            gzclose(gzfp1);
            gzclose(gzfp2);
        }
    }
    else if (inputFileName1.length() > 0 && inputFileName2.length() == 0) // 单端测序
    {
        // 输入文件流
        gzFile gzfp = gzopen(inputFileName1.c_str(), "rb");

        // 打开文件
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

            // 释放内存，关闭文件
            kseq_destroy(ks);
            gzclose(gzfp);
        }
    }

    // 输出结果
    if (outputFileName.empty()) // 如果没有指定输出文件名，则打印到标准输出
    {
        cout << "readBase:" << countOut.readBase << "\n" 
            << "readNum:" << countOut.readNum << "\n" 
            << "readLen:" << countOut.readBase/countOut.readNum 
            << endl;
    }
    else // 保存到文件
    {
        // 输出文件流
        ofstream outputFile;
        outputFile.open(outputFileName, ios::out);

        // 打开文件
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

        // 关闭文件
        outputFile.close();
    }
}