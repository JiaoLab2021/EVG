// g++ -c sample.cpp -o sample -lz -O2
#include "../include/convert.hpp"

using namespace std;

void help_convert(char** argv);

int main_convert(int argc, char** argv)
{
    // �����ļ�
    string fastaFilename;
    string outputFilename;
    
    // �������
    int c;
    while (true)
    {
        static const struct option long_options[] = {
            {"input", required_argument, 0, 'i'},
            {"out", required_argument, 0, 'o'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "i:o:h", long_options, &option_index);

        if (c == -1)
            break;

        switch (c)
        {
        case 'i':
            fastaFilename = optarg;
            break;
        case 'o':
            outputFilename = optarg;
            break;
        case 'h':
        case '?':
            help_convert(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc <= 1) {
        help_convert(argv);
        return 1;
    }

    if (fastaFilename.empty() && outputFilename.empty())
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
             << "parameter error: -i.\n";
        help_convert(argv);
        return 1;
    }

    // log
    cerr << "[" << __func__ << "::" << getTime() << "] " 
         << "Running.\n";

    convert::convert(&fastaFilename, &outputFilename);

    // log
    cerr << "[" << __func__ << "::" << getTime() << "] " 
         << "Done.\n";
    
    return 0;
}


// �����ĵ�
void help_convert(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -i -o [options]" << endl
       << "arrange the sequence of fasta files into one line" << endl
       << endl
       << "required arguments:" << endl
       << "    -i, --input   FILE    fasta file to be converted" << endl
       << "    -o, --out     FILE    output file name" << endl
       << endl
	   << "optional arguments:" << endl
       << "    -h, --help            print this help document" << endl;
}


int convert::convert(string * fastaFilename, string * outputFilename)
{
    ifstream fastaFile;
    fastaFile.open(*fastaFilename, ios::in);

    string informations;
    string outTxt;

    if(!fastaFile.is_open())
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
            << "'" << *fastaFilename << "': No such file or directory.\n";
            exit(1);
    }
    else
    {
        while(getline(fastaFile,informations))
        {
            if(informations.empty())
            {
                continue;
            }
            else
            {   
                string::size_type idx = informations.find(">");
                if ( idx != string::npos )
                {
                    outTxt += "\n" + informations + "\n";
                }
                else
                {
                    outTxt += informations;
                }
            }
        }
    }

    // �ر��ļ�
    fastaFile.close();

    ofstream outputFile;
    outputFile.open(*outputFilename);

    outputFile << strip(outTxt, '\n') << "\n";

    // �ر��ļ�
    outputFile.close();

    return 0;
}