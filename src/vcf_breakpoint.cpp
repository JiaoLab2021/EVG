// g++ vcf_breakpoint.cpp -o vcf_breakpoint -lz -O2
#include "../include/vcf_breakpoint.hpp"

using namespace std;

int main_breakpoint(int argc, char** argv)
{
    // �����ļ�
    string vcfFileName;
    
    // ����ļ�ǰ׺
    string prefix = "breakpoint";

    // �ϵ�����С
    int breakpointErrorSize = 1;

    // �������
    int c;
    while (true)
    {
        static const struct option long_options[] = {
            {"vcf", required_argument, 0, 'v'},

            {"prefix", required_argument, 0, 'p'},
            {"size", required_argument, 0, 's'},

            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "v:p:s:h", long_options, &option_index);

        if (c == -1)
            break;

        switch (c)
        {
        case 'v':
            vcfFileName = optarg;
            break;
        case 'p':
            prefix = optarg;
            break;;
        case 's':
            breakpointErrorSize = stoi(optarg);
            break;
        case 'h':
        case '?':
            help_breakpoint(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    // �������Ƿ���ȷ
    if (vcfFileName.empty())
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
             << "Parameter error.\n";
        help_breakpoint(argv);
        return 1;
    }
    
    cerr << "[" << __func__ << "::" << getTime() << "] " 
         << "Running.\n";
    
    BREAKPOINT::vcf_breakpoint(vcfFileName, prefix, breakpointErrorSize);
    
    cerr << "[" << __func__ << "::" << getTime() << "] " 
         << "Done.\n";

    return 0;
}

// �����ĵ�
void help_breakpoint(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -v [options]" << endl
       << "set error for breakpoint of variations" << endl
       << endl
       << "required arguments (Note: vcf files must be sorted):" << endl
       << "    -v, --vcf           FILE       vcf file to set breakpoint error" << endl
       << endl
	   << "optional arguments:" << endl
       << "    -p, --prefix        STRING     output prefix [breakpoint]" << endl
       << "    -s, --size          INT        breakpoint error size [1]" << endl
       << endl
       << "    -h, --help                     print this help document" << endl;
}


int BREAKPOINT::vcf_breakpoint(const string & vcfFileName, const string & prefix, const int & breakpointErrorSize)
{
    // �����ļ���
    gzFile gzfpI = gzopen(vcfFileName.c_str(), "rb");
    if(!gzfpI)
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "'" << vcfFileName << "': No such file or directory." << endl;
        exit(1);
    }

    // ��¼ת����Ľ��
    string outTxt;

    // ����ļ���
    string outFileName = prefix + ".vcf.gz";

    // ����ļ���
    gzFile outFile = gzopen(outFileName.c_str(), "wb");

    if(!outFile)
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'"
                << outFileName
                << "': No such file or directory." 
                << endl;
        exit(1);
    }

    // ��vcf�ļ�����ת��
    string information = ""; // ��ʱ�ַ���
    char line[1024]; // һ��ֻ��1024�ֽڵ�����

    while(gzgets(gzfpI, line, 1024))
    {
        information += line;
        if (information.find("\n") != string::npos) // һ�н���
        {
            if(line != nullptr && line[0] == '\0')
            {
                continue;
            }

            information = strip(information, '\n'); // ȥ�����з�

            if (information.find("#") != string::npos) // ����Ƿ���ע����
            {
                string headLine = information + "\n";
                gzwrite(outFile, headLine.c_str(), headLine.length());
            }
            else
            {
                vector<string> informationVec = split(information, "\t");

                // vcf�ĸ�����Ϣ
                string refChr = informationVec[0];
                long long int refStart = stol(informationVec[1]);

                string refSeq = informationVec[3];
                string qrySeqs = informationVec[4];
                vector<string> qrySeqsVec = split(qrySeqs, ",");
                long long int refLen = refSeq.length();
                long long int refEnd = refStart + refLen - 1;
                long long int svLen = qrySeqs.length() - refLen;


                // ����INFO��ߵ�end��Ϣ
                smatch endResult;
                regex endReg("END=(\\d+)");

                string::const_iterator iterStart = informationVec[7].begin();
                string::const_iterator iterEnd = informationVec[7].end();
                regex_search(iterStart, iterEnd, endResult, endReg);

                string endResultTmp = endResult[1];
                if (endResultTmp.empty())
                {
                    if (refEnd == stol(endResult[1]))
                    {
                        continue;
                    }
                    else
                    {
                        cerr << "[" << __func__ << "::" << getTime() << "] " <<  information << endl;
                    }
                }

                // ������
                refStart -= breakpointErrorSize;
                informationVec[1] = to_string(refStart);
                refEnd -= breakpointErrorSize;
                informationVec[7] = regex_replace(string(informationVec[7]), regex(endReg), string("END=" + to_string(refEnd)));

                // ��ӵ�����ַ�����
                for (size_t i = 0; i < informationVec.size(); i++)
                {
                    if (i > 0)
                    {
                        outTxt += "\t";
                    }
                    outTxt += informationVec[i];
                }
                outTxt += "\n";

                // ������
                if (outTxt.size() > 10000000) // ÿ10Mbд��һ��
                {
                    gzwrite(outFile, outTxt.c_str(), outTxt.length());

                    // ����ַ���
                    outTxt.clear();
                    string().swap(outTxt);
                }
            }

            // ����ַ���
            information.clear();
            string().swap(information);
        }
    }
    // �ر��ļ�
    gzclose(gzfpI);

    // ���д��һ��
    gzwrite(outFile, outTxt.c_str(), outTxt.length());

    // ����ַ���
    outTxt.clear();
    string().swap(outTxt);

    gzclose(outFile);

    return 0;
}