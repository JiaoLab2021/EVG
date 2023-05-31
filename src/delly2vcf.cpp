#include "../include/delly2vcf.hpp"

using namespace std;

int main_delly2vcf(int argc, char** argv)
{
    string dellyFilename;
    string refgenomeFilename;
    string outputFilename;

    // �������
    int c;
    while (true)
    {
        static const struct option long_options[] = 
        {
            {"vcf", required_argument, 0, 'v'},
            {"refgenome", required_argument, 0, 'r'},
            {"out", no_argument, 0, 'o'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "v:r:o:h", long_options, &option_index);

        if (c == -1)
            break;

        switch (c)
        {
        case 'v':
            dellyFilename = optarg;
            break;
        case 'r':
            refgenomeFilename = optarg;
            break;
        case 'o':
            outputFilename = optarg;
            break;
        case 'h':
        case '?':
            help_delly2vcf(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc <= 2)
    {
        help_delly2vcf(argv);
        return 1;
    }

    // ��������ļ�������
    if (outputFilename.empty())
    {
        outputFilename = split(dellyFilename, ".")[0] + ".convert.vcf";
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " 
         << "Running.\n";

    // ���������ֵ�
    map<string,string> seq_dict = DELLY2VCF::build_dict(& refgenomeFilename);

    // vcfת��
    DELLY2VCF::vcf_convert(& dellyFilename, seq_dict, & outputFilename);

    cerr << "[" << __func__ << "::" << getTime() << "] " 
         << "Done.\n";
    return 0;
}

// �����ĵ�
void help_delly2vcf(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -v -r [options]" << endl
       << "convert the format of delly to vcf" << endl
       << endl
       << "required arguments:" << endl
       << "    -v, --vcf       FILE    vcf file output by delly" << endl
       << "    -r, --refgenome FILE    refgenome genome (fasta)" << endl
       << endl
	   << "optional arguments:" << endl
       << "    -o, --out       FILE    output file name [xxx.convert.vcf]" << endl
       << endl
       << "    -h, --help              print this help document" << endl;
}


map<string,string> DELLY2VCF::build_dict(string * ref_file)
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Build index.\n";
    ifstream fasta_file;
    fasta_file.open(*ref_file, ios::in);

    string fastas;
    string fasta;

    if(!fasta_file.is_open())
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "'" << *ref_file << "': No such file or directory." << endl;
        exit(1);
    }
    else
        while(getline(fasta_file,fasta))
        {
            if(fasta.empty())
            {
                continue;
            }

            string::size_type idx;
            idx = fasta.find('>'); //���Ƿ����>���������ں�߼ӻ��з�����������ֱ�����
            if (idx == string::npos)
            {
                fastas += strip(fasta, '\n'); // ��stripɾ���������˵Ļ��з�
            }
            else
            {
                fastas += '\n' + strip(fasta, '\n') + '\t';
            }
        }

    fasta_file.close();
    fastas = strip(fastas,'\n');
    std::vector<string> fastas_split = split(fastas, "\n");

    size_t size = fastas_split.size();

    map<string,string> seq_dict;

    for (int i = 0; i < size; ++i) 
    {
        std::vector<string> fastas_split_split = split(fastas_split[i], "\t"); //��loc��seq��֣�\t
        string seq = fastas_split_split[1];
        std::vector<string> locs = split(fastas_split_split[0], " ");
        string loc = locs[0];
        loc.erase(std::remove(loc.begin(), loc.end(), '>'), loc.end()); //ɾ��loc��ͷ��>
        seq_dict[loc] = seq;
    }
    
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Successfully build index.\n";

    return seq_dict;
}

int DELLY2VCF::vcf_convert(string * dellyFilename, map<string,string> seq_dict, string * outputFilename)
{
    // ���ļ�
    ifstream dellyFile;
    dellyFile.open(*dellyFilename, ios::in);

    string informations;
    string out_txt;

    string::size_type idx;
    string::size_type idx_bnd;

    string chromosome;
    int start;
    string ref_seq;
    int end;
    string qry_seq;

    // ������ʽ��END���������滻����[]ʱ����������regex�⣬��������ɾ��
    std::regex end_reg("END=(\\d+)");
    std::regex reg1("\\[");
    std::regex reg2("\\]");

    smatch end_search; // Ѱ��informations�ֶ��е�END=
    
    srand((int)time(0));  // �����������,Ҳ���԰�0����NULL��

    if(!dellyFile.is_open())
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "'" << *dellyFilename << "': No such file or directory." << endl;
        exit(1);
    }
    else
    {
        while(getline(dellyFile,informations))
        {
            idx = informations.find("#");
            if (idx != string::npos)
            {
                std::regex reg_INFO("INFO\t.*");
                out_txt += regex_replace(string(informations), regex(reg_INFO), string("INFO")) + "\n";
            }
            else
            {
                // ɾ��������
                informations = regex_replace(string(informations), regex(reg1), string(""));
                informations = regex_replace(string(informations), regex(reg2), string(""));
                
                chromosome = split(informations, "\t")[0];
                start = stoi(split(informations, "\t")[1]);
                ref_seq = split(informations, "\t")[3];
                qry_seq = split(informations, "\t")[4];

                std::regex reg_ref_seq("\t" + ref_seq + "\t");
                std::regex reg_qry_seq("\t" + qry_seq + "\t");
                std::regex reg_PASS("PASS.*");

                
                // ��INFO�ֶ��е�ENDλ��
                if (regex_search(informations, end_search, end_reg))
                {
                    // end = stoi(end_search.format("$1")) - 1;
                    end = stoi(end_search.str(1)) - 1;
                }
                else
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: skip ->" << informations << endl;
                }

                // delly�ҵĽṹ�����еĻ�ܴ�ܴ󣬽��ñ�����С��
                if ((end - start) > 10000)
                {
                    end = start + 10000;
                }        

                // Ѱ���Ƿ���Translocation
                idx_bnd = qry_seq.find(chromosome);
                
                ref_seq = seq_dict[chromosome].substr(start-1, end-start+1);
                
                if (qry_seq == "<INV>")
                {
                    qry_seq = ref_seq;
                    reverse(qry_seq.begin(), qry_seq.end());
                }
                else if (qry_seq == "<DEL>")
                {
                    qry_seq = ref_seq[0];
                }
                else if (qry_seq == "<DUP>")
                {
                    qry_seq = ref_seq + ref_seq;
                }
                else if (qry_seq == "<INS>")
                {
                    int chromosome_len = seq_dict[chromosome].length();
                    int qry_start = rand()%(chromosome_len-10001)+0;
                    int qryLen = rand()%10000+50;
                    qry_seq = ref_seq[0] + seq_dict[chromosome].substr(qry_start, qryLen);
                }              
                else if (idx_bnd == string::npos) // Translocation.
                {
                    if (ref_seq.length() == 1)
                    {
                        cerr << "[" << __func__ << "::" << getTime() << "] " << "Skip this short variations: " << chromosome << " " << start << endl;
                        continue;
                    }
                    else
                    {
                        qry_seq = ref_seq[0];
                    }
                }
                else
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] " << "This variations could not be converted, please check code: " << informations << endl;
                }

                
                informations = regex_replace(string(informations), regex(reg_ref_seq), string("\t" + ref_seq + "\t"));
                informations = regex_replace(string(informations), regex(reg_qry_seq), string("\t" + qry_seq + "\t"));
                informations = regex_replace(string(informations), regex(reg_PASS), string("PASS\tEND=" + to_string(end)));

                out_txt += informations + "\n";
            }
        }
    }

    // ������
    ofstream outFile;

    outFile.open(* outputFilename, ios::app);
    outFile << out_txt;

    // �ͷ��ڴ�
    dellyFile.close();
    outFile.close();

    return 0;
}