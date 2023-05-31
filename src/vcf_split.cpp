// g++ vcf_split.cpp -o vcf_split -lz -O2
#include <getopt.h>
#include <cstdlib>

#include "../include/vcf_split.hpp"


using namespace std;

int main_split(int argc, char** argv)
{
    // �����ļ�
    string vcfFileName;

    // ������ȡsnp��indelλ�õ�vcf�ļ�
    string baseVcfFileName;
    
    // ��ֵ�ģʽ
    string mode = "type";

    // snp��indel��sv�ľ���
    int length = 100;

    // ����ļ�ǰ׺
    string prefix = "split";

    // �Ƿ����FILTER�н��й���
    string filterBool = "true";

    // �������
    int c;
    while (true)
    {
        static const struct option long_options[] = {
            {"vcf", required_argument, 0, 'v'},

            {"mode", required_argument, 0, 'm'},
            {"length", required_argument, 0, 'l'},
            {"basevcf", required_argument, 0, 'V'},
            {"prefix", required_argument, 0, 'p'},
            {"filter", required_argument, 0, 'f'},

            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "v:m:l:V:p:f:h", long_options, &option_index);

        if (c == -1)
            break;

        switch (c)
        {
        case 'v':
            vcfFileName = optarg;
            break;
        case 'm':
            mode = optarg;
            break;
        case 'l':
            length = stoi(optarg);
            break;
        case 'V':
            baseVcfFileName = optarg;
            break;
        case 'p':
            prefix = optarg;
            break;;
        case 'f':
            filterBool = optarg;
            break;
        case 'h':
        case '?':
            help_split(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    // �������Ƿ���ȷ
    if (vcfFileName.empty() || (mode != "type" && mode != "number"))
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error.\n";
        help_split(argv);
        return 1;
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running.\n";
    
    if (mode == "type")
    {
        VCFSPLIT::vcf_split_type(vcfFileName, prefix, filterBool);
    }
    else if (mode == "number")
    {
        VCFSPLIT::vcf_split_number(vcfFileName, baseVcfFileName, length, prefix, filterBool);
    }
    else
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error.\n";
        help_split(argv);
        return 1;
    }
    
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done.\n";

    return 0;
}

// �����ĵ�
void help_split(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -v [options]" << endl
       << "split variation by type or number of nearby SNP+Indels" << endl
       << endl
       << "required arguments (Note: vcf files must be sorted):" << endl
       << "    -v, --vcf           FILE       vcf file to be converted" << endl
       << endl
	   << "optional arguments:" << endl
       << "    -m, --mode          STRING     split standard: (type -> SNP+Indel+Ins+Del / number -> number of nearby SNP+Indel) [type]" << endl
       << "    -l, --length        INT        distance of snp and indel from sv [100]" << endl      
       << "    -V, --basevcf       INT        vcf file used to extract snp and indel locations [same as -v]" << endl      
       << "    -p, --prefix        STRING     output prefix [split]" << endl
       << "    -f, --filter        BOOL       filter using the FILTER column (true/false) [true]" << endl
       << endl
       << "    -h, --help                     print this help document" << endl;
}


// ��GT���в��
vector<string> VCFSPLIT::gt_split(const string & gtTxt)
{
    // ��ʱgt�б�
    vector<string> gtVecTmp;

    if (gtTxt.find("/") != string::npos)
    {
        gtVecTmp = split(gtTxt, "/");
    }
    else if (gtTxt.find("|") != string::npos)
    {
        gtVecTmp = split(gtTxt, "|");
    }
    else
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: GT is not separated by '/' or '|' -> " << gtTxt << endl;
        exit(1);
    }

    // ���б��������
    sort(begin(gtVecTmp), end(gtVecTmp));

    return gtVecTmp;
}

int VCFSPLIT::vcf_split_type(const string & vcfFileName, const string & prefix, const string & filterBool)
{
    // �����ļ���
    gzFile gzfpI = gzopen(vcfFileName.c_str(), "rb");
    if(!gzfpI)
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "'" << vcfFileName << "': No such file or directory." << endl;
        exit(1);
    }

    // ��¼ÿ�ֱ����ע����
    string simpleSVs, invSVs, dupSVs, otherSVs;

    // ����ļ���
    string simpleSVsFile = prefix + ".snp.indel.ins.del.vcf.gz";
    string invSVsFile = prefix + ".inv.vcf.gz";
    string dupSVsFile = prefix + ".dup.vcf.gz";
    string otherSVsFile = prefix + ".other.vcf.gz";

    // ����ļ���
    gzFile gzfpSimple = gzopen(simpleSVsFile.c_str(), "wb");
    gzFile gzfpInv = gzopen(invSVsFile.c_str(), "wb");
    gzFile gzfpDup = gzopen(dupSVsFile.c_str(), "wb");
    gzFile gzfpOther = gzopen(otherSVsFile.c_str(), "wb");
    if(!gzfpSimple || !gzfpInv || !gzfpDup || !gzfpOther)
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
            << "'" << simpleSVsFile << "' or '" << invSVsFile << "' or '" << dupSVsFile << "' or '"  << otherSVsFile
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
                gzwrite(gzfpSimple, headLine.c_str(), headLine.length());
                gzwrite(gzfpInv, headLine.c_str(), headLine.length());
                gzwrite(gzfpDup, headLine.c_str(), headLine.length());
                gzwrite(gzfpOther, headLine.c_str(), headLine.length());
            }
            else
            {
                vector<string> informationVec = split(information, "\t");


                // ���һ����Ϣ
                vector<string> lastColVec = split(informationVec[informationVec.size()-1], ":");

                // ��FORMAT�ֶ���gt��λ��
                vector<string> formatVec = split(strip(informationVec[8], '\n'), ":");
                vector<string>::iterator gtItera = find(formatVec.begin(), formatVec.end(), "GT");

                string gt;
                int gtIndex = 0;
                if (gtItera != formatVec.end()) // FORMAT����GT
                {
                    // GT��index
                    gtIndex = distance(formatVec.begin(), gtItera);
                    gt = strip(lastColVec[gtIndex], '\n');
                }
                else // û�еĻ��˳�����
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] "
                            << "Error: GT not in FORMAT column -> " 
                            << " (" << strip(information, '\n') << ")"
                            << endl;
                    exit(1);
                }

                if (gt == "0/0" || gt == "." || gt == "./." || gt == ".|." || gt == "0|0") // �����(0/0, .)��ʽ����ֱ����������������
                {
                    // ����ַ���
                    information.clear();
                    string().swap(information);
                    continue;
                }

                // ����FILTER�ֶν��й���
                string filter = informationVec[6];
                if (filter != "PASS" && filterBool == "true") // �������(PASS)��ֱ����������������ˣ��������ж�
                {
                    // ����ַ���
                    information.clear();
                    string().swap(information);
                    continue;
                }

                // ����ķ��ͽ��
                vector<string> evaluate_gt;
                if (filterBool == "true")
                {
                    evaluate_gt = gt_split(gt);
                }
                else // ��������ˣ�������Ҳ���жϣ���Ϊgramtools��GT��ʽ��1��û�зָ���
                {
                    evaluate_gt = {"1", "1"};
                }
                

                uint32_t refLen, qryLen;

                string svType = informationVec[2];
                string refSeq = informationVec[3];
                string qrySeq = informationVec[4];

                if (qrySeq.find(">") != string::npos) // graphtyper����Ľ��
                {
                    if (qrySeq.find("<INS") != string::npos) // �ַ����а���������ֶΣ�graphtyper��<INS:SVSIZE=97:BREAKPOINT1>
                    {   
                        // ���qry��������û�г�����Ϣ��û�еĻ��˳����벢����
                        if (qrySeq.find("=") == string::npos)
                        {
                            cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: no length information in the fourth column -> " << information << endl; 
                            exit(1);
                        }

                        refLen = refSeq.size();
                        qryLen = stoi(split(split(qrySeq, "=")[1], ":")[0]);
                    }
                    else if (qrySeq.find("<DEL") != string::npos) // �ַ����а���ȱʧ���ֶΣ�graphtyper�� <DEL:SVSIZE=233:COVERAGE>
                    {
                        // ���qry��������û�г�����Ϣ��û�еĻ��˳����벢����
                        if (qrySeq.find("=") == string::npos)
                        {
                            cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: no length information in the fourth column -> " << information << endl; 
                            exit(1);
                        }

                        refLen = stoi(split(split(qrySeq, "=")[1], ":")[0]);
                        qryLen = refSeq.size();
                    }
                    else if (qrySeq.find("<DUP") != string::npos) // �ַ����а����ظ����ֶΣ�graphtyper��<DUP:SVSIZE=2806:BREAKPOINT1>
                    {
                        // ���qry��������û�г�����Ϣ��û�еĻ��˳����벢����
                        if (qrySeq.find("=") == string::npos)
                        {
                            cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: no length information in the fourth column -> " << information << endl; 
                            exit(1);
                        }

                        refLen = stoi(split(split(qrySeq, "=")[1], ":")[0]);
                        qryLen = refLen * 2;
                    }
                    else // ����ʶ���ֶ�
                    {
                        cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: unknown string -> ref_seq:" << refSeq << " qey_seq:" << qrySeq << endl;
                        refLen = refSeq.size();

                        // ���qry��������û�г�����Ϣ��û�еĻ��˳����벢����
                        if (qrySeq.find("=") == string::npos)
                        {
                            cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: no length information in the fourth column -> " << information << endl; 
                            exit(1);
                        }
                        
                        qryLen = stoi(split(split(qrySeq, "=")[1], ":")[0]);
                    }
                }
                else // �����ı���
                {
                    refLen = refSeq.size();
                    
                    int maxGtNum = stoi(*max_element(evaluate_gt.begin(), evaluate_gt.end()));
                    if (maxGtNum > split(informationVec[4], ",").size())
                    {
                        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: genotyping and query sequence numbers do not match -> " << information << endl;
                        exit(1);
                    }
                    else
                    {
                        qrySeq = split(informationVec[4], ",")[maxGtNum-1]; // GT�ż�ȥ1�����е�����
                    }

                    qryLen = qrySeq.size();
                }

                // �ж�bayestyper�Ľ��Ϊduplication����ref_lenΪ1-2��qry_seqȴ�ܳ�ʱ
                if ((svType.find("Duplication") != string::npos) && (refLen <= 2)) // bayestyper���duplication��ɲ��룬�������¼��㳤��
                {
                    refLen = qryLen;
                    qryLen *= 2;
                }

                // �жϱ������ͣ������浽�ַ�����
                if ((refLen < 50 && qryLen < 50) || (refLen >= 50 && qryLen < 50) || (refLen < 50 && qryLen >= 50)) // snp+indel+del+ins
                {
                    simpleSVs += information + "\n";
                }
                else if (refLen >= 50 && qryLen >= 50)
                {
                    if (refLen*2 <= qryLen+10 && refLen*2 >= qryLen-10) // dup
                    {
                        dupSVs += information + "\n";
                    }
                    else if (refLen <= qryLen+10 && refLen >= qryLen-10) // inv
                    {
                        invSVs += information + "\n";
                    }
                    else
                    {
                        otherSVs += information + "\n";
                    }
                }
                else
                {
                    otherSVs += information + "\n";
                }

                // save result
                if (simpleSVs.size() > 10000000) // ÿ10Mbд��һ��
                {
                    gzwrite(gzfpSimple, simpleSVs.c_str(), simpleSVs.length());
                    gzwrite(gzfpInv, invSVs.c_str(), invSVs.length());
                    gzwrite(gzfpDup, dupSVs.c_str(), dupSVs.length());
                    gzwrite(gzfpOther, otherSVs.c_str(), otherSVs.length());

                    // ����ַ���
                    simpleSVs.clear();
                    invSVs.clear();
                    dupSVs.clear();
                    otherSVs.clear();
                    string().swap(simpleSVs);
                    string().swap(invSVs);
                    string().swap(dupSVs);
                    string().swap(otherSVs);
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
    gzwrite(gzfpSimple, simpleSVs.c_str(), simpleSVs.length());
    gzwrite(gzfpInv, invSVs.c_str(), invSVs.length());
    gzwrite(gzfpDup, dupSVs.c_str(), dupSVs.length());
    gzwrite(gzfpOther, otherSVs.c_str(), otherSVs.length());

    // ����ַ���
    simpleSVs.clear();
    invSVs.clear();
    dupSVs.clear();
    otherSVs.clear();
    string().swap(simpleSVs);
    string().swap(invSVs);
    string().swap(dupSVs);
    string().swap(otherSVs);


    gzclose(gzfpSimple);
    gzclose(gzfpInv);
    gzclose(gzfpDup);
    gzclose(gzfpOther);

    return 0;
}

int VCFSPLIT::vcf_split_number(
    const string & vcfFileName, 
    string & baseVcfFileName, 
    const int & length, 
    const string & prefix, 
    const string & filterBool
)
{
    if (baseVcfFileName.empty())
    {
        baseVcfFileName = vcfFileName;
    }
    
    // �����ļ���
    gzFile gzfpI1 = gzopen(baseVcfFileName.c_str(), "rb");
    if(!gzfpI1)
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'" << baseVcfFileName 
                << "': No such file or directory." << endl;
        exit(1);
    }

    // ��¼sv��ֵĽ��
    string SVs0, SVs1, SVs2, SVs4, SVs6, SVs8, SVsMore;

    // ����ļ���
    string SVs0FileName = prefix + ".0.vcf.gz";
    string SVs1FileName = prefix + ".1.vcf.gz";
    string SVs2FileName = prefix + ".2.vcf.gz";
    string SVs4FileName = prefix + ".4.vcf.gz";
    string SVs6FileName = prefix + ".6.vcf.gz";
    string SVs8FileName = prefix + ".8.vcf.gz";
    string SVsMoreFileName = prefix + ".more.vcf.gz";

    // ����ļ���
    gzFile SVs0File = gzopen(SVs0FileName.c_str(), "wb");
    gzFile SVs1File = gzopen(SVs1FileName.c_str(), "wb");
    gzFile SVs2File = gzopen(SVs2FileName.c_str(), "wb");
    gzFile SVs4File = gzopen(SVs4FileName.c_str(), "wb");
    gzFile SVs6File = gzopen(SVs6FileName.c_str(), "wb");
    gzFile SVs8File = gzopen(SVs8FileName.c_str(), "wb");
    gzFile SVsMoreFile = gzopen(SVsMoreFileName.c_str(), "wb");

    if(!SVs0File || !SVs1File || !SVs2File || !SVs4File || !SVs6File || !SVs8File || !SVsMoreFile)
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
            << "'" << SVs0FileName << "' or '" 
            << SVs1FileName << "' or '" 
            << SVs2FileName << "' or '" 
            << SVs4FileName << "' or '" 
            << SVs6FileName << "' or '" 
            << SVs8FileName << "' or '" 
            << SVsMoreFileName << "' or '" 
            << "': No such file or directory." 
            << endl;
        exit(1);
    }


    // snp+indel������
    // ��vcf�ļ���������
    string information = ""; // ��ʱ�ַ���

    // ������ϣ�� ��snp+indel��
    unordered_map<string, snpIndelInfo> snpIndelInfoMap; // unordered_map<chromosome, snpIndelInfo>

    char line[1024]; // һ��ֻ��1024�ֽڵ�����

    while(gzgets(gzfpI1, line, 1024))
    {
        information += line;
        if (information.find("\n") != string::npos) // һ�н���
        {
            if(line != nullptr && line[0] == '\0')
            {
                continue;
            }

            information = strip(information, '\n'); // ȥ�����з�

            if (information.find("#") == string::npos)
            {
                vector<string> informationVec = split(information, "\t");

                // ���һ����Ϣ
                vector<string> lastColVec = split(informationVec[informationVec.size()-1], ":");

                // ��FORMAT�ֶ���gt��λ��
                vector<string> formatVec = split(strip(informationVec[8], '\n'), ":");
                vector<string>::iterator gtItera = find(formatVec.begin(), formatVec.end(), "GT");

                string gt;
                int gtIndex = 0;
                if (gtItera != formatVec.end()) // FORMAT����GT
                {
                    // GT��index
                    gtIndex = distance(formatVec.begin(), gtItera);
                    gt = strip(lastColVec[gtIndex], '\n');
                }
                else // û�еĻ��˳�����
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "Error: GT not in FORMAT column -> " 
                        << " (" << strip(information, '\n') << ")"
                        << endl;
                    exit(1);
                }

                if (gt == "0/0" || gt == "." || gt == "./." || gt == ".|." || gt == "0|0") // �����(0/0, .)��ʽ����ֱ����������������
                {
                    // ����ַ���
                    information.clear();
                    string().swap(information);
                    continue;
                }

                // ����FILTER�ֶν��й���
                string filter = informationVec[6];
                if (filter != "PASS" && filterBool == "true") // �������(PASS)��ֱ����������������ˣ��������ж�
                {
                    // ����ַ���
                    information.clear();
                    string().swap(information);
                    continue;
                }

                // ����ķ��ͽ��
                vector<string> evaluate_gt;
                if (filterBool == "true")
                {
                    evaluate_gt = gt_split(gt);
                }
                else // ��������ˣ�������Ҳ���жϣ���Ϊgramtools��GT��ʽ��1��û�зָ���
                {
                    evaluate_gt = {"1", "1"};
                }
                
                uint32_t refLen, qryLen;

                string chromosome = informationVec[0];
                uint32_t refStart = stol(informationVec[1]);
                string svType = informationVec[2];
                string refSeq = informationVec[3];
                string qrySeq = informationVec[4];

                if (qrySeq.find(">") != string::npos) // graphtyper����Ľ��
                {
                    if (qrySeq.find("<INS") != string::npos) // �ַ����а���������ֶΣ�graphtyper��<INS:SVSIZE=97:BREAKPOINT1>
                    {   
                        // ���qry��������û�г�����Ϣ��û�еĻ��˳����벢����
                        if (qrySeq.find("=") == string::npos)
                        {
                            cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: no length information in the fourth column -> " << information << endl; 
                            exit(1);
                        }

                        refLen = refSeq.size();
                        qryLen = stoi(split(split(qrySeq, "=")[1], ":")[0]);
                    }
                    else if (qrySeq.find("<DEL") != string::npos) // �ַ����а���ȱʧ���ֶΣ�graphtyper�� <DEL:SVSIZE=233:COVERAGE>
                    {
                        // ���qry��������û�г�����Ϣ��û�еĻ��˳����벢����
                        if (qrySeq.find("=") == string::npos)
                        {
                            cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: no length information in the fourth column -> " << information << endl; 
                            exit(1);
                        }

                        refLen = stoi(split(split(qrySeq, "=")[1], ":")[0]);
                        qryLen = refSeq.size();
                    }
                    else if (qrySeq.find("<DUP") != string::npos) // �ַ����а����ظ����ֶΣ�graphtyper��<DUP:SVSIZE=2806:BREAKPOINT1>
                    {
                        // ���qry��������û�г�����Ϣ��û�еĻ��˳����벢����
                        if (qrySeq.find("=") == string::npos)
                        {
                            cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: no length information in the fourth column -> " << information << endl; 
                            exit(1);
                        }

                        refLen = stoi(split(split(qrySeq, "=")[1], ":")[0]);
                        qryLen = refLen * 2;
                    }
                    else // ����ʶ���ֶ�
                    {
                        cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: unknown string -> ref_seq:" << refSeq << " qey_seq:" << qrySeq << endl;
                        refLen = refSeq.size();

                        // ���qry��������û�г�����Ϣ��û�еĻ��˳����벢����
                        if (qrySeq.find("=") == string::npos)
                        {
                            cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: no length information in the fourth column -> " << information << endl; 
                            exit(1);
                        }
                        
                        qryLen = stoi(split(split(qrySeq, "=")[1], ":")[0]);
                    }
                }
                else // �����ı���
                {
                    refLen = refSeq.size();
                    
                    int maxGtNum = stoi(*max_element(evaluate_gt.begin(), evaluate_gt.end()));
                    if (maxGtNum > split(informationVec[4], ",").size())
                    {
                        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: genotyping and query sequence numbers do not match -> " << information << endl;
                        exit(1);
                    }
                    else
                    {
                        qrySeq = split(informationVec[4], ",")[maxGtNum-1]; // GT�ż�ȥ1�����е�����
                    }

                    qryLen = qrySeq.size();
                }

                // �ж�bayestyper�Ľ��Ϊduplication����ref_lenΪ1-2��qry_seqȴ�ܳ�ʱ
                if ((svType.find("Duplication") != string::npos) && (refLen <= 2)) // bayestyper���duplication��ɲ��룬�������¼��㳤��
                {
                    refLen = qryLen;
                    qryLen *= 2;
                }

                uint32_t refEnd;
                refEnd = refStart + refLen - 1;

                // �жϱ������ͣ������浽��ϣ����
                if ((refLen < 50 && qryLen < 50)) // ����snp+indel����ʼ����ֹλ��
                {
                    snpIndelInfoMap[chromosome].refStartVec.push_back(refStart);
                    snpIndelInfoMap[chromosome].refEndVec.push_back(refEnd);
                }
            }

            // ����ַ���
            information.clear();
            string().swap(information);
        }
    }
    // �ر��ļ�
    gzclose(gzfpI1);

    // �����ļ���
    gzFile gzfpI2 = gzopen(vcfFileName.c_str(), "rb");
    if(!gzfpI2)
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'" << vcfFileName 
                << "': No such file or directory." << endl;
        exit(1);
    }

    // sv������
    // ��vcf�ļ���������
    information = ""; // ��ʱ�ַ���

    // ������ϣ�� ���ṹ���죩
    map<string, map<uint32_t, svInfo>> svInfoMap; // map<chromosome, map<refStart, svInfo>>

    while(gzgets(gzfpI2, line, 1024))
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
                gzwrite(SVs0File, headLine.c_str(), headLine.length());
                gzwrite(SVs1File, headLine.c_str(), headLine.length());
                gzwrite(SVs2File, headLine.c_str(), headLine.length());
                gzwrite(SVs4File, headLine.c_str(), headLine.length());
                gzwrite(SVs6File, headLine.c_str(), headLine.length());
                gzwrite(SVs8File, headLine.c_str(), headLine.length());
                gzwrite(SVsMoreFile, headLine.c_str(), headLine.length());
            }
            else
            {
                vector<string> informationVec = split(information, "\t");

                // ���һ����Ϣ
                vector<string> lastColVec = split(informationVec[informationVec.size()-1], ":");

                // ��FORMAT�ֶ���gt��λ��
                vector<string> formatVec = split(strip(informationVec[8], '\n'), ":");
                vector<string>::iterator gtItera = find(formatVec.begin(), formatVec.end(), "GT");

                string gt;
                int gtIndex = 0;
                if (gtItera != formatVec.end()) // FORMAT����GT
                {
                    // GT��index
                    gtIndex = distance(formatVec.begin(), gtItera);
                    gt = strip(lastColVec[gtIndex], '\n');
                }
                else // û�еĻ��˳�����
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "Error: GT not in FORMAT column -> " 
                        << " (" << strip(information, '\n') << ")"
                        << endl;
                    exit(1);
                }

                if (gt == "0/0" || gt == "." || gt == "./." || gt == ".|." || gt == "0|0") // �����(0/0, .)��ʽ����ֱ����������������
                {
                    // ����ַ���
                    information.clear();
                    string().swap(information);
                    continue;
                }

                // ����FILTER�ֶν��й���
                string filter = informationVec[6];
                if (filter != "PASS" && filterBool == "true") // �������(PASS)��ֱ����������������ˣ��������ж�
                {
                    // ����ַ���
                    information.clear();
                    string().swap(information);
                    continue;
                }

                // ����ķ��ͽ��
                vector<string> evaluate_gt;
                if (filterBool == "true")
                {
                    evaluate_gt = gt_split(gt);
                }
                else // ��������ˣ�������Ҳ���жϣ���Ϊgramtools��GT��ʽ��1��û�зָ���
                {
                    evaluate_gt = {"1", "1"};
                }
                
                uint32_t refLen, qryLen;

                string chromosome = informationVec[0];
                uint32_t refStart = stol(informationVec[1]);
                string svType = informationVec[2];
                string refSeq = informationVec[3];
                string qrySeq = informationVec[4];

                if (qrySeq.find(">") != string::npos) // graphtyper����Ľ��
                {
                    if (qrySeq.find("<INS") != string::npos) // �ַ����а���������ֶΣ�graphtyper��<INS:SVSIZE=97:BREAKPOINT1>
                    {   
                        // ���qry��������û�г�����Ϣ��û�еĻ��˳����벢����
                        if (qrySeq.find("=") == string::npos)
                        {
                            cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: no length information in the fourth column -> " << information << endl; 
                            exit(1);
                        }

                        refLen = refSeq.size();
                        qryLen = stoi(split(split(qrySeq, "=")[1], ":")[0]);
                    }
                    else if (qrySeq.find("<DEL") != string::npos) // �ַ����а���ȱʧ���ֶΣ�graphtyper�� <DEL:SVSIZE=233:COVERAGE>
                    {
                        // ���qry��������û�г�����Ϣ��û�еĻ��˳����벢����
                        if (qrySeq.find("=") == string::npos)
                        {
                            cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: no length information in the fourth column -> " << information << endl; 
                            exit(1);
                        }

                        refLen = stoi(split(split(qrySeq, "=")[1], ":")[0]);
                        qryLen = refSeq.size();
                    }
                    else if (qrySeq.find("<DUP") != string::npos) // �ַ����а����ظ����ֶΣ�graphtyper��<DUP:SVSIZE=2806:BREAKPOINT1>
                    {
                        // ���qry��������û�г�����Ϣ��û�еĻ��˳����벢����
                        if (qrySeq.find("=") == string::npos)
                        {
                            cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: no length information in the fourth column -> " << information << endl; 
                            exit(1);
                        }

                        refLen = stoi(split(split(qrySeq, "=")[1], ":")[0]);
                        qryLen = refLen * 2;
                    }
                    else // ����ʶ���ֶ�
                    {
                        cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: unknown string -> ref_seq:" << refSeq << " qey_seq:" << qrySeq << endl;
                        refLen = refSeq.size();

                        // ���qry��������û�г�����Ϣ��û�еĻ��˳����벢����
                        if (qrySeq.find("=") == string::npos)
                        {
                            cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: no length information in the fourth column -> " << information << endl; 
                            exit(1);
                        }
                        
                        qryLen = stoi(split(split(qrySeq, "=")[1], ":")[0]);
                    }
                }
                else // �����ı���
                {
                    refLen = refSeq.size();
                    
                    int maxGtNum = stoi(*max_element(evaluate_gt.begin(), evaluate_gt.end()));
                    if (maxGtNum > split(informationVec[4], ",").size())
                    {
                        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: genotyping and query sequence numbers do not match -> " << information << endl;
                        exit(1);
                    }
                    else
                    {
                        qrySeq = split(informationVec[4], ",")[maxGtNum-1]; // GT�ż�ȥ1�����е�����
                    }

                    qryLen = qrySeq.size();
                }

                // �ж�bayestyper�Ľ��Ϊduplication����ref_lenΪ1-2��qry_seqȴ�ܳ�ʱ
                if ((svType.find("Duplication") != string::npos) && (refLen <= 2)) // bayestyper���duplication��ɲ��룬�������¼��㳤��
                {
                    refLen = qryLen;
                    qryLen *= 2;
                }

                uint32_t refEnd;
                refEnd = refStart + refLen - 1;

                // �жϱ������ͣ������浽��ϣ����
                if (refLen >= 50 || qryLen >= 50) // ����sv����ֹ�ͱ�����Ϣ
                {
                    svInfoMap[chromosome][refStart].information = information;
                    svInfoMap[chromosome][refStart].refEnd = refEnd;
                }
            }

            // ����ַ���
            information.clear();
            string().swap(information);
        }
    }
    // �ر��ļ�
    gzclose(gzfpI2);


    // ��sv������snp��indel���������в��
    // ��sv��ϣ�����ѭ��
    for (auto iter1 : svInfoMap)
    {
        string chromosome = iter1.first;

        // snp+indel��ʼ����ֹ���б�
        vector<uint32_t> snpIndelRefStartVec = snpIndelInfoMap[chromosome].refStartVec;
        vector<uint32_t> snpIndelRefEndVec = snpIndelInfoMap[chromosome].refEndVec;

        for (auto iter2 : iter1.second)
        {
            uint32_t refStart = iter2.first;
            uint32_t refEnd = iter2.second.refEnd;
            string information = iter2.second.information;

            // ���ֲ��ҷ�����ʼ����ֹ200bp�ڵ�snp��indel����
            // Ѱ����������
            int32_t startLeftIdxTmp = search_Binary_right(snpIndelRefEndVec, refStart - length);
            int32_t startRightIdxTmp = search_Binary_left(snpIndelRefEndVec, refStart);
            int32_t endLeftIdxTmp = search_Binary_right(snpIndelRefStartVec, refEnd);
            int32_t endRightIdxTmp = search_Binary_left(snpIndelRefStartVec, refEnd + length);

            // sv����snp+indel������
            int32_t leftSnpIndelNum = startRightIdxTmp - startLeftIdxTmp;
            if (leftSnpIndelNum < 0) // û��
            {
                leftSnpIndelNum = 0;
            }
            else // ��һ������
            {
                leftSnpIndelNum++;
            }
            
            int32_t rightSnpIndelNum = endRightIdxTmp - endLeftIdxTmp;
            if (rightSnpIndelNum < 0) // û��
            {
                rightSnpIndelNum = 0;
            }
            else // ��һ������
            {
                rightSnpIndelNum++;
            }

            int32_t snpIndelNum = leftSnpIndelNum + rightSnpIndelNum;
            if (snpIndelNum == 0)
            {
                SVs0 += information + "\n";
            }
            else if (snpIndelNum == 1)
            {
                SVs1 += information + "\n";
            }
            else if (snpIndelNum == 2)
            {
                SVs2 += information + "\n";
            }
            else if (snpIndelNum == 4)
            {
                SVs4 += information + "\n";
            }
            else if (snpIndelNum == 6)
            {
                SVs6 += information + "\n";
            }
            else if (snpIndelNum == 8)
            {
                SVs8 += information + "\n";
            }
            else
            {
                SVsMore += information + "\n";
            }
            
            // save result
            if (SVs0.size() > 10000000) // ÿ10Mbд��һ��
            {
                gzwrite(SVs0File, SVs0.c_str(), SVs0.length());
                gzwrite(SVs1File, SVs1.c_str(), SVs1.length());
                gzwrite(SVs2File, SVs2.c_str(), SVs2.length());
                gzwrite(SVs4File, SVs4.c_str(), SVs4.length());
                gzwrite(SVs6File, SVs6.c_str(), SVs6.length());
                gzwrite(SVs8File, SVs8.c_str(), SVs8.length());
                gzwrite(SVsMoreFile, SVsMore.c_str(), SVsMore.length());

                // ����ַ���
                SVs0.clear();
                string().swap(SVs0);
                SVs1.clear();
                string().swap(SVs1);
                SVs2.clear();
                string().swap(SVs2);
                SVs4.clear();
                string().swap(SVs4);
                SVs6.clear();
                string().swap(SVs6);
                SVs8.clear();
                string().swap(SVs8);
                SVsMore.clear();
                string().swap(SVsMore);
            }
        }
    }
    
    gzwrite(SVs0File, SVs0.c_str(), SVs0.length());
    gzwrite(SVs1File, SVs1.c_str(), SVs1.length());
    gzwrite(SVs2File, SVs2.c_str(), SVs2.length());
    gzwrite(SVs4File, SVs4.c_str(), SVs4.length());
    gzwrite(SVs6File, SVs6.c_str(), SVs6.length());
    gzwrite(SVs8File, SVs8.c_str(), SVs8.length());
    gzwrite(SVsMoreFile, SVsMore.c_str(), SVsMore.length());

    // ����ַ���
    SVs0.clear();
    string().swap(SVs0);
    SVs1.clear();
    string().swap(SVs1);
    SVs2.clear();
    string().swap(SVs2);
    SVs4.clear();
    string().swap(SVs4);
    SVs6.clear();
    string().swap(SVs6);
    SVs8.clear();
    string().swap(SVs8);
    SVsMore.clear();
    string().swap(SVsMore);


    gzclose(SVs0File);
    gzclose(SVs1File);
    gzclose(SVs2File);
    gzclose(SVs4File);
    gzclose(SVs6File);
    gzclose(SVs8File);
    gzclose(SVsMoreFile);

    return 0;
}