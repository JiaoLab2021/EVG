// g++ vcf_win_count_thread.cpp -std=c++11 -o vcf_win_count_thread -lpthread -O2
#include "../include/vcf_win.hpp"

using namespace std;

int main_win(int argc, char** argv)
{
    string vcfFilename;
    string chrLenFilename;
    string waysFilename;
    int windowSize = 200000;
    int stepSize = 100000;
    string mode = "number";

    // �������
    int c;
    while (true)
    {
        static const struct option long_options[] = {
            {"vcf", required_argument, 0, 'v'},
            {"length", required_argument, 0, 'l'},         
            {"group", required_argument, 0, 'g'},
            {"windowSize", required_argument, 0, 'w'},
            {"stepSize", required_argument, 0, 's'},
            {"mode", required_argument, 0, 'm'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "v:l:g:w:s:m:h", long_options, &option_index);

        if (c == -1)
            break;

        switch (c)
        {
        case 'v':
            vcfFilename = optarg;
            break;
        case 'l':
            chrLenFilename = optarg;
            break;
        case 'g':
            waysFilename = optarg;
            break;
        case 'w':
            windowSize = stoi(optarg);
            break;
        case 's':
            stepSize = stoi(optarg);
            break;
        case 'm':
            mode = optarg;
            break;
        case 'h':
        case '?':
            help_win(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc <= 2) {
        help_win(argv);
        return 1;
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running.\n";

    // �����ֵ�
    map<string,vector<long int>> winStartMap, winEndMap;
    tie(winStartMap, winEndMap) = VCFWIN::step_count(chrLenFilename, windowSize, stepSize);

    // ways�ķ�����Ϣ
    map<string,string> waysGroupMap = VCFWIN::ways_group(waysFilename);

    // ���ж�����Ⱦɫ��
    vector<string> chrNameVector = VCFWIN::count_chrName(chrLenFilename);
    
    // ���vcf
    vector<string> vcfFileVector = VCFWIN::split_file(vcfFilename, chrNameVector);

    // �����߳�
    string inputVcf;
    int threadNum = vcfFileVector.size();
    std::thread threads[threadNum];
    for (int i = 0; i < threadNum; i++)
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Thread start: ";
        cerr << i;
        cerr << " " + vcfFileVector[i] << endl;
        threads[i] = std::thread(VCFWIN::window_len_count, std::ref(vcfFileVector[i]), std::ref(winStartMap), std::ref(winEndMap), std::ref(windowSize), std::ref(waysGroupMap), mode);
    }
    // �߳�����
    for (auto& t: threads)
    {
        t.join();
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done.\n";
    return 0;
}

// �����ĵ�
void help_win(char** argv)

{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -v -l [options]" << endl
       << "variation statistics within a sliding window" << endl
       << endl
       << "required parameter:" << endl
       << "    -v, --vcf         FILE        vcf file merged by vcftools" << endl
       << "    -l, --length      FILE        chromosome length file (chr\tlength)" << endl
       << endl
       << "optional parameter:" << endl
       << "    -w, --windowSize  INT         window [200000]" << endl
       << "    -s, --stepSize    INT         step [100000]" << endl
       << "    -g, --group       FILE        group information for every way (way\tgroup)" << endl
       << "    -m, --mode        STRING      which indicator to count (length/number) [number]" << endl
       << endl
       << "    -h, --help                    print this help document" << endl;
}


// ͳ��Ⱦɫ������������
vector<string> VCFWIN::count_chrName(string chrLenFilename)
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Count the number and names of chromosomes.\n";
    // ���ж�����Ⱦɫ��
    string line;
    vector<string> chrNameVector;
    ifstream chrLenFile;
    chrLenFile.open(chrLenFilename, ios::in);

    if(!chrLenFile.is_open())
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
            << "'" << chrLenFilename << "': No such file or directory." << endl;
        exit(1);
    }
    else
    {
        while(getline(chrLenFile,line))
        {
            if (line.empty())
            {
                continue;
            }
            else
            {
                vector<string> lineList = split(line, "\t");
                if (lineList.size() < 2)
                {
                    lineList = split(line, " ");
                }

                if (lineList.size() < 2)
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] " 
                        << "Format error \'" << chrLenFilename << "\'" << ": chrName\tchrLen\n";
                    exit(1);
                }
                string chrName = lineList[0]; 
                chrNameVector.push_back(chrName);          
            }
        }
    }
    // �ͷ��ڴ�
    chrLenFile.close();
    return chrNameVector;
}

// ��Ⱦɫ����vcf�ļ�
vector<string> VCFWIN::split_file(string vcfFilename, vector<string> chrNameVector)
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Split vcf file.\n";
    vector<string> vcfFileVector;
    for (int i = 0; i < chrNameVector.size(); i++)
    {
        string vcfFile;  // ����ļ���
        string cmd_str;  // cmd

        if(vcfFilename.find(".gz") == string::npos && vcfFilename.find(".GZ") == string::npos)  // ��ѹ���ļ�
        {
            vcfFile = chrNameVector[i] + ".vcf";  // ���Ϊ��ͨ�ļ�
            cmd_str = "awk '{if($1==\"" + chrNameVector[i] + "\") print$0; else if($0~/#/) print$0;}' " + vcfFilename + " > " + vcfFile;
        }
        else  // ѹ���ļ�
        {
            vcfFile = chrNameVector[i] + ".vcf.gz";  // ���Ϊѹ���ļ�
            cmd_str = "zcat " + vcfFilename + " | awk '{if($1==\"" + chrNameVector[i] + "\") print$0; else if($0~/#/) print$0;}' | gzip > " + vcfFile;
        }

        vcfFileVector.push_back(vcfFile);  // �������ļ���vector��
        
        int length_cmd_str = cmd_str.length();
        char cmd_char[length_cmd_str];
        strcpy(cmd_char, cmd_str.c_str());
        system(cmd_char); 
    }
    
    return vcfFileVector;
}

// �����ֵ�
pair<map<string,vector<long int>>, map<string,vector<long int>>> VCFWIN::step_count(string chrLenFilename, int windowSize, int stepSize)
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Build step index.\n";
    // �����ֵ�
    map<string,vector<long int>> winStartMap;
    map<string,vector<long int>> winEndMap;

    // Ⱦɫ�峤�������ļ�
    string line;
    ifstream chrLenFile;
    chrLenFile.open(chrLenFilename, ios::in);

    if(!chrLenFile.is_open())
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
            << "'" << chrLenFilename << "': No such file or directory." << endl;
        exit(1);
        exit(1);
    }
    else
    {
        while(getline(chrLenFile,line))
        {
            if (line.empty())
            {
                continue;
            }
            else
            {
                vector<string> lineList = split(line, "\t");

                if (lineList.size() < 2)
                {
                    lineList = split(line, " ");
                }

                if (lineList.size() < 2)
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] " << "Format error \'" << chrLenFilename << "\': chromosome\tlength\n";
                    exit(1);
                }

                string chrName = lineList[0];
                int chrLen = stoi(lineList[1]);

                // ����
                double stepsNumber = double(chrLen)/stepSize;

                // ����Ǹ�������������ȡ��
                if (stepsNumber - int(stepsNumber) != 0)
                {
                    stepsNumber = int(stepsNumber+1);
                }

                // ���������ֵ�
                for (int i = 0; i < stepsNumber; i++)
                {
                    long int stepStart = i * stepSize + 1;
                    long int stepEnd = stepStart + windowSize - 1;
                    winStartMap[chrName].push_back(stepStart);
                    winEndMap[chrName].push_back(stepEnd);
                }
                // ���Ⱦɫ��Ľ�β
                if (winEndMap[chrName][winEndMap[chrName].size()-1]<chrLen)
                {
                    winStartMap[chrName].push_back(winEndMap[chrName][winEndMap[chrName].size()-1]+1);
                    winEndMap[chrName].push_back(chrLen);
                }
            } 
        }
    }
    // �ͷ��ڴ�
    chrLenFile.close();
    return make_pair(winStartMap, winEndMap);
}

// ways�ķ���
map<string,string> VCFWIN::ways_group(string waysFilename)
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Build group index.\n";
    // ways�����ֵ�
    map<string,string> waysGroupMap;

    if (waysFilename.empty()) // ���group�ļ�Ϊ�գ��򹹽��ٵ�group��ϣ��
    {
        waysGroupMap["0"] = "0";
        return waysGroupMap;
    }

    ifstream waysFile;
    waysFile.open(waysFilename, ios::in);

    string line;
    if (!waysFile.is_open())
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
            << "'" << waysFilename << "': No such file or directory." << endl;
        exit(1);
    }
    else
    {
        while (getline(waysFile, line))
        {
            if (line.empty())
            {
                continue;
            }
            else
            {
                line = strip(line, '\n');
                vector<string> lineList = split(line, "\t");

                if (lineList.size() < 2)
                {
                    lineList = split(line, " ");
                }

                if (lineList.size() < 2)
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] " << "Format error \'" << waysFilename << "\': way\tgroup\n";
                    exit(1);
                }

                string wayName = lineList[0];
                string wayGroup = lineList[1];
                waysGroupMap[wayName] = wayGroup;
            }
        }
    }

    // �ͷ��ڴ�
    waysFile.close();

    return waysGroupMap;
}


// ���㴰���ڱ���ĳ���
int VCFWIN::window_len_count(string vcfFilename, map<string,vector<long int>> & winStartMap, map<string,vector<long int>> & winEndMap, int windowSize, map<string,string> & waysGroupMap, string mode)
{
    // �����ļ���
    gzFile gzfp = gzopen(vcfFilename.c_str(), "rb");

    string::size_type idx;
    string::size_type idx1;
    int waysNumber = waysGroupMap.size();

    vector<string> waysVector;
    map<string,map<long int, vector<long int>>> chrWinStartLenMap;  // <chromosome, winStart, vector<length>>
    map<string,map<long int, vector<string>>> chrWinStartGroupMap;  // <chromosome, winStart, vector<group>>

    if(!gzfp)
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "'" << vcfFilename << "': No such file or directory." << endl;
        exit(1);
        
    }
    else
    {
        string information;
        char line[1024]; // һ��ֻ��1024�ֽڵ�����

        while(gzgets(gzfp, line, 1024))
        {
            information += line;

            if (information.find("\n") != string::npos) // һ�н���
            {
                idx = information.find("#");
                idx1 = information.find("#CHROM");
                if (idx1 != string::npos)
                {
                    vector<string> lineList = split(information, "\t");
                    for (int i = 0; i < lineList.size(); i++)
                    {
                        waysVector.push_back(lineList[i]);
                    }
                }
                
                if (idx == string::npos)
                {
                    // vcf�ļ�û��ע���У���ֱ���˳����롣
                    if (waysVector.size() < 1)
                    {
                        cerr << "[" << __func__ << "::" << getTime() << "] " 
                            << "Error: missing the comment line: #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  CW01    CW02......\n";
                        exit(1);
                    }

                    string chromosome{};
                    long int refStart{};
                    long int refLen{};
                    vector<int> qryLenList;

                    vector<string> lineList = split(information, "\t");

                    chromosome = lineList[0];                  
                    refStart = stoi(lineList[1]);
                    refLen = lineList[3].size();
                    long int refEnd = refStart + refLen - 1;

                    // ��qry���а�','��ֲ�ѭ�����뵽qryLenList��
                    vector<string> qry_seq_list = split(lineList[4], ",");
                    for (int i = 0; i < qry_seq_list.size(); i++)
                    {
                        qryLenList.push_back(qry_seq_list[i].size());
                    }

                    // vcf�ӵ�10�к�ʼ�ǻ�������Ϣ�����Դ�i=9��Ҳ���ǵ�10�п�ʼѭ��
                    // ��ͬ��vcf�ļ���ʽ���ܲ�ͬ����Ҫ���
                    vector<string> gt_list;
                    string gt;
                    for (int i = 9; i < lineList.size(); i++)
                    {
                        gt = lineList[i];

                        // ���gt�ǿյ�������
                        if (gt == "." || gt == "0/0" || gt =="0|0")
                        {
                            continue;
                        }
                        int gtIndex = stoi(split(lineList[i], "/")[0]) - 1;
                        int gtLen = qryLenList[gtIndex];
                        vector<string> stepsVectorKey;

                        // �ڲ����ֵ��е������� �ȼ��Ⱦɫ���ڲ����ֵ��У����ڵĻ�����
                        for (auto rit = winStartMap.begin(); rit != winStartMap.end(); ++rit)
                        {
                            stepsVectorKey.push_back(rit->first);                       
                        }

                        if (stepsVectorKey.end() == find(stepsVectorKey.begin(), stepsVectorKey.end(), chromosome))
                        {
                            cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: indexError: " << chromosome << "\t" << stepsVectorKey[0] << endl;
                            exit(1);
                        }
                        else
                        {
                            int indexRight = search_Binary_right(winEndMap[chromosome], refStart);  // ��һ����������
                            int indexLeft = search_Binary_left(winStartMap[chromosome], refStart);  // ��һ����������
                            
                            string way = waysVector[i];  // ����ϵ��way
                            string group;  // ����ϵ��group
                            if (waysGroupMap.end() == waysGroupMap.find(way))  // group
                            {
                                group = "0";
                            }
                            else
                            {
                                group = waysGroupMap[way];
                            }

                            long int winStart1 = winStartMap[chromosome][indexRight];  // ��һ������
                            long int winStart2 = winStartMap[chromosome][indexLeft];  // ������

                            // ��ref���Ⱥ�group��ӵ���ϣ����
                            if (indexLeft - indexRight == 1)  // ��һ��������Ϣ���
                            {
                                chrWinStartLenMap[chromosome][winStart1].push_back(min(winEndMap[chromosome][indexRight], refEnd)-refStart+1);  // ����
                                chrWinStartGroupMap[chromosome][winStart1].push_back(group);  // group
                            }

                            // ������
                            chrWinStartLenMap[chromosome][winStart2].push_back(min(winStart2+windowSize-1, refEnd)-refStart+1);  // ����
                            chrWinStartGroupMap[chromosome][winStart2].push_back(group);  // group
                        }
                    }
                }
                // ����ַ���
                information.clear();
                string().swap(information);
            }
        }
    }

    // �ͷ��ڴ�
    gzclose(gzfp);

    // ��¼ÿ��group��Ӧ������ <group, number>
    map<string, int> groupNumMap;
    for (auto rit = waysGroupMap.begin(); rit != waysGroupMap.end(); ++rit)
    {
        string group = rit->second;
        groupNumMap[group]++;
    }

    // save result
    ofstream outFile;
    string outputFilename = split(vcfFilename, ".")[0] + ".out";
    outFile.open(outputFilename, ios::out);

    // ��ӡ��ͷ
    outFile << "#waysNumber: " << waysNumber << endl;
    outFile << "#chr\tLocation\t"<< mode <<"\taverage";
    for (auto iter1 : groupNumMap)
    {
        outFile << "\t" << iter1.first;
    }
    outFile << endl;

    // ������
    long double sumValue;
    long double averageValue;

    for (auto iter1 : winStartMap)
    {
        int idxTmp = 0;
        string chromosome = iter1.first;

        // �̵߳�Ⱦɫ�岻��Ӧ�Ļ�����ӡ
        if (chrWinStartLenMap.find(chromosome) == chrWinStartLenMap.end())
        {
            continue;
        }

        for (auto iter2 : iter1.second)
        {
            long int winStart = iter2;
            long int winEnd = winEndMap[chromosome][idxTmp];
            map<long int, vector<long int>>::iterator findIter1 = chrWinStartLenMap[chromosome].find(winStart);
            if (findIter1 != chrWinStartLenMap[chromosome].end())
            {
                // ͳ�ƴ����ڱ���ĳ���
                if (mode == "length")
                {
                    // �ô��ڱ��쳤�ȵ��ܺ�
                    sumValue = accumulate(findIter1->second.begin(), findIter1->second.end(), 0);
                    // ���쳤��ռ���ڵı���
                    averageValue = sumValue/(winEnd-winStart+1);
                }
                // ͳ�ƴ����ڱ��������
                else
                {
                    // �ô��ڱ����������ܺ�
                    sumValue = findIter1->second.size();
                    // ƽ���ı�������
                    averageValue = sumValue/waysNumber;
                }
                
                outFile << chromosome << "\t" << winStart << "-" << winEnd << "\t" << sumValue << "\t" << averageValue;
                    
                // group������ͳ��
                for (auto iter3 : groupNumMap)
                {
                    string group = iter3.first;
                    int groupNum = iter3.second;
                    int number = count(chrWinStartGroupMap[chromosome][winStart].begin(), chrWinStartGroupMap[chromosome][winStart].end(), group);
                    outFile << "\t" << number/(float)groupNum;
                }
            }
            else
            {
                outFile << chromosome << "\t" << winStart << "-" << winEnd << "\t" << 0 << "\t" << 0;

                // group������ͳ��
                for (auto iter3 : groupNumMap)
                {
                    outFile << "\t" << 0;
                }
            }
            outFile << endl;
            idxTmp++;
        }
    }

    // �ͷ��ڴ�
    outFile.close();

    return 0;
}