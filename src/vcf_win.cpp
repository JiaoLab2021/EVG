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

    // 输入参数
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

    // 步长字典
    map<string,vector<long int>> winStartMap, winEndMap;
    tie(winStartMap, winEndMap) = VCFWIN::step_count(chrLenFilename, windowSize, stepSize);

    // ways的分组信息
    map<string,string> waysGroupMap = VCFWIN::ways_group(waysFilename);

    // 看有多少条染色体
    vector<string> chrNameVector = VCFWIN::count_chrName(chrLenFilename);
    
    // 拆分vcf
    vector<string> vcfFileVector = VCFWIN::split_file(vcfFilename, chrNameVector);

    // 开多线程
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
    // 线程阻塞
    for (auto& t: threads)
    {
        t.join();
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done.\n";
    return 0;
}

// 帮助文档
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


// 统计染色体条数和名字
vector<string> VCFWIN::count_chrName(string chrLenFilename)
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Count the number and names of chromosomes.\n";
    // 看有多少条染色体
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
    // 释放内存
    chrLenFile.close();
    return chrNameVector;
}

// 按染色体拆分vcf文件
vector<string> VCFWIN::split_file(string vcfFilename, vector<string> chrNameVector)
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Split vcf file.\n";
    vector<string> vcfFileVector;
    for (int i = 0; i < chrNameVector.size(); i++)
    {
        string vcfFile;  // 输出文件名
        string cmd_str;  // cmd

        if(vcfFilename.find(".gz") == string::npos && vcfFilename.find(".GZ") == string::npos)  // 非压缩文件
        {
            vcfFile = chrNameVector[i] + ".vcf";  // 输出为普通文件
            cmd_str = "awk '{if($1==\"" + chrNameVector[i] + "\") print$0; else if($0~/#/) print$0;}' " + vcfFilename + " > " + vcfFile;
        }
        else  // 压缩文件
        {
            vcfFile = chrNameVector[i] + ".vcf.gz";  // 输出为压缩文件
            cmd_str = "zcat " + vcfFilename + " | awk '{if($1==\"" + chrNameVector[i] + "\") print$0; else if($0~/#/) print$0;}' | gzip > " + vcfFile;
        }

        vcfFileVector.push_back(vcfFile);  // 添加输出文件到vector中
        
        int length_cmd_str = cmd_str.length();
        char cmd_char[length_cmd_str];
        strcpy(cmd_char, cmd_str.c_str());
        system(cmd_char); 
    }
    
    return vcfFileVector;
}

// 步长字典
pair<map<string,vector<long int>>, map<string,vector<long int>>> VCFWIN::step_count(string chrLenFilename, int windowSize, int stepSize)
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Build step index.\n";
    // 步长字典
    map<string,vector<long int>> winStartMap;
    map<string,vector<long int>> winEndMap;

    // 染色体长度配置文件
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

                // 步数
                double stepsNumber = double(chrLen)/stepSize;

                // 如果是浮点数，则向上取整
                if (stepsNumber - int(stepsNumber) != 0)
                {
                    stepsNumber = int(stepsNumber+1);
                }

                // 制作步长字典
                for (int i = 0; i < stepsNumber; i++)
                {
                    long int stepStart = i * stepSize + 1;
                    long int stepEnd = stepStart + windowSize - 1;
                    winStartMap[chrName].push_back(stepStart);
                    winEndMap[chrName].push_back(stepEnd);
                }
                // 最后到染色体的结尾
                if (winEndMap[chrName][winEndMap[chrName].size()-1]<chrLen)
                {
                    winStartMap[chrName].push_back(winEndMap[chrName][winEndMap[chrName].size()-1]+1);
                    winEndMap[chrName].push_back(chrLen);
                }
            } 
        }
    }
    // 释放内存
    chrLenFile.close();
    return make_pair(winStartMap, winEndMap);
}

// ways的分类
map<string,string> VCFWIN::ways_group(string waysFilename)
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Build group index.\n";
    // ways分类字典
    map<string,string> waysGroupMap;

    if (waysFilename.empty()) // 如果group文件为空，则构建假的group哈希表
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

    // 释放内存
    waysFile.close();

    return waysGroupMap;
}


// 计算窗口内变异的长度
int VCFWIN::window_len_count(string vcfFilename, map<string,vector<long int>> & winStartMap, map<string,vector<long int>> & winEndMap, int windowSize, map<string,string> & waysGroupMap, string mode)
{
    // 输入文件流
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
        char line[1024]; // 一次只读1024字节的数据

        while(gzgets(gzfp, line, 1024))
        {
            information += line;

            if (information.find("\n") != string::npos) // 一行结束
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
                    // vcf文件没有注释行，则直接退出代码。
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

                    // 把qry序列按','拆分并循环加入到qryLenList中
                    vector<string> qry_seq_list = split(lineList[4], ",");
                    for (int i = 0; i < qry_seq_list.size(); i++)
                    {
                        qryLenList.push_back(qry_seq_list[i].size());
                    }

                    // vcf从第10列后开始是基因型信息，所以从i=9，也就是第10列开始循环
                    // 不同的vcf文件格式可能不同，需要检查
                    vector<string> gt_list;
                    string gt;
                    for (int i = 9; i < lineList.size(); i++)
                    {
                        gt = lineList[i];

                        // 如果gt是空的则跳过
                        if (gt == "." || gt == "0/0" || gt =="0|0")
                        {
                            continue;
                        }
                        int gtIndex = stoi(split(lineList[i], "/")[0]) - 1;
                        int gtLen = qryLenList[gtIndex];
                        vector<string> stepsVectorKey;

                        // 在步长字典中的索引， 先检查染色体在不在字典中，不在的话报错
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
                            int indexRight = search_Binary_right(winEndMap[chromosome], refStart);  // 上一个窗口坐标
                            int indexLeft = search_Binary_left(winStartMap[chromosome], refStart);  // 下一个窗口坐标
                            
                            string way = waysVector[i];  // 该株系的way
                            string group;  // 该株系的group
                            if (waysGroupMap.end() == waysGroupMap.find(way))  // group
                            {
                                group = "0";
                            }
                            else
                            {
                                group = waysGroupMap[way];
                            }

                            long int winStart1 = winStartMap[chromosome][indexRight];  // 上一个窗口
                            long int winStart2 = winStartMap[chromosome][indexLeft];  // 本窗口

                            // 把ref长度和group添加到哈希表中
                            if (indexLeft - indexRight == 1)  // 上一个窗口信息添加
                            {
                                chrWinStartLenMap[chromosome][winStart1].push_back(min(winEndMap[chromosome][indexRight], refEnd)-refStart+1);  // 长度
                                chrWinStartGroupMap[chromosome][winStart1].push_back(group);  // group
                            }

                            // 本窗口
                            chrWinStartLenMap[chromosome][winStart2].push_back(min(winStart2+windowSize-1, refEnd)-refStart+1);  // 长度
                            chrWinStartGroupMap[chromosome][winStart2].push_back(group);  // group
                        }
                    }
                }
                // 清空字符串
                information.clear();
                string().swap(information);
            }
        }
    }

    // 释放内存
    gzclose(gzfp);

    // 记录每个group对应的数量 <group, number>
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

    // 打印表头
    outFile << "#waysNumber: " << waysNumber << endl;
    outFile << "#chr\tLocation\t"<< mode <<"\taverage";
    for (auto iter1 : groupNumMap)
    {
        outFile << "\t" << iter1.first;
    }
    outFile << endl;

    // 计算结果
    long double sumValue;
    long double averageValue;

    for (auto iter1 : winStartMap)
    {
        int idxTmp = 0;
        string chromosome = iter1.first;

        // 线程的染色体不对应的话不打印
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
                // 统计窗口内变异的长度
                if (mode == "length")
                {
                    // 该窗口变异长度的总和
                    sumValue = accumulate(findIter1->second.begin(), findIter1->second.end(), 0);
                    // 变异长度占窗口的比例
                    averageValue = sumValue/(winEnd-winStart+1);
                }
                // 统计窗口内变异的数量
                else
                {
                    // 该窗口变异数量的总和
                    sumValue = findIter1->second.size();
                    // 平均的变异数量
                    averageValue = sumValue/waysNumber;
                }
                
                outFile << chromosome << "\t" << winStart << "-" << winEnd << "\t" << sumValue << "\t" << averageValue;
                    
                // group的数量统计
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

                // group的数量统计
                for (auto iter3 : groupNumMap)
                {
                    outFile << "\t" << 0;
                }
            }
            outFile << endl;
            idxTmp++;
        }
    }

    // 释放内存
    outFile.close();

    return 0;
}