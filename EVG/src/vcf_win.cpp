// g++ vcf_win_count_thread.cpp -std=c++11 -o vcf_win_count_thread -lpthread
#include <fstream>
#include <string>
#include <iostream>
#include <numeric>
#include <getopt.h>
#include <thread>
#include "../include/get_time.hpp"
#include "../include/vcf_win.hpp"

using namespace std;

void help_win(char** argv);

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