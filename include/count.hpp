#ifndef COUNT_HPP
#define COUNT_HPP
#include <fstream>
#include <string>
#include <iostream>
#include "zlib.h"

#include "save.hpp"
#include "kseq.h"
#include "get_time.hpp"
#include <getopt.h>
#include <vector>
#include <map>
#include <regex>

using namespace std;

void help_count(char** argv);
int main_count(int argc, char** argv);

struct countStruct
{
    uint64_t readNum;
    uint64_t readBase;

    countStruct() : readNum(0), readBase(0) {}
};

class Count
{
private:
    // kseq.h 打开文件
    KSEQ_INIT(gzFile, gzread)

    string inputFileName1_;
    string inputFileName2_;
    string outputFileName_;

    countStruct countOut_;
public:
    Count(
        const string& inputFileName1, 
        const string& inputFileName2, 
        const string & outputFileName
    );

    // 打开fastq/a.gz文件
    void fastq_a_count();

    // 打开fastq/a.gz文件
    void save_result();
};

#endif