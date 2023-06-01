#ifndef COUNT_HPP
#define COUNT_HPP
#include <fstream>
#include <string>
#include <iostream>
#include "zlib.h"

#include "kseq.h"
#include "get_time.hpp"
#include <getopt.h>
#include <vector>
#include <map>
#include <regex>

using namespace std;

void help_count(char** argv);
int main_count(int argc, char** argv);

// kseq.h 打开文件
KSEQ_INIT(gzFile, gzread)

namespace count
{
    struct countStruct
    {
        long long int readNum = 0;
        long long int readBase = 0;
    };


    // 打开fastq/a.gz文件
    void fastq_a_count(
        string inputFileName1, 
        string inputFileName2, 
        countStruct & countOut, 
        const string & outputFileName
    );
}

#endif