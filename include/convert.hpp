#ifndef CONVERT_HPP
#define CONVERT_HPP

#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <getopt.h>
#include "get_time.hpp"
#include "strip_split_join.hpp"
#include "zlib.h"

#include "save.hpp"
#include "kseq.h"

using namespace std;

void help_convert(char** argv);
int main_convert(int argc, char** argv);

class Convert
{
private:
    string fastaFilename_;
    string outputFilename_;
public:
    Convert(
        const string& fastaFilename, 
        const string& outputFilename
    );

    void convert();

};

#endif