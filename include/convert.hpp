#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <getopt.h>
#include "get_time.hpp"
#include "strip_split_join.hpp"
#include "zlib.h"

#include "kseq.h"

using namespace std;

void help_convert(char** argv);
int main_convert(int argc, char** argv);

namespace convert
{
    int convert(string * fastaFilename, string * outputFilename);
}