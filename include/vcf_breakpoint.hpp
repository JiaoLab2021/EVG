#ifndef VCF_BREAKPOINT_HPP
#define VCF_BREAKPOINT_HPP

#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include "zlib.h"
#include <getopt.h>
#include <cstdlib>
#include <regex>

#include "strip_split_join.hpp"
#include "get_time.hpp"

using namespace std;

void help_breakpoint(char** argv);

int main_breakpoint(int argc, char** argv);

namespace BREAKPOINT
{
    int vcf_breakpoint(const string & vcfFileName, const string & prefix, const int & breakpointErrorSize);
}

#endif  // namespace VCF_BREAKPOINT_HPP