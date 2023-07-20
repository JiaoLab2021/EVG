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

#include "vcf_open.hpp"
#include "save.hpp"
#include "strip_split_join.hpp"
#include "get_time.hpp"

using namespace std;

void help_breakpoint(char** argv);

int main_breakpoint(int argc, char** argv);

class VCFBreakpoint
{
private:
    string vcfFileName_;
    string prefix_;
    int breakpointErrorSize_;

public:
    /**
	 * init
	 *
	 * @param vcfFileName                     input VCF  file name
     * @param prefix                          output file prefix
     * @param breakpointErrorSize             breakpoint error
     * 
	**/
    VCFBreakpoint(
        const string & vcfFileName, 
        const string & prefix, 
        const int & breakpointErrorSize
    );

    void vcf_breakpoint();
};

#endif  // namespace VCF_BREAKPOINT_HPP