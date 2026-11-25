#pragma once

#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <cstdlib>
#include <regex>

#include "vcf_open.hpp"
#include "save.hpp"
#include "strip_split_join.hpp"

using namespace std;

// Shift the start position of every variant by a fixed size.
class VCFBreakpoint
{
private:
    string vcfFileName_;
    string prefix_;
    int32_t breakpointErrorSize_;
    uint32_t threads_;

public:
    VCFBreakpoint(
        const string & vcfFileName,
        const string & prefix,
        const int32_t & breakpointErrorSize,
        uint32_t threads
    );

    void shift_breakpoints();

    void shift_record(VCFINFOSTRUCT & INFOSTRUCTTMP, string & outLine) const;
};
