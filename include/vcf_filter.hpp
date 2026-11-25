#pragma once

#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <unordered_map>
#include <cmath>
#include <cstdlib>

#include "strip_split_join.hpp"
#include "vcf_open.hpp"
#include "save.hpp"

using namespace std;

// Filter the SNPs of a VCF file by minor allele frequency and missing rate.
class VCFFilter
{
private:
    string vcfFileName_;
    string outputFileName_;
    double MAF_;
    double MISSRATE_;
    uint32_t threads_;
public:
    VCFFilter(
        const string & vcfFileName,
        const string & outputFileName,
        const double & MAF,
        const double & MISSRATE,
        uint32_t threads
    );

    void filter_vcf();

    bool filter_record(VCFINFOSTRUCT & INFOSTRUCTTMP) const;
};
