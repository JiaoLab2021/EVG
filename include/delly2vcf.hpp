#pragma once

#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <regex>
#include <unordered_map>
#include <sstream>

#include "strip_split_join.hpp"
#include "vcf_open.hpp"
#include "save.hpp"

using namespace std;

// Fill REF and ALT sequences of DELLY records from the reference.
class Delly2VCF
{
private:
    string refFileName_;
    string dellyFilename_;
    string outputFilename_;
    uint32_t threads_;

    unordered_map<string , string> FastaMap_;
public:
    Delly2VCF(
        const string& refFileName,
        const string& dellyFilename,
        const string& outputFilename,
        uint32_t threads
    );

    void build_reference_index();

    void convert_vcf();

    void convert_record(VCFINFOSTRUCT& INFOSTRUCTTMP, string& outLine) const;
};
