#pragma once

#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <map>
#include <unordered_map>
#include <sstream>

#include "strip_split_join.hpp"
#include "vcf_open.hpp"
#include "save.hpp"
#include "sort_find.hpp"

using namespace std;

// REF and ALT lengths of a data record found by the per-record worker.
struct SplitResult {
    bool keep = false;  // get_len accepted the record
    uint32_t refLen = 0;
    uint32_t qryLen = 0;
};

// Split a VCF file into several gzip VCF files.
class VCFSplit
{
private:
    struct svInfo
    {
        string information;
        uint32_t refEnd;

        svInfo() : information(""), refEnd(0) {}
    };

    struct snpIndelInfo
    {
        vector<uint32_t> refStartVec;
        vector<uint32_t> refEndVec;
    };

    string vcfFileName_;
    string baseVcfFileName_;
    int length_;
    bool filterBool_;
    string prefix_;
    uint32_t threads_;

public:
    VCFSplit(
        string vcfFileName,
        string baseVcfFileName,
        int length,
        bool filterBool,
        string prefix,
        uint32_t threads
    );

    bool get_len(
        const VCFINFOSTRUCT& INFOSTRUCTTMP,
        uint32_t& refLen,
        uint32_t& qryLen
    ) const;

    void measure_record(VCFINFOSTRUCT& INFOSTRUCTTMP, SplitResult& result) const;

    void split_by_type();

    void split_by_number();
};
