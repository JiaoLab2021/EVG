#pragma once

#include <fstream>
#include <string>
#include <iostream>
#include <regex>
#include <map>
#include <algorithm>
#include <sstream>
#include <unordered_set>
#include <set>

#include "strip_split_join.hpp"
#include "save.hpp"
#include "vcf_open.hpp"

using namespace std;

// Reference sequences, lengths and ##contig header lines.
struct refIndexStruct {
    map<string, string> chrSeqMap;
    map<string, uint32_t> chrLenMap;
    string chrLenTxt;

    refIndexStruct() : chrLenTxt("") {}
};

// What the per-record conversion worker found for one data record.
struct ConvertResult {
    uint32_t pos = 0;        // POS as read (the record's POS may be shifted by the conversion)
    bool keep = false;       // the record passed every per-record check
    bool lineReady = false;  // outLine is the record's output line; false after merge_loci changed it or if building it printed or threw
    string outLine;          // the record as written ("" if dropped), built on the worker
};

// Convert a VCF file to a form accepted by every genome graph tool.
class Convert
{
private:
    string refFileName_;
    string vcfFileName_;

    int readLen_;

    double MAF_;
    double MISSRATE_;

    refIndexStruct refIndexS_;

    string outFileName_;

    uint32_t threads_;

public:
    Convert(string refFileName, string vcfFileName, int readLen, double MAF, double MISSRATE, string outFileName, uint32_t threads);

    void build_reference_index();

    void convert_vcf();

private:
    void convert_record(VCFINFOSTRUCT& INFOSTRUCTTMP, ConvertResult& result) const;

    string format_record(VCFINFOSTRUCT& INFOSTRUCTTMP) const;

    void merge_loci(
        VCFINFOSTRUCT& preINFOSTRUCTTMP,
        VCFINFOSTRUCT& INFOSTRUCTTMP,
        const std::map<std::string, std::string>& chrSeqMap
    ) const;

    void process_sequences(VCFINFOSTRUCT& INFOSTRUCTTMP, const std::map<std::string, std::string>& chrSeqMap) const;

    void trim_update_alt_alleles(VCFINFOSTRUCT& INFOSTRUCTTMP) const;

    void check_repeated_alleles(VCFINFOSTRUCT& infoStruct) const;

    vector<string> sort_unique(vector<string> vec) const;
};
