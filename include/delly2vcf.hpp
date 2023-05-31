#ifndef delly2vcf_hpp
#define delly2vcf_hpp

#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <regex>

#include "get_time.hpp"
#include "strip_split_join.hpp"
#include "sort_find.hpp"
#include <getopt.h>

using namespace std;

void help_delly2vcf(char** argv);
int main_delly2vcf(int argc, char** argv);

namespace DELLY2VCF {
    map<string,string> build_dict(string * ref_file);

    int vcf_convert(string * dellyFilename, map<string,string> seq_dict, string * outputFilename);

}  // namespace DELLY2VCF

#endif