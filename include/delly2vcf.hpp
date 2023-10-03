#ifndef delly2vcf_hpp
#define delly2vcf_hpp

#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <regex>
#include <unordered_map>
#include "zlib.h"
#include <getopt.h>
#include <sstream>

#include "kseq.h"
#include "get_time.hpp"
#include "strip_split_join.hpp"
#include "sort_find.hpp"
#include "vcf_open.hpp"
#include "save.hpp"

using namespace std;

void help_delly2vcf(char** argv);
int main_delly2vcf(int argc, char** argv);


class Delly2VCF
{
private:
    string refFileName_;
    string dellyFilename_;
    string outputFilename_;

    KSEQ_INIT(gzFile, gzread)

    unordered_map<string , string> FastaMap_;  // map<chr, sequence>
public:
    /**
	 * init
	 *
	 * @param refFileName_
     * @param dellyFilename
     * @param outputFilename
     * 
	**/
    Delly2VCF(
        const string& refFileName,
        const string& dellyFilename, 
        const string& outputFilename
    );

    void build_fasta_index();

    void vcf_convert();
};

#endif