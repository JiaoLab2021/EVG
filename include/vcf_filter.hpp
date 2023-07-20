#ifndef vcf_filter_hpp
#define vcf_filter_hpp

#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include "zlib.h"
#include <unordered_map>
#include <cmath>
#include "strip_split_join.hpp"
#include "get_time.hpp"
#include "sort_find.hpp"
#include "vcf_open.hpp"
#include "save.hpp"
#include <getopt.h>
#include <cstdlib>


using namespace std;


void help_filter(char** argv);
int main_filter(int argc, char** argv);


class VCFFilter
{
private:
    string vcfFileName_;
    string outputFileName_;
    double MAF_;
    double MISSRATE_;
public:
    /**
	 * init
	 *
	 * @param vcfFileName         input VCF  file name
     * @param outputFileName      output file name
     * @param MAF                 MAF
     * @param MISSRATE            MISSRATE
     * 
	**/
    VCFFilter(
        const string & vcfFileName, 
        const string & outputFileName, 
        const double & MAF, 
        const double & MISSRATE
    );

    void vcf_filter();
};

#endif