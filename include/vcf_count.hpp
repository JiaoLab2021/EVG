#ifndef VCF_COUNT_HPP
#define VCF_COUNT_HPP

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
#include <getopt.h>
#include <cstdlib>


using namespace std;


void help_count(char** argv);
int main_count(int argc, char** argv);


class VCFCount
{
private:
    string vcfFileName_;
public:
    /**
	 * init
	 *
	 * @param vcfFileName       input VCF  file name
     * 
	**/
    VCFCount(const string & vcfFileName);

    void count();
};

#endif