#ifndef CHECK_VCF_SORT_HPP
#define CHECK_VCF_SORT_HPP
#include <string>
#include <iostream>
#include "zlib.h"

#include "vcf_open.hpp"
#include "get_time.hpp"
#include "strip_split_join.hpp"

int check_vcf_sort(string inputVcf);

#endif