// g++ graphvcf.cpp -c src/*.cpp -o graphvcf -lpthread -lz -O3 -march=native
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <mutex>
#include <iomanip>

#include "include/vcf_convert.hpp"
#include "include/vcf_merge.hpp"
#include "include/vcf_recall.hpp"
#include "include/vcf_split.hpp"
#include "include/vcf_breakpoint.hpp"
#include "include/vcf_count.hpp"
#include "include/vcf_filter.hpp"
#include "include/delly2vcf.hpp"
#include "include/vcf_win.hpp"

using namespace std;

// define data
#define PROGRAM_DATA "2024/06/12"
// define version
#define PROGRAM_VERSION "1.1.9"
// define author
#define PROGRAM_AUTHOR "Zezhen Du"
// define E-mail
#define PROGRAM_E_MAIL "dzz0539@gmail.com or dzz0539@163.com"

void help(char** argv);

int main(int argc, char** argv)
{
	// Print help document
	if (argc == 1) {
        help(argv);
        return 1;
    }

	cerr << "[" << __func__ << "::" << getTime() << "] " << "You are using graphvcf (v" << PROGRAM_VERSION << ")\n\n";

	// Select the subfunction
	string subcommand = argv[1];

	if (subcommand == "-h" || subcommand == "--help")
	{
		help(argv);
        return 1;
	}
	else if (subcommand == "convert")
	{
		main_convert(argc, argv);	
	}
	else if (subcommand == "merge")
	{
		main_merge(argc, argv);
	}
	else if (subcommand == "recall")
	{
		main_recall(argc, argv);
	}
	else if (subcommand == "split")
	{
		main_split(argc, argv);	
	}
	else if (subcommand == "breakpoint")
	{
		main_breakpoint(argc, argv);
	}
	else if (subcommand == "count")
	{
		main_count(argc, argv);	
	}
	else if (subcommand == "filter")
	{
		main_filter(argc, argv);	
	}
	else if (subcommand == "delly2vcf")
	{
		main_delly2vcf(argc, argv);	
	}
	else if (subcommand == "win")
	{
		main_win(argc, argv);	
	}
	else
	{
		cerr << "Error: ["<< argv[0] << "] command " << subcommand << " not found" << endl;
		help(argv);
        return 1;
	}

    return 0;
}

// Help document
void help(char** argv)
{
  cerr << "usage: " << argv[0] << " <command> [options]" << endl
       << endl
       << "data: " << PROGRAM_DATA << endl
       << "version: " << PROGRAM_VERSION << endl
       << "author: " << PROGRAM_AUTHOR << endl
	   << endl
       << "subcommands:" << endl
       << "  convert        convert VCF files to the required format for genome graph software" << endl
	   << "  merge          merge VCF files generated by genome graph software" << endl
	   << "  recall         evaluate the results of genome graph software" << endl
	   << "  split          split VCF files by type or number" << endl
	   << "  breakpoint     set error for variants breakpoints" << endl
	   << "  count          calculate the number of alleles" << endl
	   << "  filter         filter SNPs by maf and missing rate" << endl
	   << "  delly2vcf      convert delly-generated VCF file to standard format" << endl
	   << "  win            calculate the variants density in VCF file using window-based statistics" << endl
       << endl
       << "  -h, --help     print this help document" << endl
       << endl
       << "If you encounter any issues related to the code, please don't hesitate to contact us via email at " << PROGRAM_E_MAIL << "." << endl;
}