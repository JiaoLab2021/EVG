// g++ graphvcf.cpp -o graphvcf -lpthread -lz
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <mutex>
#include <iomanip>
#include "src/vcf_convert.cpp"
#include "src/vcf_split.cpp"
#include "src/vcf_recall.cpp"
#include "src/vcf_breakpoint.cpp"
#include "src/vcf_merge.cpp"
#include "src/vcf_win.cpp"
#include "src/vcf_count.cpp"
#include "src/vcf_filter.cpp"
#include "src/delly2vcf.cpp"

using namespace std;

// define version
#define PROGRAM_VERSION "0.01"

void help(char** argv);

int main(int argc, char** argv)
{
	// 打印帮助文档
	if (argc == 1) {
        help(argv);
        return 1;
    }

	// 选择的子功能
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
	else if (subcommand == "split")
	{
		main_split(argc, argv);	
	}
	else if (subcommand == "recall")
	{
		main_recall(argc, argv);
	}
	else if (subcommand == "breakpoint")
	{
		main_breakpoint(argc, argv);
	}
	else if (subcommand == "merge")
	{
		main_merge(argc, argv);
	}
	else if (subcommand == "win")
	{
		main_win(argc, argv);	
	}
	else if (subcommand == "count")
	{
		main_count(argc, argv);	
	}
	else if (subcommand == "filter")
	{
		main_filter(argc, argv);	
	}
	else if (subcommand == "DELLY2VCF")
	{
		main_delly2vcf(argc, argv);	
	}
	else
	{
		cerr << "Error: ["<< argv[0] << "] command " << subcommand << " not found" << endl;
		help(argv);
        return 1;
	}

    return 0;
}

// 帮助文档
void help(char** argv)
{
  cerr << "usage: " << argv[0] << " <command> [options]" << endl
       << endl
       << "version: " << PROGRAM_VERSION << endl
	   << endl
       << "subcommands:" << endl
       << "  -- convert     convert vcf files merged by vcftools to the format required by genome graph" << endl
	   << "  -- split       split variation by type or number" << endl
	   << "  -- recall      assessing the results of genomes graph software" << endl
	   << "  -- breakpoint  set error for breakpoint of variations" << endl
	   << "  -- merge       merge the output of genomes graph software" << endl
	   << "  -- win         variation statistics within a sliding window" << endl
	   << "  -- count       count the number of alleles" << endl
	   << "  -- filter      filter SNPs by maf and missing rate" << endl
	   << "  -- DELLY2VCF   convert the format of delly to vcf" << endl
       << endl
       << "  -h, --help     print this help document" << endl;
}