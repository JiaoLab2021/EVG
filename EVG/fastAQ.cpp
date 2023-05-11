// g++ fastAQ.cpp -o fastAQ -lpthread -lz
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <mutex>
#include <iomanip>
#include "src/count.cpp"
#include "src/sample_fast.cpp"
#include "src/convert.cpp"

using namespace std;
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
	else if (subcommand == "count")
	{
		main_count(argc, argv);	
	}
	else if (subcommand == "sample")
	{
		main_sample(argc, argv);	
	}
	else if (subcommand == "convert")
	{
		main_convert(argc, argv);	
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
       << "subcommands:" << endl
       << "  -- count       calculate the number of bases and reads of fasta/q" << endl
	   << "  -- sample      sample sequences by frac" << endl
	   << "  -- convert     arrange fasta sequence into one line" << endl
       << endl
       << "  -h, --help     print this help document" << endl;
}