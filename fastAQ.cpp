// g++ fastAQ.cpp -o fastAQ -lpthread -lz -O3
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <mutex>
#include <iomanip>

#include "include/count.hpp"
#include "src/sample_fast.cpp"
#include "include/convert.hpp"

// define data
#define PROGRAM_DATA "2024/03/18"
// define version
#define PROGRAM_VERSION "1.1.5"
// define author
#define PROGRAM_AUTHOR "Zezhen Du"
// define E-mail
#define PROGRAM_E_MAIL "dzz0539@gmail.com or dzz0539@163.com"

using namespace std;
void help(char** argv);

int main(int argc, char** argv)
{
	// Print help document
	if (argc == 1) {
        help(argv);
        return 1;
    }

	cerr << "[" << __func__ << "::" << getTime() << "] " << "You are using fastAQ (v" << PROGRAM_VERSION << ")\n\n";

	// Select the subfunction
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
       << "  count          calculate the number of bases and reads of fasta/q" << endl
	   << "  sample         sample sequences by frac" << endl
	   << "  convert        arrange fasta sequence into one line" << endl
       << endl
       << "  -h, --help     print this help document" << endl
       << endl
       << "If you encounter any issues related to the code, please don't hesitate to contact us via email at " << PROGRAM_E_MAIL << "." << endl;
}