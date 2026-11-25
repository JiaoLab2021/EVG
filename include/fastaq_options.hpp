#pragma once

#include <cstdint>
#include <string>

namespace fastaq {

// fastAQ subcommands, in help order.
enum class Command { Count, Sample, Convert };

// ------ Options ------
// Options of fastAQ count; an empty output file means stdout.
struct CountOptions {
    std::string input_file1;
    std::string input_file2;
    std::string output_file;
};

// Options of fastAQ sample with their defaults; the help prints these defaults. An empty output file means stdout.
struct SampleOptions {
    std::string input_file;
    std::string output_file;
    double frac = 0;
    long int seed = 11;
    std::string mode = "2";
};

// Options of fastAQ convert; an empty output file means stdout.
struct ConvertOptions {
    std::string input_file;
    std::string output_file;
};

// The selected subcommand, the thread count and the options of every subcommand.
struct FastaqConfig {
    Command command = Command::Count;
    uint32_t threads = 4;
    CountOptions count;
    SampleOptions sample;
    ConvertOptions convert;
};

// Parsed configuration; run is false after help or a parameter error, status is then the exit code.
struct ParseResult {
    FastaqConfig config;
    bool run = false;
    int status = 0;
};

// Parse the subcommand in argv[1] and its options.
ParseResult parse_command(int argc, char** argv);

// Print the commands and program information to stderr.
void print_help(const char* executable);

}
