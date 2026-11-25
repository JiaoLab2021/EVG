#pragma once

#include <cstdint>
#include <string>

namespace graphvcf {

enum class Command { Convert, Merge, Recall, Split, Breakpoint, Count, Filter, Delly2vcf, Win };

// ------ Options (defaults are the single source for parsing and help) ------
struct ConvertOptions {
    std::string reference_file;
    std::string vcf_file;
    int read_length = 350;
    double maf = 0;
    double missing_rate = 1.0;
    std::string output_file = "vcfConvert.out.vcf.gz";
};

struct MergeOptions {
    std::string vcf_file;
    std::string paragraph_file;
    std::string graphtyper2_file;
    std::string bayestyper_file;
    std::string vg_map_file;
    std::string vg_giraffe_file;
    std::string graphaligner_file;
    std::string pangenie_file;
    std::string sample_name = "";
    std::string mode = "all";
    std::string output_file = "";
};

struct RecallOptions {
    std::string truth_file;
    std::string evaluate_file;
    std::string model = "recall";
    std::string score_field = "";
    std::string roc_type = "genotype";
};

struct SplitOptions {
    std::string vcf_file;
    std::string base_vcf_file;
    std::string mode = "type";
    int length = 100;
    std::string prefix = "split";
    bool filter = false;
};

struct BreakpointOptions {
    std::string vcf_file;
    std::string prefix = "breakpoint";
    int32_t error_size = 1;
};

struct CountOptions {
    std::string vcf_file;
};

struct FilterOptions {
    std::string vcf_file;
    double maf = 0.01;
    double missing_rate = 0.1;
    std::string output_file = "";
};

struct Delly2vcfOptions {
    std::string vcf_file;
    std::string reference_file;
    std::string output_file;
};

struct WinOptions {
    std::string vcf_file;
    std::string length_file;
    std::string group_file;
    std::string output_file;
    uint32_t window_size = 200000;
    uint32_t step_size = 100000;
    std::string mode = "number";
};

// The subcommand to run, the thread count and the options of every subcommand.
struct GraphvcfConfig {
    Command command = Command::Convert;
    uint32_t threads = 4;
    ConvertOptions convert;
    MergeOptions merge;
    RecallOptions recall;
    SplitOptions split;
    BreakpointOptions breakpoint;
    CountOptions count;
    FilterOptions filter;
    Delly2vcfOptions delly2vcf;
    WinOptions win;
};

// Parsed command line; run is false when only help was printed or the options were rejected.
struct ParseResult {
    GraphvcfConfig config;
    bool run = false;
    int status = 0;
};

// Parse the subcommand and its options; print help or the error on -h, bad options and parameter errors.
ParseResult parse_command(int argc, char** argv);

// Print the commands and program information.
void print_help(const char* executable);

}
