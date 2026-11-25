#include "../include/fastaq_options.hpp"

#include <cmath>
#include <getopt.h>
#include <stdexcept>
#include <string>
#include <sys/stat.h>

#include "../include/HelpPrinter.hpp"
#include "../include/ProgramMetadata.hpp"
#include "../include/logger.hpp"

namespace fastaq {

namespace {

// Parse a thread count: digits only, 1 to 256; throw std::logic_error otherwise.
uint32_t parse_threads(const char* text) {
    const std::string value(text);
    if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos) throw std::invalid_argument("threads");
    if (value.size() > 9) throw std::out_of_range("threads");
    const unsigned long threads = std::stoul(value);
    if (threads < 1 || threads > 256) throw std::out_of_range("threads");
    return static_cast<uint32_t>(threads);
}

// True when both names exist and are the same file (same device and inode: same name, symlink or hard link).
bool same_file(const std::string& name1, const std::string& name2) {
    struct stat stat1;
    struct stat stat2;
    return stat(name1.c_str(), &stat1) == 0 && stat(name2.c_str(), &stat2) == 0 && stat1.st_dev == stat2.st_dev
           && stat1.st_ino == stat2.st_ino;
}

// Add the -t line with the default thread count.
void add_threads_help(HelpPrinter& help) {
    help.line("-t, --threads", "INT", "number of threads [" + std::to_string(FastaqConfig().threads) + "]");
}

// Print the count help.
void help_count(const char* executable) {
    diagnostic_stream() << "Usage: " << executable << " count -i FILE [options]\n\n"
                        << "count bases and reads in FASTA/FASTQ files\n\n";

    HelpPrinter help;
    help.section("Input/output");
    help.line("-i, --input", "FILE", "input FASTA/FASTQ, plain or gzip; give -i twice for paired reads, one per mate");
    help.line("-o, --out", "FILE", "output file, gzip if the name contains .gz [stdout]");
    help.blank();
    help.section("General options");
    add_threads_help(help);
    help.line("-h, --help", "", "show this help");
    help.print(diagnostic_stream());
}

// Print the sample help using the option defaults.
void help_sample(const char* executable) {
    const SampleOptions defaults;
    diagnostic_stream() << "Usage: " << executable << " sample -i FILE -f FLOAT [options]\n\n"
                        << "randomly sample sequences by fraction or number\n\n";

    HelpPrinter help;
    help.section("Input/output");
    help.line("-i, --input", "FILE", "input FASTA/FASTQ, plain or gzip");
    help.line("-o, --out", "FILE", "output file, gzip if the name contains .gz [stdout]");
    help.blank();
    help.section("Sampling");
    help.line("-f, --frac", "FLOAT", "fraction (0-1) or number (>1) of sequences to sample");
    help.line("-s, --seed", "INT", "random seed [" + std::to_string(defaults.seed) + "]");
    help.line("-m, --mode", "INT", "for a number: 1 = one pass, samples kept in memory; 2 = two passes, less memory [" + defaults.mode + "]");
    help.blank();
    help.section("General options");
    add_threads_help(help);
    help.line("-h, --help", "", "show this help");
    help.print(diagnostic_stream());
}

// Print the convert help.
void help_convert(const char* executable) {
    diagnostic_stream() << "Usage: " << executable << " convert -i FILE [options]\n\n"
                        << "write each FASTA sequence on one line\n\n";

    HelpPrinter help;
    help.section("Input/output");
    help.line("-i, --input", "FILE", "input FASTA, plain text");
    help.line("-o, --out", "FILE", "output file, gzip if the name contains .gz [stdout]");
    help.blank();
    help.section("General options");
    add_threads_help(help);
    help.line("-h, --help", "", "show this help");
    help.print(diagnostic_stream());
}

// Parse count options; the first -i is read 1 and any later -i replaces read 2.
ParseResult parse_count(int argc, char** argv) {
    ParseResult result;
    result.config.command = Command::Count;
    CountOptions& options = result.config.count;

    static const struct option long_options[] = {
        {"input", required_argument, 0, 'i'},
        {"out", required_argument, 0, 'o'},
        {"threads", required_argument, 0, 't'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "i:o:t:h", long_options, nullptr)) != -1) {
        try {
            switch (c) {
            case 'i':
                if (options.input_file1.empty()) {
                    options.input_file1 = optarg;
                } else {
                    options.input_file2 = optarg;
                }
                break;
            case 'o':
                options.output_file = optarg;
                break;
            case 't':
                result.config.threads = parse_threads(optarg);
                break;
            case 'h':
                help_count(argv[0]);
                return result;
            case '?':
                help_count(argv[0]);
                result.status = 1;
                return result;
            }
        } catch (const std::logic_error&) {  // invalid or out-of-range number
            error_stream() << "Parameter error: -" << static_cast<char>(c)
                           << ". Invalid value '" << optarg << "'.\n\n";
            help_count(argv[0]);
            result.status = 1;
            return result;
        }
    }

    if (options.input_file1.empty()) {
        error_stream() << "Parameter error: -i.\n\n";
        help_count(argv[0]);
        result.status = 1;
        return result;
    }

    result.run = true;
    return result;
}

// Parse sample options; -f is required and must be a fraction in (0, 1) or a count in (1, 2^63).
ParseResult parse_sample(int argc, char** argv) {
    ParseResult result;
    result.config.command = Command::Sample;
    SampleOptions& options = result.config.sample;
    std::string frac_text;

    static const struct option long_options[] = {
        {"input", required_argument, 0, 'i'},
        {"frac", required_argument, 0, 'f'},
        {"seed", required_argument, 0, 's'},
        {"mode", required_argument, 0, 'm'},
        {"out", required_argument, 0, 'o'},
        {"threads", required_argument, 0, 't'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "i:f:s:m:o:t:h", long_options, nullptr)) != -1) {
        try {
            switch (c) {
            case 'i':
                options.input_file = optarg;
                break;
            case 'o':
                options.output_file = optarg;
                break;
            case 't':
                result.config.threads = parse_threads(optarg);
                break;
            case 'f':
                options.frac = std::stod(optarg);
                frac_text = optarg;
                break;
            case 's':
                options.seed = std::stol(optarg);
                break;
            case 'm':
                options.mode = optarg;
                break;
            case 'h':
                help_sample(argv[0]);
                return result;
            case '?':
                help_sample(argv[0]);
                result.status = 1;
                return result;
            }
        } catch (const std::logic_error&) {  // invalid or out-of-range number
            error_stream() << "Parameter error: -" << static_cast<char>(c)
                           << ". Invalid value '" << optarg << "'.\n\n";
            help_sample(argv[0]);
            result.status = 1;
            return result;
        }
    }

    if (options.input_file.empty()) {
        error_stream() << "Parameter error: -i.\n\n";
        help_sample(argv[0]);
        result.status = 1;
        return result;
    }

    if (frac_text.empty()) {
        error_stream() << "Parameter error: -f is required.\n\n";
        help_sample(argv[0]);
        result.status = 1;
        return result;
    }

    if (!std::isfinite(options.frac) || options.frac <= 0 || options.frac == 1
        || options.frac >= std::ldexp(1.0, 63)) {
        error_stream() << "Parameter error: -f. Invalid value '" << frac_text << "'.\n\n";
        help_sample(argv[0]);
        result.status = 1;
        return result;
    }

    if (!options.output_file.empty() && same_file(options.input_file, options.output_file)) {
        error_stream() << "Parameter error: -o. '" << options.output_file << "' is the input file.\n\n";
        help_sample(argv[0]);
        result.status = 1;
        return result;
    }

    result.run = true;
    return result;
}

// Parse convert options.
ParseResult parse_convert(int argc, char** argv) {
    ParseResult result;
    result.config.command = Command::Convert;
    ConvertOptions& options = result.config.convert;

    static const struct option long_options[] = {
        {"input", required_argument, 0, 'i'},
        {"out", required_argument, 0, 'o'},
        {"threads", required_argument, 0, 't'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "i:o:t:h", long_options, nullptr)) != -1) {
        try {
            switch (c) {
            case 'i':
                options.input_file = optarg;
                break;
            case 'o':
                options.output_file = optarg;
                break;
            case 't':
                result.config.threads = parse_threads(optarg);
                break;
            case 'h':
                help_convert(argv[0]);
                return result;
            case '?':
                help_convert(argv[0]);
                result.status = 1;
                return result;
            }
        } catch (const std::logic_error&) {  // invalid or out-of-range number
            error_stream() << "Parameter error: -" << static_cast<char>(c)
                           << ". Invalid value '" << optarg << "'.\n\n";
            help_convert(argv[0]);
            result.status = 1;
            return result;
        }
    }

    if (options.input_file.empty()) {
        error_stream() << "Parameter error: -i.\n\n";
        help_convert(argv[0]);
        result.status = 1;
        return result;
    }

    if (!options.output_file.empty() && same_file(options.input_file, options.output_file)) {
        error_stream() << "Parameter error: -o. '" << options.output_file << "' is the input file.\n\n";
        help_convert(argv[0]);
        result.status = 1;
        return result;
    }

    result.run = true;
    return result;
}

}

// Print the commands and program information.
void print_help(const char* executable) {
    diagnostic_stream() << "Usage: " << executable << " <command> [options]\n\n"
                        << "Count, sample and convert FASTA/FASTQ files.\n\n"
                        << "Version: " << program::version << '\n'
                        << "Date: " << program::build_date << '\n'
                        << "Author: " << program::author << '\n'
                        << "Contact: " << program::email << "\n\n";

    HelpPrinter help;
    help.section("Commands");
    help.line("count", "", "count bases and reads in FASTA/FASTQ files");
    help.line("sample", "", "randomly sample sequences by fraction or number");
    help.line("convert", "", "write each FASTA sequence on one line");
    help.blank();
    help.section("General options");
    help.line("-h, --help", "", "show this help");
    help.line("--version", "", "show the program version");
    help.print(diagnostic_stream());
}

// Parse the subcommand in argv[1] and its options; run is false after help or a parameter error.
ParseResult parse_command(int argc, char** argv) {
    const std::string command = argv[1];
    if (command == "count") return parse_count(argc, argv);
    if (command == "sample") return parse_sample(argc, argv);
    if (command == "convert") return parse_convert(argc, argv);

    error_stream() << command << ": command not found\n\n";
    print_help(argv[0]);
    ParseResult result;
    result.status = 1;
    return result;
}

}
