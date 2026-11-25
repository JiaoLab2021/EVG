#include "../include/graphvcf_options.hpp"

#include <cctype>
#include <getopt.h>
#include <sstream>
#include <stdexcept>
#include <string>

#include "../include/HelpPrinter.hpp"
#include "../include/ProgramMetadata.hpp"
#include "../include/logger.hpp"

using namespace std;

namespace graphvcf {

// Print the commands and program information.
void print_help(const char* executable) {
    diagnostic_stream() << "Usage: " << executable << " <command> [options]\n\n"
                        << "Convert, merge, evaluate and summarise VCF files for genome graph genotyping.\n\n"
                        << "Version: " << program::version << '\n'
                        << "Date: " << program::build_date << '\n'
                        << "Author: " << program::author << '\n'
                        << "Contact: " << program::email << "\n\n";

    HelpPrinter help;
    help.section("Commands");
    help.line("convert", "", "convert a VCF file to the input format of the genome graph tools");
    help.line("merge", "", "merge the VCF files of one sample genotyped by several genome graph tools");
    help.line("recall", "", "evaluate genotyping results against a truth VCF file");
    help.line("split", "", "split a VCF file by variant type or by the number of nearby SNPs and indels");
    help.line("breakpoint", "", "shift variant start positions by a fixed number of bases");
    help.line("count", "", "count alleles by variant type and length");
    help.line("filter", "", "filter SNPs by minor allele frequency and missing rate");
    help.line("delly2vcf", "", "convert a DELLY VCF file to a VCF file with REF and ALT sequences");
    help.line("win", "", "count variants in sliding windows along each chromosome");
    help.blank();
    help.section("General options");
    help.line("-h, --help", "", "show this help");
    help.line("--version", "", "show the program version");
    help.print(diagnostic_stream());
}

namespace {

// Return an option as shown in the help: "-x" for a letter code, "--name" for a long-only (digit) code.
string option_name(const struct option* long_options, int code) {
    if (isalpha(code)) return string("-") + static_cast<char>(code);
    for (const struct option* entry = long_options; entry->name != nullptr; ++entry) {
        if (entry->val == code) return string("--") + entry->name;
    }
    return string(1, static_cast<char>(code));
}

// Parse a thread count: digits only, 1 to 256; throw std::logic_error otherwise.
uint32_t parse_threads(const char* text) {
    const string value(text);
    if (value.empty() || value.find_first_not_of("0123456789") != string::npos) throw invalid_argument("threads");
    if (value.size() > 9) throw out_of_range("threads");
    const unsigned long threads = stoul(value);
    if (threads < 1 || threads > 256) throw out_of_range("threads");
    return static_cast<uint32_t>(threads);
}

// Add the -t line (or --threads when -t is taken) with the default thread count.
void add_threads_help(HelpPrinter& help, const string& option = "-t, --threads") {
    help.line(option, "INT", "number of threads [" + to_string(GraphvcfConfig().threads) + "]");
}

// Print the convert help.
void help_convert(const char* executable) {
    const ConvertOptions defaults;
    ostringstream maf, geno;
    maf << defaults.maf;
    geno << defaults.missing_rate;

    diagnostic_stream() << "Usage: " << executable << " convert -r FILE -v FILE [options]\n\n"
                        << "convert a VCF file to the input format of the genome graph tools\n\n";

    HelpPrinter help;
    help.section("Input/output");
    help.line("-r, --reference", "FILE", "reference FASTA, plain or gzip");
    help.line("-v, --vcf", "FILE", "input VCF");
    help.line("-o, --out", "FILE", "output VCF, gzip if the name contains .gz; CHROMOSOME.NAME is also written to the current directory [" + defaults.output_file + "]");
    help.blank();
    help.section("Filtering");
    help.line("-l, --length", "INT", "read length; skip variants starting before this position, values above 1000 turn this off [" + to_string(defaults.read_length) + "]");
    help.line("--maf", "FLOAT", "drop SNPs whose minor allele frequency among called alleles is below FLOAT; the SNP filter runs only when --maf > 0 and --geno < 1 [" + maf.str() + "]");
    help.line("--geno", "FLOAT", "drop SNPs with missing rate above FLOAT; when the SNP filter runs, SNPs with fewer than 2 called samples are dropped too [" + geno.str() + "]");
    help.blank();
    help.section("General options");
    add_threads_help(help);
    help.line("-h, --help", "", "show this help");
    help.print(diagnostic_stream());
}

// Print the merge help.
void help_merge(const char* executable) {
    const MergeOptions defaults;

    diagnostic_stream() << "Usage: " << executable << " merge -v FILE -n STR [options]\n\n"
                        << "merge the VCF files of one sample genotyped by several genome graph tools\n\n";

    HelpPrinter help;
    help.section("Input/output");
    help.line("-v, --vcf", "FILE", "population VCF used for genotyping; every input VCF must be sorted");
    help.line("--VG-MAP", "FILE", "VG-MAP output VCF");
    help.line("--VG-Giraffe", "FILE", "VG-Giraffe output VCF");
    help.line("--GraphAligner", "FILE", "GraphAligner output VCF");
    help.line("--Paragraph", "FILE", "Paragraph output VCF");
    help.line("--BayesTyper", "FILE", "BayesTyper output VCF");
    help.line("--GraphTyper2", "FILE", "GraphTyper2 output VCF");
    help.line("--PanGenie", "FILE", "PanGenie output VCF");
    help.line("-n, --name", "STR", "sample name, a sample column of the tool VCFs");
    help.line("-o, --out", "FILE", "output VCF, gzip if the name contains .gz [" + (defaults.output_file.empty() ? string("stdout") : defaults.output_file) + "]");
    help.blank();
    help.section("Merging");
    help.line("-m, --mode", "STR", "all (population VCF; also writes the high-confidence h.<out>) or specific (sample-specific VCF) [" + defaults.mode + "]");
    help.blank();
    help.section("General options");
    add_threads_help(help);
    help.line("-h, --help", "", "show this help");
    help.print(diagnostic_stream());
}

// Print the recall help.
void help_recall(const char* executable) {
    const RecallOptions defaults;

    diagnostic_stream() << "Usage: " << executable << " recall -t FILE -e FILE [options]\n\n"
                        << "evaluate genotyping results against a truth VCF file\n\n";

    HelpPrinter help;
    help.section("Input/output");
    help.line("-t, --vcf1", "FILE", "truth VCF, sorted");
    help.line("-e, --vcf2", "FILE", "VCF to evaluate, sorted; results are written to the current directory");
    help.blank();
    help.section("Evaluation");
    help.line("-m, --model", "STR", "evaluation mode: recall, genotype or gramtools (convert gramtools output, then genotype) [" + defaults.model + "]");
    help.line("-s, --score-field", "STR", "VCF field used as the ROC score: FORMAT.<name> or INFO.<name> [" + (defaults.score_field.empty() ? string("none") : defaults.score_field) + "]");
    help.line("-r, --roc", "STR", "calls counted as true in the ROC: recall (matched site) or genotype (matched genotype) [" + defaults.roc_type + "]");
    help.blank();
    help.section("General options");
    add_threads_help(help, "--threads");
    help.line("-h, --help", "", "show this help");
    help.print(diagnostic_stream());
}

// Print the split help.
void help_split(const char* executable) {
    const SplitOptions defaults;

    diagnostic_stream() << "Usage: " << executable << " split -v FILE [options]\n\n"
                        << "split a VCF file by variant type or by the number of nearby SNPs and indels\n\n";

    HelpPrinter help;
    help.section("Input/output");
    help.line("-v, --vcf", "FILE", "input VCF, sorted");
    help.line("-V, --basevcf", "FILE", "VCF giving the SNP and indel positions for -m number [same as -v]");
    help.line("-p, --prefix", "STR", "output prefix [" + defaults.prefix + "]");
    help.blank();
    help.section("Splitting");
    help.line("-m, --mode", "STR", "type (SNP+indel+ins+del, inv, dup, other) or number (count of SNPs and indels near each SV) [" + defaults.mode + "]");
    help.line("-l, --length", "INT", "distance from the SV breakpoints searched for SNPs and indels in -m number [" + to_string(defaults.length) + "]");
    help.line("-f, --filter", "", "skip records whose FILTER is not PASS; take the ALT allele from the genotype instead of the first ALT");
    help.blank();
    help.section("General options");
    add_threads_help(help);
    help.line("-h, --help", "", "show this help");
    help.print(diagnostic_stream());
}

// Print the breakpoint help.
void help_breakpoint(const char* executable) {
    const BreakpointOptions defaults;

    diagnostic_stream() << "Usage: " << executable << " breakpoint -v FILE [options]\n\n"
                        << "shift variant start positions by a fixed number of bases\n\n";

    HelpPrinter help;
    help.section("Input/output");
    help.line("-v, --vcf", "FILE", "input VCF");
    help.line("-p, --prefix", "STR", "output prefix; writes <prefix>.vcf.gz [" + defaults.prefix + "]");
    help.blank();
    help.section("Breakpoints");
    help.line("-s, --size", "INT", "bases subtracted from POS (not below 1); END is recomputed [" + to_string(defaults.error_size) + "]");
    help.blank();
    help.section("General options");
    add_threads_help(help);
    help.line("-h, --help", "", "show this help");
    help.print(diagnostic_stream());
}

// Print the count help.
void help_count(const char* executable) {
    diagnostic_stream() << "Usage: " << executable << " count -v FILE [options]\n\n"
                        << "count alleles by variant type and length\n\n";

    HelpPrinter help;
    help.section("Input/output");
    help.line("-v, --vcf", "FILE", "input VCF; the tables are written to stdout");
    help.blank();
    help.section("General options");
    add_threads_help(help);
    help.line("-h, --help", "", "show this help");
    help.print(diagnostic_stream());
}

// Print the filter help.
void help_filter(const char* executable) {
    const FilterOptions defaults;
    ostringstream maf, geno;
    maf << defaults.maf;
    geno << defaults.missing_rate;

    diagnostic_stream() << "Usage: " << executable << " filter -v FILE [options]\n\n"
                        << "filter SNPs by minor allele frequency and missing rate\n\n";

    HelpPrinter help;
    help.section("Input/output");
    help.line("-v, --vcf", "FILE", "input VCF");
    help.line("-o, --out", "FILE", "output VCF, gzip if the name contains .gz [" + (defaults.output_file.empty() ? string("stdout") : defaults.output_file) + "]");
    help.blank();
    help.section("Filtering");
    help.line("--maf", "FLOAT", "drop SNPs whose minor allele frequency among called alleles is below FLOAT [" + maf.str() + "]");
    help.line("--geno", "FLOAT", "drop SNPs with missing rate above FLOAT; SNPs with fewer than 2 called samples are always dropped [" + geno.str() + "]");
    help.blank();
    help.section("General options");
    add_threads_help(help);
    help.line("-h, --help", "", "show this help");
    help.print(diagnostic_stream());
}

// Print the delly2vcf help.
void help_delly2vcf(const char* executable) {
    diagnostic_stream() << "Usage: " << executable << " delly2vcf -r FILE -v FILE [options]\n\n"
                        << "convert a DELLY VCF file to a VCF file with REF and ALT sequences\n\n";

    HelpPrinter help;
    help.section("Input/output");
    help.line("-r, --reference", "FILE", "reference FASTA, plain or gzip");
    help.line("-v, --vcf", "FILE", "VCF from DELLY");
    help.line("-o, --out", "FILE", "output VCF, gzip if the name contains .gz [-v up to its first '.', plus .convert.vcf]");
    help.blank();
    help.section("General options");
    add_threads_help(help);
    help.line("-h, --help", "", "show this help");
    help.print(diagnostic_stream());
}

// Print the win help.
void help_win(const char* executable) {
    const WinOptions defaults;

    diagnostic_stream() << "Usage: " << executable << " win -v FILE -l FILE [options]\n\n"
                        << "count variants in sliding windows along each chromosome\n\n";

    HelpPrinter help;
    help.section("Input/output");
    help.line("-v, --vcf", "FILE", "multi-sample VCF");
    help.line("-l, --length", "FILE", "chromosome lengths, one 'chromosome length' pair per line");
    help.line("-g, --group", "FILE", "sample groups, one 'sample group' pair per line [all samples in group 0]");
    help.line("-o, --out", "FILE", "output table, gzip if the name contains .gz [" + (defaults.output_file.empty() ? string("stdout") : defaults.output_file) + "]");
    help.blank();
    help.section("Windows");
    help.line("-w, --windowSize", "INT", "window size in bp [" + to_string(defaults.window_size) + "]");
    help.line("-s, --stepSize", "INT", "step between window starts in bp [" + to_string(defaults.step_size) + "]");
    help.line("-m, --mode", "STR", "value per window: number (non-reference genotypes of all samples) or length (their bases) [" + defaults.mode + "]");
    help.blank();
    help.section("General options");
    add_threads_help(help);
    help.line("-h, --help", "", "show this help");
    help.print(diagnostic_stream());
}

// Parse the convert options; -r and -v are required.
ParseResult parse_convert(int argc, char** argv) {
    ParseResult result;
    result.config.command = Command::Convert;
    ConvertOptions& options = result.config.convert;

    static const struct option long_options[] = {
        {"reference", required_argument, 0, 'r'},
        {"vcf", required_argument, 0, 'v'},

        {"length", required_argument, 0, 'l'},
        {"maf", required_argument, 0, '1'},
        {"geno", required_argument, 0, '2'},
        {"out", required_argument, 0, 'o'},
        {"threads", required_argument, 0, 't'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "r:v:l:1:2:o:t:h", long_options, nullptr)) != -1) {
        try {
            switch (c) {
            case 'r':
                options.reference_file = optarg;
                break;
            case 'v':
                options.vcf_file = optarg;
                break;
            case 'l':
                options.read_length = stoi(optarg);
                break;
            case '1':
                options.maf = stod(optarg);
                if (!(options.maf >= 0 && options.maf <= 1)) throw std::out_of_range("maf");
                break;
            case '2':
                options.missing_rate = stod(optarg);
                if (!(options.missing_rate >= 0 && options.missing_rate <= 1)) throw std::out_of_range("geno");
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
        } catch (const std::logic_error&) {
            error_stream() << "Parameter error: " << option_name(long_options, c) << ". Invalid value '" << optarg << "'.\n\n";
            help_convert(argv[0]);
            result.status = 1;
            return result;
        }
    }

    if (argc <= 2) {
        help_convert(argv[0]);
        result.status = 1;
        return result;
    }

    if (options.reference_file.empty() || options.vcf_file.empty()) {
        error_stream() << "Parameter error: -r and -v are required.\n\n";
        help_convert(argv[0]);
        result.status = 1;
        return result;
    }

    result.run = true;
    return result;
}

// Parse the merge options; -v, -n and a valid -m are required.
ParseResult parse_merge(int argc, char** argv) {
    ParseResult result;
    result.config.command = Command::Merge;
    MergeOptions& options = result.config.merge;

    static const struct option long_options[] = {
        {"vcf", required_argument, 0, 'v'},
        {"Paragraph", required_argument, 0, '1'},
        {"GraphTyper2", required_argument, 0, '2'},
        {"BayesTyper", required_argument, 0, '3'},
        {"VG-MAP", required_argument, 0, '4'},
        {"VG-Giraffe", required_argument, 0, '5'},
        {"GraphAligner", required_argument, 0, '6'},
        {"PanGenie", required_argument, 0, '7'},
        {"name", required_argument, 0, 'n'},

        {"mode", required_argument, 0, 'm'},
        {"out", required_argument, 0, 'o'},
        {"threads", required_argument, 0, 't'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "v:1:2:3:4:5:6:7:n:m:o:t:h", long_options, nullptr)) != -1) {
        try {
            switch (c) {
            case 'v':
                options.vcf_file = optarg;
                break;
            case '1':
                options.paragraph_file = optarg;
                break;
            case '2':
                options.graphtyper2_file = optarg;
                break;
            case '3':
                options.bayestyper_file = optarg;
                break;
            case '4':
                options.vg_map_file = optarg;
                break;
            case '5':
                options.vg_giraffe_file = optarg;
                break;
            case '6':
                options.graphaligner_file = optarg;
                break;
            case '7':
                options.pangenie_file = optarg;
                break;
            case 'n':
                options.sample_name = optarg;
                break;
            case 'm':
                options.mode = optarg;
                break;
            case 'o':
                options.output_file = optarg;
                break;
            case 't':
                result.config.threads = parse_threads(optarg);
                break;
            case 'h':
                help_merge(argv[0]);
                return result;
            case '?':
                help_merge(argv[0]);
                result.status = 1;
                return result;
            }
        } catch (const std::logic_error&) {
            error_stream() << "Parameter error: " << option_name(long_options, c) << ". Invalid value '" << optarg << "'.\n\n";
            help_merge(argv[0]);
            result.status = 1;
            return result;
        }
    }

    if (argc <= 2) {
        help_merge(argv[0]);
        result.status = 1;
        return result;
    }

    if (options.vcf_file.empty()) {
        error_stream() << "Parameter error: -v is required.\n\n";
        help_merge(argv[0]);
        result.status = 1;
        return result;
    }

    if (options.mode != "specific" && options.mode != "all") {
        error_stream() << "Parameter error: -m. Use all or specific.\n\n";
        help_merge(argv[0]);
        result.status = 1;
        return result;
    }

    if (options.sample_name.size() == 0) {
        error_stream() << "Parameter error: -n is required.\n\n";
        help_merge(argv[0]);
        result.status = 1;
        return result;
    }

    result.run = true;
    return result;
}

// Parse the recall options; -t, -e, a valid -m and a valid -r are required.
ParseResult parse_recall(int argc, char** argv) {
    ParseResult result;
    result.config.command = Command::Recall;
    RecallOptions& options = result.config.recall;

    static const struct option long_options[] = {
        {"vcf1", required_argument, 0, 't'},
        {"vcf2", required_argument, 0, 'e'},
        {"model", required_argument, 0, 'm'},

        {"score-field", required_argument, 0, 's'},
        {"roc", required_argument, 0, 'r'},
        {"threads", required_argument, 0, '1'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "t:e:m:s:r:h", long_options, nullptr)) != -1) {
        try {
            switch (c) {
            case 't':
                options.truth_file = optarg;
                break;
            case 'e':
                options.evaluate_file = optarg;
                break;
            case 'm':
                options.model = optarg;
                break;
            case 's':
                options.score_field = optarg;
                break;
            case 'r':
                options.roc_type = optarg;
                break;
            case '1':
                result.config.threads = parse_threads(optarg);
                break;
            case 'h':
                help_recall(argv[0]);
                return result;
            case '?':
                help_recall(argv[0]);
                result.status = 1;
                return result;
            }
        } catch (const std::logic_error&) {
            error_stream() << "Parameter error: " << option_name(long_options, c) << ". Invalid value '" << optarg << "'.\n\n";
            help_recall(argv[0]);
            result.status = 1;
            return result;
        }
    }

    if (options.truth_file.empty() || options.evaluate_file.empty()) {
        error_stream() << "Parameter error: -t and -e are required.\n\n";
        help_recall(argv[0]);
        result.status = 1;
        return result;
    }

    if (options.roc_type != "recall" && options.roc_type != "genotype") {
        error_stream() << "Parameter error: -r. Use recall or genotype.\n\n";
        help_recall(argv[0]);
        result.status = 1;
        return result;
    }

    if (options.model != "recall" && options.model != "genotype" && options.model != "gramtools") {
        error_stream() << "Parameter error: -m. Use recall, genotype or gramtools.\n\n";
        help_recall(argv[0]);
        result.status = 1;
        return result;
    }

    result.run = true;
    return result;
}

// Parse the split options; -v and a valid -m are required.
ParseResult parse_split(int argc, char** argv) {
    ParseResult result;
    result.config.command = Command::Split;
    SplitOptions& options = result.config.split;

    static const struct option long_options[] = {
        {"vcf", required_argument, 0, 'v'},

        {"mode", required_argument, 0, 'm'},
        {"length", required_argument, 0, 'l'},
        {"basevcf", required_argument, 0, 'V'},
        {"prefix", required_argument, 0, 'p'},
        {"filter", no_argument, 0, 'f'},

        {"threads", required_argument, 0, 't'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "v:m:l:V:p:ft:h", long_options, nullptr)) != -1) {
        try {
            switch (c) {
            case 'v':
                options.vcf_file = optarg;
                break;
            case 'm':
                options.mode = optarg;
                break;
            case 'l':
                options.length = stoi(optarg);
                break;
            case 'V':
                options.base_vcf_file = optarg;
                break;
            case 'p':
                options.prefix = optarg;
                break;
            case 'f':
                options.filter = true;
                break;
            case 't':
                result.config.threads = parse_threads(optarg);
                break;
            case 'h':
                help_split(argv[0]);
                return result;
            case '?':
                help_split(argv[0]);
                result.status = 1;
                return result;
            }
        } catch (const std::logic_error&) {
            error_stream() << "Parameter error: " << option_name(long_options, c) << ". Invalid value '" << optarg << "'.\n\n";
            help_split(argv[0]);
            result.status = 1;
            return result;
        }
    }

    if (options.vcf_file.empty()) {
        error_stream() << "Parameter error: -v is required.\n\n";
        help_split(argv[0]);
        result.status = 1;
        return result;
    }

    if (options.mode != "type" && options.mode != "number") {
        error_stream() << "Parameter error: -m. Use type or number.\n\n";
        help_split(argv[0]);
        result.status = 1;
        return result;
    }

    result.run = true;
    return result;
}

// Parse the breakpoint options; -v is required.
ParseResult parse_breakpoint(int argc, char** argv) {
    ParseResult result;
    result.config.command = Command::Breakpoint;
    BreakpointOptions& options = result.config.breakpoint;

    static const struct option long_options[] = {
        {"vcf", required_argument, 0, 'v'},

        {"prefix", required_argument, 0, 'p'},
        {"size", required_argument, 0, 's'},

        {"threads", required_argument, 0, 't'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "v:p:s:t:h", long_options, nullptr)) != -1) {
        try {
            switch (c) {
            case 'v':
                options.vcf_file = optarg;
                break;
            case 'p':
                options.prefix = optarg;
                break;
            case 's':
                options.error_size = stol(optarg);
                break;
            case 't':
                result.config.threads = parse_threads(optarg);
                break;
            case 'h':
                help_breakpoint(argv[0]);
                return result;
            case '?':
                help_breakpoint(argv[0]);
                result.status = 1;
                return result;
            }
        } catch (const std::logic_error&) {
            error_stream() << "Parameter error: " << option_name(long_options, c) << ". Invalid value '" << optarg << "'.\n\n";
            help_breakpoint(argv[0]);
            result.status = 1;
            return result;
        }
    }

    if (options.vcf_file.empty()) {
        error_stream() << "Parameter error: -v is required.\n\n";
        help_breakpoint(argv[0]);
        result.status = 1;
        return result;
    }

    result.run = true;
    return result;
}

// Parse the count options; -v is required.
ParseResult parse_count(int argc, char** argv) {
    ParseResult result;
    result.config.command = Command::Count;
    CountOptions& options = result.config.count;

    static const struct option long_options[] = {
        {"vcf", required_argument, 0, 'v'},

        {"threads", required_argument, 0, 't'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "v:t:h", long_options, nullptr)) != -1) {
        try {
            switch (c) {
            case 'v':
                options.vcf_file = optarg;
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
        } catch (const std::logic_error&) {
            error_stream() << "Parameter error: " << option_name(long_options, c) << ". Invalid value '" << optarg << "'.\n\n";
            help_count(argv[0]);
            result.status = 1;
            return result;
        }
    }

    if (options.vcf_file.empty()) {
        error_stream() << "Parameter error: -v is required.\n\n";
        help_count(argv[0]);
        result.status = 1;
        return result;
    }

    result.run = true;
    return result;
}

// Parse the filter options; -v is required.
ParseResult parse_filter(int argc, char** argv) {
    ParseResult result;
    result.config.command = Command::Filter;
    FilterOptions& options = result.config.filter;

    static const struct option long_options[] = {
        {"vcf", required_argument, 0, 'v'},

        {"maf", required_argument, 0, '1'},
        {"geno", required_argument, 0, '2'},
        {"out", required_argument, 0, 'o'},
        {"threads", required_argument, 0, 't'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "v:1:2:o:t:h", long_options, nullptr)) != -1) {
        try {
            switch (c) {
            case 'v':
                options.vcf_file = optarg;
                break;
            case '1':
                options.maf = stod(optarg);
                if (!(options.maf >= 0 && options.maf <= 1)) throw std::out_of_range("maf");
                break;
            case '2':
                options.missing_rate = stod(optarg);
                if (!(options.missing_rate >= 0 && options.missing_rate <= 1)) throw std::out_of_range("geno");
                break;
            case 'o':
                options.output_file = optarg;
                break;
            case 't':
                result.config.threads = parse_threads(optarg);
                break;
            case 'h':
                help_filter(argv[0]);
                return result;
            case '?':
                help_filter(argv[0]);
                result.status = 1;
                return result;
            }
        } catch (const std::logic_error&) {
            error_stream() << "Parameter error: " << option_name(long_options, c) << ". Invalid value '" << optarg << "'.\n\n";
            help_filter(argv[0]);
            result.status = 1;
            return result;
        }
    }

    if (options.vcf_file.empty()) {
        error_stream() << "Parameter error: -v is required.\n\n";
        help_filter(argv[0]);
        result.status = 1;
        return result;
    }

    result.run = true;
    return result;
}

// Parse the delly2vcf options; -r and -v are required.
ParseResult parse_delly2vcf(int argc, char** argv) {
    ParseResult result;
    result.config.command = Command::Delly2vcf;
    Delly2vcfOptions& options = result.config.delly2vcf;

    static const struct option long_options[] =
    {
        {"vcf", required_argument, 0, 'v'},
        {"reference", required_argument, 0, 'r'},
        {"out", required_argument, 0, 'o'},
        {"threads", required_argument, 0, 't'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "v:r:o:t:h", long_options, nullptr)) != -1) {
        try {
            switch (c) {
            case 'v':
                options.vcf_file = optarg;
                break;
            case 'r':
                options.reference_file = optarg;
                break;
            case 'o':
                options.output_file = optarg;
                break;
            case 't':
                result.config.threads = parse_threads(optarg);
                break;
            case 'h':
                help_delly2vcf(argv[0]);
                return result;
            case '?':
                help_delly2vcf(argv[0]);
                result.status = 1;
                return result;
            }
        } catch (const std::logic_error&) {
            error_stream() << "Parameter error: " << option_name(long_options, c) << ". Invalid value '" << optarg << "'.\n\n";
            help_delly2vcf(argv[0]);
            result.status = 1;
            return result;
        }
    }

    if (argc <= 2)
    {
        help_delly2vcf(argv[0]);
        result.status = 1;
        return result;
    }

    if (options.reference_file.empty() || options.vcf_file.empty()) {
        error_stream() << "Parameter error: -r and -v are required.\n\n";
        help_delly2vcf(argv[0]);
        result.status = 1;
        return result;
    }

    result.run = true;
    return result;
}

// Parse the win options; -v and -l are required.
ParseResult parse_win(int argc, char** argv) {
    ParseResult result;
    result.config.command = Command::Win;
    WinOptions& options = result.config.win;

    static const struct option long_options[] = {
        {"vcf", required_argument, 0, 'v'},
        {"length", required_argument, 0, 'l'},
        {"group", required_argument, 0, 'g'},
        {"out", required_argument, 0, 'o'},

        {"windowSize", required_argument, 0, 'w'},
        {"stepSize", required_argument, 0, 's'},
        {"mode", required_argument, 0, 'm'},
        {"threads", required_argument, 0, 't'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "v:l:g:o:w:s:m:t:h", long_options, nullptr)) != -1) {
        try {
            switch (c) {
            case 'v':
                options.vcf_file = optarg;
                break;
            case 'l':
                options.length_file = optarg;
                break;
            case 'g':
                options.group_file = optarg;
                break;
            case 'o':
                options.output_file = optarg;
                break;
            case 'w':
                options.window_size = stoi(optarg);
                break;
            case 's':
                options.step_size = stoi(optarg);
                break;
            case 'm':
                options.mode = optarg;
                break;
            case 't':
                result.config.threads = parse_threads(optarg);
                break;
            case 'h':
                help_win(argv[0]);
                return result;
            case '?':
                help_win(argv[0]);
                result.status = 1;
                return result;
            }
        } catch (const std::logic_error&) {
            error_stream() << "Parameter error: " << option_name(long_options, c) << ". Invalid value '" << optarg << "'.\n\n";
            help_win(argv[0]);
            result.status = 1;
            return result;
        }
    }

    if (argc <= 2) {
        help_win(argv[0]);
        result.status = 1;
        return result;
    }

    if (options.vcf_file.empty() || options.length_file.empty()) {
        error_stream() << "Parameter error: -v and -l are required.\n\n";
        help_win(argv[0]);
        result.status = 1;
        return result;
    }

    result.run = true;
    return result;
}

}

// Parse the subcommand and its options; print help or the error on -h, bad options and parameter errors.
ParseResult parse_command(int argc, char** argv) {
    const string command = argv[1];
    if (command == "convert") return parse_convert(argc, argv);
    if (command == "merge") return parse_merge(argc, argv);
    if (command == "recall") return parse_recall(argc, argv);
    if (command == "split") return parse_split(argc, argv);
    if (command == "breakpoint") return parse_breakpoint(argc, argv);
    if (command == "count") return parse_count(argc, argv);
    if (command == "filter") return parse_filter(argc, argv);
    if (command == "delly2vcf") return parse_delly2vcf(argc, argv);
    if (command == "win") return parse_win(argc, argv);

    error_stream() << command << ": command not found\n\n";
    print_help(argv[0]);
    ParseResult result;
    result.status = 1;
    return result;
}

}
