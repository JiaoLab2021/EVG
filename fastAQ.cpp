#include <exception>
#include <iomanip>
#include <string>

#include <htslib/hts_log.h>

#include "include/ProgramMetadata.hpp"
#include "include/fastaq_convert.hpp"
#include "include/fastaq_count.hpp"
#include "include/fastaq_options.hpp"
#include "include/fastaq_sample.hpp"
#include "include/get_time.hpp"
#include "include/logger.hpp"
#include "include/sys.hpp"

namespace {

// Run the parsed subcommand and return its exit status.
int run_command(const fastaq::FastaqConfig& config) {
    switch (config.command) {
    case fastaq::Command::Count: {
        const fastaq::CountOptions& options = config.count;
        Count counter(options.input_file1, options.input_file2, options.output_file, config.threads);
        counter.count_reads();
        counter.save_result();
        return 0;
    }
    case fastaq::Command::Sample: {
        const fastaq::SampleOptions& options = config.sample;
        stk_sample(options.input_file, options.output_file, options.frac, options.seed, options.mode, config.threads);
        return 0;
    }
    case fastaq::Command::Convert: {
        const fastaq::ConvertOptions& options = config.convert;
        FastaConvert converter(options.input_file, options.output_file, config.threads);
        converter.convert();
        return 0;
    }
    }
    return 1;
}

}

// Run one fastAQ subcommand and return its exit status.
int main(int argc, char** argv) {
    hts_set_log_level(HTS_LOG_OFF);  // htslib errors reach the user as EVG errors, in order
    if (argc < 2) {
        fastaq::print_help(argv[0]);
        return 1;
    }
    const std::string command = argv[1];
    if (command == "-h" || command == "--help") {
        fastaq::print_help(argv[0]);
        return 0;
    }
    if (command == "--version") {
        output_stream() << program::version << '\n';
        return 0;
    }
    try {
        const fastaq::ParseResult parsed = fastaq::parse_command(argc, argv);
        if (!parsed.run) return parsed.status;

        const double start_time = realtime();
        log_stream() << "Version: " << program::version << '\n';
        log_stream() << "Date: " << getTime() << '\n';
        const int status = run_command(parsed.config);
        if (status != 0) return status;
        stage_stream() << "Real time: " << std::fixed << std::setprecision(3) << (realtime() - start_time)
                       << " sec; CPU: " << cputime() << " sec; Peak RSS: "
                       << peakrss() / 1024.0 / 1024.0 / 1024.0 << " GB\n";
        return 0;
    } catch (const std::exception& error) {
        error_stream() << error.what() << '\n';
        return 1;
    }
}
