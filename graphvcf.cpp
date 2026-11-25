#include <exception>
#include <iomanip>
#include <string>

#include <htslib/hts_log.h>

#include "include/ProgramMetadata.hpp"
#include "include/delly2vcf.hpp"
#include "include/get_time.hpp"
#include "include/graphvcf_options.hpp"
#include "include/logger.hpp"
#include "include/strip_split_join.hpp"
#include "include/sys.hpp"
#include "include/vcf_breakpoint.hpp"
#include "include/vcf_convert.hpp"
#include "include/vcf_count.hpp"
#include "include/vcf_filter.hpp"
#include "include/vcf_merge.hpp"
#include "include/vcf_recall.hpp"
#include "include/vcf_split.hpp"
#include "include/vcf_win.hpp"

namespace {

// Run the parsed subcommand: build its module object from the options, call the steps in order and return the exit status.
int run_command(const graphvcf::GraphvcfConfig& config) {
    switch (config.command) {
    case graphvcf::Command::Convert: {
        const graphvcf::ConvertOptions& options = config.convert;
        const int read_length = (options.read_length > 1000) ? 0 : options.read_length;
        Convert ConvertClass(options.reference_file, options.vcf_file, read_length, options.maf, options.missing_rate, options.output_file, config.threads);
        ConvertClass.build_reference_index();
        ConvertClass.convert_vcf();
        return 0;
    }
    case graphvcf::Command::Merge: {
        const graphvcf::MergeOptions& options = config.merge;
        VCFMerge VCFMergeClass(
            options.vcf_file,
            options.paragraph_file,
            options.graphtyper2_file,
            options.bayestyper_file,
            options.vg_map_file,
            options.vg_giraffe_file,
            options.graphaligner_file,
            options.pangenie_file,
            options.sample_name,
            options.mode,
            options.output_file,
            config.threads
        );
        VCFMergeClass.run_index_merge();
        const bool similar_sample_found = VCFMergeClass.vcf_merge_filter();
        VCFMergeClass.result_save();

        if (options.mode == "all" && similar_sample_found) {
            log_stream() << "Generating high-confidence VCF file in progress." << std::endl;
            VCFMergeClass.mMode = "specific";

            std::string directory;
            std::string fileName;
            VCFMergeClass.extract_directory_FileName(options.output_file, directory, fileName);
            VCFMergeClass.mOutputFileName = directory + "h." + fileName;

            VCFMergeClass.vcf_merge_filter(true);
            VCFMergeClass.result_save();
        }
        return 0;
    }
    case graphvcf::Command::Recall: {
        const graphvcf::RecallOptions& options = config.recall;
        VCFRecall VCFRecallClass(options.truth_file, options.evaluate_file, options.score_field, options.roc_type, config.threads);
        VCFRecallClass.build_true_index();
        if (options.model == "recall") {
            VCFRecallClass.evaluate_recall();
        } else if (options.model == "genotype") {
            VCFRecallClass.evaluate_genotype();
        } else {
            VCFRecallClass.convert_gramtools();
            VCFRecallClass.evaluate_genotype();
        }
        return 0;
    }
    case graphvcf::Command::Split: {
        const graphvcf::SplitOptions& options = config.split;
        VCFSplit VCFSplitClass(options.vcf_file, options.base_vcf_file, options.length, options.filter, options.prefix, config.threads);
        if (options.mode == "type") {
            VCFSplitClass.split_by_type();
        } else {
            VCFSplitClass.split_by_number();
        }
        return 0;
    }
    case graphvcf::Command::Breakpoint: {
        const graphvcf::BreakpointOptions& options = config.breakpoint;
        VCFBreakpoint VCFBreakpointClass(options.vcf_file, options.prefix, options.error_size, config.threads);
        VCFBreakpointClass.shift_breakpoints();
        return 0;
    }
    case graphvcf::Command::Count: {
        VCFCount VCFCountClass(config.count.vcf_file, config.threads);
        VCFCountClass.count_variants();
        VCFCountClass.count_length_bins();
        VCFCountClass.print_result();
        return 0;
    }
    case graphvcf::Command::Filter: {
        const graphvcf::FilterOptions& options = config.filter;
        VCFFilter VCFFilterClass(options.vcf_file, options.output_file, options.maf, options.missing_rate, config.threads);
        VCFFilterClass.filter_vcf();
        return 0;
    }
    case graphvcf::Command::Delly2vcf: {
        const graphvcf::Delly2vcfOptions& options = config.delly2vcf;
        const std::string output_file = options.output_file.empty()
            ? split(options.vcf_file, ".")[0] + ".convert.vcf"
            : options.output_file;
        Delly2VCF Delly2VCFClass(options.reference_file, options.vcf_file, output_file, config.threads);
        Delly2VCFClass.build_reference_index();
        Delly2VCFClass.convert_vcf();
        return 0;
    }
    case graphvcf::Command::Win: {
        const graphvcf::WinOptions& options = config.win;
        VCFWin VCFWinClass(
            options.window_size,
            options.step_size,
            options.length_file,
            options.vcf_file,
            options.group_file,
            options.output_file,
            options.mode,
            config.threads
        );
        VCFWinClass.load_chromosome_lengths();
        VCFWinClass.build_windows();
        VCFWinClass.load_groups();
        VCFWinClass.count_window_variants();
        VCFWinClass.save_result();
        return 0;
    }
    }
    return 1;  // not reached: every Command is handled above
}

}

// Parse the command line, run one graphvcf subcommand and return its exit status.
int main(int argc, char** argv) {
    hts_set_log_level(HTS_LOG_OFF);  // htslib errors reach the user as EVG errors, in order
    if (argc < 2) {
        graphvcf::print_help(argv[0]);
        return 1;
    }
    const std::string command = argv[1];
    if (command == "-h" || command == "--help") {
        graphvcf::print_help(argv[0]);
        return 0;
    }
    if (command == "--version") {
        output_stream() << program::version << '\n';
        return 0;
    }
    try {
        const graphvcf::ParseResult parsed = graphvcf::parse_command(argc, argv);
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
