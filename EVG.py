#!/usr/bin/env python3

import argparse
import errno
import fcntl
import logging
import os
import re
import resource
import shutil
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "src"))
from evg_pipeline import bayestyper, bwa, graphaligner, graphtyper2, merge, pangenie, paragraph, prepare, vg
from evg_pipeline.common import CODE_DIR, ENV, get_version, parse_version, run_jobs

SOFTWARE = ("VG-MAP", "VG-Giraffe", "GraphAligner", "Paragraph", "BayesTyper", "GraphTyper2", "PanGenie")

logger = logging.getLogger("EVG")


def get_cpu_time():
    """Return the CPU seconds used by EVG and its finished child processes."""
    self_usage = resource.getrusage(resource.RUSAGE_SELF)
    child_usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    return self_usage.ru_utime + self_usage.ru_stime + child_usage.ru_utime + child_usage.ru_stime


START_TIME = time.perf_counter()
START_CPU = get_cpu_time()


def format_number(number, width):
    """Format a number as the C++ logger does: six decimals cut to `width` characters, keeping the integer part."""
    text = f"{number:.6f}"
    return text[:max(width, text.find("."))]


class LogFormatter(logging.Formatter):
    """Format records as "[L::function::wall*cpu] message" like graphvcf and fastAQ; stage titles follow a blank line."""

    def format(self, record):
        wall = time.perf_counter() - START_TIME
        level = {logging.WARNING: "W", logging.ERROR: "E"}.get(record.levelno, "M")
        prefix = f"[{level}::{record.funcName[:10].ljust(10, '.')}::{format_number(wall, 6)}*" \
                 f"{format_number((get_cpu_time() - START_CPU) / wall if wall > 0 else 0.0, 4)}] "
        return ("\n" if getattr(record, "stage", False) else "") + prefix + record.getMessage()


def log_stage(title):
    """Log a stage title after a blank line."""
    logger.info(title, extra={"stage": True}, stacklevel=2)


def read_metadata():
    """Read the program name, version, date, author and contact from include/ProgramMetadata.hpp."""
    with open(os.path.join(CODE_DIR, "include", "ProgramMetadata.hpp")) as f:
        return dict(re.findall(r'inline constexpr const char\* (\w+)\s*=\s*"([^"]*)";', f.read()))


def parse_args(metadata):
    """Parse the command line and make the input paths absolute."""
    parser = argparse.ArgumentParser(
        description=f"{metadata['description']}\n\n"
                    f"Version: {metadata['version']}\n"
                    f"Date: {metadata['build_date']}\n"
                    f"Author: {metadata['author']}\n"
                    f"Contact: {metadata['email']}",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--version", action="version", version=metadata["version"])
    input_option = parser.add_argument_group(title="required input")
    input_option.add_argument("-r", dest="reference", metavar="FILE", help="input FASTA reference", type=str, required=True)
    input_option.add_argument("-v", dest="vcf", metavar="FILE", help="population VCF (plain or .gz)", type=str, required=True)
    input_option.add_argument("-s", dest="samples", metavar="FILE", help="samples file; each line: NAME READ1 [READ2]", type=str, required=True)
    depth_option = parser.add_argument_group(title="depth")
    depth_option.add_argument("--depth", dest="depth", metavar="FLOAT", help="read depth for genotyping [%(default)s]", type=float, default=15)
    algorithm_option = parser.add_argument_group(title="algorithm")
    algorithm_option.add_argument("--software", dest="software", metavar="STR", help=f"genotyping software: {', '.join(SOFTWARE)} [auto]", type=str, nargs="+", choices=SOFTWARE, default=None)
    algorithm_option.add_argument("--mode", dest="mode", metavar="STR", help="software mode: precise, fast [%(default)s]", type=str, choices=("precise", "fast"), default="precise")
    parallel_option = parser.add_argument_group(title="parallel")
    parallel_option.add_argument("--threads", dest="threads", metavar="INT", help="threads per job [%(default)s]", type=int, default=10)
    parallel_option.add_argument("--jobs", dest="jobs", metavar="INT", help="jobs run in parallel [%(default)s]", type=int, default=3)
    parallel_option.add_argument("--number", dest="number", metavar="INT", help="samples per BayesTyper/Paragraph batch [%(default)s]", type=int, default=5)
    resume_option = parser.add_argument_group(title="resume").add_mutually_exclusive_group()
    resume_option.add_argument("--force", action="store_true", help="clear existing output directories")
    resume_option.add_argument("--resume", action="store_true", help="continue a stopped run; skip commands that finished with the same options and unchanged files")
    args = parser.parse_args()
    args.reference = os.path.abspath(args.reference)
    args.vcf = os.path.abspath(args.vcf)
    args.samples = os.path.abspath(args.samples)
    return args


def print_config(args):
    """Log the run settings."""
    log_stage("Configuration ...")
    logger.info(f"  - Reference: {args.reference}")
    logger.info(f"  - VCF: {args.vcf}")
    logger.info(f"  - Samples: {args.samples}")
    logger.info(f"  - Software: {' '.join(args.software) if args.software else 'auto'}")
    logger.info(f"  - Mode: {args.mode}")
    logger.info(f"  - Depth: {args.depth:g}x")
    logger.info(f"  - Threads: {args.threads}")
    logger.info(f"  - Jobs: {args.jobs}")
    logger.info(f"  - Batch size: {args.number}")
    if args.force:
        logger.info("  - Force: yes")
    if args.resume:
        logger.info("  - Resume: yes")


def check_inputs(args):
    """Raise if --threads is out of range or an input file or an external tool is missing."""
    if not 1 <= args.threads <= 256:
        raise ValueError(f"--threads {args.threads}: must be 1 to 256.")
    for path in (args.reference, args.vcf, args.samples):
        if not os.path.exists(path):
            raise RuntimeError(f"'{path}' No such file or directory.")
    for tool in ("graphvcf", "fastAQ"):
        if not os.access(os.path.join(CODE_DIR, tool), os.X_OK):
            raise RuntimeError(f"{os.path.join(CODE_DIR, tool)}: not found next to EVG.py.")
    for tool in ("tabix", "bwa", "samtools", "vg", "GraphAligner", "paragraph", "bayesTyper", "graphtyper", "PanGenie"):
        if shutil.which(tool, path=ENV["PATH"]) is None:
            raise RuntimeError(f"{tool}: command not found in PATH.")


def make_dir(path, args):
    """Create a step directory; clear it with --force, keep it with --resume, otherwise fail if it exists."""
    if os.path.isdir(path):
        if args.force:
            shutil.rmtree(path)
            os.makedirs(path)
            logger.info(f"  - Cleared: {path}")
        elif args.resume:
            logger.info(f"  - Reusing: {path}")
        else:
            raise RuntimeError(f"'{path}' already exists. Use '--force' to overwrite or '--resume' to continue the run.")
    else:
        os.makedirs(path)
    return path


def read_samples(samples_file):
    """Read the samples file into {name: (read1, read2)} in file order; relative read paths follow the launch directory."""
    try:
        with open(samples_file, "r") as f:
            lines = f.readlines()
    except IOError as error:
        raise ValueError(f"Cannot read the samples file: {error}")
    if not lines:
        raise ValueError(f"'{samples_file}': empty samples file.")
    samples = {}
    for line in lines:
        if not line.strip() or "#" in line or "read1_path" in line:
            continue
        fields = line.split()
        if len(fields) < 2:
            raise ValueError(f"Invalid entry in samples file. Expected 2 or 3 fields, got: {line.strip()}")
        samples[fields[0]] = (os.path.abspath(fields[1]), os.path.abspath(fields[2]) if len(fields) > 2 else "")
    for name in samples:
        if not re.fullmatch(r"[\w+][\w.+-]*", name):
            raise ValueError(f"Sample name '{name}' may only contain letters, digits, '_', '.', '+' and '-', and must not start with '.' or '-'.")
        if f"h.{name}" in samples:
            raise ValueError(f"Sample names '{name}' and 'h.{name}' would write the same file in merge/.")
    return samples


def split_batches(names, number):
    """Split sorted sample names into batches of `number` names (one batch if number <= 0)."""
    batches = []
    batch = []
    for name in sorted(names):
        batch.append(name)
        if len(batch) == number:
            batches.append(batch)
            batch = []
    if batch:
        batches.append(batch)
    return batches


def downsample_samples(args, samples, ref_base, convert_dir):
    """Count and downsample every sample (--jobs samples at once, --threads each); return {name: {read1, read2, depth, length}}."""
    tasks = [(prepare.downsample_reads, sample_name, read1, read2, ref_base, args.depth, args.threads, convert_dir)
             for sample_name, (read1, read2) in samples.items()]
    return dict(zip(samples, run_jobs(tasks, args.jobs)))


def select_software(args, reads, ref_base, convert_dir):
    """Return --software, or pick tools from read length, depth and genome size; drop PanGenie older than 3.0.0 or of unknown version."""
    if args.software:
        software = list(dict.fromkeys(args.software))
    else:
        depth_min = min([1000] + [read["depth"] for read in reads.values()])
        read_len_min = min([100000] + [read["length"] for read in reads.values()])
        if read_len_min >= 10000:
            software = ["GraphAligner"]
        else:
            software = ["Paragraph", "BayesTyper"]
            software.append("VG-Giraffe" if ref_base > 200 * 1024 * 1024 else "VG-MAP")
            software.append("GraphTyper2" if read_len_min > 130 and depth_min > 5 else "PanGenie")
    if "PanGenie" in software:
        pangenie_version = get_version("PanGenie", convert_dir)
        if not pangenie_version:
            software.remove("PanGenie")
            logger.warning("Could not read the PanGenie version. It will be skipped. Please check that PanGenie 3.0.0 or later is installed.")
        elif parse_version(pangenie_version) < (3, 0, 0):
            software.remove("PanGenie")
            logger.warning(f"PanGenie version ({pangenie_version}) is below 3.0.0. It will be skipped. Please upgrade PanGenie if you wish to use it.")
    logger.info(f"  - Software: {' '.join(software)}")
    return software


def build_indexes(args, software, reference, inputs, convert_dir):
    """Build the BWA index and the vg, GraphAligner and PanGenie indexes the software needs; return the index dirs."""
    indexes = {
        "VG-MAP": os.path.join(convert_dir, "map"),
        "VG-Giraffe": os.path.join(convert_dir, "giraffe"),
        "GraphAligner": os.path.join(convert_dir, "GraphAligner"),
        "PanGenie": os.path.join(convert_dir, "PanGenie"),
    }
    tasks = [(bwa.build_bwa_index, reference, convert_dir)]
    for tool, workflow in (("VG-MAP", "map"), ("VG-Giraffe", "giraffe")):
        if tool in software:
            make_dir(indexes[tool], args)
            tasks.append((vg.build_vg_index, reference, inputs["vcf"], workflow, args.threads, indexes[tool]))
    if "GraphAligner" in software:
        make_dir(indexes["GraphAligner"], args)
        tasks.append((graphaligner.build_graphaligner_index, reference, inputs["vcf"], args.threads, indexes["GraphAligner"]))
    if "PanGenie" in software:
        make_dir(indexes["PanGenie"], args)
        tasks.append((pangenie.build_pangenie_index, reference, inputs["vcf"], args.threads, indexes["PanGenie"]))
    run_jobs(tasks, args.jobs)
    return indexes


def map_samples(args, reads, reference, work_dir):
    """Map every sample with BWA MEM (--jobs jobs); return {name: BAM}."""
    bwa_dir = make_dir(os.path.join(work_dir, "bwa"), args)
    tasks = [(bwa.map_reads, sample_name, read["read1"], read["read2"], reference, args.threads, bwa_dir)
             for sample_name, read in reads.items()]
    return dict(zip(reads, run_jobs(tasks, args.jobs)))


def genotype_samples(args, software, reads, bams, reference, inputs, indexes, work_dir):
    """Run the selected genotyping software; return {sample: [(software, VCF), ...]} in software order."""
    genotype_dir = make_dir(os.path.join(work_dir, "genotype"), args)

    # ------ Sample lists ------
    batches = split_batches(reads, args.number)
    if "GraphTyper2" in software:
        graphtyper_list = graphtyper2.write_graphtyper_list(bams, genotype_dir)
    if "BayesTyper" in software:
        bayestyper_lists = bayestyper.write_bayestyper_lists(batches, bams, genotype_dir)
    if "Paragraph" in software:
        paragraph_lists = paragraph.write_paragraph_lists(batches, reads, bams, genotype_dir)

    # ------ Jobs ------
    tasks = []
    task_samples = []
    for tool in software:
        tool_dir = make_dir(os.path.join(genotype_dir, tool), args)
        if tool in ("VG-MAP", "VG-Giraffe"):
            for sample_name, read in reads.items():
                tasks.append((vg.genotype_vg, tool, sample_name, read["read1"], read["read2"], read["depth"],
                              args.threads, indexes[tool], tool_dir))
                task_samples.append((tool, [sample_name]))
        elif tool == "GraphAligner":
            for sample_name, read in reads.items():
                tasks.append((graphaligner.genotype_graphaligner, sample_name, read["read1"], read["depth"],
                              args.threads, indexes[tool], tool_dir))
                task_samples.append((tool, [sample_name]))
        elif tool == "PanGenie":
            for sample_name, read in reads.items():
                tasks.append((pangenie.genotype_pangenie, sample_name, read["read1"], read["read2"],
                              args.threads, indexes[tool], tool_dir))
                task_samples.append((tool, [sample_name]))
        elif tool == "GraphTyper2":
            tasks.append((graphtyper2.genotype_graphtyper2, reference, inputs["vcf"], graphtyper_list, list(bams.values()),
                          inputs["region"], args.threads, tool_dir))
            task_samples.append((tool, list(reads)))
        elif tool == "Paragraph":
            for number, list_file in enumerate(paragraph_lists):
                batch_dir = make_dir(os.path.join(tool_dir, str(number)), args)
                tasks.append((paragraph.genotype_paragraph, reference, inputs["sv"], list_file,
                              [bams[sample_name] for sample_name in batches[number]], args.threads, batch_dir))
                task_samples.append((tool, batches[number]))
        elif tool == "BayesTyper":
            for number, list_file in enumerate(bayestyper_lists):
                batch_dir = make_dir(os.path.join(tool_dir, str(number)), args)
                tasks.append((bayestyper.genotype_bayestyper, reference, inputs["vcf"], list_file, bams,
                              args.threads, batch_dir))
                task_samples.append((tool, batches[number]))

    # ------ Results ------
    calls = {sample_name: [] for sample_name in reads}
    for (tool, sample_names), vcf_file in zip(task_samples, run_jobs(tasks, args.jobs)):
        for sample_name in sample_names:
            calls[sample_name].append((tool, vcf_file))
    return calls


def merge_samples(args, calls, vcf, work_dir):
    """Merge each sample's calls with graphvcf merge (--jobs jobs); return {sample: merged VCF} in sorted order."""
    merge_dir = make_dir(os.path.join(work_dir, "merge"), args)
    sample_names = sorted(calls)
    tasks = [(merge.merge_sample, vcf, sample_name, calls[sample_name], args.threads, merge_dir)
             for sample_name in sample_names]
    merged = dict(zip(sample_names, run_jobs(tasks, args.jobs)))
    for sample_name, merge_vcf_file in merged.items():
        logger.info(f"  - {sample_name}: {merge_vcf_file}")
    return merged


def main():
    """Run the EVG workflow: prepare inputs, build indexes, map, genotype and merge."""
    handler = logging.StreamHandler()
    handler.setFormatter(LogFormatter())
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

    metadata = read_metadata()
    args = parse_args(metadata)
    logger.info(f"Version: {metadata['version']}")
    logger.info(f"Date: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Command: {' '.join(sys.argv)}")
    print_config(args)
    try:
        check_inputs(args)
        samples = read_samples(args.samples)
        work_dir = os.getcwd()
        lock_file = open(os.path.join(work_dir, ".evg.lock"), "w")
        try:
            fcntl.flock(lock_file, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as error:
            if error.errno in (errno.EAGAIN, errno.EACCES):
                raise RuntimeError(f"Another EVG run is using {work_dir}.")
            raise
        convert_dir = make_dir(os.path.join(work_dir, "convert"), args)

        log_stage("Counting reference bases ...")
        ref_base = prepare.count_fasta(args.reference, args.threads, convert_dir)

        log_stage(f"Counting and downsampling reads to {args.depth:g}x ...")
        reads = downsample_samples(args, samples, ref_base, convert_dir)

        log_stage("Converting the reference genome and VCF ...")
        reference = prepare.convert_reference(args.reference, convert_dir)
        ref_base = prepare.count_fasta(reference, args.threads, convert_dir)
        inputs = prepare.convert_vcf(reference, args.vcf, args.mode, args.threads, convert_dir)

        log_stage("Selecting genotyping software ...")
        software = select_software(args, reads, ref_base, convert_dir)

        log_stage("Building indexes ...")
        indexes = build_indexes(args, software, reference, inputs, convert_dir)

        bams = {}
        if any(tool in software for tool in ("GraphTyper2", "BayesTyper", "Paragraph")):
            log_stage("Mapping reads with BWA MEM ...")
            bams = map_samples(args, reads, reference, work_dir)

        log_stage("Genotyping ...")
        calls = genotype_samples(args, software, reads, bams, reference, inputs, indexes, work_dir)

        log_stage("Merging the genotypes of each sample ...")
        merge_samples(args, calls, inputs["vcf"], work_dir)
    except KeyboardInterrupt:
        logger.error("Interrupted; continue later with --resume.")
        sys.exit(130)
    except Exception as error:
        message = str(error) if isinstance(error, (RuntimeError, ValueError)) else f"{type(error).__name__}: {error}"
        for line in message.split("\n"):
            logger.error(line)
        sys.exit(1)

    # ------ Resource usage ------
    peak_rss = max(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    log_stage(f"Real time: {time.perf_counter() - START_TIME:.3f} sec; CPU: {get_cpu_time():.3f} sec; "
              f"Peak RSS: {peak_rss / 1024 / 1024:.3f} GB")


if __name__ == "__main__":
    main()
