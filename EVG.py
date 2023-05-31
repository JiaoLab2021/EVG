#!/usr/bin/python3

# -*- coding: utf-8 -*-


__data__ = "2023/05/30"
__version__ = "1.0.1"
__author__ = "Zezhen Du"
__email__ = "dzz0539@gmail.com or dzz0539@163.com"

import os
import argparse
import sys
import shutil
from multiprocessing import Pool
import concurrent.futures
import time
import logging

# code path
code_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(code_dir, "src/"))
import convert, select_software, merge, bwa, GraphTyper2, \
    GraphAligner, BayesTyper, Paragraph, VG_MAP_Giraffe, PanGenie

# environment variable
env_path = {'PATH': os.environ.get('PATH')}

# global parameters
force = False
restart = False
threads = 10
jobs_num = 3
mode = "precise"
need_depth = 15.0
merge_mode = "specific"


# log
logger = logging.getLogger('SynDiv')
formatter = logging.Formatter('[%(asctime)s] %(message)s')
handler = logging.StreamHandler()  # output to the console
handler.setFormatter(formatter)
logger.addHandler(handler)


# Print the result of merge
def print_merge(
    vcf_merge_files_map: dict
):
    """
    :param vcf_merge_files_map: Merged vcf file, map<sample_name, vcf_file>
    :return:
    """
    for key, value in vcf_merge_files_map.items():
        log = "{}: {}\n".format(key, value)
        logger.error(log)

    return 0


# Create a directory
def makedir(
    path_dir: str,
    force_tmp: bool,
    restart_tmp: bool
):
    """
    :param path_dir:     The folder path to be created
    :param force_tmp:    Whether to forcibly delete existing directories
    :param restart_tmp   Whether to continue the program based on the status of the generated files. If this parameter is specified and force is not specified, the software will skip the existing directory and no error will be reported
    :return: 0
    """
    if os.path.isdir(path_dir):
        if force_tmp:  # If forced to delete, empty the directory and create a new one
            shutil.rmtree(path_dir)
            os.makedirs(path_dir)
            log = '[EVG.makedir] \'{}\' already exists, clear and recreate.'.format(path_dir)
            logger.error(log)
        elif restart_tmp:
            log = '[EVG.makedir] \'{}\' already exists, restart is used, ' \
                  'skipping emptying folders.\n'.format(path_dir)
            logger.error(log)
        else:  # Print a warning and exit code if not forced to delete
            log = '[EVG.makedir] \'{}\' already exists. used --force to overwrite or --restart to restart workflow.\n'.format(path_dir)
            logger.error(log)
            raise SystemExit(1)
    else:
        os.makedirs(path_dir)


# Parse the samples file
def get_samples_path(
        samples_file: str
):
    """
    :param samples_file: file containing sample names for sequencing data, 'name\tread1_path\tread2_path'
    :return samples_map  map<sample_name, list<fastq_path>>
    """
    samples_map = {}
    with open(samples_file, 'r') as f:
        samples_list = f.readlines()

    # Judge the length of samples_list
    if len(samples_list) == 0:
        log = '[EVG.get_parser] Error: empty file -> {}.\n'.format(samples_file)
        logger.error(log)
        # os.killpg(os.getpgid(os.getpid()), signal.SIGKILL)
        raise SystemExit(1)

    for sample in samples_list:
        if not sample or "#" in sample or "read1_path" in sample:
            continue

        sample_list = sample.split()
        if len(sample_list) > 1:
            name = sample_list[0]
            read1_path = sample_list[1]
            if len(sample_list) > 2:  # paired-end sequencing
                read2_path = sample_list[2]
                samples_map[name] = [os.path.abspath(read1_path), os.path.abspath(read2_path)]
            else:  # single-end sequencing
                samples_map[name] = [os.path.abspath(read1_path), ""]
    return samples_map


def get_parser():
    """
    :return: parser_map = {
        "reference_file": reference_file,
        "vcf_file": vcf_file,
        "samples_file": samples_file,
        "vcf_out_name": vcf_out_name,
        "path": path,
        "software_list": software_list,
        "line_vcf_file": line_vcf_file
    }
    """
    # log
    logger = logging.getLogger('getParser')

    logger.error(f"data: {__data__}")
    logger.error(f"version: {__version__}")
    logger.error(f"author: {__author__}")
    logger.error(f"\nIf you encounter any issues related to the code, please don't hesitate to contact us via email at {__email__}.\n")

    # parameter
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    input_option = parser.add_argument_group(title="required input")
    input_option.add_argument("-r", dest="reference", help="input FASTA reference",
                              type=argparse.FileType('r'), required=True)
    input_option.add_argument("-v", dest="vcf", help="input the merged vcf file",
                              type=argparse.FileType('r'), required=True)
    # Parse the samples file
    input_option.add_argument("-s", dest="samples", help="samples file (name read1_path read2_path)",
                              type=argparse.FileType('r'), required=True)

    optional_option = parser.add_argument_group(title="optional input")
    optional_option.add_argument("-V", dest="VCF", help="input line-specific vcf file [-v]",
                                 type=argparse.FileType('r'))

    # Set how many x data to take to type
    depth_option = parser.add_argument_group(title="depth")
    depth_option.add_argument("--depth", dest="depth", help="read depth for genotyping [15]",
                              type=float, default=15)

    # user-defined software
    algorithm_option = parser.add_argument_group(
        title="algorithm"
    )
    algorithm_option.add_argument(
        "--software", dest="software", help="genotyping software [auto]", type=str, nargs='+',
        choices=['VG-MAP', 'VG-Giraffe', 'GraphAligner', 'Paragraph', 'BayesTyper', 'GraphTyper2', 'PanGenie'],
        default="auto"
    )
    # process
    algorithm_option.add_argument(
        "--mode", dest="mode", help="software mode [precise]", type=str, choices=['precise', 'fast'], default="precise"
    )

    # parallel parameter
    parallel_option = parser.add_argument_group(
        title="parallel"
    )
    # thread
    parallel_option.add_argument(
        "--threads", dest="threads", help="number of threads [10]", type=int, default=10
    )
    # number of processes
    parallel_option.add_argument(
        "--jobs", dest="jobs", help="run n jobs in parallel [3]", type=int, default=3
    )

    # How to deal with existing directories and files
    resume_option = parser.add_argument_group(
        title="resume"
    )
    resume_option.add_argument(
        "--force", action="store_true", help="erase already existing output directory [False]"
    )
    resume_option.add_argument(
        "--restart", action="store_true", help="attempt to restart existing workflow based "
                                               "on the state of output files [False]"
    )

    args = parser.parse_args()

    # parsing parameters
    reference_file = args.reference.name
    vcf_file = args.vcf.name
    try:
        line_vcf_file = args.VCF.name
    except AttributeError:
        line_vcf_file = vcf_file

    # Sequencing data files
    samples_file = args.samples.name

    vcf_out_name = "convert." + os.path.basename(vcf_file.replace(".gz", ""))
    software_list = args.software

    # Modify global variables
    global force, restart, threads, jobs_num, mode, need_depth, merge_mode
    force = args.force
    restart = args.restart
    threads = args.threads
    jobs_num = args.jobs
    mode = args.mode
    need_depth = args.depth
    if line_vcf_file == vcf_file:  # Merge process algorithm, if the same, use all algorithm
        merge_mode = "all"

    # Complete the path
    reference_file = os.path.abspath(reference_file)
    vcf_file = os.path.abspath(vcf_file)
    line_vcf_file = os.path.abspath(line_vcf_file)
    samples_file = os.path.abspath(samples_file)

    # current path
    path = os.getcwd()

    # dictionary of output parameters
    parser_map = {
        "reference_file": reference_file,
        "vcf_file": vcf_file,
        "samples_file": samples_file,
        "vcf_out_name": vcf_out_name,
        "path": path,
        "software_list": software_list,
        "line_vcf_file": line_vcf_file
    }

    return parser_map


# bwa
def run_bwa(
    reference_file: str,
    sample_name: str,
    fastq_file1: str,
    fastq_file2: str,
    real_depth: int,
    read_len: float,
    work_path: str
):
    """
    :param reference_file:   Transformed reference genome
    :param sample_name:      sample name
    :param fastq_file1:      read1
    :param fastq_file2:      read2
    :param real_depth:       sequencing dapth
    :param read_len:         read length
    :param work_path:        work path
    :return: stdout, stderr, log_out, bam_info_map = {
            "sample_name": sample_name,
            "fastq_file1": fastq_file1,
            "fastq_file2": fastq_file2,
            "real_depth": real_depth,
            "read_len": read_len,
            "bam_file": bam_file
        }
    """

    os.chdir(work_path)
    stdout, stderr, log_out, bam_file = bwa.mem(
        reference_file,
        threads,
        sample_name, 
        env_path, 
        restart,
        fastq_file1,
        fastq_file2
    )

    # Report an error if there is a problem with the exit code
    if log_out:
        return stdout, stderr, log_out, {}

    # fastq path, sequencing depth, read length
    bam_info_map = {
        "sample_name": sample_name,
        "fastq_file1": fastq_file1,
        "fastq_file2": fastq_file2,
        "real_depth": real_depth,
        "read_len": read_len,
        "bam_file": bam_file
    }

    return stdout, stderr, log_out, bam_info_map


# GraphTyper2
def run_graphtyper(
    reference_file: str,
    vcf_file: str,
    region_file: str,
    bam2graphtyper_file: str,
    work_path: str
):
    """
    :param reference_file:        reference genome
    :param vcf_file:              vcf
    :param region_file:           configuration file of region
    :param bam2graphtyper_file:   configuration file of BAM
    :param work_path:             work path
    :return: "GraphTyper2", sample_name, graphtyper_vcf_file
    """

    # Create folder and switch paths
    os.chdir(work_path)
    graphtyper_path = os.path.join(work_path, "GraphTyper2")
    makedir(
        graphtyper_path,
        force_tmp=force,
        restart_tmp=restart
    )
    os.chdir(graphtyper_path)

    # genotype
    stdout, stderr, log_out, graphtyper_vcf_file = GraphTyper2.main(
        reference_file,
        vcf_file,
        bam2graphtyper_file,
        region_file, 
        env_path, 
        threads,
        restart
    )

    return stdout, stderr, log_out, "GraphTyper2", "", graphtyper_vcf_file


# Paragraph
def run_paragraph(
        reference_file: str,
        vcf_file: str,
        bam2paragraph_file: str,
        work_path: str
):
    """
    :param reference_file:      reference genome
    :param vcf_file:            vcf file
    :param bam2paragraph_file:  configuration file of BAM
    :param work_path:           work path
    :return: "Paragraph", sample_name, paragraph_vcf_file
    """

    # Create folder and switch paths
    os.chdir(work_path)
    paragraph_path = os.path.join(work_path, "Paragraph")
    makedir(
        paragraph_path,
        force_tmp=force,
        restart_tmp=restart
    )
    os.chdir(paragraph_path)

    # genotype
    stdout, stderr, log_out, paragraph_vcf_file = Paragraph.main(
        reference_file,
        vcf_file,
        bam2paragraph_file, 
        env_path, 
        threads,
        restart
    )

    return stdout, stderr, log_out, "Paragraph", "", paragraph_vcf_file


# BayesTyper
def run_bayestyper(
        reference_file: str,
        vcf_file: str,
        bam2bayestyper_file: str,
        work_path: str,
        bam_infos_map: str
):
    """
    :param reference_file:        reference genome
    :param vcf_file:              vcf file
    :param bam2bayestyper_file:   configuration file
    :param work_path:             work path
    :param bam_infos_map:         the informations of BAM
    :return:
    """

    # Create folder and switch paths
    os.chdir(work_path)
    bayestyper_path = os.path.join(work_path, "BayesTyper")
    makedir(
        bayestyper_path,
        force_tmp=force,
        restart_tmp=restart
    )
    os.chdir(bayestyper_path)

    # genotype
    stdout, stderr, log_out, bayestyper_vcf_file = BayesTyper.main(
        reference_file,
        vcf_file,
        bam2bayestyper_file,
        bam_infos_map, 
        env_path, 
        threads,
        restart
    )

    return stdout, stderr, log_out, "BayesTyper", "", bayestyper_vcf_file


# vg_map_giraffe
def run_vg_map_giraffe(
    sample_name: str,
    fastq_file1: str,
    fastq_file2: str,
    real_depth: int,
    software: str,
    software_work_path: str,
    index_dir: str
):
    """
    :param sample_name:          sample name
    :param fastq_file1:          read1
    :param fastq_file2:          read2
    :param real_depth:           read depth
    :param software:             software
    :param software_work_path:   work path
    :param index_dir:            the path of index
    :return: software, sample_name, map_vcf_file
    """

    # Create folder and switch paths
    os.chdir(software_work_path)
    map_path = os.path.join(software_work_path, sample_name)
    makedir(
        map_path,
        force_tmp=force,
        restart_tmp=restart
    )
    os.chdir(map_path)

    # genotype
    stdout, stderr, log_out, map_vcf_file = VG_MAP_Giraffe.main(
        sample_name,
        fastq_file1,
        fastq_file2, 
        env_path, 
        threads,
        real_depth,
        software,
        index_dir,
        restart
    )

    return stdout, stderr, log_out, software, sample_name, map_vcf_file


# GraphAligner
def run_graphaligner(
    sample_name: str,
    fastq_file1: str,
    real_depth: int,
    software_work_path: str,
    index_dir: str
):
    """
    :param sample_name:          sample name
    :param fastq_file1:          read
    :param real_depth:           read depth
    :param software_work_path:   work path
    :param index_dir:            the path of index
    :return: "GraphAligner", sample_name, graphaligner_vcf_file
    """

    # genotype
    os.chdir(software_work_path)
    graphaligner_path = os.path.join(software_work_path, sample_name)
    makedir(
        graphaligner_path,
        force_tmp=force,
        restart_tmp=restart
    )
    os.chdir(graphaligner_path)

    stdout, stderr, log_out, graphaligner_vcf_file = GraphAligner.main(
        sample_name,
        fastq_file1,
        threads,
        real_depth,
        index_dir,
        env_path, 
        restart
    )

    return stdout, stderr, log_out, "GraphAligner", sample_name, graphaligner_vcf_file


# PanGenie
def run_pangenie(
    sample_name: str,
    reference_file: str,
    vcf_file: str,
    fastq_file1: str,
    fastq_file2: str,
    software_work_path: str
):
    """

    :param sample_name:        sample name
    :param reference_file:     reference genome
    :param vcf_file:           vcf file
    :param fastq_file1:        read1
    :param fastq_file2:        read2
    :param software_work_path: work path
    :return: "PanGenie", sample_name, pangenie_vcf_file
    """

    # Create folder and switch paths
    os.chdir(software_work_path)
    pangenie_path = os.path.join(software_work_path, sample_name)
    makedir(
        pangenie_path,
        force_tmp=force,
        restart_tmp=restart
    )
    os.chdir(pangenie_path)

    # genotype
    stdout, stderr, log_out, pangenie_vcf_file = PanGenie.main(
        reference_file,
        vcf_file,
        sample_name,
        fastq_file1,
        fastq_file2,
        env_path, 
        threads,
        restart
    )

    return stdout, stderr, log_out, "PanGenie", sample_name, pangenie_vcf_file


# Conversion of reference genome and vcf files
def ref_vcf_convert(
        parser_map,
        work_path
):
    """
    :param parser_map: get_parser(), Returned parameter dictionary
    :param work_path:  work path
    :return: stdout, stderr, log_out, convert_out_map = {
                "reference_file": reference_file,
                "vcf_out_name": vcf_out_name,
                "vcf_sv_out_name": vcf_sv_out_name,
                "region_file": region_file,
                "fasta_base": fasta_base,
                "map_index_dir": map_index_dir,
                "giraffe_index_dir": giraffe_index_dir
            }
    """
    # switch paths
    os.chdir(work_path)

    # ################################### reference to convert ###################################
    stdout, stderr, log_out, reference_file = convert.convert_reference(
        parser_map["reference_file"],
        env_path, 
        restart
    )
    # Report an error if there is a problem with the exit code
    if log_out:
        return stdout, stderr, log_out, {}
    # Evaluate the fasta file
    stdout, stderr, log_out, fasta_base = convert.fasta_count(code_dir, env_path, reference_file)
    # Report an error if there is a problem with the exit code
    if log_out:
        return stdout, stderr, log_out, {}

    # ################################### Build an index on the reference ###################################
    stdout, stderr, log_out = bwa.index(
        reference_file, 
        env_path, 
        restart
    )
    # Report an error if there is a problem with the exit code
    if log_out:
        return stdout, stderr, log_out, {}

    # ################################### Convert vcf files ###################################
    stdout, stderr, log_out, vcf_out_name, vcf_sv_out_name, region_file = convert.vcf_convert(
        code_dir,
        reference_file,
        parser_map["vcf_file"],
        350,
        parser_map["vcf_out_name"],
        mode, 
        env_path, 
        restart
    )
    # Report an error if there is a problem with the exit code
    if log_out:
        return stdout, stderr, log_out, {}

    # ################################### Compress vcf files ###################################
    if vcf_sv_out_name != vcf_out_name:  # If the two files are the same, vcf_sv_out_name will not be compressed, otherwise an error will be reported
        stdout, stderr, log_out, vcf_sv_out_name = convert.bgzip_vcf(
            vcf_sv_out_name,
            env_path, 
            restart
        )
        # Report an error if there is a problem with the exit code
        if log_out:
            return stdout, stderr, log_out, {}
        stdout, stderr, log_out, vcf_out_name = convert.bgzip_vcf(
            vcf_out_name,
            env_path, 
            restart
        )
        # Report an error if there is a problem with the exit code
        if log_out:
            return stdout, stderr, log_out, {}
    else:
        stdout, stderr, log_out, vcf_out_name = convert.bgzip_vcf(
            vcf_out_name,
            env_path, 
            restart
        )
        vcf_sv_out_name = vcf_out_name
        # Report an error if there is a problem with the exit code
        if log_out:
            return stdout, stderr, log_out, {}

    # ################################### Build vg graph and index ###################################
    # First select the vg software
    select_software_list = []
    if isinstance(parser_map["software_list"], str):  # Program Automatic Judgment Software
        if fasta_base > 200000000:  # If the genome is larger than 200Mb, use giraffe to run
            select_software_list.append("giraffe")
        else:  # Otherwise run with vg_map
            select_software_list.append("map")
    else:  # user-defined software
        if "VG-MAP" in parser_map["software_list"]:
            select_software_list.append("map")
        if "VG-Giraffe" in parser_map["software_list"]:
            select_software_list.append("giraffe")

    # multi-process process pool
    if len(select_software_list) > 0:  # If the number of files is greater than 1, then build the index
        # multi-Progress
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
            to_do = []
            # Map and giraffe build indexes
            for software in select_software_list:
                # switch paths
                index_dir = os.path.join(work_path, software)
                makedir(
                    index_dir,
                    force_tmp=force,
                    restart_tmp=restart
                )
                os.chdir(index_dir)

                # submit task
                future = executor.submit(
                    VG_MAP_Giraffe.vg_autoindex,
                    reference_file,
                    vcf_out_name,
                    software, 
                    env_path, 
                    threads,
                    index_dir,
                    restart
                )
                # save result
                to_do.append(future)

                time.sleep(0.05)

            # Back to the main working path
            os.chdir(work_path)

            # GraphAligner build index
            if "GraphAligner" in parser_map["software_list"]:
                # switch paths
                index_dir = os.path.join(work_path, "GraphAligner")
                makedir(
                    index_dir,
                    force_tmp=force,
                    restart_tmp=restart
                )
                os.chdir(index_dir)

                # submit task
                future = executor.submit(
                    GraphAligner.vg_index,
                    reference_file,
                    vcf_out_name,
                    threads,
                    index_dir, 
                    env_path, 
                    restart
                )
                # save result
                to_do.append(future)

            # return to main path
            os.chdir(work_path)

            # Obtain the return value of multithreading. If log_out exists, it indicates that the operation failed and the exit code
            for future in concurrent.futures.as_completed(to_do):  # concurrent execution
                stdout, stderr, log_out = future.result()  # check return value
                if log_out:
                    return stdout, stderr, log_out, {}

    # path to the index file
    map_index_dir = os.path.join(work_path, "map")
    giraffe_index_dir = os.path.join(work_path, "giraffe")
    graphaligner_index_dir = os.path.join(work_path, "GraphAligner")

    # The path to store the converted file
    convert_out_map = {
        "reference_file": reference_file,
        "vcf_out_name": vcf_out_name,
        "vcf_sv_out_name": vcf_sv_out_name,
        "region_file": region_file,
        "fasta_base": fasta_base,
        "map_index_dir": map_index_dir,
        "giraffe_index_dir": giraffe_index_dir,
        "GraphAligner_index_dir": graphaligner_index_dir
    }

    return stdout, stderr, log_out, convert_out_map


# Convert sequencing files
def read_convert(
    work_path: str,
    sample_name: str,
    fastq_file1: str,
    fastq_file2: str,
    fasta_base: float
):
    """
    :param work_path:    work path
    :param sample_name:  sample name
    :param fastq_file1:  read1
    :param fastq_file2:  read2
    :param fasta_base:   the base of reference genome
    :return: stdout, stderr, log_out, fastq_info_map = {
                "sample_name": sample_name,
                "fastq_file1": fastq_out_file1,
                "fastq_file2": fastq_out_file2,
                "real_depth": real_depth,
                "read_len": read_len
            }
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # switch paths
    os.chdir(work_path)

    # ################################### Transform the second-generation data ###################################
    # Evaluate fastq files
    if fastq_file1:
        stdout, stderr, log_out, fastq_base, read_num, read_len = convert.fastq_count(
            code_dir, 
            env_path, 
            fastq_file1,
            fastq_file2
        )
        # Report an error if there is a problem with the exit code
        if log_out:
            return stdout, stderr, log_out, {}

        # downsample
        stdout, stderr, log_out, fastq_out_file1, fastq_out_file2 = convert.downsample(
            code_dir,
            fastq_file1,
            fastq_file2,
            fastq_base,
            fasta_base,
            need_depth, 
            env_path, 
            restart
        )
        # Report an error if there is a problem with the exit code
        if log_out:
            return stdout, stderr, log_out, {}

        # actual sequencing depth
        real_depth = min(int(fastq_base) / int(fasta_base), need_depth)
    else:
        read_len = 0
        fastq_out_file1, fastq_out_file2 = "", ""
        real_depth = 0

    # Determine whether the sequencing depth is 0
    if real_depth == 0:
        log_out = '[EVG.run_convert] sequencing depth is 0, please check the input file file.'
        return stdout, stderr, log_out, {}

    # fastq path, sequencing depth, read length
    fastq_info_map = {
        "sample_name": sample_name,
        "fastq_file1": fastq_out_file1,
        "fastq_file2": fastq_out_file2,
        "real_depth": real_depth,
        "read_len": read_len
    }

    return stdout, stderr, log_out, fastq_info_map


# Each line runs independently
def run_genotype(
    work_path: str,
    parser_map,
    convert_out_map,
    bam_infos_map
):
    """
    :param work_path:        work path
    :param parser_map:       the return of get_parser
    :param convert_out_map:  Hash table after conversion of vcf, reference and sequencing files
    :param bam_infos_map:    run_bwa return value
    :return: stdout, stderr, log_out, vcf_merge_files_map
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # switch paths
    os.chdir(work_path)

    # Minimum sequencing depth and sequencing length for obtaining sequencing data
    depth_min = 1000
    read_len_mix = 100000
    for key, value in bam_infos_map.items():
        depth_min = min(depth_min, value["real_depth"])
        read_len_mix = min(read_len_mix, value["read_len"])

    # run select_software
    if isinstance(parser_map["software_list"], str):  # Program Automatic Judgment Software
        select_software_list = select_software.main(depth_min, read_len_mix, convert_out_map["fasta_base"])
    else:  # user-defined software
        select_software_list = parser_map["software_list"]

    # ################################### config ###################################
    # Generate configuration files for GraphTyper2, BayesTyper and Paragraph
    bam2graphtyper_file = ""
    bam2bayestyper_file = ""
    bam2paragraph_file = ""
    if "GraphTyper2" in select_software_list:
        bam2graphtyper_file = bwa.bam2graphtyper(
            bam_infos_map
        )
    if "BayesTyper" in select_software_list:
        bam2bayestyper_file = bwa.bam2bayestyper(
            bam_infos_map
        )
    if "Paragraph" in select_software_list:
        bam2paragraph_file = bwa.bam2paragraph(
            bam_infos_map
        )

    # ################################### genotype ###################################
    # run genotype
    # multi-process process pool
    pool = Pool(processes=jobs_num)

    # thread pool return value
    pool_out_list = []  # Get results from asynchronously submitted tasks

    # Save the location of typing results
    genotype_outs_map = {}  # map<sample_name, map<software, vcf_path>>

    for software in select_software_list:
        # back to working path
        os.chdir(work_path)

        if software in ["VG-MAP", "VG-Giraffe", "GraphAligner", "PanGenie"]:
            # Create folder and switch paths
            os.chdir(work_path)
            software_work_path = os.path.join(work_path, software)
            makedir(
                software_work_path,
                force_tmp=force,
                restart_tmp=restart
            )
            os.chdir(software_work_path)

            for key, value in bam_infos_map.items():
                # back to working path
                os.chdir(work_path)

                sample_name = key
                fastq_file1 = value["fastq_file1"]
                fastq_file2 = value["fastq_file2"]
                real_depth = value["real_depth"]

                if "VG-MAP" == software:
                    # Multi-process submission tasks
                    pool_out = pool.apply_async(run_vg_map_giraffe, args=(
                        sample_name,
                        fastq_file1,
                        fastq_file2,
                        real_depth,
                        software,
                        software_work_path,
                        convert_out_map["map_index_dir"],
                    ), error_callback=throw_exception)
                    # Save the return value of multiple threads
                    pool_out_list.append(pool_out)
                elif "VG-Giraffe" == software:
                    # Multi-process submission tasks
                    pool_out = pool.apply_async(run_vg_map_giraffe, args=(
                        sample_name,
                        fastq_file1,
                        fastq_file2,
                        real_depth,
                        software,
                        software_work_path,
                        convert_out_map["giraffe_index_dir"],
                    ), error_callback=throw_exception)
                    # Save the return value of multiple threads
                    pool_out_list.append(pool_out)
                elif "GraphAligner" == software:
                    # Multi-process submission tasks
                    pool_out = pool.apply_async(run_graphaligner, args=(
                        sample_name,
                        fastq_file1,
                        real_depth,
                        software_work_path,
                        convert_out_map["GraphAligner_index_dir"],
                    ), error_callback=throw_exception)
                    # Save the return value of multiple threads
                    pool_out_list.append(pool_out)
                else:
                    # Multi-process submission tasks
                    pool_out = pool.apply_async(run_pangenie, args=(
                        sample_name,
                        convert_out_map["reference_file"],
                        convert_out_map["vcf_out_name"],
                        fastq_file1,
                        fastq_file2,
                        software_work_path,
                    ), error_callback=throw_exception)
                    # Save the return value of multiple threads
                    pool_out_list.append(pool_out)
                # thread interval
                time.sleep(0.05)
                # back to working path
                os.chdir(work_path)
        elif software == "Paragraph":
            # Multi-process submission tasks
            pool_out = pool.apply_async(run_paragraph, args=(
                convert_out_map["reference_file"],
                convert_out_map["vcf_sv_out_name"],
                bam2paragraph_file,
                work_path,
            ), error_callback=throw_exception)
            # Save the return value of multiple threads
            pool_out_list.append(pool_out)
        elif software == "BayesTyper":
            # Multi-process submission tasks
            pool_out = pool.apply_async(run_bayestyper, args=(
                convert_out_map["reference_file"],
                convert_out_map["vcf_out_name"],
                bam2bayestyper_file,
                work_path,
                bam_infos_map
            ), error_callback=throw_exception)
            # Save the return value of multiple threads
            pool_out_list.append(pool_out)
        elif software == "GraphTyper2":
            # Multi-process submission tasks
            pool_out = pool.apply_async(run_graphtyper, args=(
                convert_out_map["reference_file"],
                convert_out_map["vcf_out_name"],
                convert_out_map["region_file"],
                bam2graphtyper_file,
                work_path,
            ), error_callback=throw_exception)
            # Save the return value of multiple threads
            pool_out_list.append(pool_out)

        # thread interval
        time.sleep(0.05)
        # back to working path
        os.chdir(work_path)

    # Get multithreaded return value
    for pool_out in pool_out_list:
        stdout, stderr, log_out, software, sample_name, vcf_file = pool_out.get()

        # Report an error if there is a problem with the exit code
        if log_out:
            return stdout, stderr, log_out, {}

        if software in ["VG-MAP", "VG-Giraffe", "GraphAligner", "PanGenie"]:
            if software not in genotype_outs_map:
                genotype_outs_map[software] = {}
            genotype_outs_map[software][sample_name] = vcf_file
        else:
            genotype_outs_map[software] = vcf_file

        # back to working path
        os.chdir(work_path)

    # Close the thread pool
    pool.close()
    pool.join()

    # ################################################ graphvcf merge ################################################
    # multi-process process pool
    pool = Pool(processes=jobs_num*threads)  # The merge step uses less resources, so all cores are used

    # thread pool return value
    pool_out_list = []  # Get results from asynchronously submitted tasks

    # Store the merged result
    vcf_merge_files_map = {}  # map<sample_name, vcf_file>

    # merged result
    sample_name_tmp_list = bam_infos_map.keys()  # Get all samples_name first
    for sample_name_tmp in sample_name_tmp_list:
        vcf_out_tmp_list = []  # Temporary list for submitting tasks
        software_tmp_list = []  # Temporary list for submitting tasks
        for software in select_software_list:  # software not in results file
            if software not in genotype_outs_map.keys():
                log = '[EVG.genotype] {}: {} results are missing, skipped.'.format(sample_name_tmp, software)
                sys.stdout.write(log)
                continue
            software_tmp_list.append(software)
            if isinstance(genotype_outs_map[software], dict):
                if sample_name_tmp not in genotype_outs_map[software].keys():  # sample_name is not in the result file
                    log = '[EVG.genotype] {}: {} results are missing, skipped.'.format(sample_name_tmp, software)
                    sys.stdout.write(log)
                    continue
                vcf_out_tmp_list.append(genotype_outs_map[software][sample_name_tmp])
            else:
                vcf_out_tmp_list.append(genotype_outs_map[software])

        # sample_name is not in the result file, graphvcf merge
        pool_out = pool.apply_async(
            merge.main, args=(
                code_dir,
                work_path,
                parser_map["line_vcf_file"],
                sample_name_tmp,
                merge_mode,
                vcf_out_tmp_list,
                software_tmp_list,
                env_path, 
                restart,
            ), error_callback=throw_exception
        )

        # Save the return value of multiple threads
        pool_out_list.append(pool_out)

        # Multi-thread interval
        time.sleep(0.05)

    for pool_out in pool_out_list:  # The return value of read
        stdout, stderr, log_out, sample_name, merge_vcf_file = pool_out.get()

        # Report an error if there is a problem with the exit code
        if log_out:
            return stdout, stderr, log_out, {}

        # Assigned to the total hash table
        vcf_merge_files_map[sample_name] = merge_vcf_file

        # Close the thread pool
    pool.close()
    pool.join()

    return stdout, stderr, log_out, vcf_merge_files_map


# capture subprocess exit status
def throw_exception(name):
    # print log
    log = '[EVG.genotype] %s' % name.__cause__
    logger.error(log)
    # os.killpg(os.getpgid(os.getpid()), signal.SIGKILL)
    raise SystemExit(1)


def main():
    # parsing parameters
    parser_map = get_parser()

    # ################################################ EVG convert ################################################
    # Return to the main working path
    os.chdir(parser_map["path"])
    # multi-process process pool
    pool = Pool(processes=threads)

    # thread pool return value
    pool_out_list = []  # Get results from asynchronously submitted tasks

    # Store multi-threaded results, fastq path, sequencing depth, read length
    fastq_infos_map = {}

    # Get sequencing file path and output pathname
    samples_map = get_samples_path(parser_map["samples_file"])

    # reference and vcf, Convert storage paths and create folders
    convert_dir = os.path.join(parser_map["path"], "convert")  # Conversion file storage path
    makedir(
        convert_dir,
        force_tmp=force,
        restart_tmp=restart
    )

    # Return to the main working path
    os.chdir(parser_map["path"])
    # Evaluate the fasta file
    stdout, stderr, log_out, fasta_base = convert.fasta_count(
        code_dir, 
        env_path, 
        parser_map["reference_file"]
    )
    # Report an error if there is a problem with the exit code
    if log_out:
        logger.error(log_out.strip())
        exit(1)

    # Convert reference and vcf
    # Multi-thread submission
    pool_convert_out = pool.apply_async(
        ref_vcf_convert, args=(
            parser_map,
            convert_dir,
        ), error_callback=throw_exception
    )

    # convert read
    # Multi-thread submission
    for key, value in samples_map.items():
        # Return to the main working path
        os.chdir(parser_map["path"])

        sample_name = key
        fastq_file1 = value[0]
        fastq_file2 = value[1]

        # Multi-threaded conversion of sequencing files
        pool_out = pool.apply_async(read_convert, args=(
            convert_dir,
            sample_name,
            fastq_file1,
            fastq_file2,
            fasta_base,
        ), error_callback=throw_exception)

        # Save the return value of multiple threads
        pool_out_list.append(pool_out)

        # Multi-thread interval
        time.sleep(0.05)

    # Get multithreaded return value
    stdout, stderr, log_out, convert_out_map = pool_convert_out.get()  # vcf和reference的返回�?
    # Report an error if there is a problem with the exit code
    if log_out:
        logger.error(log_out.strip())
        exit(1)
    for pool_out in pool_out_list:  # The return value of read
        stdout, stderr, log_out, fastq_info_map = pool_out.get()
        # Report an error if there is a problem with the exit code
        if log_out:
            logger.error(log_out.strip())
            exit(1)
        # Assigned to the total hash table
        fastq_infos_map[fastq_info_map["sample_name"]] = {
            "fastq_file1": fastq_info_map["fastq_file1"],
            "fastq_file2": fastq_info_map["fastq_file2"],
            "real_depth": fastq_info_map["real_depth"],
            "read_len": fastq_info_map["read_len"]
        }

    # Close the thread pool
    pool.close()
    pool.join()

    # ################################################ bwa mem ################################################
    # Return to the main working path
    os.chdir(parser_map["path"])
    # Reinitialize the thread pool
    pool = Pool(processes=jobs_num)

    # thread pool return value
    pool_out_list = []  # Get results from asynchronously submitted tasks

    # Store multi-threaded results, bam path, sequencing depth, read length
    bam_infos_map = {}

    # Sequence comparison storage path
    bwa_dir = os.path.join(parser_map["path"], "bwa")  # the path of bwa mem

    # Create a directory
    makedir(
        bwa_dir,
        force_tmp=force,
        restart_tmp=restart
    )

    # Multi-threaded sequence alignment
    for key, value in fastq_infos_map.items():
        sample_name = key
        fastq_file1 = value["fastq_file1"]
        fastq_file2 = value["fastq_file2"]
        real_depth = value["real_depth"]
        read_len = value["read_len"]

        # bwa mem
        # Multi-threaded conversion of sequencing files
        pool_out = pool.apply_async(run_bwa, args=(
            convert_out_map["reference_file"],
            sample_name,
            fastq_file1,
            fastq_file2,
            real_depth,
            read_len,
            bwa_dir
        ), error_callback=throw_exception)

        # Save the return value of multiple threads
        pool_out_list.append(pool_out)

        # Multi-thread interval
        time.sleep(0.05)

    # Get multithreaded return value
    for pool_out in pool_out_list:
        stdout, stderr, log_out, bam_info_map = pool_out.get()
        # Report an error if there is a problem with the exit code
        if log_out:
            logger.error(log_out.strip())
            exit(1)
        # Assigned to the total hash table
        bam_infos_map[bam_info_map["sample_name"]] = {
            "fastq_file1": bam_info_map["fastq_file1"],
            "fastq_file2": bam_info_map["fastq_file2"],
            "real_depth": bam_info_map["real_depth"],
            "read_len": bam_info_map["read_len"],
            "bam_file": bam_info_map["bam_file"]
        }

    # Close the thread pool
    pool.close()
    pool.join()

    # ################################################ genotype ################################################
    # Return to the main working path
    os.chdir(parser_map["path"])

    # Typing result path
    genotype_dir = os.path.join(parser_map["path"], "genotype")  # the path of bwa mem
    # Create a directory
    makedir(
        genotype_dir,
        force_tmp=force,
        restart_tmp=restart
    )
    # genotyping
    stdout, stderr, log_out, vcf_merge_files_map = run_genotype(
        genotype_dir,
        parser_map,
        convert_out_map,
        bam_infos_map
    )
    # Report an error if there is a problem with the exit code
    if log_out:
        logger.error(log_out.strip())
        exit(1)

    # printout result
    print_merge(vcf_merge_files_map)

    return 0


if __name__ == '__main__':
    main()
