#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import os
import run_cmd
from getsize import getsize


# bwa index
def index(
        reference_file: str, 
        env_path, 
        restart: bool
):
    """
    :param reference_file:  reference genome
    :param env_path:        env path
    :param restart:         Whether to check if the file exists and skip this step
    :return: stdout, stderr, log_out
    """

    print(reference_file)

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # bwa index building and comparison
    # bam file path
    cmd = "bwa index {}".format(reference_file)

    # Check if the file exists
    if restart:
        # <= 0
        if getsize("{}.amb".format(reference_file)) <= 0 or \
                getsize("{}.ann".format(reference_file)) <= 0 or \
                getsize("{}.bwt".format(reference_file)) <= 0 or \
                getsize("{}.fai".format(reference_file)) <= 0 or \
                getsize("{}.pac".format(reference_file)) <= 0 or \
                getsize("{}.sa".format(reference_file)) <= 0:
            # submit task
            stdout, stderr, log_out = run_cmd.run(cmd, "bwa.index", env_path)
    else:  # If restart is not specified, run directly
        # submit task
        stdout, stderr, log_out = run_cmd.run(cmd, "bwa.index", env_path)

    return stdout, stderr, log_out


# bwa mem
def mem(
        reference_file: str,
        threads: int,
        sample_name: str,
        env_path, 
        restart: bool,
        fastq_file1: str,
        fastq_file2: str = ""
):
    """
    :param reference_file: reference genome
    :param threads:        Threads
    :param sample_name:    sample name
    :param env_path:       env path
    :param restart: Whether to check if the file exists and skip this step
    :param fastq_file1:    read1
    :param fastq_file2:    read2
    :return:
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # bwa index building and comparison
    # bam file path
    bam_file = os.path.abspath(sample_name + ".bam")
    cmd = "bwa mem -R '@RG\\tID:foo\\tSM:{}\\tLB:library1' -t {} {} {} {} | samtools view -b -S | samtools " \
          "sort -@ {} -o {} && " \
          "samtools index {}".\
        format(sample_name, threads, reference_file, fastq_file1,
               fastq_file2, threads, bam_file, bam_file)

    # Check if the file exists
    if restart:
        # check file
        file_size = getsize(
            bam_file
        )
        # <= 0
        if file_size <= 0:
            # submit task
            stdout, stderr, log_out = run_cmd.run(cmd, "bwa.mem", env_path)
    else:  # If restart is not specified, run directly
        # submit task
        stdout, stderr, log_out = run_cmd.run(cmd, "bwa.mem", env_path)

    return stdout, stderr, log_out, bam_file


# save result for GraphTyper2
def bam2graphtyper(
        bam_infos_map
):
    out_txt = ""
    for key, value in bam_infos_map.items():
        out_txt += value["bam_file"] + "\n"

    # remove special characters
    out_txt = out_txt.strip()

    with open("bam_for_GraphTyper2.txt", "w") as f:
        f.write(out_txt)

    return os.path.abspath("bam_for_GraphTyper2.txt")


# save result for BayesTyper
def bam2bayestyper(
        bam_infos_map
):
    out_txt = ""

    for key, value in bam_infos_map.items():
        out_txt += key + "\tF\t" + os.path.basename(value["bam_file"]) + "\n"

    # remove special characters
    out_txt = out_txt.strip()

    with open("bam_for_BayesTyper.tsv", "w") as f:
        f.write(out_txt)

    return os.path.abspath("bam_for_BayesTyper.tsv")


# save result for ParaGraph
def bam2paragraph(
        bam_infos_map
):
    out_txt = "id\tpath\tdepth\tread length\n"

    for key, value in bam_infos_map.items():
        out_txt += key + \
                   "\t" + \
                   value["bam_file"] + \
                   "\t" + \
                   str(value["real_depth"]) + \
                   "\t" + \
                   str(value["read_len"]) + \
                   "\n"

    # remove special characters
    out_txt = out_txt.strip()

    with open("bam_for_Paragraph.txt", "w") as f:
        f.write(out_txt)

    return os.path.abspath("bam_for_Paragraph.txt")
