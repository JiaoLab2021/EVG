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

    with open("bam_for_GraphTyper2.txt", "w") as f:
        f.write(out_txt)

    return os.path.abspath("bam_for_GraphTyper2.txt")


# save result for BayesTyper
def bam2bayestyper(
    bam_infos_map
):
    out_samplename_list = []
    out_txt_list = []
    out_txt = ""
    
    # Only 30 samples can be processed each time
    num = 0
    for key, value in bam_infos_map.items():
        num += 1
        out_txt += key + "\tF\t" + os.path.basename(value["bam_file"]) + "\n"
        out_samplename_list.append(key)

        if num == 30:
            out_txt_list.append(out_txt)
            out_txt = ""
            num = 0

    # last save
    if num > 0:
        out_txt_list.append(out_txt)
        out_txt = ""
        num = 0

    out_file_list = []
    num = 0
    for out_txt in out_txt_list:
        # save the file name
        out_file = os.path.abspath(f"bam_for_BayesTyper_{num}.tsv")
        out_file_list.append(out_file)
        num += 1
        # save to file
        with open(out_file, "w") as f:
            f.write(out_txt)

    return out_samplename_list, out_file_list


# save result for ParaGraph
def bam2paragraph(
    bam_infos_map
):
    out_samplename_list = []
    out_txt_list = []
    head_txt = "id\tpath\tdepth\tread length\n"

    out_txt = ""

    # Only 10 samples can be processed each time
    num = 0
    for key, value in bam_infos_map.items():
        num += 1
        out_txt += key + \
                   "\t" + \
                   value["bam_file"] + \
                   "\t" + \
                   str(value["real_depth"]) + \
                   "\t" + \
                   str(value["read_len"]) + \
                   "\n"
        out_samplename_list.append(key)

        if num == 10:
            out_txt_list.append(out_txt)
            out_txt = ""
            num = 0

    # last save
    if num > 0:
        out_txt_list.append(out_txt)
        out_txt = ""
        num = 0

    out_file_list = []
    num = 0
    for out_txt in out_txt_list:
        # save the file name
        out_file = os.path.abspath(f"bam_for_Paragraph_{num}.txt")
        out_file_list.append(out_file)
        num += 1
        # save to file
        with open(out_file, "w") as f:
            f.write(head_txt + out_txt)

    return out_samplename_list, out_file_list
