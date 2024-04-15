#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import os
import run_cmd
from getsize import getsize


# bwa index
def index(
    software_work_path, 
    reference_file: str, 
    env_path, 
    restart: bool
):
    stdout = stderr = log_out = ""

    os.chdir(software_work_path)

    cmd = f"bwa index {reference_file} -p {os.path.join(software_work_path, reference_file)}"

    # Check if the file exists
    if restart:
        if getsize(f"{reference_file}.amb") <= 0 or \
                getsize(f"{reference_file}.ann") <= 0 or \
                getsize(f"{reference_file}.bwt") <= 0 or \
                getsize(f"{reference_file}.fai") <= 0 or \
                getsize(f"{reference_file}.pac") <= 0 or \
                getsize(f"{reference_file}.sa") <= 0:
            stdout, stderr, log_out = run_cmd.run(cmd, env_path)
    else:  # If restart is not specified, run directly
        stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    return stdout, stderr, log_out

# bwa mem
def mem(
    bwa_dir, 
    reference_file: str,
    threads: int,
    sample_name: str,
    env_path, 
    restart: bool,
    read1: str,
    read2: str = ""
):
    stdout = stderr = log_out = ""
    
    # change working path
    os.chdir(bwa_dir)

    # bam file path
    bam_file = os.path.join(bwa_dir, sample_name + ".bam")
    cmd = f"bwa mem -R '@RG\\tID:foo\\tSM:{sample_name}\\tLB:library1' -t {threads} {reference_file} {read1} {read2} | samtools view -b -S | samtools sort -@ {threads} -o {bam_file} && samtools index {bam_file}"

    # Check if the file exists
    if restart:
        file_size = getsize(bam_file)
        if file_size <= 0:
            stdout, stderr, log_out = run_cmd.run(cmd, env_path)
    else:  # If restart is not specified, run directly
        stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    return stdout, stderr, log_out, sample_name, bam_file

# save result for GraphTyper2
def bam2graphtyper(read_infos_map):
    out_txt = ""
    for _, value in read_infos_map.items():
        out_txt += value["bam"] + "\n"
    with open("bam_for_GraphTyper2.txt", "w") as f:
        f.write(out_txt)
    return os.path.abspath("bam_for_GraphTyper2.txt")

# save result for BayesTyper
def bam2bayestyper(read_infos_map, number: int):
    work_dir = os.getcwd()

    out_samplename_list = []
    out_txt_list = []
    out_txt = ""
    
    # 'number' samples can be processed each time
    num = 0
    num_tmp = 0
    for sample_name in sorted(read_infos_map.keys()):
        value = read_infos_map[sample_name]
        num += 1
        out_txt += sample_name + "\tF\t" + os.path.join(work_dir, "BayesTyper", str(num_tmp), os.path.basename(value["bam"])) + "\n"
        out_samplename_list.append(sample_name)
        if num == number:
            out_txt_list.append(out_txt)
            out_txt = ""
            num = 0
            num_tmp += 1

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
def bam2paragraph(read_infos_map, number: int):
    out_samplename_list = []
    out_txt_list = []
    head_txt = "id\tpath\tdepth\tread length\n"

    out_txt = ""

    # 'number' samples can be processed each time
    num = 0
    for sample_name in sorted(read_infos_map.keys()):
        value = read_infos_map[sample_name]
        num += 1
        out_txt += sample_name + "\t" + value["bam"] + "\t" + str(value["depth"]) + "\t" + str(value["length"]) + "\n"
        out_samplename_list.append(sample_name)

        if num == number:
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
