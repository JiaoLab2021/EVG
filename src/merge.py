#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import os
import run_cmd
from getsize import getsize


# 合并结果
def main(
    code_dir: str,
    work_path: str,
    line_vcf_file: str,
    sample_name: str,
    mode: str,
    genotype_vcf_file_list,
    select_software_list,
    env_path, 
    restart: bool
):
    """
    :param code_dir: code directory
    :param work_path: work path
    :param line_vcf_file: line的vcf
    :param sample_name: sample name
    :param mode: mode
    :param genotype_vcf_file_list:
    :param select_software_list:
    :param env_path: environment variable
    :param restart: Whether to check if the file exists and skip this step
    :return: stdout, stderr, log_out, sample_name, merge_vcf_file
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # output file
    merge_vcf_file = "{}.vcf.gz".format(sample_name)

    # back to working directory
    os.chdir(work_path)

    # code path
    code_path = os.path.join(code_dir, "graphvcf merge")

    cmd = code_path + " -v " + line_vcf_file + " "
    for index in range(len(genotype_vcf_file_list)):
        file = genotype_vcf_file_list[index]
        software = select_software_list[index]
        cmd += "--" + software + " " + file + " "
    cmd += "-n {} -m {} -o {}".format(sample_name, mode, merge_vcf_file)

    # Check if the file exists
    if restart:
        # check file
        file_size = getsize(
            "{}.vcf.gz".format(sample_name)
        )
        # <= 0
        if file_size <= 0:
            # submit task
            stdout, stderr, log_out = run_cmd.run(cmd, "graphvcf merge", env_path)
    else:  # If restart is not specified, run directly
        # submit task
        stdout, stderr, log_out = run_cmd.run(cmd, "graphvcf merge", env_path)

    merge_vcf_file = os.path.abspath(merge_vcf_file)

    return stdout, stderr, log_out, sample_name, merge_vcf_file
