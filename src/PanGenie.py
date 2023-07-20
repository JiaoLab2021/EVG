#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import os
import run_cmd
from getsize import getsize


# check file
def check(
    vcf_file: str,
    sample_name: str,
    fastq_file1: str,
    fastq_file2: str,
    env_path, 
    restart: bool
):
    """
    :param vcf_file: the path of vcf file
    :param sample_name sample name
    :param fastq_file1: read1
    :param fastq_file2: read2
    :param env_path: environment variable
    :param restart: Whether to check if the file exists and skip this step
    :return: stdout, stderr, log_out, vcf_out_file, fastq_out_file
    """

    #  initial log
    stdout = ""
    stderr = ""
    log_out = ""

    # output file path
    vcf_out_file = os.path.abspath(sample_name + ".vcf")
    fastq_out_file = os.path.abspath("{}.fq".format(sample_name))

    if '.gz' in vcf_file or ".GZ" in vcf_file:  # If it is a compressed file, decompress it
        cmd = f"gunzip -c {vcf_file} 1>{vcf_out_file}"

        # Check if the file exists
        if restart:
            # check file
            file_size = getsize(
                vcf_out_file
            )
            # <= 0
            if file_size <= 0:
                # submit task
                stdout, stderr, log_out = run_cmd.run(cmd, "PanGenie.gunzip", env_path)
        else:  # If restart is not specified, run directly
            # submit task
            stdout, stderr, log_out = run_cmd.run(cmd, "PanGenie.gunzip", env_path)

        # Report an error if there is a problem with the exit code
        if log_out:
            return stdout, stderr, log_out, "", ""
    else:  # Otherwise assign the same path
        vcf_out_file = vcf_file

    cmd = ""
    if '.gz' in fastq_file1 or ".GZ" in fastq_file1:  # fastq1
        cmd = f"gunzip -c {fastq_file1} 1>{fastq_out_file}"
        if fastq_file2:  # fastq2
            if '.gz' in fastq_file2 or ".GZ" in fastq_file2:
                cmd += f" && gunzip -c {fastq_file2} 1>>{fastq_out_file}"
            else:
                cmd += f" && cat {fastq_file2} 1>>{fastq_out_file}"
    else:  # not a compressed file
        cmd = f"cat {fastq_file1} 1>{fastq_out_file}"
        if fastq_file2:  # fastq2
            if '.gz' in fastq_file2 or ".GZ" in fastq_file2:
                cmd += f" && gunzip -c {fastq_file2} 1>>{fastq_out_file}"
            else:
                cmd += f" && cat {fastq_file2} 1>>{fastq_out_file}"
    # Check if the file exists
    if restart:
        # check file
        file_size = getsize(
            fastq_out_file
        )
        # <= 0
        if file_size <= 0:
            # submit task
            stdout, stderr, log_out = run_cmd.run(cmd, "PanGenie.gunzip", env_path)
    else:  # If restart is not specified, run directly
        # submit task
        stdout, stderr, log_out = run_cmd.run(cmd, "PanGenie.gunzip", env_path)

    # Report an error if there is a problem with the exit code
    if log_out:
        return stdout, stderr, log_out, "", ""

    return stdout, stderr, log_out, vcf_out_file, fastq_out_file


# genotype
def main(
    reference_file: str,
    vcf_file: str,
    sample_name: str,
    fastq_file1: str,
    fastq_file2: str,
    env_path, 
    threads: int,
    restart: bool
):
    """
    :param reference_file: reference genome
    :param vcf_file: vcf file
    :param sample_name sample name
    :param fastq_file1: read1
    :param fastq_file2: read2
    :param env_path: environment variable
    :param threads: Threads
    :param restart: Whether to check if the file exists and skip this step
    :return:
    """

    # Check if files are compressed and merge
    stdout, stderr, log_out, vcf_file, fastq_file = check(
        vcf_file,
        sample_name,
        fastq_file1,
        fastq_file2,
        env_path, 
        restart
    )
    # Report an error if there is a problem with the exit code
    if log_out:
        return stdout, stderr, log_out, ""

    cmd = f"PanGenie -s {sample_name} -i {fastq_file} -r {reference_file} -v {vcf_file} -t {threads} -j {threads} -o {sample_name}"

    # output vcf path
    out_vcf_file = os.path.abspath(f"{sample_name}_genotyping.vcf")

    # Check if the file exists
    if restart:
        # check file
        file_size = getsize(
            out_vcf_file
        )
        # <= 0
        if file_size <= 0:
            # submit task
            stdout, stderr, log_out = run_cmd.run(cmd, "PanGenie", env_path)
    else:  # If restart is not specified, run directly
        # submit task
        stdout, stderr, log_out = run_cmd.run(cmd, "PanGenie", env_path)

    return stdout, stderr, log_out, os.path.abspath(out_vcf_file)
