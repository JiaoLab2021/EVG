#!/usr/bin/python3

# -*- coding: utf-8 -*-
import os
import run_cmd
from getsize import getsize


# merge vcf
def merge_vcf(
    env_path, 
    restart: bool
):
    """
    :param env_path: environment variable
    :param restart: Whether to check if the file exists and skip this step
    :return: stdout, stderr, log_out, vcf_file
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    cmd = "for i in $(ls */ | grep 'vcf.gz' | grep -v 'tbi' | head -n1); do zcat */$i | grep '#' > graphtyper.vcf; done && " \
          "zcat */*.vcf.gz >> graphtyper.vcf"

    # Check if the file exists
    if restart:
        # check file
        file_size = getsize(
            "graphtyper.vcf"
        )
        # <= 0
        if file_size <= 0:
            # submit task
            stdout, stderr, log_out = run_cmd.run(cmd, "GraphTyper2.genotype_merge", env_path)
    else:  # If restart is not specified, run directly
        # submit task
        stdout, stderr, log_out = run_cmd.run(cmd, "GraphTyper2.genotype_merge", env_path)

    vcf_file = os.path.abspath("graphtyper.vcf")

    return stdout, stderr, log_out, vcf_file


# genotype
def main(
    reference_file: str,
    vcf_file: str,
    bam2graphtyper_file: str,
    region_file: str, 
    env_path, 
    threads: int,
    restart: bool
):
    """
    :param reference_file: reference genome
    :param vcf_file: vcf file
    :param bam2graphtyper_file: configure file
    :param region_file: configure file
    :param env_path: environment variable
    :param threads: Threads
    :param restart: Whether to check if the file exists and skip this step
    :return: stdout, stderr, log_out, vcf_out_file
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    cmd = "graphtyper genotype_sv {} {} --sams={} --region_file={} --output=./ " \
          "--threads {}".\
        format(reference_file, vcf_file, bam2graphtyper_file, region_file, threads)

    # Check if the file exists
    if restart:
        # check file
        file_size = getsize(
            "graphtyper.vcf"
        )
        # <= 0
        if file_size <= 0:
            # submit task
            stdout, stderr, log_out = run_cmd.run(cmd, "GraphTyper2.genotype_sv", env_path)
    else:  # If restart is not specified, run directly
        # submit task
        stdout, stderr, log_out = run_cmd.run(cmd, "GraphTyper2.genotype_sv", env_path)

    # Report an error if there is a problem with the exit code
    if log_out:
        return stdout, stderr, log_out, ""

    # Merge vcf files
    stdout, stderr, log_out, vcf_out_file = merge_vcf(
        env_path, 
        restart
    )

    return stdout, stderr, log_out, vcf_out_file
