#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import os
import run_cmd
from getsize import getsize


# genotype
def main(
    reference_file: str,
    vcf_file: str,
    bam2paragraph_file: str,
    env_path, 
    threads: int,
    restart: bool
):
    """
    :param reference_file: reference genome
    :param vcf_file: vcf file
    :param bam2paragraph_file: configure file
    :param env_path: environment variable
    :param threads: Threads
    :param restart: Whether to check if the file exists and skip this step
    :return: stdout, stderr, log_out, vcf_out_file
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # run
    cmd = '''tmpn=`mktemp -u paragraph_XXXXX` && 
        tmpd="/tmp/${tmpn}" && 
        mkdir ${tmpd} && 
        echo ${tmpd} && 
        export TMP=${tmpd} &&'''\
          + " multigrmpy.py -i " + vcf_file \
          + " -m " + bam2paragraph_file \
          + " -r " + reference_file \
          + " --threads " + str(threads) \
          + ''' --scratch-dir ${tmpd}''' \
          + " -M 1000 -o ./ && rm -r ${tmpd}"

    vcf_out_file = os.path.join(os.getcwd(), "genotypes.vcf.gz")

    # Check if the file exists
    if restart:
        # check file
        file_size = getsize(
            vcf_out_file
        )
        # <= 0
        if file_size <= 0:
            # submit task
            stdout, stderr, log_out = run_cmd.run(cmd, "Paragraph.multigrmpy.py", env_path)
    else:  # If restart is not specified, run directly
        # submit task
        stdout, stderr, log_out = run_cmd.run(cmd, "Paragraph.multigrmpy.py", env_path)

    return stdout, stderr, log_out, vcf_out_file
