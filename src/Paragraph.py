#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import os
import run_cmd
from getsize import getsize


# genotype
def main(
    paragraph_sample_path, 
    reference_file: str,
    vcf_file: str,
    bam2paragraph_file: str,
    threads: int,
    env_path, 
    restart: bool
):
    os.chdir(paragraph_sample_path)

    # log
    stdout = stderr = log_out = ""

    # run
    cmd = '''tmpn=`mktemp -u paragraph_XXXXX` && tmpd="/tmp/${tmpn}" && mkdir ${tmpd} && echo ${tmpd} && export TMP=${tmpd} &&'''\
          + f" multigrmpy.py -i {vcf_file} -m {bam2paragraph_file} -r {reference_file} --threads {threads}" + ''' --scratch-dir ${tmpd}''' \
          + f" -M 1000 -o {paragraph_sample_path}" + " && rm -r ${tmpd}"

    vcf_out_file = os.path.join(paragraph_sample_path, "genotypes.vcf.gz")

    # Check if the file exists
    if restart:
        file_size = getsize(vcf_out_file)
        if file_size <= 0:
            stdout, stderr, log_out = run_cmd.run(cmd, env_path)
    else:  # If restart is not specified, run directly
        stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    return stdout, stderr, log_out, vcf_out_file
