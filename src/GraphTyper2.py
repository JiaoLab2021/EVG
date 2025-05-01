#!/usr/bin/env python3

# -*- coding: utf-8 -*-
import os
import run_cmd
from getsize import getsize


# merge vcf
def merge_vcf(
    software_work_path, 
    env_path, 
    restart: bool
):
    # log
    stdout = stderr = log_out = ""

    output_file = os.path.join(software_work_path, "graphtyper.vcf")
    
    cmd = f"for i in $(ls {software_work_path}/*/ | grep 'vcf.gz' | grep -v 'tbi' | head -n1); do zcat {software_work_path}/*/$i | grep '#' > {output_file}; done && zcat {software_work_path}/*/*.vcf.gz >> {output_file}"

    # Check if restart is True and file is empty or non-existent, or restart is not specified
    if (restart and getsize("graphtyper.vcf") <= 0) or (not restart):
        stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    return stdout, stderr, log_out, output_file


# genotype
def main(
    software_work_path, 
    reference_file: str,
    vcf_file: str,
    bam2graphtyper_file: str,
    region_file: str, 
    threads: int,
    env_path, 
    restart: bool
):
    os.chdir(software_work_path)

    # log
    stdout = stderr = log_out = ""

    cmd = f"graphtyper genotype_sv {reference_file} {vcf_file} --sams={bam2graphtyper_file} --region_file={region_file} --output={software_work_path} --threads {threads}"

    # Check if restart is True and file is empty or non-existent, or restart is not specified
    if (restart and getsize("graphtyper.vcf") <= 0) or (not restart):
        stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    # Report an error if there is a problem with the exit code
    if log_out:
        return stdout, stderr, log_out, ""

    # Merge vcf files
    stdout, stderr, log_out, vcf_out_file = merge_vcf(software_work_path, env_path, restart)

    return stdout, stderr, log_out, vcf_out_file
