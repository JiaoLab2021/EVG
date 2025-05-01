#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import os
import run_cmd
from getsize import getsize


# check file
def check_vcf(
    index_dir, 
    vcf_file: str,
    env_path, 
    restart: bool
):
    stdout = stderr = log_out = ""

    # output file path
    vcf_out_file = os.path.join(index_dir, "convert.vcf")

    # awk cmd
    awk_cmd = r"""awk 'BEGIN{FS="\t"; OFS="\t"; start=0; end=0} /^#/ {print; next} $1!=chr {chr=$1; start=0; end=0} {if(start==0) {start=$2; end=$2+length($4)-1} else if($2>end) {start=$2; end=$2+length($4)-1} else{next}} {print}'"""
    if vcf_file.endswith(".gz") or vcf_file.endswith(".GZ"):
        cmd = f"zcat {vcf_file} | {awk_cmd} 1>{vcf_out_file}"
    else:
        cmd = f"cat {vcf_file} | {awk_cmd} 1>{vcf_out_file}"

    # Check if restart is True and file is empty or non-existent, or restart is not specified
    if (restart and getsize(vcf_out_file) <= 0) or (not restart):
        stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    return stdout, stderr, log_out, vcf_out_file


# check file
def check_read(
    software_work_path, 
    sample_name: str,
    read1: str,
    read2: str,
    env_path, 
    restart: bool
):
    # initial log
    stdout = stderr = log_out = ""

    # output file path
    read_out = os.path.join(software_work_path, f"{sample_name}.fq")

    if read1.endswith(".gz") or read1.endswith(".GZ"):
        cmd = f"zcat {read1} 1>{read_out}"
    else:
        cmd = f"cat {read1} 1>{read_out}"
    if read2:  # fastq2
        if read2.endswith(".gz") or read2.endswith(".GZ"):
            cmd += f" && zcat {read2} 1>>{read_out}"
        else:
            cmd += f" && cat {read2} 1>>{read_out}"

    # Check if restart is True and file is empty or non-existent, or restart is not specified
    if (restart and getsize(read_out) <= 0) or (not restart):
        stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    return stdout, stderr, log_out, read_out

# index
def run_index(
    reference_file: str, 
    vcf_file: str,
    index_dir: str,
    env_path, 
    threads: int,
    restart: bool
):
    # change working path
    os.chdir(index_dir)

    # output path
    output_prefix = os.path.join(index_dir, "out")

    stdout, stderr, log_out, vcf_out_file = check_vcf(index_dir, vcf_file, env_path, restart)

    # Report an error if there is a problem with the exit code
    if log_out:
        return stdout, stderr, log_out, ""
    
    cmd = f"PanGenie-index -r {reference_file} -v {vcf_out_file} -t {threads} -o {output_prefix}"

    # submit task
    stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    return stdout, stderr, log_out

# genotype
def main(
    software_work_path, 
    sample_name: str,
    read1: str,
    read2: str,
    index_path: str,
    threads: int,
    env_path, 
    restart: bool
):
    os.chdir(software_work_path)

    # Check if files are compressed and merge
    stdout, stderr, log_out, read = check_read(software_work_path, sample_name, read1, read2, env_path, restart)
    # Report an error if there is a problem with the exit code
    if log_out:
        return stdout, stderr, log_out, ""

    cmd = f"PanGenie -s {sample_name} -i {read} -f {index_path}/out -t {threads} -j {threads} -o {os.path.join(software_work_path, sample_name)}"

    # output vcf path
    out_vcf_file = os.path.join(software_work_path, f"{sample_name}_genotyping.vcf")

    # Check if restart is True and file is empty or non-existent, or restart is not specified
    if (restart and getsize(out_vcf_file) <= 0) or (not restart):
        stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    os.remove(read)

    return stdout, stderr, log_out, out_vcf_file
