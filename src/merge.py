#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import os
import run_cmd
from getsize import getsize


# merge result
def main(
    code_dir: str,
    work_path: str,
    line_vcf_file: str,
    sample_name: str,
    genotype_vcf_file_list: list,
    select_software_list: list, 
    bam2bayestyper_samplename_list: list, 
    bam2paragraph_samplename_list: list, 
    number: int, 
    env_path, 
    restart: bool
):
    stdout = stderr = log_out = ""
    
    # back to working directory
    os.chdir(work_path)

    # output file
    merge_vcf_file = os.path.join(work_path, f"{sample_name}.vcf.gz")

    # code path
    code_path = os.path.join(code_dir, "graphvcf merge")

    cmd = f"{code_path} -v {line_vcf_file} "

    for index in range(len(genotype_vcf_file_list)):

        file = genotype_vcf_file_list[index]
        software = select_software_list[index]

        if software == "BayesTyper":
            try:
                fileIndex = bam2bayestyper_samplename_list.index(sample_name)
                fileTmp = file[fileIndex//number]
                cmd += f"--{software} {fileTmp} "
            except IndexError:
                continue
        elif software == "Paragraph":
            try:
                fileIndex = bam2paragraph_samplename_list.index(sample_name)
                fileTmp = file[fileIndex//number]
                cmd += f"--{software} {fileTmp} "
            except IndexError:
                continue
        else:
            cmd += f"--{software} {file} "
            
    cmd += f"-n {sample_name} -o {merge_vcf_file}"

    # Check if the file exists
    if restart:
        file_size = getsize(merge_vcf_file)
        if file_size <= 0:
            stdout, stderr, log_out = run_cmd.run(cmd, env_path)
    else:  # If restart is not specified, run directly
        stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    return stdout, stderr, log_out, os.path.join(work_path, merge_vcf_file)
