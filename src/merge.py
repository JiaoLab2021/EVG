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

        cmd += f"--{software} {file} "
            
    cmd += f"-n {sample_name} -o {merge_vcf_file}"

    # Check if restart is True and file is empty or non-existent, or restart is not specified
    if (restart and getsize(merge_vcf_file) <= 0) or (not restart):
        stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    return stdout, stderr, log_out, os.path.join(work_path, merge_vcf_file)
