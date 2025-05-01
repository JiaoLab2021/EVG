#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import os
import concurrent.futures
import run_cmd
from getsize import getsize


# log
import logging
logger = logging.getLogger('convert')
formatter = logging.Formatter('[%(asctime)s] %(message)s')
handler = logging.StreamHandler()  # output to the console
handler.setFormatter(formatter)
logger.addHandler(handler)


# Convert non-ATGCNatgcn in reference to N
def convert_reference(
    reference_file: str,
    env_path, 
    restart: bool
):
    stdout = stderr = log_out = ""

    # converted filename
    out_reference_file = "convert." + os.path.basename(reference_file)

    # ############################# awk #############################
    # check fasta
    cmd = '''awk '{if ($1~/>/) {print $0} else {$0=toupper($0); gsub(/[^ATGCNatgcn]/,"N"); print $0}}' ''' + reference_file + " 1>" + out_reference_file

    # Check if restart is True and file is empty or non-existent, or restart is not specified
    if (restart and getsize(out_reference_file) <= 0) or (not restart):
        stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    # Check whether the log is normal, and exit early if it is not normal
    if log_out:
        return stdout, stderr, log_out, out_reference_file

    # ############################# samtools faidx #############################
    # build index
    cmd = f"samtools faidx {out_reference_file}"

    # submit task
    stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    return stdout, stderr, log_out, os.path.abspath(out_reference_file)


# Sort vcf files
def bgzip_vcf(
    vcf_file: str,
    env_path, 
    threads, 
    restart: bool
):
    stdout = stderr = log_out = ""
    
    vcf_file_gz = vcf_file + ".gz"

    cmd = f"bgzip -@ {threads} -f {vcf_file} && tabix -f {vcf_file_gz}"

    # Check if restart is True and file is empty or non-existent, or restart is not specified
    if (restart and (getsize(vcf_file_gz) <= 28 or getsize(vcf_file_gz + ".tbi") <= 72)) or (not restart):
        stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    return stdout, stderr, log_out, os.path.abspath(vcf_file_gz)


# Evaluate fastq files
def fastq_count(
    code_dir: str,
    env_path, 
    read1: str,
    read2: str = ""
):
    # code path
    code_path = os.path.join(code_dir, "fastAQ count")

    # Evaluate fastq size
    if read2:  # Next-generation paired-end sequencing
        cmd = f"{code_path} -i {read1} -i {read2}"
    else:  # Third-generation sequencing data or second-generation single-end sequencing
        cmd = f"{code_path} -i {read1}"

    # submit task
    stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    fastq_base = 0
    read_num = 0
    read_len = 0
    for i in stdout.split("\n"):
        if "readBase" in i:
            fastq_base = int(i.strip().split(":")[1])
        if "readNum" in i:
            read_num = int(i.strip().split(":")[1])
        if "readLen" in i:
            read_len = int(i.strip().split(":")[1])

    # Is the result of the inspection correct?
    if fastq_base == 0 or read_num == 0 or read_len == 0:
        log = '[EVG.fastAQ count] The fastq file is wrong, please check the parameters.\n'
        logger.error(log)

    return stdout, stderr, log_out, fastq_base, read_num, read_len


# Estimating Genome Size
def fasta_count(
    code_dir: str,
    env_path, 
    reference_file: str
):
    # code path
    code_path = os.path.join(code_dir, "fastAQ count")

    # Assess fasta size
    cmd = f"{code_path} -i {reference_file}"

    # submit task
    stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    fasta_base = 0
    for i in stdout.split("\n"):
        if "readBase" in i:
            fasta_base = int(i.strip().split(":")[1])

    # Is the result of the inspection correct?
    if fasta_base == 0:
        log_out = f'[EVG.fastAQ count] Error: The size of the reference genome is 0. Please check the file: {reference_file}.\n'

    return stdout, stderr, log_out, fasta_base


# Downsample fastq files
def downsample(
    code_dir: str,
    read1: str,
    read2: str,
    fastq_base: int,
    fasta_base: int,
    need_depth: float,
    env_path, 
    restart: bool
):
    # Initialize log
    stdout = stderr = log_out = ""

    # code path
    code_path = os.path.join(code_dir, "fastAQ sample")

    # downsampling
    need_base = fasta_base * need_depth
    need_ratio = round(need_base / fastq_base, 3)
    read_depth = fastq_base / fasta_base

    if need_ratio >= 1:  # If the sequencing data is less than the set value, it will be skipped and no downsampling will be performed
        log = "Insufficient sequencing data ({:.2f}x/{}x), skip downsampling step: {}, {}".format(read_depth, need_depth, os.path.basename(read1), os.path.basename(read2) if read2 else "")
        logger.error(log)
        read1_out = read1
        read2_out = read2
    else:
        read1_out = f"sample.{need_ratio}.{os.path.basename(read1)}"
        if not read1_out.endswith(('.gz', '.GZ')):
            read1_out += ".gz"
        read1_out = os.path.abspath(read1_out)

        if read2:
            read2_out = f"sample.{need_ratio}.{os.path.basename(read2)}"
            if not read2_out.endswith(('.gz', '.GZ')):
                read2_out += ".gz"
            read2_out = os.path.abspath(read2_out)

            # Check if the file exists
            if restart and getsize(read1_out) > 22 and getsize(read2_out) > 22:
                return "", "", log_out, os.path.abspath(read1_out), os.path.abspath(read2_out)

            cmd1 = f"{code_path} -i {read1} -f {need_ratio} 2>/dev/stderr | gzip > {read1_out}"
            cmd2 = f"{code_path} -i {read2} -f {need_ratio} 2>/dev/stderr | gzip > {read2_out}"
            with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
                to_do = []
                future1 = executor.submit(run_cmd.run, cmd1, env_path)
                future2 = executor.submit(run_cmd.run, cmd2, env_path)
                to_do.append(future1)
                to_do.append(future2)

                for future in concurrent.futures.as_completed(to_do):  # concurrent execution
                    stdout, stderr, log_out_tmp = future.result()  # check return value
                    if log_out_tmp:
                        log_out = log_out_tmp
        else:  # Third-generation sequencing data or second-generation single-end sequencing
            cmd = f"{code_path} -i {read1} -f {need_ratio} 2>/dev/stderr | gzip > {read1_out}"
            read2_out = ""

            # Check if the file exists
            if restart and getsize(read1_out) > 22:
                return stdout, stderr, log_out, os.path.abspath(read1_out), ""

            # submit task
            stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    return stdout, stderr, log_out, os.path.abspath(read1_out), os.path.abspath(read2_out) if read2_out else ""


# vcf convert
def vcf_convert(
    code_dir: str,
    reference_file: str,
    vcf_file: str,
    read_len: int,
    mode: str,
    env_path, 
    restart: bool
):
    # log
    stdout = stderr = log_out = ""

    # code path
    code_path = os.path.join(code_dir, "graphvcf convert")

    # output filename of convert result
    vcf_out_name = f"""convert.{os.path.basename(vcf_file.replace(".gz", ""))}"""

    # path completion
    region_file = os.path.abspath("CHROMOSOME.NAME")

    # ################################### graphvcf convert ###################################
    # graphvcf convert
    cmd = f"{code_path} -r {reference_file} -v {vcf_file} -l {read_len} -o {vcf_out_name}"

    # Check if restart is True and file is empty or non-existent, or restart is not specified
    if (restart and (getsize(vcf_out_name) <= 0 or getsize(region_file) <= 0)) or (not restart):
        stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    # Report an error if there is a problem with the exit code
    if log_out:
        return stdout, stderr, log_out, vcf_out_name, "", ""

    # ################################### sort ###################################
    # vcfsort
    cmd = f"""grep '#' {vcf_out_name} > {vcf_out_name+".sort"} && grep -v '#' {vcf_out_name} | sort -k 1,1 -k 2,2n -t $'\t' | awk -F "\t" '!a[$1,$2]++' 1>>{vcf_out_name+".sort"} && mv {vcf_out_name+".sort"} {vcf_out_name}"""

    # Check if restart is True and file is empty or non-existent, or restart is not specified
    if (restart and getsize(vcf_out_name) > 0) or (not restart):
        stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    # Report an error if there is a problem with the exit code
    if log_out:
        return stdout, stderr, log_out, vcf_out_name, "", region_file

    # vcf is divided by category
    if mode == "fast":
        sv_out_name = "sv." + vcf_out_name
        cmd = f"cat {vcf_out_name}" + ''' | awk -F '\t' 'BEGIN{FS=OFS="\t"} {if($0~/#/) {print $0} else if(length($4)<50 && length($5)<50) {pass} else {print $0}}' > ''' + sv_out_name
        # Check if the file exists
        if restart:
            # If vcf_out_name does not exist but vcf_out_name+gz does, indicating it's a file compressed using bgzip, then use zcat to open the file.
            if getsize(vcf_out_name) <= 0 and getsize(vcf_out_name + ".gz") > 28:
                cmd = "z" + cmd

            if getsize(sv_out_name) <= 0 and (getsize(vcf_out_name) > 0 or getsize(vcf_out_name + ".gz") > 28):
                stdout, stderr, log_out = run_cmd.run(cmd, env_path)
            else:
                sv_out_name = vcf_out_name
        else:  # If restart is not specified, run directly
            stdout, stderr, log_out = run_cmd.run(cmd, env_path)
    else:
        sv_out_name = vcf_out_name

    return stdout, stderr, log_out, os.path.abspath(vcf_out_name), os.path.abspath(sv_out_name), region_file
