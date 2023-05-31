#!/usr/bin/python3

# -*- coding: utf-8 -*-

import os
import concurrent.futures
import run_cmd
from getsize import getsize


# log
import logging
logger = logging.getLogger('SynDiv')
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
    """
    :param reference_file:    original reference genome
    :param env_path:          environment variable
    :param restart:           Whether to check if the file exists and skip this step
    :return: stdout, stderr, log_out, out_reference_file
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # converted filename
    out_reference_file = "convert." + os.path.basename(reference_file)

    # ############################# awk #############################
    # check fasta
    cmd = '''awk '{if ($1~/>/) {print $0} else {$0=toupper($0); gsub(/[^ATGCNatgcn]/,"N"); print $0}}' ''' + \
          reference_file + \
          " 1>" + \
          out_reference_file

    # Check if the file exists
    if restart:
        # check file
        file_size = getsize(
            out_reference_file
        )
        # <= 0
        if file_size <= 0:
            # submit task
            stdout, stderr, log_out = run_cmd.run(cmd, "convert_reference.awk", env_path)
    else:  # If restart is not specified, run directly
        # submit task
        stdout, stderr, log_out = run_cmd.run(cmd, "convert_reference.awk", env_path)

    # Check whether the log is normal, and exit early if it is not normal
    if log_out:
        return stdout, stderr, log_out, out_reference_file

    # ############################# samtools faidx #############################
    # build index
    cmd = "samtools faidx " + out_reference_file

    # submit task
    stdout, stderr, log_out = run_cmd.run(cmd, "convert_reference.faidx", env_path)

    # Check whether the log is normal, and exit early if it is not normal
    if log_out:
        return stdout, stderr, log_out, out_reference_file

    # Complete the path
    out_reference_file = os.path.abspath(out_reference_file)

    return stdout, stderr, log_out, out_reference_file


# Sort vcf files
def bgzip_vcf(
    vcf_file: str,
    env_path, 
    restart: bool
):
    """
    :param vcf_file:   vcf file
    :param env_path:   environment variable
    :param restart:    Whether to check if the file exists and skip this step
    :return: stdout, stderr, log_out, out_vcf_file
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    out_vcf_file = vcf_file + ".gz"

    cmd = "bgzip -f {} && tabix -f {}".format(vcf_file, out_vcf_file)

    # Check if the file exists
    if restart:
        # <= 0
        if getsize(out_vcf_file) <= 28 and getsize(out_vcf_file + ".tbi") <= 72:
            # submit task
            stdout, stderr, log_out = run_cmd.run(cmd, "bgzip", env_path)
    else:  # If restart is not specified, run directly
        # submit task
        stdout, stderr, log_out = run_cmd.run(cmd, "bgzip", env_path)

    # Complete the path
    out_vcf_file = os.path.abspath(out_vcf_file)

    return stdout, stderr, log_out, out_vcf_file


# Evaluate fastq files
def fastq_count(
    code_dir: str,
    env_path, 
    fastq_file1: str,
    fastq_file2: str = ""
):
    """
    :param code_dir:    code path
    :param env_path:    environment variable
    :param fastq_file1: read1
    :param fastq_file2: read2
    :return: stdout, stderr, log_out, fastq_base, read_num, read_len
    """
    # code path
    code_path = os.path.join(code_dir, "fastAQ count")

    # Evaluate fastq size
    if fastq_file2:  # Next-generation paired-end sequencing
        cmd = '{} -i {} -i {}'.format(code_path, fastq_file1, fastq_file2)
    else:  # Third-generation sequencing data or second-generation single-end sequencing
        cmd = '{} -i {}'.format(code_path, fastq_file1)

    # submit task
    stdout, stderr, log_out = run_cmd.run(cmd, "fastAQ count", env_path)

    fastq_base = 0
    read_num = 0
    read_len = 0
    for i in stdout.decode().split("\n"):
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
        raise Exception(log)

    return stdout, stderr, log_out, fastq_base, read_num, read_len


# Estimating Genome Size
def fasta_count(
    code_dir: str,
    env_path, 
    reference_file: str
):
    """
    :param code_dir:       code path
    :param env_path:       environment variable
    :param reference_file: reference genome
    :return: stdout, stderr, log_out, fasta_base
    """
    # code path
    code_path = os.path.join(code_dir, "fastAQ count")

    # Assess fasta size
    cmd = '{} -i {}'.format(code_path, reference_file)

    # submit task
    stdout, stderr, log_out = run_cmd.run(cmd, "fasta_count", env_path)

    fasta_base = 0
    for i in stdout.decode().split("\n"):
        if "readBase" in i:
            fasta_base = int(i.strip().split(":")[1])

    # Is the result of the inspection correct?
    if fasta_base == 0:
        log_out = '[EVG.fasta_count] The fasta file is wrong, please check the parameters.\n'

    return stdout, stderr, log_out, fasta_base


# Downsample fastq files
def downsample(
        code_dir: str,
        fastq_file1: str,
        fastq_file2: str,
        fastq_base: int,
        fasta_base: int,
        need_depth: float,
        env_path, 
        restart: bool
):
    """
    :param code_dir:       code path
    :param fastq_file1:    read1
    :param fastq_file2:    read2
    :param fastq_base:     the base of fastq
    :param fasta_base:     the base of fasta
    :param need_depth:     need depth
    :param env_path:       environment variable
    :param restart: Whether to check if the file exists and skip this step
    :return:
    """
    # Initialize log
    stdout = ""
    stderr = ""
    log_out = ""

    # code path
    code_path = os.path.join(code_dir, "fastAQ sample")

    # downsampling
    need_base = fasta_base * need_depth
    need_ratio = round(need_base / fastq_base, 3)
    read_depth = fastq_base / fasta_base

    if need_ratio >= 1:  # If the sequencing data is less than the set value, it will be skipped and no downsampling will be performed
        log = '[EVG.fastAQ sample] Insufficient sequencing data ({:.2f}×/{}×), skip downsampling step.\n'. \
                  format(read_depth, need_depth)
        logger.error(log)
        fastq_out_file1 = fastq_file1
        fastq_out_file2 = fastq_file2
    else:
        fastq_out_file1 = "sample." + str(need_ratio) + "." + os.path.basename(fastq_file1)
        fastq_out_file1 = fastq_out_file1.replace(".gz", "").replace(".GZ", "")
        fastq_out_file1 = os.path.abspath(fastq_out_file1)  # path completion

        if fastq_file2:  # Next-generation paired-end sequencing
            # Output the filename and complete the path
            fastq_out_file2 = "sample." + str(need_ratio) + "." + os.path.basename(fastq_file2)
            fastq_out_file2 = fastq_out_file2.replace(".gz", "").replace(".GZ", "")
            fastq_out_file2 = os.path.abspath(fastq_out_file2)  # path completion

            # Check if the file exists
            if restart:
                # check file
                file_size1 = getsize(
                    fastq_out_file1
                )
                file_size2 = getsize(
                    fastq_out_file2
                )
                # Both files exist, skip this step
                if file_size1 > 0 and file_size2 > 0:
                    return stdout, stderr, log_out, os.path.abspath(fastq_out_file1), os.path.abspath(fastq_out_file2)

            cmd1 = '{} -i {} -f {} 1>{}'.format(code_path, fastq_file1, need_ratio, fastq_out_file1)
            cmd2 = '{} -i {} -f {} 1>{}'.format(code_path, fastq_file2, need_ratio, fastq_out_file2)
            with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
                to_do = []
                future1 = executor.submit(run_cmd.run, cmd1, "fastAQ sample", env_path)
                future2 = executor.submit(run_cmd.run, cmd2, "fastAQ sample", env_path)
                to_do.append(future1)
                to_do.append(future2)

                for future in concurrent.futures.as_completed(to_do):  # concurrent execution
                    stdout, stderr, log_out_tmp = future.result()  # check return value
                    if log_out_tmp:
                        log_out = log_out_tmp
        else:  # Third-generation sequencing data or second-generation single-end sequencing
            cmd = '{} -i {} -f {} 1>{}'.format(code_path, fastq_file1, need_ratio, fastq_out_file1)
            fastq_out_file2 = ""

            # Check if the file exists
            if restart:
                # check file
                file_size1 = getsize(
                    fastq_out_file1
                )
                # file exists, skip this step
                if file_size1 > 0:
                    return stdout, stderr, log_out, os.path.abspath(fastq_out_file1), ""

            # submit task
            stdout, stderr, log_out = run_cmd.run(cmd, "fastAQ sample", env_path)

    return stdout, stderr, log_out, os.path.abspath(fastq_out_file1), os.path.abspath(fastq_out_file2)


# vcf convert
def vcf_convert(
    code_dir: str,
    reference_file: str,
    vcf_file: str,
    read_len: int,
    out_name: str,
    mode: str,
    env_path, 
    restart: bool
):
    """
    :param code_dir: code path
    :param reference_file: reference genome
    :param vcf_file: vcf file
    :param read_len: read length
    :param out_name: the vcf filename of output
    :param mode:     mode
    :param env_path: environment variable
    :param restart: Whether to check if the file exists and skip this step
    :return: stdout, stderr, log_out, out_name, sv_out_name, region_file
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # code path
    code_path = os.path.join(code_dir, "graphvcf convert")

    # ################################### graphvcf convert ###################################
    # graphvcf convert
    cmd = '{} -r {} -v {} -l {} -o {}'.format(
        code_path,
        reference_file,
        vcf_file,
        read_len,
        out_name)

    # Check if the file exists
    if restart:
        # <= 0
        if getsize(out_name) <= 0 and getsize(out_name + ".gz") <= 0:  # If the file compressed by vcf or bgzip exists, it will be skipped
            # submit task
            stdout, stderr, log_out = run_cmd.run(cmd, "graphvcf convert", env_path)
    else:  # If restart is not specified, run directly
        # submit task
        stdout, stderr, log_out = run_cmd.run(cmd, "graphvcf convert", env_path)

    # Report an error if there is a problem with the exit code
    if log_out:
        return stdout, stderr, log_out, out_name, "", ""

    # path completion
    region_file = os.path.abspath("CHROMOSOME.NAME")

    # ################################### sort ###################################
    # vcfsort
    cmd = '''grep '#' {} > {} && grep -v '#' {} | sort --parallel=10 -k 1,1 -k 2,2n -t $'\t' | '''.format(
        out_name,
        out_name + ".sort",
        out_name
    )
    cmd += '''awk -F "\t" '!a[$1,$2]++' 1>>{} && mv {} {}'''.format(
        out_name+".sort",
        out_name+".sort",
        out_name
    )

    # Check if the file exists
    if restart:
        # <= 0
        if getsize(out_name) > 0:  # If the file compressed by vcf or bgzip exists, it will be skipped
            # submit task
            stdout, stderr, log_out = run_cmd.run(cmd, "sort", env_path)
    else:  # If restart is not specified, run directly
        # submit task
        stdout, stderr, log_out = run_cmd.run(cmd, "sort", env_path)

    # Report an error if there is a problem with the exit code
    if log_out:
        return stdout, stderr, log_out, out_name, "", region_file

    # vcf is divided by category
    if mode == "fast":
        sv_out_name = "sv." + out_name
        cmd = ""
        if getsize(out_name) > 0:
            cmd = '''cat ''' + out_name + ''' | awk -F '\t' 'BEGIN{FS=OFS="\t"} {if($0~/#/) 
            {print $0} else if(length($4)<50 && length($5)<50) {pass} else {print $0}}' > ''' + sv_out_name
        elif getsize(out_name + ".gz") > 28:
            cmd = '''zcat ''' + out_name + ".gz" + ''' | awk -F '\t' 'BEGIN{FS=OFS="\t"} {if($0~/#/) 
                {print $0} else if(length($4)<50 && length($5)<50) {pass} else {print $0}}' > ''' + sv_out_name
        # Check if the file exists
        if restart:
            # check file
            file_size = getsize(
                sv_out_name
            )
            # <= 0
            if file_size <= 0:
                # submit task
                stdout, stderr, log_out = run_cmd.run(cmd, "split", env_path)
        else:  # If restart is not specified, run directly
            # submit task
            stdout, stderr, log_out = run_cmd.run(cmd, "split", env_path)
    else:
        sv_out_name = out_name

    out_name = os.path.abspath(out_name)
    sv_out_name = os.path.abspath(sv_out_name)

    return stdout, stderr, log_out, out_name, sv_out_name, region_file
