#!/usr/bin/python3

# -*- coding: utf-8 -*-

import os
import run_cmd
from getsize import getsize


# autoindex
def vg_autoindex(
    reference_file: str,
    vcf_file: str,
    software: str,
    env_path, 
    threads: int,
    index_dir: str,
    restart: bool
):
    """
    :param reference_file: reference genome
    :param vcf_file: vcf file
    :param software: map/giraffe
    :param env_path: environment variable
    :param threads: Threads
    :param index_dir: Threads
    :param restart: Whether to check if the file exists and skip this step
    :return: stdout, stderr, log_out
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # change working path
    os.chdir(index_dir)

    # file prefix
    cmd = "mkdir temp && vg autoindex -T ./temp/ -t {} -R XG --workflow {} -r {} -v {} -p out && " \
          "vg snarls -t 10 out.xg > out.snarls && rm -rf temp".format(
            threads,
            software,
            reference_file,
            vcf_file)

    # Check if the file exists
    if restart:
        # map
        if software == "map":
            # <= 0
            if getsize("out.gcsa") <= 0 or \
                    getsize("out.gcsa.lcp") <= 0 or \
                    getsize("out.snarls") <= 0 or \
                    getsize("out.xg") <= 0:
                # submit task
                stdout, stderr, log_out = run_cmd.run(cmd, "vg.autoindex", env_path)
        # giraffe
        else:
            # <= 0
            if getsize("out.dist") <= 0 or \
                    getsize("out.giraffe.gbz") <= 0 or \
                    getsize("out.min") <= 0 or \
                    getsize("out.snarls") <= 0 or \
                    getsize("out.xg") <= 0:
                # submit task
                stdout, stderr, log_out = run_cmd.run(cmd, "vg.autoindex", env_path)
    else:  # If restart is not specified, run directly
        # submit task
        stdout, stderr, log_out = run_cmd.run(cmd, "vg.autoindex", env_path)

    return stdout, stderr, log_out


# vg map
def vg_map(
    sample_name: str,
    fastq_file1: str,
    fastq_file2: str,
    env_path, 
    threads: int,
    index_dir: str,
    restart: bool
):
    """
    :param sample_name: sample name
    :param fastq_file1: read1
    :param fastq_file2: read2
    :param env_path: environment variable
    :param threads: Threads
    :param index_dir: index path
    :param restart: Whether to check if the file exists and skip this step
    :return:
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # Index and output file paths
    gcsa_file = os.path.join(index_dir, "out.gcsa")
    xg_file = os.path.join(index_dir, "out.xg")
    gam_file = os.path.abspath("{}.gam".format(sample_name))

    # map
    if fastq_file2:
        cmd = 'vg map ' \
              '-t {} ' \
              '-g {} ' \
              '-x {} ' \
              '-f {} ' \
              '-f {} ' \
              '1>{}'.format(threads, gcsa_file, xg_file, fastq_file1, fastq_file2, gam_file)
    else:
        cmd = 'vg map ' \
              '-t {} ' \
              '-g {} ' \
              '-x {} ' \
              '-f {} ' \
              '1>{}'.format(threads, gcsa_file, xg_file, fastq_file1, gam_file)

    # Check if the file exists
    if restart:
        # check file
        file_size = getsize(
            gam_file
        )
        # <= 0
        if file_size <= 0:
            # submit task
            stdout, stderr, log_out = run_cmd.run(cmd, "vg.map", env_path)
    else:  # If restart is not specified, run directly
        # submit task
        stdout, stderr, log_out = run_cmd.run(cmd, "vg.map", env_path)

    return stdout, stderr, log_out


# vg map
def giraffe(
        sample_name: str,
        fastq_file1: str,
        fastq_file2: str,
        env_path,
        threads: float,
        index_dir: str,
        restart: bool
):
    """
    :param sample_name: sample name
    :param fastq_file1: read1
    :param fastq_file2: read2
    :param env_path: environment variable
    :param threads: Threads
    :param index_dir: index path
    :param restart: Whether to check if the file exists and skip this step
    :return: stdout, stderr, log_out
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # Index and output file paths
    xg_file = os.path.join(index_dir, "out.xg")
    gbz_file = os.path.join(index_dir, "out.giraffe.gbz")
    minimizer_file = os.path.join(index_dir, "out.min")
    dist_file = os.path.join(index_dir, "out.dist")
    gam_file = os.path.abspath("{}.gam".format(sample_name))

    # map
    if fastq_file2:
        cmd = 'vg giraffe ' \
              '-t {} ' \
              '-x {} ' \
              '-Z {} ' \
              '-m {} ' \
              '-d {} ' \
              '-f {} ' \
              '-f {} ' \
              '1>{}'.format(threads,
                            xg_file,
                            gbz_file,
                            minimizer_file,
                            dist_file,
                            fastq_file1,
                            fastq_file2,
                            gam_file)
    else:
        cmd = 'vg giraffe ' \
              '-t {} ' \
              '-x {} ' \
              '-Z {} ' \
              '-m {} ' \
              '-d {} ' \
              '-f {} ' \
              '1>{}'.format(threads,
                            xg_file,
                            gbz_file,
                            minimizer_file,
                            dist_file,
                            fastq_file1,
                            gam_file)

    # Check if the file exists
    if restart:
        # check file
        file_size = getsize(
            gam_file
        )
        # <= 0
        if file_size <= 0:
            # submit task
            stdout, stderr, log_out = run_cmd.run(cmd, "vg.giraffe", env_path)
    else:  # If restart is not specified, run directly
        # submit task
        stdout, stderr, log_out = run_cmd.run(cmd, "vg.giraffe", env_path)

    return stdout, stderr, log_out


# call
def vg_call(
        sample_name: str,
        env_path, 
        threads: int,
        depth: float,
        index_dir: str,
        restart: bool
):
    """
    :param sample_name: sample name
    :param env_path: environment variable
    :param threads: Threads
    :param depth: depth
    :param index_dir: index path
    :param restart: Whether to check if the file exists and skip this step
    :return: stdout, stderr, log_out, vcf_file
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # output file path
    xg_file = os.path.join(index_dir, "out.xg")
    gam_file = os.path.abspath("{}.gam".format(sample_name))
    pack_file = os.path.abspath("{}.pack".format(sample_name))
    snarls_file = os.path.join(index_dir, "out.snarls")
    vcf_file = os.path.abspath(sample_name + ".vcf")

    # vg pack
    filter_depth = min(3, int(depth / 2))  # Filter by comparison depth greater than 3
    cmd = 'vg pack ' \
          '-t {} ' \
          '-Q {} ' \
          '-x {} ' \
          '-g {} ' \
          '-o {}'.format(threads,
                         filter_depth,
                         xg_file,
                         gam_file,
                         pack_file)

    # Check if the file exists
    if restart:
        # check file
        file_size = getsize(
            pack_file
        )
        # <= 0
        if file_size <= 0:
            # submit task
            stdout, stderr, log_out = run_cmd.run(cmd, "vg.pack", env_path)
    else:  # If restart is not specified, run directly
        # submit task
        stdout, stderr, log_out = run_cmd.run(cmd, "vg.pack", env_path)

    # Check whether the log is normal, and exit early if it is not normal
    if log_out:
        return stdout, stderr, log_out, vcf_file

    # vg call
    cmd = 'vg call ' \
          '-t {} ' \
          '-s {} ' \
          '{} ' \
          '-k {} ' \
          '-r {} ' \
          '1>{}'.format(threads,
                        sample_name,
                        xg_file,
                        pack_file,
                        snarls_file,
                        vcf_file)

    # Check if the file exists
    if restart:
        # check file
        file_size = getsize(
            vcf_file
        )
        # <= 0
        if file_size <= 0:
            # submit task
            stdout, stderr, log_out = run_cmd.run(cmd, "vg.call")
    else:  # If restart is not specified, run directly
        # submit task
        stdout, stderr, log_out = run_cmd.run(cmd, "vg.call")

    return stdout, stderr, log_out, vcf_file


def main(
    sample_name: str,
    fastq_file1: str,
    fastq_file2: str, 
    env_path, 
    threads: int,
    depth: float,
    software: str,
    index_dir: str,
    restart: bool
):
    """
    :param sample_name: sample name
    :param fastq_file1: read1
    :param fastq_file2: read2
    :param env_path: environment variable
    :param threads: Threads
    :param depth: depth
    :param software: software
    :param index_dir: index path
    :param restart: Whether to check if the file exists and skip this step
    :return:
    """
    if software == "VG-MAP":
        stdout, stderr, log_out = vg_map(
            sample_name,
            fastq_file1,
            fastq_file2, 
            env_path, 
            threads,
            index_dir,
            restart
        )
        # Report an error if there is a problem with the exit code
        if log_out:
            return stdout, stderr, log_out, ""
    elif software == "VG-Giraffe":
        stdout, stderr, log_out = giraffe(
            sample_name,
            fastq_file1,
            fastq_file2, 
            env_path, 
            threads,
            index_dir,
            restart
        )
        # Report an error if there is a problem with the exit code
        if log_out:
            return stdout, stderr, log_out, ""
    else:
        # print log
        log_out = '[vg] software arguments error. (VG-MAP/VG-Giraffe).\n'
        return "", "", log_out, ""

    stdout, stderr, log_out, vcf_out_file = vg_call(
        sample_name,
        env_path, 
        threads,
        depth,
        index_dir,
        restart,
    )

    return stdout, stderr, log_out, vcf_out_file
