#!/usr/bin/python3

# -*- coding: utf-8 -*-

import os
import run_cmd
from getsize import getsize


# index
def vg_index(
    reference_file: str,
    vcf_file: str,
    threads: int,
    index_dir: str,
    env_path, 
    restart: bool
):
    """
    :param reference_file: reference genome
    :param vcf_file: vcf file
    :param threads: Threads
    :param index_dir: Threads
    :param env_path: environment variable
    :param restart:
    :return: Whether to check if the file exists and skip this step
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # change working path
    os.chdir(index_dir)

    # file prefix
    cmd = "vg construct -t {} -r {} -v {} 1>out.vg && " \
          "mkdir temp && vg index -t {} -b temp/ -x out.xg out.vg && rm -rf temp && " \
          "vg snarls -t {} out.xg 1>out.snarls".format(threads, reference_file, vcf_file, threads, threads)

    # Check if the file exists
    if restart:
        # <= 0
        if getsize("out.snarls") <= 0 or \
                getsize("out.vg") <= 0 or \
                getsize("out.xg") <= 0:
            # submit task
            stdout, stderr, log_out = run_cmd.run(cmd, "GraphAligner.vg.index", env_path)
    else:  # If restart is not specified, run directly
        # submit task
        stdout, stderr, log_out = run_cmd.run(cmd, "GraphAligner.vg.index", env_path)

    return stdout, stderr, log_out


# vg map
def graphaligner(
    sample_name: str,
    fastq_file: str,
    threads: int,
    index_dir: str,
    env_path, 
    restart: bool
):
    """
    :param sample_name: sample name
    :param fastq_file: sequencing read
    :param threads: Threads
    :param index_dir: index path
    :param env_path: environment variable
    :param restart: Whether to check if the file exists and skip this step
    :return:
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # output file path
    vg_file = os.path.join(index_dir, "out.vg")
    gam_file = os.path.abspath("{}.gam".format(sample_name))

    # map
    cmd = 'GraphAligner ' \
          '-t {} ' \
          '-g {} ' \
          '-x vg ' \
          '-f {} ' \
          '-a {}'.format(threads, vg_file, fastq_file, gam_file)

    # Check if the file exists
    if restart:
        # check file
        file_size = getsize(
            gam_file
        )
        # <= 0
        if file_size <= 0:
            # submit task
            stdout, stderr, log_out = run_cmd.run(cmd, "GraphAligner.graphaligner", env_path)
    else:  # If restart is not specified, run directly
        # submit task
        stdout, stderr, log_out = run_cmd.run(cmd, "GraphAligner.graphaligner", env_path)

    return stdout, stderr, log_out


# call
def vg_call(
    sample_name: str,
    threads: int,
    depth: float,
    index_dir: str,
    env_path, 
    restart: bool
):
    """
    :param sample_name: sample name
    :param threads: Threads
    :param depth: depth
    :param index_dir: index path
    :param env_path: environment variable
    :param restart: Whether to check if the file exists and skip this step
    :return:
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # output file path
    xg_file = os.path.join(index_dir, "out.xg")
    snarls_file = os.path.join(index_dir, "out.snarls")
    gam_file = os.path.abspath("{}.gam".format(sample_name))
    pack_file = os.path.abspath("{}.pack".format(sample_name))
    vcf_file = os.path.abspath(sample_name + ".vcf")

    # vg pack
    filter_depth = min(0, int(depth/2))  # Filter by comparison depth greater than 0
    cmd = 'vg pack ' \
          '-t {} ' \
          '-Q {} ' \
          '-x {} ' \
          '-g {} ' \
          '-o {}'.format(threads, filter_depth, xg_file, gam_file, pack_file)

    # Check if the file exists
    if restart:
        # check file
        file_size = getsize(
            pack_file
        )
        # <= 0
        if file_size <= 0:
            # submit task
            stdout, stderr, log_out = run_cmd.run(cmd, "GraphAligner.vg.pack", env_path)
    else:  # If restart is not specified, run directly
        # submit task
        stdout, stderr, log_out = run_cmd.run(cmd, "GraphAligner.vg.pack", env_path)

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
          '1>{}'.format(threads, sample_name, xg_file, pack_file, snarls_file, vcf_file)

    # Check if the file exists
    if restart:
        # check file
        file_size = getsize(
            vcf_file
        )
        # <= 0
        if file_size <= 0:
            # submit task
            stdout, stderr, log_out = run_cmd.run(cmd, "GraphAligner.vg.call", env_path)
    else:  # If restart is not specified, run directly
        # submit task
        stdout, stderr, log_out = run_cmd.run(cmd, "GraphAligner.vg.call", env_path)

    return stdout, stderr, log_out, vcf_file


def main(
    sample_name: str,
    fastq_file: str,
    threads: int,
    depth: float,
    index_dir: str,
    env_path, 
    restart: bool
):
    """
    :param sample_name: sample name
    :param fastq_file: sequencing read
    :param threads: Threads
    :param depth: depth
    :param index_dir: index path
    :param env_path: environment variable
    :param restart: Whether to check if the file exists and skip this step
    :return:
    """

    stdout, stderr, log_out = graphaligner(
        sample_name,
        fastq_file,
        threads,
        index_dir,
        env_path, 
        restart
    )
    # Report an error if there is a problem with the exit code
    if log_out:
        return stdout, stderr, log_out, ""

    stdout, stderr, log_out, vcf_out_file = vg_call(
        sample_name,
        threads,
        depth,
        index_dir,
        env_path, 
        restart
    )

    return stdout, stderr, log_out, vcf_out_file
