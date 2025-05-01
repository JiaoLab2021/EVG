#!/usr/bin/env python3

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
    # change working path
    os.chdir(index_dir)

    # log
    stdout = stderr = log_out = ""

    # file prefix
    temp_dir = os.path.join(index_dir, "temp")
    prefix = os.path.join(index_dir, "out")
    cmd = f"vg construct -t {threads} -r {reference_file} -v {vcf_file} 1>{prefix + '.vg'} && mkdir {temp_dir} && vg index -t {threads} -b {temp_dir} -x {prefix + '.xg'} {prefix + '.vg'} && rm -rf {temp_dir} && vg snarls -t {threads} {prefix + '.xg'} 1>{prefix + '.snarls'}"

    # Check if restart is True and file is empty or non-existent, or restart is not specified
    if (
        restart and (
            getsize("out.snarls") <= 0 or \
            getsize("out.vg") <= 0 or \
            getsize("out.xg") <= 0
        )
    ) or (not restart):
        stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    return stdout, stderr, log_out


# vg map
def graphaligner(
    software_work_path, 
    sample_name: str,
    read: str,
    threads: int,
    index_dir: str,
    env_path, 
    restart: bool
):
    # log
    stdout = stderr = log_out = ""

    # output file path
    vg_file = os.path.join(index_dir, "out.vg")
    gam_file = os.path.join(software_work_path, f"{sample_name}.gam")

    # map
    cmd = f"GraphAligner -t {threads} -g {vg_file} -x vg -f {read} -a {gam_file}"

    # Check if restart is True and file is empty or non-existent, or restart is not specified
    if (restart and getsize(gam_file) <= 0) or (not restart):
        stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    return stdout, stderr, log_out


# call
def vg_call(
    software_work_path, 
    sample_name: str,
    threads: int,
    depth: float,
    index_dir: str,
    env_path, 
    restart: bool
):
    # log
    stdout = stderr = log_out = ""

    # output file path
    xg_file = os.path.join(index_dir, "out.xg")
    snarls_file = os.path.join(index_dir, "out.snarls")
    gam_file = os.path.join(software_work_path, f"{sample_name}.gam")
    pack_file = os.path.join(software_work_path, f"{sample_name}.pack")
    vcf_file = os.path.join(software_work_path, f"{sample_name}.vcf")

    # vg pack
    filter_depth = min(0, int(depth/2))  # Filter by comparison depth greater than 0
    cmd = f"vg pack -t {threads} -Q {filter_depth} -x {xg_file} -g {gam_file} -o {pack_file}"

    # Check if restart is True and file is empty or non-existent, or restart is not specified
    if (restart and getsize(pack_file) <= 0) or (not restart):
        stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    # Check whether the log is normal, and exit early if it is not normal
    if log_out:
        return stdout, stderr, log_out, vcf_file

    # vg call
    cmd = f"vg call -t {threads} -s {sample_name} {xg_file} -k {pack_file} -r {snarls_file} 1>{vcf_file}"

    # Check if restart is True and file is empty or non-existent, or restart is not specified
    if (restart and getsize(vcf_file) <= 0) or (not restart):
        stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    return stdout, stderr, log_out, vcf_file


def main(
    software_work_path, 
    sample_name: str,
    read: str,
    depth: float,
    index_dir: str,
    threads: int,
    env_path, 
    restart: bool
):
    os.chdir(software_work_path)

    stdout, stderr, log_out = graphaligner(
        software_work_path, 
        sample_name,
        read,
        threads,
        index_dir,
        env_path, 
        restart
    )
    # Report an error if there is a problem with the exit code
    if log_out:
        return stdout, stderr, log_out, ""

    stdout, stderr, log_out, vcf_out_file = vg_call(
        software_work_path, 
        sample_name,
        threads,
        depth,
        index_dir,
        env_path, 
        restart
    )

    return stdout, stderr, log_out, vcf_out_file
