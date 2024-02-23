#!/usr/bin/env python3

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
    stdout = stderr = log_out = ""

    # change working path
    os.chdir(index_dir)

    # file prefix
    temp_dir = os.path.join(index_dir, "temp")
    prefix = os.path.join(index_dir, "out")
    cmd = f"mkdir {temp_dir} && vg autoindex -T {temp_dir} -t {threads} -R XG --workflow {software} -r {reference_file} -v {vcf_file} -p {prefix} && vg snarls -t 10 {prefix + '.xg'} > {prefix + '.snarls'} && rm -rf {temp_dir}"

    # Check if the file exists
    if restart:
        # map
        if software == "map":
            if getsize("out.gcsa") <= 0 or \
                    getsize("out.gcsa.lcp") <= 0 or \
                    getsize("out.snarls") <= 0 or \
                    getsize("out.xg") <= 0:
                stdout, stderr, log_out = run_cmd.run(cmd, env_path)
        # giraffe
        else:
            if getsize("out.dist") <= 0 or \
                    getsize("out.giraffe.gbz") <= 0 or \
                    getsize("out.min") <= 0 or \
                    getsize("out.snarls") <= 0 or \
                    getsize("out.xg") <= 0:
                stdout, stderr, log_out = run_cmd.run(cmd, env_path)
    else:  # If restart is not specified, run directly
        stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    return stdout, stderr, log_out


# vg map
def vg_map(
    software_work_path, 
    sample_name: str,
    read1: str,
    read2: str,
    env_path, 
    threads: int,
    index_dir: str,
    restart: bool
):
    stdout = stderr = log_out = ""

    # Index and output file paths
    gcsa_file = os.path.join(index_dir, "out.gcsa")
    xg_file = os.path.join(index_dir, "out.xg")
    gam_file = os.path.join(software_work_path, f"{sample_name}.gam")

    # map
    if read2:
        cmd = f"vg map -t {threads} -g {gcsa_file} -x {xg_file} -f {read1} -f {read2} 1>{gam_file}"
    else:
        cmd = f"vg map -t {threads} -g {gcsa_file} -x {xg_file} -f {read1} 1>{gam_file}"

    # Check if the file exists
    if restart:
        file_size = getsize(gam_file)
        if file_size <= 0:
            stdout, stderr, log_out = run_cmd.run(cmd, env_path)
    else:  # If restart is not specified, run directly
        stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    return stdout, stderr, log_out


# vg map
def giraffe(
    software_work_path, 
    sample_name: str,
    read1: str,
    read2: str,
    env_path,
    threads: float,
    index_dir: str,
    restart: bool
):
    stdout = stderr = log_out = ""

    # Index and output file paths
    xg_file = os.path.join(index_dir, "out.xg")
    gbz_file = os.path.join(index_dir, "out.giraffe.gbz")
    minimizer_file = os.path.join(index_dir, "out.min")
    dist_file = os.path.join(index_dir, "out.dist")
    gam_file = os.path.join(software_work_path, f"{sample_name}.gam")

    # map
    if read2:
        cmd = f"vg giraffe -t {threads} -x {xg_file} -Z {gbz_file} -m {minimizer_file} -d {dist_file} -f {read1} -f {read2} 1>{gam_file}"
    else:
        cmd = f"vg giraffe -t {threads} -x {xg_file} -Z {gbz_file} -m {minimizer_file} -d {dist_file} -f {read1} 1>{gam_file}"

    # Check if the file exists
    if restart:
        file_size = getsize(gam_file)
        if file_size <= 0:
            stdout, stderr, log_out = run_cmd.run(cmd, env_path)
    else:  # If restart is not specified, run directly
        stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    return stdout, stderr, log_out


# call
def vg_call(
    software_work_path, 
    sample_name: str,
    env_path, 
    threads: int,
    depth: float,
    index_dir: str,
    restart: bool
):
    stdout = stderr = log_out = ""

    # output file path
    xg_file = os.path.join(index_dir, "out.xg")
    gam_file = os.path.join(software_work_path, f"{sample_name}.gam")
    pack_file = os.path.join(software_work_path, f"{sample_name}.pack")
    snarls_file = os.path.join(index_dir, "out.snarls")
    vcf_file = os.path.join(software_work_path, f"{sample_name}.vcf")

    # vg pack
    filter_depth = min(3, int(depth / 2))  # Filter by comparison depth greater than 3
    cmd = f"vg pack -t {threads} -Q {filter_depth} -x {xg_file} -g {gam_file} -o {pack_file}"

    # Check if the file exists
    if restart:
        file_size = getsize(pack_file)
        if file_size <= 0:
            stdout, stderr, log_out = run_cmd.run(cmd, env_path)
    else:  # If restart is not specified, run directly
        stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    # Check whether the log is normal, and exit early if it is not normal
    if log_out:
        return stdout, stderr, log_out, vcf_file

    # vg call
    cmd = f"vg call -t {threads} -s {sample_name} {xg_file} -k {pack_file} -r {snarls_file} 1>{vcf_file}"

    # Check if the file exists
    if restart:
        file_size = getsize(vcf_file)
        if file_size <= 0:
            stdout, stderr, log_out = run_cmd.run(cmd, env_path)
    else:  # If restart is not specified, run directly
        stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    return stdout, stderr, log_out, vcf_file


def main(
    software_work_path, 
    sample_name: str,
    read1: str,
    read2: str, 
    depth: float,
    software: str,
    index_dir: str,
    threads: int,
    env_path, 
    restart: bool
):
    stdout = stderr = log_out = ""
    
    os.chdir(software_work_path)

    if software == "VG-MAP":
        stdout, stderr, log_out = vg_map(
            software_work_path, 
            sample_name,
            read1,
            read2, 
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
            software_work_path, 
            sample_name,
            read1,
            read2, 
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
        software_work_path, 
        sample_name,
        env_path, 
        threads,
        depth,
        index_dir,
        restart,
    )

    return stdout, stderr, log_out, vcf_out_file
