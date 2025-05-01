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
    # Initialize output variables
    stdout = stderr = log_out = ""

    # change working path
    os.chdir(index_dir)

    # file prefix
    temp_dir = os.path.join(index_dir, "temp")
    prefix = os.path.join(index_dir, "out")
    
    # Prepare the command to run based on software type
    cmd = f"mkdir {temp_dir} && vg autoindex -T {temp_dir} -t {threads} -R XG --workflow {software} -r {reference_file} -v {vcf_file} -p {prefix} && vg snarls -t 10 {prefix + '.xg'} > {prefix + '.snarls'} && rm -rf {temp_dir}"

    # Define the files to check based on software and VG version
    required_files = {
        "map": ["out.gcsa", "out.gcsa.lcp", "out.snarls", "out.xg"],
        "giraffe": ["out.dist", "out.giraffe.gbz", "out.snarls", "out.xg"]
    }

    # Check vg version (vg giraffe)
    vg_version = run_cmd.get_version("vg", env_path)
    # Set the minimizer file and zipcodes file based on the VG version
    minimizer_file = os.path.join(index_dir, "out.shortread.withzip.min" if vg_version >= "1.63.0" else "out.min")
    zipcodes_file = os.path.join(index_dir, "out.shortread.zipcodes" if vg_version >= "1.63.0" else "")
    # Add version-specific checks for Giraffe
    if vg_version >= "1.63.0":
        required_files["giraffe"].extend([minimizer_file, zipcodes_file])
    else:
        required_files["giraffe"].append(minimizer_file)
    
    # Check if the necessary files exist and are non-empty
    # Function to check if all required files are non-empty
    def check_files(software_type):
        return all(getsize(file) > 0 for file in required_files[software_type])

    # Check if restart is True and file is empty or non-existent, or restart is not specified
    # Possible cases:
    # 1. (True and not True) → False: restart is True, but the file exists and is non-empty, so the task will not run.
    # 2. (True and not False) → True: restart is True, and the file is empty or non-existent, so the task will run.
    # 3. (False and not True) → True: restart is False, so the task will run regardless of file existence.
    # 4. (False and not False) → True: restart is False, so the task will run regardless of file existence.
    if (restart and not check_files(software)) or (not restart):
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

    # Prepare the base command
    cmd = f"vg map -t {threads} -g {gcsa_file} -x {xg_file} -f {read1}"

    # Add second read if it exists
    if read2:
        cmd += f" -f {read2}"

    # Redirect output to gam file
    cmd += f" 1>{gam_file}"

    # Check if restart is True and gam_file is empty or non-existent, or restart is not specified
    if (restart and getsize(gam_file) <= 0) or (not restart):
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
    gam_file = os.path.join(software_work_path, f"{sample_name}.gam")
    xg_file = os.path.join(index_dir, "out.xg")
    gbz_file = os.path.join(index_dir, "out.giraffe.gbz")
    dist_file = os.path.join(index_dir, "out.dist")
    minimizer_file = os.path.join(index_dir, "out.min")
    zipcodes_file = ""

    # Check vg version
    vg_version = run_cmd.get_version("vg", env_path)

    # From version 1.63.0 onward, the min file is renamed to out.shortread.withzip.min and a zipcodes file is generated.
    if vg_version >= "1.63.0":
        minimizer_file = os.path.join(index_dir, "out.shortread.withzip.min")
        zipcodes_file = os.path.join(index_dir, "out.shortread.zipcodes")

    # Prepare base command
    cmd = f"vg giraffe -t {threads} -x {xg_file} -Z {gbz_file} -m {minimizer_file} -d {dist_file} -f {read1}"

    # Add second read if it exists
    if read2:
        cmd += f" -f {read2}"

    # Add zipcodes_file if it exists
    if zipcodes_file:
        cmd += f" -z {zipcodes_file}"

    # Redirect output to gam file
    cmd += f" 1>{gam_file}"

    # Check if restart is True and file is empty or non-existent, or restart is not specified
    if (restart and getsize(gam_file) <= 0) or (not restart):
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

    # Check if restart is True and is empty or non-existent, or restart is not specified
    if (restart and getsize(pack_file) <= 0) or (not restart):
        stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    # Check whether the log is normal, and exit early if it is not normal
    if log_out:
        return stdout, stderr, log_out, vcf_file

    # vg call
    cmd = f"vg call -t {threads} -s {sample_name} {xg_file} -k {pack_file} -r {snarls_file} 1>{vcf_file}"

    # Check if restart is True and is empty or non-existent, or restart is not specified
    if (restart and getsize(vcf_file) <= 0) or (not restart):
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

    # Map software to function
    software_functions = {
        "VG-MAP": vg_map,
        "VG-Giraffe": giraffe
    }

    if software in software_functions:
        stdout, stderr, log_out = software_functions[software](
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
        # Invalid software argument
        return "", "", '[vg] software arguments error. (VG-MAP/VG-Giraffe).\n', ""

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
