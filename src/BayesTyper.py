#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import os
import run_cmd
from getsize import getsize


# genotype
def main(
    bayestyper_sample_path, 
    reference_file: str,
    vcf_file: str,
    bam2bayestyper_file: str,
    bam_infos_map, 
    threads: int,
    env_path, 
    restart: bool
):
    os.chdir(bayestyper_sample_path)

    # log
    stdout = stderr = log_out = ""

    # ----------------------------------- bayesTyper k-mers ----------------------------------- #
    try:
        with open(bam2bayestyper_file, "r") as f:
            for information in f.readlines():
                informations_split = information.strip().split()
                sample_name = informations_split[0]  # Sample name
                prefix = informations_split[2]  # Makefile prefix
                bam_file = informations_split[2]

                # ----------------------------------- ln -sf ----------------------------------- #
                cmd = f"ln -sf {bam_infos_map[sample_name]['bam']} {bam_file}"
                # submit task
                stdout, stderr, log_out = run_cmd.run(cmd, env_path)
                # Report an error if there is a problem with the exit code
                if log_out:
                    return stdout, stderr, log_out, ""

                # kmc
                cmd = f"kmc -k55 -ci1 -t1 -fbam {bam_file} {prefix} {bayestyper_sample_path}"

                # Check if the file exists
                if restart:
                    if getsize(prefix + ".kmc_pre") <= 0 or getsize(prefix + ".kmc_suf") <= 0:
                        stdout, stderr, log_out = run_cmd.run(cmd, env_path)
                else:  # If restart is not specified, run directly
                    stdout, stderr, log_out = run_cmd.run(cmd, env_path)

                # Report an error if there is a problem with the exit code
                if log_out:
                    return stdout, stderr, log_out, ""

                # bayesTyperTools makeBloom
                cmd = f"bayesTyperTools makeBloom -k {prefix} -p {threads}"

                # Check if the file exists
                if restart:
                    if getsize(prefix + ".bloomData") <= 0 or getsize(prefix + ".bloomMeta") <= 0:
                        stdout, stderr, log_out = run_cmd.run(cmd, env_path)
                else:  # If restart is not specified, run directly
                    stdout, stderr, log_out = run_cmd.run(cmd, env_path)

                # Report an error if there is a problem with the exit code
                if log_out:
                    return stdout, stderr, log_out, ""
    except FileNotFoundError:
        log_out = "FileNotFoundError: [Errno 2] No such file or directory: '{bam2bayestyper_file}'.\n"
        return "", "", log_out, ""

    # ----------------------------------- bayesTyper cluster ----------------------------------- #
    # bayesTyper cluster
    cmd = f"bayesTyper cluster -v {vcf_file} -s {bam2bayestyper_file} -g {reference_file} -p {threads} -o {os.path.join(bayestyper_sample_path, 'bayestyper')}"

    # Check if the file exists
    if restart:
        if getsize("bayestyper_cluster_data/intercluster_regions.txt.gz") <= 0 or \
                getsize("bayestyper_cluster_data/multigroup_kmers.bloomData") <= 0 or \
                getsize("bayestyper_cluster_data/multigroup_kmers.bloomMeta") <= 0 or \
                getsize("bayestyper_cluster_data/parameter_kmers.fa.gz") <= 0:
            stdout, stderr, log_out = run_cmd.run(cmd, env_path)
    else:  # If restart is not specified, run directly
        stdout, stderr, log_out = run_cmd.run(cmd, env_path)

    # Report an error if there is a problem with the exit code
    if log_out:
        return stdout, stderr, log_out, ""

    # ----------------------------------- bayesTyper genotype ----------------------------------- #
    # List all directories under the current directory starting with "bayestyper_unit_"
    unit_folders = [os.path.join(bayestyper_sample_path, folder) for folder in os.listdir(bayestyper_sample_path) if folder.startswith("bayestyper_unit_")]

    # If the list is empty, return null
    if len(unit_folders) == 0:
        log_out = "Warning: BayesTyper does not generate any folders starting with 'bayestyper_unit_'.\n"
        return "", "", log_out, ""

    # Iterate through all the unit folders and run the genotyping command
    for unit_folder in unit_folders:
        # Define the command to run BayesTyper genotype
        cmd = f"bayesTyper genotype -v {os.path.join(unit_folder, 'variant_clusters.bin')} -c {os.path.join(bayestyper_sample_path, 'bayestyper_cluster_data')} -s {bam2bayestyper_file} -g {reference_file} -o {os.path.join(unit_folder, 'bayestyper')} -z -p {threads} --noise-genotyping"

        # Check if the file exists
        if restart:
            vcf_gz_file = os.path.join(unit_folder, "bayestyper.vcf.gz")
            file_size = getsize(vcf_gz_file)
            if file_size <= 0:
                stdout, stderr, log_out = run_cmd.run(cmd, env_path)
        else:  # If restart is not specified, run the command directly
            stdout, stderr, log_out = run_cmd.run(cmd, env_path)

        # Report an error if there is a problem with the exit code
        if log_out:
            return stdout, stderr, log_out, ""    

    # ----------------------------------- merge ----------------------------------- #
    # Merge the output VCF files using a shell script
    vcf_files = [os.path.join(unit_folder, "bayestyper.vcf.gz") for unit_folder in unit_folders]
    vcf_out_file = os.path.join(bayestyper_sample_path, "bayestyper.vcf.gz")

    # head
    head_cmd = f"zcat {vcf_files[0]} | grep '^#' | gzip -c > {vcf_out_file}"
    stdout, stderr, log_out = run_cmd.run(head_cmd, env_path)
    # Report an error if there is a problem with the exit code
    if log_out:
        return stdout, stderr, log_out, ""

    # merge
    merge_cmd = f"zcat {' '.join(vcf_files)} | grep -v '^#' | sort -k 1,1 -k 2,2n -t $'\t' | gzip -c >> {vcf_out_file}"
    stdout, stderr, log_out = run_cmd.run(merge_cmd, env_path)
    # Report an error if there is a problem with the exit code
    if log_out:
        return stdout, stderr, log_out, ""

    return stdout, stderr, log_out, vcf_out_file
