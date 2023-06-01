#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import os
import run_cmd
from getsize import getsize


# genotype
def main(
    reference_file: str,
    vcf_file: str,
    bam2bayestyper_file: str,
    bam_infos_map, 
    env_path, 
    threads: int,
    restart: bool
):
    """
    :param reference_file:       reference genome
    :param vcf_file:             vcf file
    :param bam2bayestyper_file:  configure fule
    :param bam_infos_map:        information of BAM
    :param env_path:             env path
    :param threads:              thread
    :param restart:              Whether to check if the file exists and skip this step
    :return:
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # First traverse the bam file and soft-link it to the current working path
    for key, value in bam_infos_map.items():
        bam_file = value["bam_file"]
        # ln
        cmd = "ln -sf {} .".format(bam_file)
        # submit task
        stdout, stderr, log_out = run_cmd.run(cmd, "BayesTyper.ln", env_path)
        # Report an error if there is a problem with the exit code
        if log_out:
            return stdout, stderr, log_out, ""

    # Traversing bam file hashes
    work_dir = os.getcwd()
    try:
        with open(bam2bayestyper_file, "r") as f:
            for information in f.readlines():
                informations_split = information.strip().split()
                prefix = informations_split[0]
                prefix = prefix + ".bam"  # Makefile prefix
                bam_file = informations_split[2]

                # kmc
                cmd = "kmc -k55 -ci1 -t1 -fbam {} {} {}".format(bam_file, prefix, work_dir)

                # Check if the file exists
                if restart:
                    # <= 0
                    if getsize(prefix + ".kmc_pre") <= 0 or getsize(prefix + ".kmc_suf") <= 0:
                        # submit task
                        stdout, stderr, log_out = run_cmd.run(cmd, "BayesTyper.kmc", env_path)
                else:  # If restart is not specified, run directly
                    # submit task
                    stdout, stderr, log_out = run_cmd.run(cmd, "BayesTyper.kmc", env_path)

                # Report an error if there is a problem with the exit code
                if log_out:
                    return stdout, stderr, log_out, ""

                # bayesTyperTools makeBloom
                cmd = "bayesTyperTools makeBloom -k {} -p {}".format(prefix, threads)

                # Check if the file exists
                if restart:
                    # <= 0
                    if getsize(prefix + ".bloomData") <= 0 or getsize(prefix + ".bloomMeta") <= 0:
                        # submit task
                        stdout, stderr, log_out = run_cmd.run(cmd, "BayesTyper.makeBloom", env_path)
                else:  # If restart is not specified, run directly
                    # submit task
                    stdout, stderr, log_out = run_cmd.run(cmd, "BayesTyper.makeBloom", env_path)

                # Report an error if there is a problem with the exit code
                if log_out:
                    return stdout, stderr, log_out, ""
    except FileNotFoundError:
        log_out = "[EVG.{}] FileNotFoundError: [Errno 2] No such file or directory: '{}'.\n".format(
                        "BayesTyper",
                        bam2bayestyper_file
                    )
        return "", "", log_out, ""

    # bayesTyper cluster
    cmd = "bayesTyper cluster -v {} -s {} -g {} -p {}".format(
        vcf_file,
        bam2bayestyper_file,
        reference_file,
        threads
    )

    # Check if the file exists
    if restart:
        # <= 0
        if getsize("bayestyper_cluster_data/intercluster_regions.txt.gz") <= 0 or \
                getsize("bayestyper_cluster_data/multigroup_kmers.bloomData") <= 0 or \
                getsize("bayestyper_cluster_data/multigroup_kmers.bloomMeta") <= 0 or \
                getsize("bayestyper_cluster_data/parameter_kmers.fa.gz") <= 0:
            # submit task
            stdout, stderr, log_out = run_cmd.run(cmd, "BayesTyper.cluster", env_path)
    else:  # If restart is not specified, run directly
        # submit task
        stdout, stderr, log_out = run_cmd.run(cmd, "BayesTyper.cluster", env_path)

    # Report an error if there is a problem with the exit code
    if log_out:
        return stdout, stderr, log_out, ""

    # bayesTyper genotype
    cmd = "bayesTyper genotype -v bayestyper_unit_1/variant_clusters.bin -c bayestyper_" \
          "cluster_data -s {} -g {} -o bayestyper_unit_1/bayestyper -z -p {} --noise-genotyping".\
        format(bam2bayestyper_file, reference_file, threads)

    # Check if the file exists
    if restart:
        # check file
        file_size = getsize(
            "bayestyper_unit_1/bayestyper.vcf.gz"
        )
        # <= 0
        if file_size <= 0:
            # submit task
            stdout, stderr, log_out = run_cmd.run(cmd, "BayesTyper.genotype", env_path)
    else:  # If restart is not specified, run directly
        # submit task
        stdout, stderr, log_out = run_cmd.run(cmd, "BayesTyper.genotype", env_path)

    vcf_out_file = os.path.join(os.getcwd(), "bayestyper_unit_1/bayestyper.vcf.gz")

    return stdout, stderr, log_out, vcf_out_file
