#!/usr/bin/python3
# coding=gb2312

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
    :param reference_file: 参考基因组
    :param vcf_file: vcf文件
    :param bam2bayestyper_file: 配置文件
    :param bam_infos_map: bam文件信息
    :param env_path: 环境变量
    :param threads: 线程�?
    :param restart: 是否检查文件是否存在，并跳过该步骤
    :return:
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # 先遍历bam文件，将其软连接到当前工作路�?
    for key, value in bam_infos_map.items():
        bam_file = value["bam_file"]
        # 软连�?
        cmd = "ln -sf {} .".format(bam_file)
        # 提交任务
        stdout, stderr, log_out = run_cmd.run(cmd, "BayesTyper.ln", env_path)
        # 如果退出代码有问题，则报错
        if log_out:
            return stdout, stderr, log_out, ""

    # 遍历bam文件哈希�?
    work_dir = os.getcwd()
    try:
        with open(bam2bayestyper_file, "r") as f:
            for information in f.readlines():
                informations_split = information.strip().split()
                prefix = informations_split[0]
                prefix = prefix + ".bam"  # 生成文件的前缀
                bam_file = informations_split[2]

                # kmc
                cmd = "kmc -k55 -ci1 -t1 -fbam {} {} {}".format(bam_file, prefix, work_dir)

                # 检查文件是否存�?
                if restart:
                    # 如果小于等于0
                    if getsize(prefix + ".kmc_pre") <= 0 or getsize(prefix + ".kmc_suf") <= 0:
                        # 提交任务
                        stdout, stderr, log_out = run_cmd.run(cmd, "BayesTyper.kmc", env_path)
                else:  # 如果没有指定restart，直接运�?
                    # 提交任务
                    stdout, stderr, log_out = run_cmd.run(cmd, "BayesTyper.kmc", env_path)

                # 如果退出代码有问题，则报错
                if log_out:
                    return stdout, stderr, log_out, ""

                # bayesTyperTools makeBloom
                cmd = "bayesTyperTools makeBloom -k {} -p {}".format(prefix, threads)

                # 检查文件是否存�?
                if restart:
                    # 如果小于等于0
                    if getsize(prefix + ".bloomData") <= 0 or getsize(prefix + ".bloomMeta") <= 0:
                        # 提交任务
                        stdout, stderr, log_out = run_cmd.run(cmd, "BayesTyper.makeBloom", env_path)
                else:  # 如果没有指定restart，直接运�?
                    # 提交任务
                    stdout, stderr, log_out = run_cmd.run(cmd, "BayesTyper.makeBloom", env_path)

                # 如果退出代码有问题，则报错
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

    # 检查文件是否存�?
    if restart:
        # 如果小于等于0
        if getsize("bayestyper_cluster_data/intercluster_regions.txt.gz") <= 0 or \
                getsize("bayestyper_cluster_data/multigroup_kmers.bloomData") <= 0 or \
                getsize("bayestyper_cluster_data/multigroup_kmers.bloomMeta") <= 0 or \
                getsize("bayestyper_cluster_data/parameter_kmers.fa.gz") <= 0:
            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "BayesTyper.cluster", env_path)
    else:  # 如果没有指定restart，直接运�?
        # 提交任务
        stdout, stderr, log_out = run_cmd.run(cmd, "BayesTyper.cluster", env_path)

    # 如果退出代码有问题，则报错
    if log_out:
        return stdout, stderr, log_out, ""

    # bayesTyper genotype
    cmd = "bayesTyper genotype -v bayestyper_unit_1/variant_clusters.bin -c bayestyper_" \
          "cluster_data -s {} -g {} -o bayestyper_unit_1/bayestyper -z -p {} --noise-genotyping".\
        format(bam2bayestyper_file, reference_file, threads)

    # 检查文件是否存�?
    if restart:
        # 检查文�?
        file_size = getsize(
            "bayestyper_unit_1/bayestyper.vcf.gz"
        )
        # 如果小于等于0
        if file_size <= 0:
            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "BayesTyper.genotype", env_path)
    else:  # 如果没有指定restart，直接运�?
        # 提交任务
        stdout, stderr, log_out = run_cmd.run(cmd, "BayesTyper.genotype", env_path)

    vcf_out_file = os.path.join(os.getcwd(), "bayestyper_unit_1/bayestyper.vcf.gz")

    return stdout, stderr, log_out, vcf_out_file
