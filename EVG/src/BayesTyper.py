#!/usr/bin/python3
# Created on 2022/5/11
# @author: du
# email: dzz0539@163.com

import os
import datetime
import run_cmd
from getsize import getsize


# genotype
def main(
        reference_file: str,
        vcf_file: str,
        bam2bayestyper_file: str,
        bam_infos_map,
        threads: int,
        restart: bool
):
    """
    :param reference_file: 参考基因组
    :param vcf_file: vcf文件
    :param bam2bayestyper_file: 配置文件
    :param bam_infos_map: bam文件信息
    :param threads: 线程数
    :param restart: 是否检查文件是否存在，并跳过该步骤
    :return:
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # 先遍历bam文件，将其软连接到当前工作路径
    for key, value in bam_infos_map.items():
        bam_file = value["bam_file"]
        # 软连接
        cmd = "ln -sf {} .".format(bam_file)
        # 提交任务
        stdout, stderr, log_out = run_cmd.run(cmd, "BayesTyper.ln")
        # 如果退出代码有问题，则报错
        if log_out:
            return stdout, stderr, log_out, ""

    # 遍历bam文件哈希表
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

                # 检查文件是否存在
                if restart:
                    # 如果小于等于0
                    if getsize(prefix + ".kmc_pre") <= 0 or getsize(prefix + ".kmc_suf") <= 0:
                        # 提交任务
                        stdout, stderr, log_out = run_cmd.run(cmd, "BayesTyper.kmc")
                else:  # 如果没有指定restart，直接运行
                    # 提交任务
                    stdout, stderr, log_out = run_cmd.run(cmd, "BayesTyper.kmc")

                # 如果退出代码有问题，则报错
                if log_out:
                    return stdout, stderr, log_out, ""

                # bayesTyperTools makeBloom
                cmd = "bayesTyperTools makeBloom -k {} -p {}".format(prefix, threads)

                # 检查文件是否存在
                if restart:
                    # 如果小于等于0
                    if getsize(prefix + ".bloomData") <= 0 or getsize(prefix + ".bloomMeta") <= 0:
                        # 提交任务
                        stdout, stderr, log_out = run_cmd.run(cmd, "BayesTyper.makeBloom")
                else:  # 如果没有指定restart，直接运行
                    # 提交任务
                    stdout, stderr, log_out = run_cmd.run(cmd, "BayesTyper.makeBloom")

                # 如果退出代码有问题，则报错
                if log_out:
                    return stdout, stderr, log_out, ""
    except FileNotFoundError:
        log_out = "[" + str(datetime.datetime.now()).split(".")[0] + \
                    "] [EVG.{}] FileNotFoundError: [Errno 2] No such file or directory: '{}'.\n".format(
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

    # 检查文件是否存在
    if restart:
        # 如果小于等于0
        if getsize("bayestyper_cluster_data/intercluster_regions.txt.gz") <= 0 or \
                getsize("bayestyper_cluster_data/multigroup_kmers.bloomData") <= 0 or \
                getsize("bayestyper_cluster_data/multigroup_kmers.bloomMeta") <= 0 or \
                getsize("bayestyper_cluster_data/parameter_kmers.fa.gz") <= 0:
            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "BayesTyper.cluster")
    else:  # 如果没有指定restart，直接运行
        # 提交任务
        stdout, stderr, log_out = run_cmd.run(cmd, "BayesTyper.cluster")

    # 如果退出代码有问题，则报错
    if log_out:
        return stdout, stderr, log_out, ""

    # bayesTyper genotype
    cmd = "bayesTyper genotype -v bayestyper_unit_1/variant_clusters.bin -c bayestyper_" \
          "cluster_data -s {} -g {} -o bayestyper_unit_1/bayestyper -z -p {} --noise-genotyping".\
        format(bam2bayestyper_file, reference_file, threads)

    # 检查文件是否存在
    if restart:
        # 检查文件
        file_size = getsize(
            "bayestyper_unit_1/bayestyper.vcf.gz"
        )
        # 如果小于等于0
        if file_size <= 0:
            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "BayesTyper.genotype")
    else:  # 如果没有指定restart，直接运行
        # 提交任务
        stdout, stderr, log_out = run_cmd.run(cmd, "BayesTyper.genotype")

    vcf_out_file = os.path.join(os.getcwd(), "bayestyper_unit_1/bayestyper.vcf.gz")

    return stdout, stderr, log_out, vcf_out_file
