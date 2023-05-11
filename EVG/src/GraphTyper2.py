#!/usr/bin/python3
# Created on 2022/5/11
# @author: du
# email: dzz0539@163.com

import os
import run_cmd
from getsize import getsize


# merge vcf
def merge_vcf(
        restart: bool
):
    """
    :param restart: 是否检查文件是否存在，并跳过该步骤
    :return: stdout, stderr, log_out, vcf_file
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    cmd = "cat */*.vcf.gz > graphtyper.vcf.gz && " \
          "gunzip graphtyper.vcf.gz && zcat */000000001-001000000.vcf.gz | grep '#' > out.vcf && " \
          "grep -v '#' graphtyper.vcf >> out.vcf && " \
          "mv out.vcf graphtyper.vcf"

    # 检查文件是否存在
    if restart:
        # 检查文件
        file_size = getsize(
            "graphtyper.vcf"
        )
        # 如果小于等于0
        if file_size <= 0:
            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "GraphTyper2.genotype_sv")
    else:  # 如果没有指定restart，直接运行
        # 提交任务
        stdout, stderr, log_out = run_cmd.run(cmd, "GraphTyper2.genotype_sv")

    vcf_file = os.path.abspath("graphtyper.vcf")

    return stdout, stderr, log_out, vcf_file


# genotype
def main(
        reference_file: str,
        vcf_file: str,
        bam2graphtyper_file: str,
        region_file: str,
        threads: int,
        restart: bool
):
    """
    :param reference_file: 参考基因组
    :param vcf_file: vcf文件
    :param bam2graphtyper_file: 配置文件
    :param region_file: 配置文件
    :param threads: 线程数
    :param restart: 是否检查文件是否存在，并跳过该步骤
    :return: stdout, stderr, log_out, vcf_out_file
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    cmd = "graphtyper genotype_sv {} {} --sams={} --region_file={} --output=./ " \
          "--threads {}".\
        format(reference_file, vcf_file, bam2graphtyper_file, region_file, threads)

    # 检查文件是否存在
    if restart:
        # 检查文件
        file_size = getsize(
            "graphtyper.vcf"
        )
        # 如果小于等于0
        if file_size <= 0:
            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "GraphTyper2.genotype_sv")
    else:  # 如果没有指定restart，直接运行
        # 提交任务
        stdout, stderr, log_out = run_cmd.run(cmd, "GraphTyper2.genotype_sv")

    # 如果退出代码有问题，则报错
    if log_out:
        return stdout, stderr, log_out, ""

    # 合并vcf文件
    stdout, stderr, log_out, vcf_out_file = merge_vcf(
        restart
    )

    return stdout, stderr, log_out, vcf_out_file
