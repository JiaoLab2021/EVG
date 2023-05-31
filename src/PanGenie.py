#!/usr/bin/python3
# coding=gb2312

import os
import run_cmd
from getsize import getsize


# check file
def check(
        vcf_file: str,
        sample_name: str,
        fastq_file1: str,
        fastq_file2: str,
        env_path, 
        restart: bool
):
    """
    :param vcf_file: vcf文件路径
    :param sample_name 样品名
    :param fastq_file1: 测序文件1
    :param fastq_file2: 测序文件2
    :param env_path: 环境变量
    :param restart: 是否检查文件是否存在，并跳过该步骤
    :return: stdout, stderr, log_out, os.path.abspath(vcf_out_file), os.path.abspath(fastq_out_file)
    """

    #  初始log
    stdout = ""
    stderr = ""
    log_out = ""

    # 输出文件路径
    vcf_out_file = sample_name + ".vcf"
    fastq_out_file = "{}.fq".format(sample_name)

    if '.gz' in vcf_file or ".GZ" in vcf_file:  # 如果是压缩文件，解压
        cmd = "gunzip -c {} 1>{}".format(vcf_file, vcf_out_file)

        # 检查文件是否存在
        if restart:
            # 检查文件
            file_size = getsize(
                vcf_out_file
            )
            # 如果小于等于0
            if file_size <= 0:
                # 提交任务
                stdout, stderr, log_out = run_cmd.run(cmd, "PanGenie.gunzip", env_path)
        else:  # 如果没有指定restart，直接运行
            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "PanGenie.gunzip", env_path)

        # 如果退出代码有问题，则报错
        if log_out:
            return stdout, stderr, log_out, "", "", ""
    else:  # 否则赋值一样的路径
        vcf_out_file = vcf_file

    cmd = ""
    if '.gz' in fastq_file1 or ".GZ" in fastq_file1:  # fastq1
        cmd = "gunzip -c {} 1>{}".format(fastq_file1, fastq_out_file)
        if fastq_file2:  # fastq2
            if '.gz' in fastq_file2 or ".GZ" in fastq_file2:
                cmd += " && gunzip -c {} 1>>{}".format(fastq_file2, fastq_out_file)
            else:
                cmd += " && cat {} 1>>{}".format(fastq_file2, fastq_out_file)
    else:  # 不是压缩文件
        cmd = "cat {} 1>{}".format(fastq_file1, fastq_out_file)
        if fastq_file2:  # fastq2
            if '.gz' in fastq_file2 or ".GZ" in fastq_file2:
                cmd += " && gunzip -c {} 1>>{}".format(fastq_file2, fastq_out_file)
            else:
                cmd += " && cat {} 1>>{}".format(fastq_file2, fastq_out_file)
    # 检查文件是否存在
    if restart:
        # 检查文件
        file_size = getsize(
            fastq_out_file
        )
        # 如果小于等于0
        if file_size <= 0:
            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "PanGenie.gunzip", env_path)
    else:  # 如果没有指定restart，直接运行
        # 提交任务
        stdout, stderr, log_out = run_cmd.run(cmd, "PanGenie.gunzip", env_path)

    # 如果退出代码有问题，则报错
    if log_out:
        return stdout, stderr, log_out, "", "", ""

    return stdout, stderr, log_out, os.path.abspath(vcf_out_file), os.path.abspath(fastq_out_file)


# merge reads
def merge(
    sample_name: str,
    fastq_file1: str,
    fastq_file2: str,
    env_path, 
    restart: bool
):
    """
    :param fastq_file1: 测序文件1
    :param sample_name: 样品名
    :param fastq_file2: 测序文件2
    :param env_path: 环境变量
    :param restart: 是否检查文件是否存在，并跳过该步骤
    :return:
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # 输出文件路径
    fastq_out_file = "{}.fq".format(sample_name)

    if fastq_file2:
        cmd = "cat {} {} {}".format(fastq_file1, fastq_file2, fastq_out_file)

        # 检查文件是否存在
        if restart:
            # 检查文件
            file_size = getsize(
                fastq_out_file
            )
            # 如果小于等于0
            if file_size <= 0:
                # 提交任务
                stdout, stderr, log_out = run_cmd.run(cmd, "PanGenie.gunzip", env_path)
        else:  # 如果没有指定restart，直接运行
            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "PanGenie.gunzip", env_path)

        return stdout, stderr, log_out, os.path.abspath(fastq_out_file)
    else:
        return stdout, stderr, log_out, os.path.abspath(fastq_file1)


# genotype
def main(
    reference_file: str,
    vcf_file: str,
    sample_name: str,
    fastq_file1: str,
    fastq_file2: str,
    env_path, 
    threads: int,
    restart: bool
):
    """
    :param reference_file: 参考基因组
    :param vcf_file: vcf文件
    :param sample_name 样品名
    :param fastq_file1: 测序文件1
    :param fastq_file2: 测序文件2
    :param env_path: 环境变量
    :param threads: 线程数
    :param restart: 是否检查文件是否存在，并跳过该步骤
    :return:
    """

    # 检查文件是否为压缩文件并合并
    stdout, stderr, log_out, vcf_file, fastq_file = check(
        vcf_file,
        sample_name,
        fastq_file1,
        fastq_file2,
        env_path, 
        restart
    )
    # 如果退出代码有问题，则报错
    if log_out:
        return stdout, stderr, log_out, ""

    cmd = "PanGenie -s {} -i {} -r {} -v {} -t {} -j {} -o {}".\
        format(sample_name, fastq_file, reference_file, vcf_file, threads, threads, sample_name)

    # 输出的vcf路径
    out_vcf_file = os.path.abspath("{}_genotyping.vcf").format(sample_name)

    # 检查文件是否存在
    if restart:
        # 检查文件
        file_size = getsize(
            out_vcf_file
        )
        # 如果小于等于0
        if file_size <= 0:
            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "PanGenie", env_path)
    else:  # 如果没有指定restart，直接运行
        # 提交任务
        stdout, stderr, log_out = run_cmd.run(cmd, "PanGenie", env_path)

    return stdout, stderr, log_out, os.path.abspath(out_vcf_file)
