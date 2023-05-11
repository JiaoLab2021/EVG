#!/usr/bin/python3
# Created on 2022/5/11
# @author: du
# email: dzz0539@163.com

import os
import run_cmd
from getsize import getsize


# bwa index
def index(
        reference_file: str,
        restart: bool
):
    """
    :param reference_file: 参考基因组
    :param restart: 是否检查文件是否存在，并跳过该步骤
    :return: stdout, stderr, log_out
    """

    print(reference_file)

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # bwa构建索引以及比对
    # bam文件路径
    cmd = "bwa index {} 2>log.index.txt".format(reference_file)

    # 检查文件是否存在
    if restart:
        # 如果小于等于0
        if getsize("{}.amb".format(reference_file)) <= 0 or \
                getsize("{}.ann".format(reference_file)) <= 0 or \
                getsize("{}.bwt".format(reference_file)) <= 0 or \
                getsize("{}.fai".format(reference_file)) <= 0 or \
                getsize("{}.pac".format(reference_file)) <= 0 or \
                getsize("{}.sa".format(reference_file)) <= 0:
            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "bwa.index")
    else:  # 如果没有指定restart，直接运行
        # 提交任务
        stdout, stderr, log_out = run_cmd.run(cmd, "bwa.index")

    return stdout, stderr, log_out


# bwa mem
def mem(
        reference_file: str,
        threads: int,
        sample_name: str,
        restart: bool,
        fastq_file1: str,
        fastq_file2: str = ""
):
    """
    :param reference_file: 参考基因组
    :param threads: 线程数
    :param sample_name: 样品名称
    :param restart: 是否检查文件是否存在，并跳过该步骤
    :param fastq_file1: 测序文件1
    :param fastq_file2: 测序文件2
    :return:
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # bwa构建索引以及比对
    # bam文件路径
    bam_file = os.path.abspath(sample_name + ".bam")
    cmd = "bwa mem -R '@RG\\tID:foo\\tSM:{}\\tLB:library1' -t {} {} {} {} | samtools view -b -S | samtools " \
          "sort -@ {} -o {} && " \
          "samtools index {}".\
        format(sample_name, threads, reference_file, fastq_file1,
               fastq_file2, threads, bam_file, bam_file)

    # 检查文件是否存在
    if restart:
        # 检查文件
        file_size = getsize(
            bam_file
        )
        # 如果小于等于0
        if file_size <= 0:
            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "bwa.mem")
    else:  # 如果没有指定restart，直接运行
        # 提交任务
        stdout, stderr, log_out = run_cmd.run(cmd, "bwa.mem")

    return stdout, stderr, log_out, bam_file


# save result for GraphTyper2
def bam2graphtyper(
        bam_infos_map
):
    out_txt = ""
    for key, value in bam_infos_map.items():
        out_txt += value["bam_file"] + "\n"

    # 去掉特殊字符
    out_txt = out_txt.strip()

    with open("bam_for_GraphTyper2.txt", "w") as f:
        f.write(out_txt)

    return os.path.abspath("bam_for_GraphTyper2.txt")


# save result for BayesTyper
def bam2bayestyper(
        bam_infos_map
):
    out_txt = ""

    for key, value in bam_infos_map.items():
        out_txt += key + "\tF\t" + os.path.basename(value["bam_file"]) + "\n"

    # 去掉特殊字符
    out_txt = out_txt.strip()

    with open("bam_for_BayesTyper.tsv", "w") as f:
        f.write(out_txt)

    return os.path.abspath("bam_for_BayesTyper.tsv")


# save result for ParaGraph
def bam2paragraph(
        bam_infos_map
):
    out_txt = "id\tpath\tdepth\tread length\n"

    for key, value in bam_infos_map.items():
        out_txt += key + \
                   "\t" + \
                   value["bam_file"] + \
                   "\t" + \
                   str(value["real_depth"]) + \
                   "\t" + \
                   str(value["read_len"]) + \
                   "\n"

    # 去掉特殊字符
    out_txt = out_txt.strip()

    with open("bam_for_Paragraph.txt", "w") as f:
        f.write(out_txt)

    return os.path.abspath("bam_for_Paragraph.txt")
