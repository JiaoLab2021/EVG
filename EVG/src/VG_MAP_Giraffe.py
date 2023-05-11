#!/usr/bin/python3
# Created on 2022/5/11
# @author: du
# email: dzz0539@163.com

import os
import datetime
import sys
import run_cmd
from getsize import getsize


# autoindex
def vg_autoindex(
        reference_file: str,
        vcf_file: str,
        software: str,
        threads: int,
        index_dir: str,
        restart: bool
):
    """
    :param reference_file: 参考基因组
    :param vcf_file: vcf文件
    :param software: map/giraffe
    :param threads: 线程数
    :param index_dir: 索引路径
    :param restart: 是否检查文件是否存在，并跳过该步骤
    :return: stdout, stderr, log_out
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # 更改工作路径
    os.chdir(index_dir)

    # 文件前缀
    cmd = "mkdir temp && vg autoindex -T ./temp/ -t {} -R XG --workflow {} -r {} -v {} -p out && " \
          "vg snarls -t 10 out.xg > out.snarls && rm -rf temp".format(
            threads,
            software,
            reference_file,
            vcf_file)

    # 检查文件是否存在
    if restart:
        # map
        if software == "map":
            # 如果小于等于0
            if getsize("out.gcsa") <= 0 or \
                    getsize("out.gcsa.lcp") <= 0 or \
                    getsize("out.snarls") <= 0 or \
                    getsize("out.xg") <= 0:
                # 提交任务
                stdout, stderr, log_out = run_cmd.run(cmd, "vg.autoindex")
        # giraffe
        else:
            # 如果小于等于0
            if getsize("out.dist") <= 0 or \
                    getsize("out.giraffe.gbz") <= 0 or \
                    getsize("out.min") <= 0 or \
                    getsize("out.snarls") <= 0 or \
                    getsize("out.xg") <= 0:
                # 提交任务
                stdout, stderr, log_out = run_cmd.run(cmd, "vg.autoindex")
    else:  # 如果没有指定restart，直接运行
        # 提交任务
        stdout, stderr, log_out = run_cmd.run(cmd, "vg.autoindex")

    return stdout, stderr, log_out


# vg map
def vg_map(
        sample_name: str,
        fastq_file1: str,
        fastq_file2: str,
        threads: int,
        index_dir: str,
        restart: bool
):
    """
    :param sample_name: 样品名
    :param fastq_file1: 测序文件1
    :param fastq_file2: 测序文件2
    :param threads: 线程数
    :param index_dir: index的路径
    :param restart: 是否检查文件是否存在，并跳过该步骤
    :return:
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # 索引和输出文件路径
    gcsa_file = os.path.join(index_dir, "out.gcsa")
    xg_file = os.path.join(index_dir, "out.xg")
    gam_file = os.path.abspath("{}.gam".format(sample_name))

    # map
    if fastq_file2:
        cmd = 'vg map ' \
              '-t {} ' \
              '-g {} ' \
              '-x {} ' \
              '-f {} ' \
              '-f {} ' \
              '1>{}'.format(threads, gcsa_file, xg_file, fastq_file1, fastq_file2, gam_file)
    else:
        cmd = 'vg map ' \
              '-t {} ' \
              '-g {} ' \
              '-x {} ' \
              '-f {} ' \
              '1>{}'.format(threads, gcsa_file, xg_file, fastq_file1, gam_file)

    # 检查文件是否存在
    if restart:
        # 检查文件
        file_size = getsize(
            gam_file
        )
        # 如果小于等于0
        if file_size <= 0:
            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "vg.map")
    else:  # 如果没有指定restart，直接运行
        # 提交任务
        stdout, stderr, log_out = run_cmd.run(cmd, "vg.map")

    return stdout, stderr, log_out


# vg map
def giraffe(
        sample_name: str,
        fastq_file1: str,
        fastq_file2: str,
        threads: float,
        index_dir: str,
        restart: bool
):
    """
    :param sample_name: 样品名
    :param fastq_file1: 测序文件1
    :param fastq_file2: 测序文件2
    :param threads: 线程数
    :param index_dir: 索引的路径
    :param restart: 是否检查文件是否存在，并跳过该步骤
    :return: stdout, stderr, log_out
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # 索引和输出文件路径
    xg_file = os.path.join(index_dir, "out.xg")
    gbz_file = os.path.join(index_dir, "out.giraffe.gbz")
    minimizer_file = os.path.join(index_dir, "out.min")
    dist_file = os.path.join(index_dir, "out.dist")
    gam_file = os.path.abspath("{}.gam".format(sample_name))

    # map
    if fastq_file2:
        cmd = 'vg giraffe ' \
              '-t {} ' \
              '-x {} ' \
              '-Z {} ' \
              '-m {} ' \
              '-d {} ' \
              '-f {} ' \
              '-f {} ' \
              '1>{}'.format(threads,
                            xg_file,
                            gbz_file,
                            minimizer_file,
                            dist_file,
                            fastq_file1,
                            fastq_file2,
                            gam_file)
    else:
        cmd = 'vg giraffe ' \
              '-t {} ' \
              '-x {} ' \
              '-Z {} ' \
              '-m {} ' \
              '-d {} ' \
              '-f {} ' \
              '1>{}'.format(threads,
                            xg_file,
                            gbz_file,
                            minimizer_file,
                            dist_file,
                            fastq_file1,
                            gam_file)

    # 检查文件是否存在
    if restart:
        # 检查文件
        file_size = getsize(
            gam_file
        )
        # 如果小于等于0
        if file_size <= 0:
            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "vg.giraffe")
    else:  # 如果没有指定restart，直接运行
        # 提交任务
        stdout, stderr, log_out = run_cmd.run(cmd, "vg.giraffe")

    return stdout, stderr, log_out


# call
def vg_call(
        sample_name: str,
        threads: int,
        depth: float,
        index_dir: str,
        restart: bool
):
    """
    :param sample_name: 样本名
    :param threads: 线程数
    :param depth: 深度
    :param index_dir: 索引的路径
    :param restart: 是否检查文件是否存在，并跳过该步骤
    :return: stdout, stderr, log_out, vcf_file
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # 输出文件路径
    xg_file = os.path.join(index_dir, "out.xg")
    gam_file = os.path.abspath("{}.gam".format(sample_name))
    pack_file = os.path.abspath("{}.pack".format(sample_name))
    snarls_file = os.path.join(index_dir, "out.snarls")
    vcf_file = os.path.abspath(sample_name + ".vcf")

    # vg pack
    filter_depth = min(3, int(depth / 2))  # 按比对深度大于3来过滤
    cmd = 'vg pack ' \
          '-t {} ' \
          '-Q {} ' \
          '-x {} ' \
          '-g {} ' \
          '-o {}'.format(threads,
                         filter_depth,
                         xg_file,
                         gam_file,
                         pack_file)

    # 检查文件是否存在
    if restart:
        # 检查文件
        file_size = getsize(
            pack_file
        )
        # 如果小于等于0
        if file_size <= 0:
            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "vg.pack")
    else:  # 如果没有指定restart，直接运行
        # 提交任务
        stdout, stderr, log_out = run_cmd.run(cmd, "vg.pack")

    # 检查log是否正常，不正常提前退出
    if log_out:
        return stdout, stderr, log_out, vcf_file

    # vg call
    cmd = 'vg call ' \
          '-t {} ' \
          '-s {} ' \
          '{} ' \
          '-k {} ' \
          '-r {} ' \
          '1>{}'.format(threads,
                        sample_name,
                        xg_file,
                        pack_file,
                        snarls_file,
                        vcf_file)

    # 检查文件是否存在
    if restart:
        # 检查文件
        file_size = getsize(
            vcf_file
        )
        # 如果小于等于0
        if file_size <= 0:
            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "vg.call")
    else:  # 如果没有指定restart，直接运行
        # 提交任务
        stdout, stderr, log_out = run_cmd.run(cmd, "vg.call")

    return stdout, stderr, log_out, vcf_file


def main(
        sample_name: str,
        fastq_file1: str,
        fastq_file2: str,
        threads: int,
        depth: float,
        software: str,
        index_dir: str,
        restart: bool
):
    """
    :param sample_name: 样本名
    :param fastq_file1: 测序文件1
    :param fastq_file2: 测序文件2
    :param threads: 线程数
    :param depth: 深度
    :param software: 软件
    :param index_dir: 索引的路径
    :param restart: 是否检查文件是否存在，并跳过该步骤
    :return:
    """
    if software == "VG-MAP":
        stdout, stderr, log_out = vg_map(
            sample_name,
            fastq_file1,
            fastq_file2,
            threads,
            index_dir,
            restart
        )
        # 如果退出代码有问题，则报错
        if log_out:
            return stdout, stderr, log_out, ""
    elif software == "VG-Giraffe":
        stdout, stderr, log_out = giraffe(
            sample_name,
            fastq_file1,
            fastq_file2,
            threads,
            index_dir,
            restart
        )
        # 如果退出代码有问题，则报错
        if log_out:
            return stdout, stderr, log_out, ""
    else:
        # print log
        log_out = "[" + str(datetime.datetime.now()).split(".")[0] + \
                  '] [vg] software arguments error. (VG-MAP/VG-Giraffe).\n'
        return "", "", log_out, ""

    stdout, stderr, log_out, vcf_out_file = vg_call(
        sample_name,
        threads,
        depth,
        index_dir,
        restart,
    )

    return stdout, stderr, log_out, vcf_out_file
