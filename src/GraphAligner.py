#!/usr/bin/python3
# coding=gb2312

import os
import run_cmd
from getsize import getsize


# index
def vg_index(
    reference_file: str,
    vcf_file: str,
    threads: int,
    index_dir: str,
    env_path, 
    restart: bool
):
    """
    :param reference_file: 参考基因组
    :param vcf_file: vcf文件
    :param threads: 线程数
    :param index_dir: 索引路径
    :param env_path: 环境变量
    :param restart:
    :return: 是否检查文件是否存在，并跳过该步骤
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # 更改工作路径
    os.chdir(index_dir)

    # 文件前缀
    cmd = "vg construct -t {} -r {} -v {} 1>out.vg && " \
          "mkdir temp && vg index -t {} -b temp/ -x out.xg out.vg && rm -rf temp && " \
          "vg snarls -t {} out.xg 1>out.snarls".format(threads, reference_file, vcf_file, threads, threads)

    # 检查文件是否存在
    if restart:
        # 如果小于等于0
        if getsize("out.snarls") <= 0 or \
                getsize("out.vg") <= 0 or \
                getsize("out.xg") <= 0:
            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "GraphAligner.vg.index", env_path)
    else:  # 如果没有指定restart，直接运行
        # 提交任务
        stdout, stderr, log_out = run_cmd.run(cmd, "GraphAligner.vg.index", env_path)

    return stdout, stderr, log_out


# vg map
def graphaligner(
    sample_name: str,
    fastq_file: str,
    threads: int,
    index_dir: str,
    env_path, 
    restart: bool
):
    """
    :param sample_name: 样品名
    :param fastq_file: 测序文件
    :param threads: 线程数
    :param index_dir: 索引的路径
    :param env_path: 环境变量
    :param restart: 是否检查文件是否存在，并跳过该步骤
    :return:
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # 输出文件路径
    vg_file = os.path.join(index_dir, "out.vg")
    gam_file = os.path.abspath("{}.gam".format(sample_name))

    # map
    cmd = 'GraphAligner ' \
          '-t {} ' \
          '-g {} ' \
          '-x vg ' \
          '-f {} ' \
          '-a {}'.format(threads, vg_file, fastq_file, gam_file)

    # 检查文件是否存在
    if restart:
        # 检查文件
        file_size = getsize(
            gam_file
        )
        # 如果小于等于0
        if file_size <= 0:
            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "GraphAligner.graphaligner", env_path)
    else:  # 如果没有指定restart，直接运行
        # 提交任务
        stdout, stderr, log_out = run_cmd.run(cmd, "GraphAligner.graphaligner", env_path)

    return stdout, stderr, log_out


# call
def vg_call(
    sample_name: str,
    threads: int,
    depth: float,
    index_dir: str,
    env_path, 
    restart: bool
):
    """
    :param sample_name: 样本名
    :param threads: 线程数
    :param depth: 深度
    :param index_dir: 索引的路径
    :param env_path: 环境变量
    :param restart: 是否检查文件是否存在，并跳过该步骤
    :return:
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # 输出文件路径
    xg_file = os.path.join(index_dir, "out.xg")
    snarls_file = os.path.join(index_dir, "out.snarls")
    gam_file = os.path.abspath("{}.gam".format(sample_name))
    pack_file = os.path.abspath("{}.pack".format(sample_name))
    vcf_file = os.path.abspath(sample_name + ".vcf")

    # vg pack
    filter_depth = min(0, int(depth/2))  # 按比对深度大于0来过滤
    cmd = 'vg pack ' \
          '-t {} ' \
          '-Q {} ' \
          '-x {} ' \
          '-g {} ' \
          '-o {}'.format(threads, filter_depth, xg_file, gam_file, pack_file)

    # 检查文件是否存在
    if restart:
        # 检查文件
        file_size = getsize(
            pack_file
        )
        # 如果小于等于0
        if file_size <= 0:
            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "GraphAligner.vg.pack", env_path)
    else:  # 如果没有指定restart，直接运行
        # 提交任务
        stdout, stderr, log_out = run_cmd.run(cmd, "GraphAligner.vg.pack", env_path)

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
          '1>{}'.format(threads, sample_name, xg_file, pack_file, snarls_file, vcf_file)

    # 检查文件是否存在
    if restart:
        # 检查文件
        file_size = getsize(
            vcf_file
        )
        # 如果小于等于0
        if file_size <= 0:
            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "GraphAligner.vg.call", env_path)
    else:  # 如果没有指定restart，直接运行
        # 提交任务
        stdout, stderr, log_out = run_cmd.run(cmd, "GraphAligner.vg.call", env_path)

    return stdout, stderr, log_out, vcf_file


def main(
    sample_name: str,
    fastq_file: str,
    threads: int,
    depth: float,
    index_dir: str,
    env_path, 
    restart: bool
):
    """
    :param sample_name: 样本名
    :param fastq_file: 测序文件
    :param threads: 线程数
    :param depth: 深度
    :param index_dir: 索引的路径
    :param env_path: 环境变量
    :param restart: 是否检查文件是否存在，并跳过该步骤
    :return:
    """

    stdout, stderr, log_out = graphaligner(
        sample_name,
        fastq_file,
        threads,
        index_dir,
        env_path, 
        restart
    )
    # 如果退出代码有问题，则报错
    if log_out:
        return stdout, stderr, log_out, ""

    stdout, stderr, log_out, vcf_out_file = vg_call(
        sample_name,
        threads,
        depth,
        index_dir,
        env_path, 
        restart
    )

    return stdout, stderr, log_out, vcf_out_file
