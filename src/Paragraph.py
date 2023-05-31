#!/usr/bin/python3
# coding=gb2312

import os
import run_cmd
from getsize import getsize


# genotype
def main(
    reference_file: str,
    vcf_file: str,
    bam2paragraph_file: str,
    env_path, 
    threads: int,
    restart: bool
):
    """
    :param reference_file: 参考基因组
    :param vcf_file: vcf文件
    :param bam2paragraph_file: 配置文件
    :param env_path: 环境变量
    :param threads: 线程数
    :param restart: 是否检查文件是否存在，并跳过该步骤
    :return: stdout, stderr, log_out, vcf_out_file
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # run
    cmd = '''tmpn=`mktemp -u paragraph_XXXXX` && 
        tmpd="/tmp/${tmpn}" && 
        mkdir ${tmpd} && 
        echo ${tmpd} && 
        export TMP=${tmpd} &&'''\
          + " multigrmpy.py -i " + vcf_file \
          + " -m " + bam2paragraph_file \
          + " -r " + reference_file \
          + " --threads " + str(threads) \
          + ''' --scratch-dir ${tmpd}''' \
          + " -M 1000 -o ./ && rm -r ${tmpd}"

    vcf_out_file = os.path.join(os.getcwd(), "genotypes.vcf.gz")

    # 检查文件是否存在
    if restart:
        # 检查文件
        file_size = getsize(
            vcf_out_file
        )
        # 如果小于等于0
        if file_size <= 0:
            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "Paragraph.multigrmpy.py", env_path)
    else:  # 如果没有指定restart，直接运行
        # 提交任务
        stdout, stderr, log_out = run_cmd.run(cmd, "Paragraph.multigrmpy.py", env_path)

    return stdout, stderr, log_out, vcf_out_file
