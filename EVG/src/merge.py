#!/usr/bin/python3
# Created on 2022/5/5
# @author: du
# email: dzz0539@163.com

import os
import run_cmd
from getsize import getsize


# 合并结果
def main(
        code_dir: str,
        work_path: str,
        line_vcf_file: str,
        sample_name: str,
        mode: str,
        genotype_vcf_file_list,
        select_software_list,
        restart: bool
):
    """
    :param code_dir: 代码目录
    :param work_path: 工作路径
    :param line_vcf_file: line的vcf
    :param sample_name: 品种名
    :param mode: 模式
    :param genotype_vcf_file_list:
    :param select_software_list:
    :param restart: 是否检查文件是否存在，并跳过该步骤
    :return: stdout, stderr, log_out, sample_name, merge_vcf_file
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # 输出文件
    merge_vcf_file = "{}.vcf.gz".format(sample_name)

    # 回到工作目录
    os.chdir(work_path)

    # 代码的路径
    code_path = os.path.join(code_dir, "src", "graphvcf merge")

    cmd = code_path + " -v " + line_vcf_file + " "
    for index in range(len(genotype_vcf_file_list)):
        file = genotype_vcf_file_list[index]
        software = select_software_list[index]
        cmd += "--" + software + " " + file + " "
    cmd += "-n {} -m {} -o {}".format(sample_name, mode, merge_vcf_file)

    # 检查文件是否存在
    if restart:
        # 检查文件
        file_size = getsize(
            "{}.vcf.gz".format(sample_name)
        )
        # 如果小于等于0
        if file_size <= 0:
            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "graphvcf merge")
    else:  # 如果没有指定restart，直接运行
        # 提交任务
        stdout, stderr, log_out = run_cmd.run(cmd, "graphvcf merge")

    merge_vcf_file = os.path.abspath(merge_vcf_file)

    return stdout, stderr, log_out, sample_name, merge_vcf_file
