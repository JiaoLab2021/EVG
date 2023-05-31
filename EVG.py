#!/usr/bin/python3
# coding=gb2312

__data__ = "2023/05/30"
__version__ = "1.0.1"
__author__ = "Zezhen Du"
__email__ = "dzz0539@gmail.com or dzz0539@163.com"

import os
import argparse
import sys
import shutil
from multiprocessing import Pool
import concurrent.futures
import time
import logging

# 代码路径
code_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(code_dir, "src/"))
import convert, select_software, merge, bwa, GraphTyper2, \
    GraphAligner, BayesTyper, Paragraph, VG_MAP_Giraffe, PanGenie

# 环境变量
env_path = {'PATH': os.environ.get('PATH')}

# 全局参数
force = False
restart = False
threads = 10
jobs_num = 3
mode = "precise"
need_depth = 15.0
merge_mode = "specific"


# log
logger = logging.getLogger('SynDiv')
formatter = logging.Formatter('[%(asctime)s] %(message)s')
handler = logging.StreamHandler()  # 输出到控制台
handler.setFormatter(formatter)
logger.addHandler(handler)


# 打印merge的结果
def print_merge(
    vcf_merge_files_map: dict
):
    """
    :param vcf_merge_files_map: 合并后的vcf文件 map<sample_name, vcf_file>
    :return:
    """
    for key, value in vcf_merge_files_map.items():
        log = "{}: {}\n".format(key, value)
        logger.error(log)

    return 0


# 创建目录
def makedir(
        path_dir: str,
        force_tmp: bool,
        restart_tmp: bool
):
    """
    :param path_dir: 需要创建的文件夹路径
    :param force_tmp: 是否强制删除已经存在的目录
    :param restart_tmp 是否根据生成文件的状态继续程序。如果指定该参数，且force没被指定，软件会跳过已经存在的目录，不再报错
    :return: 0
    """
    if os.path.isdir(path_dir):
        if force_tmp:  # 如果强制删除，则清空目录后重新创建一个
            shutil.rmtree(path_dir)
            os.makedirs(path_dir)
            log = '[EVG.makedir] \'{}\' already exists, clear and recreate.'.format(path_dir)
            logger.error(log)
        elif restart_tmp:
            log = '[EVG.makedir] \'{}\' already exists, restart is used, ' \
                  'skipping emptying folders.\n'.format(path_dir)
            logger.error(log)
        else:  # 如果不强制删除，打印警告并退出代码
            log = '[EVG.makedir] \'{}\' already exists. used --force to overwrite or --restart to restart workflow.\n'.format(path_dir)
            logger.error(log)
            raise SystemExit(1)
    else:
        os.makedirs(path_dir)


# 解析samples文件
def get_samples_path(
        samples_file: str
):
    """
    :param samples_file: 包含测序数据的样品名称的文件 'name\tread1_path\tread2_path'
    :return samples_map  map<sample_name, list<fastq_path>>
    """
    samples_map = {}
    with open(samples_file, 'r') as f:
        samples_list = f.readlines()

    # 判断samples_list长度
    if len(samples_list) == 0:
        log = '[EVG.get_parser] Error: empty file -> {}.\n'.format(samples_file)
        logger.error(log)
        # os.killpg(os.getpgid(os.getpid()), signal.SIGKILL)
        raise SystemExit(1)

    for sample in samples_list:
        if not sample or "#" in sample or "read1_path" in sample:
            continue

        sample_list = sample.split()
        if len(sample_list) > 1:
            name = sample_list[0]
            read1_path = sample_list[1]
            if len(sample_list) > 2:  # 双端测序
                read2_path = sample_list[2]
                samples_map[name] = [os.path.abspath(read1_path), os.path.abspath(read2_path)]
            else:  # 单端测序
                samples_map[name] = [os.path.abspath(read1_path), ""]
    return samples_map


def get_parser():
    """
    :return: parser_map = {
        "reference_file": reference_file,
        "vcf_file": vcf_file,
        "samples_file": samples_file,
        "vcf_out_name": vcf_out_name,
        "path": path,
        "software_list": software_list,
        "line_vcf_file": line_vcf_file
    }
    """
    # log
    logger = logging.getLogger('getParser')

    logger.error(f"data: {__data__}")
    logger.error(f"version: {__version__}")
    logger.error(f"author: {__author__}")
    logger.error(f"\nIf you encounter any issues related to the code, please don't hesitate to contact us via email at {__email__}.\n")

    # 参数
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    input_option = parser.add_argument_group(title="required input")
    input_option.add_argument("-r", dest="reference", help="input FASTA reference",
                              type=argparse.FileType('r'), required=True)
    input_option.add_argument("-v", dest="vcf", help="input the merged vcf file",
                              type=argparse.FileType('r'), required=True)
    # 解析samples文件
    input_option.add_argument("-s", dest="samples", help="samples file (name read1_path read2_path)",
                              type=argparse.FileType('r'), required=True)

    optional_option = parser.add_argument_group(title="optional input")
    optional_option.add_argument("-V", dest="VCF", help="input line-specific vcf file [-v]",
                                 type=argparse.FileType('r'))

    # 设置取多少×的数据来分型
    depth_option = parser.add_argument_group(title="depth")
    depth_option.add_argument("--depth", dest="depth", help="read depth for genotyping [15]",
                              type=float, default=15)

    # 用户自定义软件
    algorithm_option = parser.add_argument_group(
        title="algorithm"
    )
    algorithm_option.add_argument(
        "--software", dest="software", help="genotyping software [auto]", type=str, nargs='+',
        choices=['VG-MAP', 'VG-Giraffe', 'GraphAligner', 'Paragraph', 'BayesTyper', 'GraphTyper2', 'PanGenie'],
        default="auto"
    )
    # 进程
    algorithm_option.add_argument(
        "--mode", dest="mode", help="software mode [precise]", type=str, choices=['precise', 'fast'], default="precise"
    )

    # 并行参数
    parallel_option = parser.add_argument_group(
        title="parallel"
    )
    # 线程
    parallel_option.add_argument(
        "--threads", dest="threads", help="number of threads [10]", type=int, default=10
    )
    # 进程数
    parallel_option.add_argument(
        "--jobs", dest="jobs", help="run n jobs in parallel [3]", type=int, default=3
    )

    # 处理已经存在的目录和文件的方式
    resume_option = parser.add_argument_group(
        title="resume"
    )
    resume_option.add_argument(
        "--force", action="store_true", help="erase already existing output directory [False]"
    )
    resume_option.add_argument(
        "--restart", action="store_true", help="attempt to restart existing workflow based "
                                               "on the state of output files [False]"
    )

    args = parser.parse_args()

    # 解析参数
    reference_file = args.reference.name
    vcf_file = args.vcf.name
    try:
        line_vcf_file = args.VCF.name
    except AttributeError:
        line_vcf_file = vcf_file

    # 测序数据文件
    samples_file = args.samples.name

    vcf_out_name = "convert." + os.path.basename(vcf_file.replace(".gz", ""))
    software_list = args.software

    # 修改全局变量
    global force, restart, threads, jobs_num, mode, need_depth, merge_mode
    force = args.force
    restart = args.restart
    threads = args.threads
    jobs_num = args.jobs
    mode = args.mode
    need_depth = args.depth
    if line_vcf_file == vcf_file:  # 合并过程算法，如果一样，用all算法
        merge_mode = "all"

    # 将路径补全
    reference_file = os.path.abspath(reference_file)
    vcf_file = os.path.abspath(vcf_file)
    line_vcf_file = os.path.abspath(line_vcf_file)
    samples_file = os.path.abspath(samples_file)

    # 当前路径
    path = os.getcwd()

    # 输出参数字典
    parser_map = {
        "reference_file": reference_file,
        "vcf_file": vcf_file,
        "samples_file": samples_file,
        "vcf_out_name": vcf_out_name,
        "path": path,
        "software_list": software_list,
        "line_vcf_file": line_vcf_file
    }

    return parser_map


# bwa
def run_bwa(
        reference_file: str,
        sample_name: str,
        fastq_file1: str,
        fastq_file2: str,
        real_depth: int,
        read_len: float,
        work_path: str
):
    """
    :param reference_file: 转换后的参考基因组
    :param sample_name: 样本名字
    :param fastq_file1: 测序文件1
    :param fastq_file2: 测序文件2
    :param real_depth: 测序数据的深度
    :param read_len: 数据读长
    :param work_path: 工作路径
    :return: stdout, stderr, log_out, bam_info_map = {
            "sample_name": sample_name,
            "fastq_file1": fastq_file1,
            "fastq_file2": fastq_file2,
            "real_depth": real_depth,
            "read_len": read_len,
            "bam_file": bam_file
        }
    """

    os.chdir(work_path)
    stdout, stderr, log_out, bam_file = bwa.mem(
        reference_file,
        threads,
        sample_name, 
        env_path, 
        restart,
        fastq_file1,
        fastq_file2
    )

    # 如果退出代码有问题，则报错
    if log_out:
        return stdout, stderr, log_out, {}

    # fastq路径、测序深度、读长
    bam_info_map = {
        "sample_name": sample_name,
        "fastq_file1": fastq_file1,
        "fastq_file2": fastq_file2,
        "real_depth": real_depth,
        "read_len": read_len,
        "bam_file": bam_file
    }

    return stdout, stderr, log_out, bam_info_map


# GraphTyper2
def run_graphtyper(
        reference_file: str,
        vcf_file: str,
        region_file: str,
        bam2graphtyper_file: str,
        work_path: str
):
    """
    :param reference_file: 参考基因组
    :param vcf_file: 变异文件
    :param region_file: 配置文件
    :param bam2graphtyper_file: bam配置文件路径
    :param work_path: genotype工作路径
    :return: "GraphTyper2", sample_name, graphtyper_vcf_file
    """

    # 创建文件夹并切换路径
    os.chdir(work_path)
    graphtyper_path = os.path.join(work_path, "GraphTyper2")
    makedir(
        graphtyper_path,
        force_tmp=force,
        restart_tmp=restart
    )
    os.chdir(graphtyper_path)

    # genotype
    stdout, stderr, log_out, graphtyper_vcf_file = GraphTyper2.main(
        reference_file,
        vcf_file,
        bam2graphtyper_file,
        region_file, 
        env_path, 
        threads,
        restart
    )

    return stdout, stderr, log_out, "GraphTyper2", "", graphtyper_vcf_file


# Paragraph
def run_paragraph(
        reference_file: str,
        vcf_file: str,
        bam2paragraph_file: str,
        work_path: str
):
    """
    :param reference_file: 参考基因组
    :param vcf_file: 编译文件
    :param bam2paragraph_file: bam配置文件路径
    :param work_path: genotype工作路径
    :return: "Paragraph", sample_name, paragraph_vcf_file
    """

    # 创建文件夹并切换路径
    os.chdir(work_path)
    paragraph_path = os.path.join(work_path, "Paragraph")
    makedir(
        paragraph_path,
        force_tmp=force,
        restart_tmp=restart
    )
    os.chdir(paragraph_path)

    # genotype
    stdout, stderr, log_out, paragraph_vcf_file = Paragraph.main(
        reference_file,
        vcf_file,
        bam2paragraph_file, 
        env_path, 
        threads,
        restart
    )

    return stdout, stderr, log_out, "Paragraph", "", paragraph_vcf_file


# BayesTyper
def run_bayestyper(
        reference_file: str,
        vcf_file: str,
        bam2bayestyper_file: str,
        work_path: str,
        bam_infos_map: str
):
    """
    :param reference_file: 参考基因组
    :param vcf_file: 变异文件
    :param bam2bayestyper_file: 配置文件路径
    :param work_path: genotype工作路径
    :param bam_infos_map: bam所有信息
    :return:
    """

    # 创建文件夹并切换路径
    os.chdir(work_path)
    bayestyper_path = os.path.join(work_path, "BayesTyper")
    makedir(
        bayestyper_path,
        force_tmp=force,
        restart_tmp=restart
    )
    os.chdir(bayestyper_path)

    # genotype
    stdout, stderr, log_out, bayestyper_vcf_file = BayesTyper.main(
        reference_file,
        vcf_file,
        bam2bayestyper_file,
        bam_infos_map, 
        env_path, 
        threads,
        restart
    )

    return stdout, stderr, log_out, "BayesTyper", "", bayestyper_vcf_file


# vg_map_giraffe
def run_vg_map_giraffe(
        sample_name: str,
        fastq_file1: str,
        fastq_file2: str,
        real_depth: int,
        software: str,
        software_work_path: str,
        index_dir: str
):
    """
    :param sample_name: 品种名称
    :param fastq_file1: 测序文件1
    :param fastq_file2: 测序文件2
    :param real_depth: 测序深度
    :param software: 软件
    :param software_work_path: 软件工作路径
    :param index_dir: 索引路径
    :return: software, sample_name, map_vcf_file
    """

    # 创建文件夹并切换路径
    os.chdir(software_work_path)
    map_path = os.path.join(software_work_path, sample_name)
    makedir(
        map_path,
        force_tmp=force,
        restart_tmp=restart
    )
    os.chdir(map_path)

    # genotype
    stdout, stderr, log_out, map_vcf_file = VG_MAP_Giraffe.main(
        sample_name,
        fastq_file1,
        fastq_file2, 
        env_path, 
        threads,
        real_depth,
        software,
        index_dir,
        restart
    )

    return stdout, stderr, log_out, software, sample_name, map_vcf_file


# GraphAligner
def run_graphaligner(
        sample_name: str,
        fastq_file1: str,
        real_depth: int,
        software_work_path: str,
        index_dir: str
):
    """
    :param sample_name: 品种名称
    :param fastq_file1: 测序文件1
    :param real_depth: 测序深度
    :param software_work_path: 软件工作路径
    :param index_dir: 索引路径
    :return: "GraphAligner", sample_name, graphaligner_vcf_file
    """

    # genotype
    os.chdir(software_work_path)
    graphaligner_path = os.path.join(software_work_path, sample_name)
    makedir(
        graphaligner_path,
        force_tmp=force,
        restart_tmp=restart
    )
    os.chdir(graphaligner_path)

    stdout, stderr, log_out, graphaligner_vcf_file = GraphAligner.main(
        sample_name,
        fastq_file1,
        threads,
        real_depth,
        index_dir,
        env_path, 
        restart
    )

    return stdout, stderr, log_out, "GraphAligner", sample_name, graphaligner_vcf_file


# PanGenie
def run_pangenie(
    sample_name: str,
    reference_file: str,
    vcf_file: str,
    fastq_file1: str,
    fastq_file2: str,
    software_work_path: str
):
    """

    :param sample_name: 品种名
    :param reference_file: 参考基因组
    :param vcf_file: 变异文件
    :param fastq_file1: 测序文件1
    :param fastq_file2: 测序文件2
    :param software_work_path: 软件工作路径
    :return: "PanGenie", sample_name, pangenie_vcf_file
    """

    # 创建文件夹并切换路径
    os.chdir(software_work_path)
    pangenie_path = os.path.join(software_work_path, sample_name)
    makedir(
        pangenie_path,
        force_tmp=force,
        restart_tmp=restart
    )
    os.chdir(pangenie_path)

    # genotype
    stdout, stderr, log_out, pangenie_vcf_file = PanGenie.main(
        reference_file,
        vcf_file,
        sample_name,
        fastq_file1,
        fastq_file2,
        env_path, 
        threads,
        restart
    )

    return stdout, stderr, log_out, "PanGenie", sample_name, pangenie_vcf_file


# 对参考基因组和vcf文件转换
def ref_vcf_convert(
        parser_map,
        work_path
):
    """
    :param parser_map: get_parser()返回的参数字典
    :param work_path: 工作路径
    :return: stdout, stderr, log_out, convert_out_map = {
                "reference_file": reference_file,
                "vcf_out_name": vcf_out_name,
                "vcf_sv_out_name": vcf_sv_out_name,
                "region_file": region_file,
                "fasta_base": fasta_base,
                "map_index_dir": map_index_dir,
                "giraffe_index_dir": giraffe_index_dir
            }
    """
    # 更改工作路径
    os.chdir(work_path)

    # ################################### reference进行转换 ###################################
    stdout, stderr, log_out, reference_file = convert.convert_reference(
        parser_map["reference_file"],
        env_path, 
        restart
    )
    # 如果退出代码有问题，则报错
    if log_out:
        return stdout, stderr, log_out, {}
    # 对fasta文件进行评估
    stdout, stderr, log_out, fasta_base = convert.fasta_count(code_dir, env_path, reference_file)
    # 如果退出代码有问题，则报错
    if log_out:
        return stdout, stderr, log_out, {}

    # ################################### 对reference构建索引 ###################################
    stdout, stderr, log_out = bwa.index(
        reference_file, 
        env_path, 
        restart
    )
    # 如果退出代码有问题，则报错
    if log_out:
        return stdout, stderr, log_out, {}

    # ################################### 对vcf文件进行转换 ###################################
    stdout, stderr, log_out, vcf_out_name, vcf_sv_out_name, region_file = convert.vcf_convert(
        code_dir,
        reference_file,
        parser_map["vcf_file"],
        350,
        parser_map["vcf_out_name"],
        mode, 
        env_path, 
        restart
    )
    # 如果退出代码有问题，则报错
    if log_out:
        return stdout, stderr, log_out, {}

    # ################################### 对vcf文件压缩 ###################################
    if vcf_sv_out_name != vcf_out_name:  # 如果两个文件一样，则不对vcf_sv_out_name进行压缩，否则会报错
        stdout, stderr, log_out, vcf_sv_out_name = convert.bgzip_vcf(
            vcf_sv_out_name,
            env_path, 
            restart
        )
        # 如果退出代码有问题，则报错
        if log_out:
            return stdout, stderr, log_out, {}
        stdout, stderr, log_out, vcf_out_name = convert.bgzip_vcf(
            vcf_out_name,
            env_path, 
            restart
        )
        # 如果退出代码有问题，则报错
        if log_out:
            return stdout, stderr, log_out, {}
    else:
        stdout, stderr, log_out, vcf_out_name = convert.bgzip_vcf(
            vcf_out_name,
            env_path, 
            restart
        )
        vcf_sv_out_name = vcf_out_name
        # 如果退出代码有问题，则报错
        if log_out:
            return stdout, stderr, log_out, {}

    # ################################### 构建vg图和索引 ###################################
    # 先选择vg的软件
    select_software_list = []
    if isinstance(parser_map["software_list"], str):  # 程序自动判定软件
        if fasta_base > 200000000:  # 如果基因组大于200Mb，用giraffe跑
            select_software_list.append("giraffe")
        else:  # 否则用vg_map跑
            select_software_list.append("map")
    else:  # 用户自定义软件
        if "VG-MAP" in parser_map["software_list"]:
            select_software_list.append("map")
        if "VG-Giraffe" in parser_map["software_list"]:
            select_software_list.append("giraffe")

    # 多进程进程池
    if len(select_software_list) > 0:  # 如果文件数量大于1了再构建索引
        # 多进程
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
            to_do = []
            # map和giraffe构建索引
            for software in select_software_list:
                # 更改工作路径
                index_dir = os.path.join(work_path, software)
                makedir(
                    index_dir,
                    force_tmp=force,
                    restart_tmp=restart
                )
                os.chdir(index_dir)

                # 提交任务
                future = executor.submit(
                    VG_MAP_Giraffe.vg_autoindex,
                    reference_file,
                    vcf_out_name,
                    software, 
                    env_path, 
                    threads,
                    index_dir,
                    restart
                )
                # 保存结果
                to_do.append(future)

                time.sleep(0.05)

            # 回到主工作路径
            os.chdir(work_path)

            # GraphAligner构建索引
            if "GraphAligner" in parser_map["software_list"]:
                # 更改工作路径
                index_dir = os.path.join(work_path, "GraphAligner")
                makedir(
                    index_dir,
                    force_tmp=force,
                    restart_tmp=restart
                )
                os.chdir(index_dir)

                # 提交任务
                future = executor.submit(
                    GraphAligner.vg_index,
                    reference_file,
                    vcf_out_name,
                    threads,
                    index_dir, 
                    env_path, 
                    restart
                )
                # 保存结果
                to_do.append(future)

            # 返回主路径
            os.chdir(work_path)

            # 获取多线程返回值，如果log_out存在表明运行失败，退出代码
            for future in concurrent.futures.as_completed(to_do):  # 并发执行
                stdout, stderr, log_out = future.result()  # 检查返回值
                if log_out:
                    return stdout, stderr, log_out, {}

    # 索引文件的路径
    map_index_dir = os.path.join(work_path, "map")
    giraffe_index_dir = os.path.join(work_path, "giraffe")
    graphaligner_index_dir = os.path.join(work_path, "GraphAligner")

    # 存储转化后文件的路径
    convert_out_map = {
        "reference_file": reference_file,
        "vcf_out_name": vcf_out_name,
        "vcf_sv_out_name": vcf_sv_out_name,
        "region_file": region_file,
        "fasta_base": fasta_base,
        "map_index_dir": map_index_dir,
        "giraffe_index_dir": giraffe_index_dir,
        "GraphAligner_index_dir": graphaligner_index_dir
    }

    return stdout, stderr, log_out, convert_out_map


# 对测序文件进行转换
def read_convert(
        work_path: str,
        sample_name: str,
        fastq_file1: str,
        fastq_file2: str,
        fasta_base: float
):
    """
    :param work_path: 工作路径
    :param sample_name: 样本名
    :param fastq_file1: 测序文件1
    :param fastq_file2: 测序文件2
    :param fasta_base: 参考基因组base
    :return: stdout, stderr, log_out, fastq_info_map = {
                "sample_name": sample_name,
                "fastq_file1": fastq_out_file1,
                "fastq_file2": fastq_out_file2,
                "real_depth": real_depth,
                "read_len": read_len
            }
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # 更改工作路径
    os.chdir(work_path)

    # ################################### 对二代数据进行转换 ###################################
    # 对fastq文件进行评估
    if fastq_file1:
        stdout, stderr, log_out, fastq_base, read_num, read_len = convert.fastq_count(
            code_dir, 
            env_path, 
            fastq_file1,
            fastq_file2
        )
        # 如果退出代码有问题，则报错
        if log_out:
            return stdout, stderr, log_out, {}

        # downsample
        stdout, stderr, log_out, fastq_out_file1, fastq_out_file2 = convert.downsample(
            code_dir,
            fastq_file1,
            fastq_file2,
            fastq_base,
            fasta_base,
            need_depth, 
            env_path, 
            restart
        )
        # 如果退出代码有问题，则报错
        if log_out:
            return stdout, stderr, log_out, {}

        # 实际的测序深度
        real_depth = min(int(fastq_base) / int(fasta_base), need_depth)
    else:
        read_len = 0
        fastq_out_file1, fastq_out_file2 = "", ""
        real_depth = 0

    # 判断测序深度是不是0
    if real_depth == 0:
        log_out = '[EVG.run_convert] sequencing depth is 0, please check the input file file.'
        return stdout, stderr, log_out, {}

    # fastq路径、测序深度、读长
    fastq_info_map = {
        "sample_name": sample_name,
        "fastq_file1": fastq_out_file1,
        "fastq_file2": fastq_out_file2,
        "real_depth": real_depth,
        "read_len": read_len
    }

    return stdout, stderr, log_out, fastq_info_map


# 每个line单独跑
def run_genotype(
        work_path: str,
        parser_map,
        convert_out_map,
        bam_infos_map
):
    """
    :param work_path: genotype工作路径
    :param parser_map: get_parser的返回值
    :param convert_out_map: vcf、reference和测序文件转换后的哈希表
    :param bam_infos_map: run_bwa返回值
    :return: stdout, stderr, log_out, vcf_merge_files_map
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # 更改工作路径
    os.chdir(work_path)

    # 获取测序数据的最低测序深度和测序长度
    depth_min = 1000
    read_len_mix = 100000
    for key, value in bam_infos_map.items():
        depth_min = min(depth_min, value["real_depth"])
        read_len_mix = min(read_len_mix, value["read_len"])

    # run select_software
    if isinstance(parser_map["software_list"], str):  # 程序自动判定软件
        select_software_list = select_software.main(depth_min, read_len_mix, convert_out_map["fasta_base"])
    else:  # 用户自定义软件
        select_software_list = parser_map["software_list"]

    # ################################### config ###################################
    # 生成GraphTyper2、BayesTyper和Paragraph的配置文件
    bam2graphtyper_file = ""
    bam2bayestyper_file = ""
    bam2paragraph_file = ""
    if "GraphTyper2" in select_software_list:
        bam2graphtyper_file = bwa.bam2graphtyper(
            bam_infos_map
        )
    if "BayesTyper" in select_software_list:
        bam2bayestyper_file = bwa.bam2bayestyper(
            bam_infos_map
        )
    if "Paragraph" in select_software_list:
        bam2paragraph_file = bwa.bam2paragraph(
            bam_infos_map
        )

    # ################################### genotype ###################################
    # run genotype
    # 多进程进程池
    pool = Pool(processes=jobs_num)

    # 线程池返回值
    pool_out_list = []  # 从异步提交任务获取结果

    # 保存分型结果位置
    genotype_outs_map = {}  # map<sample_name, map<software, vcf_path>>

    for software in select_software_list:
        # 回到工作路径
        os.chdir(work_path)

        if software in ["VG-MAP", "VG-Giraffe", "GraphAligner", "PanGenie"]:
            # 创建文件夹并切换路径
            os.chdir(work_path)
            software_work_path = os.path.join(work_path, software)
            makedir(
                software_work_path,
                force_tmp=force,
                restart_tmp=restart
            )
            os.chdir(software_work_path)

            for key, value in bam_infos_map.items():
                # 回到工作路径
                os.chdir(work_path)

                sample_name = key
                fastq_file1 = value["fastq_file1"]
                fastq_file2 = value["fastq_file2"]
                real_depth = value["real_depth"]

                if "VG-MAP" == software:
                    # 多进程提交任务
                    pool_out = pool.apply_async(run_vg_map_giraffe, args=(
                        sample_name,
                        fastq_file1,
                        fastq_file2,
                        real_depth,
                        software,
                        software_work_path,
                        convert_out_map["map_index_dir"],
                    ), error_callback=throw_exception)
                    # 保存多线程的返回值
                    pool_out_list.append(pool_out)
                elif "VG-Giraffe" == software:
                    # 多进程提交任务
                    pool_out = pool.apply_async(run_vg_map_giraffe, args=(
                        sample_name,
                        fastq_file1,
                        fastq_file2,
                        real_depth,
                        software,
                        software_work_path,
                        convert_out_map["giraffe_index_dir"],
                    ), error_callback=throw_exception)
                    # 保存多线程的返回值
                    pool_out_list.append(pool_out)
                elif "GraphAligner" == software:
                    # 多进程提交任务
                    pool_out = pool.apply_async(run_graphaligner, args=(
                        sample_name,
                        fastq_file1,
                        real_depth,
                        software_work_path,
                        convert_out_map["GraphAligner_index_dir"],
                    ), error_callback=throw_exception)
                    # 保存多线程的返回值
                    pool_out_list.append(pool_out)
                else:
                    # 多进程提交任务
                    pool_out = pool.apply_async(run_pangenie, args=(
                        sample_name,
                        convert_out_map["reference_file"],
                        convert_out_map["vcf_out_name"],
                        fastq_file1,
                        fastq_file2,
                        software_work_path,
                    ), error_callback=throw_exception)
                    # 保存多线程的返回值
                    pool_out_list.append(pool_out)
                # 线程间隔
                time.sleep(0.05)
                # 回到工作路径
                os.chdir(work_path)
        elif software == "Paragraph":
            # 多进程提交任务
            pool_out = pool.apply_async(run_paragraph, args=(
                convert_out_map["reference_file"],
                convert_out_map["vcf_sv_out_name"],
                bam2paragraph_file,
                work_path,
            ), error_callback=throw_exception)
            # 保存多线程的返回值
            pool_out_list.append(pool_out)
        elif software == "BayesTyper":
            # 多进程提交任务
            pool_out = pool.apply_async(run_bayestyper, args=(
                convert_out_map["reference_file"],
                convert_out_map["vcf_out_name"],
                bam2bayestyper_file,
                work_path,
                bam_infos_map
            ), error_callback=throw_exception)
            # 保存多线程的返回值
            pool_out_list.append(pool_out)
        elif software == "GraphTyper2":
            # 多进程提交任务
            pool_out = pool.apply_async(run_graphtyper, args=(
                convert_out_map["reference_file"],
                convert_out_map["vcf_out_name"],
                convert_out_map["region_file"],
                bam2graphtyper_file,
                work_path,
            ), error_callback=throw_exception)
            # 保存多线程的返回值
            pool_out_list.append(pool_out)

        # 线程间隔
        time.sleep(0.05)
        # 回到工作路径
        os.chdir(work_path)

    # 获取多线程返回值
    for pool_out in pool_out_list:
        stdout, stderr, log_out, software, sample_name, vcf_file = pool_out.get()

        # 如果退出代码有问题，则报错
        if log_out:
            return stdout, stderr, log_out, {}

        if software in ["VG-MAP", "VG-Giraffe", "GraphAligner", "PanGenie"]:
            if software not in genotype_outs_map:
                genotype_outs_map[software] = {}
            genotype_outs_map[software][sample_name] = vcf_file
        else:
            genotype_outs_map[software] = vcf_file

        # 回到工作路径
        os.chdir(work_path)

    # 关闭线程池
    pool.close()
    pool.join()

    # ################################################ graphvcf merge ################################################
    # 多进程进程池
    pool = Pool(processes=jobs_num*threads)  # 合并步骤用资源少，因此把所有内核用上

    # 线程池返回值
    pool_out_list = []  # 从异步提交任务获取结果

    # 存储合并后的结果
    vcf_merge_files_map = {}  # map<sample_name, vcf_file>

    # 合并结果
    sample_name_tmp_list = bam_infos_map.keys()  # 先获取所有的samples_name
    for sample_name_tmp in sample_name_tmp_list:
        vcf_out_tmp_list = []  # 用于提交任务的临时列表
        software_tmp_list = []  # 用于提交任务的临时列表
        for software in select_software_list:  # 软件没在结果文件中
            if software not in genotype_outs_map.keys():
                log = '[EVG.genotype] {}: {} results are missing, skipped.'.format(sample_name_tmp, software)
                sys.stdout.write(log)
                continue
            software_tmp_list.append(software)
            if isinstance(genotype_outs_map[software], dict):
                if sample_name_tmp not in genotype_outs_map[software].keys():  # sample_name不在结果文件中
                    log = '[EVG.genotype] {}: {} results are missing, skipped.'.format(sample_name_tmp, software)
                    sys.stdout.write(log)
                    continue
                vcf_out_tmp_list.append(genotype_outs_map[software][sample_name_tmp])
            else:
                vcf_out_tmp_list.append(genotype_outs_map[software])

        # 多线程提交 graphvcf merge
        pool_out = pool.apply_async(
            merge.main, args=(
                code_dir,
                work_path,
                parser_map["line_vcf_file"],
                sample_name_tmp,
                merge_mode,
                vcf_out_tmp_list,
                software_tmp_list,
                env_path, 
                restart,
            ), error_callback=throw_exception
        )

        # 保存多线程的返回值
        pool_out_list.append(pool_out)

        # 多线程间隔
        time.sleep(0.05)

    for pool_out in pool_out_list:  # read的返回值
        stdout, stderr, log_out, sample_name, merge_vcf_file = pool_out.get()

        # 如果退出代码有问题，则报错
        if log_out:
            return stdout, stderr, log_out, {}

        # 赋值给总的哈希表
        vcf_merge_files_map[sample_name] = merge_vcf_file

        # 关闭线程池
    pool.close()
    pool.join()

    return stdout, stderr, log_out, vcf_merge_files_map


# 捕获子进程退出状态
def throw_exception(name):
    # print log
    log = '[EVG.genotype] %s' % name.__cause__
    logger.error(log)
    # os.killpg(os.getpgid(os.getpid()), signal.SIGKILL)
    raise SystemExit(1)


def main():
    # 解析参数
    parser_map = get_parser()

    # ################################################ EVG convert ################################################
    # 返回主工作路径
    os.chdir(parser_map["path"])
    # 多进程进程池
    pool = Pool(processes=threads)

    # 线程池返回值
    pool_out_list = []  # 从异步提交任务获取结果

    # 存储多线程结果，fastq路径、测序深度、读长
    fastq_infos_map = {}

    # 获取测序文件路径和输出路径名
    samples_map = get_samples_path(parser_map["samples_file"])

    # reference and vcf 转换存储路径并创建文件夹
    convert_dir = os.path.join(parser_map["path"], "convert")  # 转换文件存放路径
    makedir(
        convert_dir,
        force_tmp=force,
        restart_tmp=restart
    )

    # 返回主工作路径
    os.chdir(parser_map["path"])
    # 对fasta文件进行评估
    stdout, stderr, log_out, fasta_base = convert.fasta_count(
        code_dir, 
        env_path, 
        parser_map["reference_file"]
    )
    # 如果退出代码有问题，则报错
    if log_out:
        logger.error(log_out.strip())
        exit(1)

    # 对reference和vcf进行转换
    # 多线程提交
    pool_convert_out = pool.apply_async(
        ref_vcf_convert, args=(
            parser_map,
            convert_dir,
        ), error_callback=throw_exception
    )

    # 转换read
    # 多线程提交
    for key, value in samples_map.items():
        # 返回主工作路径
        os.chdir(parser_map["path"])

        sample_name = key
        fastq_file1 = value[0]
        fastq_file2 = value[1]

        # 多线程对测序文件进行转换
        pool_out = pool.apply_async(read_convert, args=(
            convert_dir,
            sample_name,
            fastq_file1,
            fastq_file2,
            fasta_base,
        ), error_callback=throw_exception)

        # 保存多线程的返回值
        pool_out_list.append(pool_out)

        # 多线程间隔
        time.sleep(0.05)

    # 获取多线程返回值
    stdout, stderr, log_out, convert_out_map = pool_convert_out.get()  # vcf和reference的返回值
    # 如果退出代码有问题，则报错
    if log_out:
        logger.error(log_out.strip())
        exit(1)
    for pool_out in pool_out_list:  # read的返回值
        stdout, stderr, log_out, fastq_info_map = pool_out.get()
        # 如果退出代码有问题，则报错
        if log_out:
            logger.error(log_out.strip())
            exit(1)
        # 赋值给总的哈希表
        fastq_infos_map[fastq_info_map["sample_name"]] = {
            "fastq_file1": fastq_info_map["fastq_file1"],
            "fastq_file2": fastq_info_map["fastq_file2"],
            "real_depth": fastq_info_map["real_depth"],
            "read_len": fastq_info_map["read_len"]
        }

    # 关闭线程池
    pool.close()
    pool.join()

    # ################################################ bwa mem ################################################
    # 返回主工作路径
    os.chdir(parser_map["path"])
    # 重新初始化线程池
    pool = Pool(processes=jobs_num)

    # 线程池返回值
    pool_out_list = []  # 从异步提交任务获取结果

    # 存储多线程的结果，bam路径、测序深度、读长
    bam_infos_map = {}

    # 序列比对存储路径
    bwa_dir = os.path.join(parser_map["path"], "bwa")  # bwa mem 比对路径

    # 创建目录
    makedir(
        bwa_dir,
        force_tmp=force,
        restart_tmp=restart
    )

    # 多线程进行序列比对
    for key, value in fastq_infos_map.items():
        sample_name = key
        fastq_file1 = value["fastq_file1"]
        fastq_file2 = value["fastq_file2"]
        real_depth = value["real_depth"]
        read_len = value["read_len"]

        # bwa mem
        # 多线程对测序文件进行转换
        pool_out = pool.apply_async(run_bwa, args=(
            convert_out_map["reference_file"],
            sample_name,
            fastq_file1,
            fastq_file2,
            real_depth,
            read_len,
            bwa_dir
        ), error_callback=throw_exception)

        # 保存多线程的返回值
        pool_out_list.append(pool_out)

        # 多线程间隔
        time.sleep(0.05)

    # 获取多线程返回值
    for pool_out in pool_out_list:
        stdout, stderr, log_out, bam_info_map = pool_out.get()
        # 如果退出代码有问题，则报错
        if log_out:
            logger.error(log_out.strip())
            exit(1)
        # 赋值给总的哈希表
        bam_infos_map[bam_info_map["sample_name"]] = {
            "fastq_file1": bam_info_map["fastq_file1"],
            "fastq_file2": bam_info_map["fastq_file2"],
            "real_depth": bam_info_map["real_depth"],
            "read_len": bam_info_map["read_len"],
            "bam_file": bam_info_map["bam_file"]
        }

    # 关闭线程池
    pool.close()
    pool.join()

    # ################################################ genotype ################################################
    # 返回主工作路径
    os.chdir(parser_map["path"])

    # 分型结果路径
    genotype_dir = os.path.join(parser_map["path"], "genotype")  # bwa mem 比对路径
    # 创建目录
    makedir(
        genotype_dir,
        force_tmp=force,
        restart_tmp=restart
    )
    # 分型
    stdout, stderr, log_out, vcf_merge_files_map = run_genotype(
        genotype_dir,
        parser_map,
        convert_out_map,
        bam_infos_map
    )
    # 如果退出代码有问题，则报错
    if log_out:
        logger.error(log_out.strip())
        exit(1)

    # 打印输出结果
    print_merge(vcf_merge_files_map)

    return 0


if __name__ == '__main__':
    main()
