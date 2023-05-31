#!/usr/bin/python3
# coding=gb2312

import os
import concurrent.futures
import run_cmd
from getsize import getsize


# log
import logging
logger = logging.getLogger('SynDiv')
formatter = logging.Formatter('[%(asctime)s] %(message)s')
handler = logging.StreamHandler()  # 输出到控制台
handler.setFormatter(formatter)
logger.addHandler(handler)


# 将reference中的非ATGCNatgcn转为N
def convert_reference(
        reference_file: str,
        env_path, 
        restart: bool
):
    """
    :param reference_file: 原始的参考基因组
    :param env_path: 环境变量
    :param restart: 是否检查文件是否存在，并跳过该步骤
    :return: stdout, stderr, log_out, out_reference_file
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # 转换后的文件名
    out_reference_file = "convert." + os.path.basename(reference_file)

    # ############################# awk #############################
    # 检查fasta
    cmd = '''awk '{if ($1~/>/) {print $0} else {$0=toupper($0); gsub(/[^ATGCNatgcn]/,"N"); print $0}}' ''' + \
          reference_file + \
          " 1>" + \
          out_reference_file

    # 检查文件是否存在
    if restart:
        # 检查文件
        file_size = getsize(
            out_reference_file
        )
        # 如果小于等于0
        if file_size <= 0:
            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "convert_reference.awk", env_path)
    else:  # 如果没有指定restart，直接运行
        # 提交任务
        stdout, stderr, log_out = run_cmd.run(cmd, "convert_reference.awk", env_path)

    # 检查log是否正常，不正常提前退出
    if log_out:
        return stdout, stderr, log_out, out_reference_file

    # ############################# samtools faidx #############################
    # 构建索引
    cmd = "samtools faidx " + out_reference_file

    # 提交任务
    stdout, stderr, log_out = run_cmd.run(cmd, "convert_reference.faidx", env_path)

    # 检查log是否正常，不正常提前退出
    if log_out:
        return stdout, stderr, log_out, out_reference_file

    # 将路径补全
    out_reference_file = os.path.abspath(out_reference_file)

    return stdout, stderr, log_out, out_reference_file


# 对vcf文件进行排序
def bgzip_vcf(
    vcf_file: str,
    env_path, 
    restart: bool
):
    """
    :param vcf_file: vcf文件
    :param env_path: 环境变量
    :param restart: 是否检查文件是否存在，并跳过该步骤
    :return: stdout, stderr, log_out, out_vcf_file
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    out_vcf_file = vcf_file + ".gz"

    cmd = "bgzip -f {} && tabix -f {}".format(vcf_file, out_vcf_file)

    # 检查文件是否存在
    if restart:
        # 如果小于等于0
        if getsize(out_vcf_file) <= 28 and getsize(out_vcf_file + ".tbi") <= 72:
            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "bgzip", env_path)
    else:  # 如果没有指定restart，直接运行
        # 提交任务
        stdout, stderr, log_out = run_cmd.run(cmd, "bgzip", env_path)

    # 将路径补全
    out_vcf_file = os.path.abspath(out_vcf_file)

    return stdout, stderr, log_out, out_vcf_file


# 对fastq文件进行评估
def fastq_count(
    code_dir: str,
    env_path, 
    fastq_file1: str,
    fastq_file2: str = ""
):
    """
    :param code_dir: 代码路径
    :param env_path: 环境变量
    :param fastq_file1: 测序文件1
    :param fastq_file2: 测序文件2
    :return: stdout, stderr, log_out, fastq_base, read_num, read_len
    """
    # 代码的路径
    code_path = os.path.join(code_dir, "fastAQ count")

    # 评估fastq大小
    if fastq_file2:  # 二代双端测序
        cmd = '{} -i {} -i {}'.format(code_path, fastq_file1, fastq_file2)
    else:  # 三代测序数据或二代单端测序
        cmd = '{} -i {}'.format(code_path, fastq_file1)

    # 提交任务
    stdout, stderr, log_out = run_cmd.run(cmd, "fastAQ count", env_path)

    fastq_base = 0
    read_num = 0
    read_len = 0
    for i in stdout.decode().split("\n"):
        if "readBase" in i:
            fastq_base = int(i.strip().split(":")[1])
        if "readNum" in i:
            read_num = int(i.strip().split(":")[1])
        if "readLen" in i:
            read_len = int(i.strip().split(":")[1])

    # 检查结果对不对
    if fastq_base == 0 or read_num == 0 or read_len == 0:
        log = '[EVG.fastAQ count] The fastq file is wrong, please check the parameters.\n'
        logger.error(log)
        raise Exception(log)

    return stdout, stderr, log_out, fastq_base, read_num, read_len


# 对基因组大小进行评估
def fasta_count(
    code_dir: str,
    env_path, 
    reference_file: str
):
    """
    :param code_dir: 代码路径
    :param env_path: 环境变量
    :param reference_file: 参考基因组
    :return: stdout, stderr, log_out, fasta_base
    """
    # 代码的路径
    code_path = os.path.join(code_dir, "fastAQ count")

    # 评估fasta大小
    cmd = '{} -i {}'.format(code_path, reference_file)

    # 提交任务
    stdout, stderr, log_out = run_cmd.run(cmd, "fasta_count", env_path)

    fasta_base = 0
    for i in stdout.decode().split("\n"):
        if "readBase" in i:
            fasta_base = int(i.strip().split(":")[1])

    # 检查结果对不对
    if fasta_base == 0:
        log_out = '[EVG.fasta_count] The fasta file is wrong, please check the parameters.\n'

    return stdout, stderr, log_out, fasta_base


# 对fastq文件进行下采样
def downsample(
        code_dir: str,
        fastq_file1: str,
        fastq_file2: str,
        fastq_base: int,
        fasta_base: int,
        need_depth: float,
        env_path, 
        restart: bool
):
    """
    :param code_dir: 代码路径
    :param fastq_file1: 测序文件1
    :param fastq_file2: 测序文件2
    :param fastq_base: 测序文件碱基数
    :param fasta_base: 参考基因组碱基数
    :param need_depth: 需要的深度
    :param env_path: 环境变量
    :param restart: 是否检查文件是否存在，并跳过该步骤
    :return:
    """
    # 初始化log
    stdout = ""
    stderr = ""
    log_out = ""

    # 代码的路径
    code_path = os.path.join(code_dir, "fastAQ sample")

    # downsampling
    need_base = fasta_base * need_depth
    need_ratio = round(need_base / fastq_base, 3)
    read_depth = fastq_base / fasta_base

    if need_ratio >= 1:  # 测序数据小于设定值，跳过，不进行下采样
        log = '[EVG.fastAQ sample] Insufficient sequencing data ({:.2f}×/{}×), skip downsampling step.\n'. \
                  format(read_depth, need_depth)
        logger.error(log)
        fastq_out_file1 = fastq_file1
        fastq_out_file2 = fastq_file2
    else:
        fastq_out_file1 = "sample." + str(need_ratio) + "." + os.path.basename(fastq_file1)
        fastq_out_file1 = fastq_out_file1.replace(".gz", "").replace(".GZ", "")
        fastq_out_file1 = os.path.abspath(fastq_out_file1)  # 补全路径

        if fastq_file2:  # 二代双端测序
            # 输出文件名并补全路径
            fastq_out_file2 = "sample." + str(need_ratio) + "." + os.path.basename(fastq_file2)
            fastq_out_file2 = fastq_out_file2.replace(".gz", "").replace(".GZ", "")
            fastq_out_file2 = os.path.abspath(fastq_out_file2)  # 补全路径

            # 检查文件是否存在
            if restart:
                # 检查文件
                file_size1 = getsize(
                    fastq_out_file1
                )
                file_size2 = getsize(
                    fastq_out_file2
                )
                # 两个文件都存在，跳过该步骤
                if file_size1 > 0 and file_size2 > 0:
                    return stdout, stderr, log_out, os.path.abspath(fastq_out_file1), os.path.abspath(fastq_out_file2)

            cmd1 = '{} -i {} -f {} 1>{}'.format(code_path, fastq_file1, need_ratio, fastq_out_file1)
            cmd2 = '{} -i {} -f {} 1>{}'.format(code_path, fastq_file2, need_ratio, fastq_out_file2)
            with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
                to_do = []
                future1 = executor.submit(run_cmd.run, cmd1, "fastAQ sample", env_path)
                future2 = executor.submit(run_cmd.run, cmd2, "fastAQ sample", env_path)
                to_do.append(future1)
                to_do.append(future2)

                for future in concurrent.futures.as_completed(to_do):  # 并发执行
                    stdout, stderr, log_out_tmp = future.result()  # 检查返回值
                    if log_out_tmp:
                        log_out = log_out_tmp
        else:  # 三代测序数据或二代单端测序
            cmd = '{} -i {} -f {} 1>{}'.format(code_path, fastq_file1, need_ratio, fastq_out_file1)
            fastq_out_file2 = ""

            # 检查文件是否存在
            if restart:
                # 检查文件
                file_size1 = getsize(
                    fastq_out_file1
                )
                # 文件存在，跳过该步骤
                if file_size1 > 0:
                    return stdout, stderr, log_out, os.path.abspath(fastq_out_file1), ""

            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "fastAQ sample", env_path)

    return stdout, stderr, log_out, os.path.abspath(fastq_out_file1), os.path.abspath(fastq_out_file2)


# vcf转换
def vcf_convert(
    code_dir: str,
    reference_file: str,
    vcf_file: str,
    read_len: int,
    out_name: str,
    mode: str,
    env_path, 
    restart: bool
):
    """
    :param code_dir: 代码路径
    :param reference_file: 参考基因组
    :param vcf_file: vcf文件
    :param read_len: 读长
    :param out_name: 输出vcf文件名
    :param mode: 模式
    :param env_path: 环境变量
    :param restart: 是否检查文件是否存在，并跳过该步骤
    :return: stdout, stderr, log_out, out_name, sv_out_name, region_file
    """

    # log
    stdout = ""
    stderr = ""
    log_out = ""

    # 代码的路径
    code_path = os.path.join(code_dir, "graphvcf convert")

    # ################################### graphvcf convert ###################################
    # graphvcf convert
    cmd = '{} -r {} -v {} -l {} -o {}'.format(
        code_path,
        reference_file,
        vcf_file,
        read_len,
        out_name)

    # 检查文件是否存在
    if restart:
        # 如果小于等于0
        if getsize(out_name) <= 0 and getsize(out_name + ".gz") <= 0:  # vcf或者bgzip压缩后的文件存在则跳过
            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "graphvcf convert", env_path)
    else:  # 如果没有指定restart，直接运行
        # 提交任务
        stdout, stderr, log_out = run_cmd.run(cmd, "graphvcf convert", env_path)

    # 如果退出代码有问题，则报错
    if log_out:
        return stdout, stderr, log_out, out_name, "", ""

    # 路径补全
    region_file = os.path.abspath("CHROMOSOME.NAME")

    # ################################### sort ###################################
    # vcf排序
    cmd = '''grep '#' {} > {} && grep -v '#' {} | sort --parallel=10 -k 1,1 -k 2,2n -t $'\t' | '''.format(
        out_name,
        out_name + ".sort",
        out_name
    )
    cmd += '''awk -F "\t" '!a[$1,$2]++' 1>>{} && mv {} {}'''.format(
        out_name+".sort",
        out_name+".sort",
        out_name
    )

    # 检查文件是否存在
    if restart:
        # 如果小于等于0
        if getsize(out_name) > 0:  # vcf或者bgzip压缩后的文件存在则跳过
            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "sort", env_path)
    else:  # 如果没有指定restart，直接运行
        # 提交任务
        stdout, stderr, log_out = run_cmd.run(cmd, "sort", env_path)

    # 如果退出代码有问题，则报错
    if log_out:
        return stdout, stderr, log_out, out_name, "", region_file

    # vcf按类别进行划分
    if mode == "fast":
        sv_out_name = "sv." + out_name
        cmd = ""
        if getsize(out_name) > 0:
            cmd = '''cat ''' + out_name + ''' | awk -F '\t' 'BEGIN{FS=OFS="\t"} {if($0~/#/) 
            {print $0} else if(length($4)<50 && length($5)<50) {pass} else {print $0}}' > ''' + sv_out_name
        elif getsize(out_name + ".gz") > 28:
            cmd = '''zcat ''' + out_name + ".gz" + ''' | awk -F '\t' 'BEGIN{FS=OFS="\t"} {if($0~/#/) 
                {print $0} else if(length($4)<50 && length($5)<50) {pass} else {print $0}}' > ''' + sv_out_name
        # 检查文件是否存在
        if restart:
            # 检查文件
            file_size = getsize(
                sv_out_name
            )
            # 如果小于等于0
            if file_size <= 0:
                # 提交任务
                stdout, stderr, log_out = run_cmd.run(cmd, "split", env_path)
        else:  # 如果没有指定restart，直接运行
            # 提交任务
            stdout, stderr, log_out = run_cmd.run(cmd, "split", env_path)
    else:
        sv_out_name = out_name

    out_name = os.path.abspath(out_name)
    sv_out_name = os.path.abspath(sv_out_name)

    return stdout, stderr, log_out, out_name, sv_out_name, region_file
