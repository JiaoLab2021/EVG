#!/usr/bin/python3
# coding=gb2312

import subprocess
import tempfile
import os
import logging


# log
logger = logging.getLogger('SynDiv')
formatter = logging.Formatter('[%(asctime)s] %(message)s')
handler = logging.StreamHandler()  # 输出到控制台
handler.setFormatter(formatter)
logger.addHandler(handler)


# run
def run(command, subcommands, envPath):
    # 创建输出文件
    stdout_file = tempfile.NamedTemporaryFile(delete=False)
    stderr_file = tempfile.NamedTemporaryFile(delete=False)

    # 提交任务
    proc = subprocess.Popen(command, shell=True, stdout=stdout_file, stderr=stderr_file, env=envPath)

    log = f'[EVG.{subcommands}] CMD: {command}\n'
    logger.error(log)

    # 重置log用于判断是否正常退出
    log = ""

    # 等待命令执行完成
    proc.wait()

    # 读取输出文件并关闭文件句柄
    with open(stdout_file.name, 'rb') as f:
        stdout_data = f.read()
    os.unlink(stdout_file.name)

    with open(stderr_file.name, 'rb') as f:
        stderr_data = f.read()
    os.unlink(stderr_file.name)

    exit_state = proc.returncode
    stdout = stdout_data.strip().decode('utf-8')
    stderr = stderr_data.strip().decode('utf-8')

    # 标准输出和错误输出
    if exit_state != 0 or "FileNotFoundError" in stderr or \
            "command not found" in stderr or \
            "error" in stderr or \
            "Error" in stderr:
        log = f'[EVG.{subcommands}] Status: {exit_state} -> {stderr}.'

    return stdout, stderr, log
