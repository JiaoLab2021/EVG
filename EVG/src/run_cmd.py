#!/usr/bin/python3
# Created on 2022/11/7
# @author: du
# email: dzz0539@163.com

import subprocess
import datetime
import sys
import logging


# run
def run(cmd, subcommands):
    # log
    logger = logging.getLogger('run_cmd')

    # 提交任务
    cmd_out = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    log = "[" + str(datetime.datetime.now()).split(".")[0] + \
         '] [EVG.{}] CMD: {}\n'.format(subcommands, cmd)
    logger.error(log)

    # 标准输出和错误输出
    stdout, stderr = cmd_out.communicate()

    # 重置log用于判断是否正常退出
    log = ""

    # 检查程序退出状态并打印log
    exit_state = cmd_out.returncode
    if exit_state != 0 or "FileNotFoundError" in str(stderr) or \
            "command not found" in str(stderr) or \
            "error" in str(stderr) or \
            "Error" in str(stderr):
        log = "[" + str(datetime.datetime.now()).split(".")[0] + \
              '] [EVG.{}] Status: {} -> {}.\n'.format(subcommands, exit_state, stderr.decode().strip())

    return stdout, stderr, log
