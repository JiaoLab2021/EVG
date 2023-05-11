#!/usr/bin/python3
# Created on 2023/3/10
# @author: du
# email: dzz0539@163.com

import datetime
import sys
import os


# 判断文件是否存在，若存在是否为空，返回为文件的大小 -1->不存在  0->大小为0  1->跳过该步骤
def getsize(
        file_path: str
):
    """
    :param file_path:    需要判断的文件路径
    :return: filesize    -1->不存在  0->大小为0  1->跳过该步骤
    """
    if os.path.exists(file_path):
        size = os.path.getsize(file_path)
        if not size:
            return 0
        else:
            log = "[" + str(datetime.datetime.now()).split(".")[0] + \
                  '] [EVG.getsize] file exists, skip: {}\n'.format(file_path)
            sys.stderr.write(log)
            return size
    else:
        return -1
