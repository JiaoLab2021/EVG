#!/usr/bin/python3
# coding=gb2312

import os
import logging


# log
logger = logging.getLogger('SynDiv')
formatter = logging.Formatter('[%(asctime)s] %(message)s')
handler = logging.StreamHandler()  # 输出到控制台
handler.setFormatter(formatter)
logger.addHandler(handler)


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
            log = '[EVG.getsize] file exists, skip: {}\n'.format(file_path)
            logger.error(log)
            return size
    else:
        return -1
