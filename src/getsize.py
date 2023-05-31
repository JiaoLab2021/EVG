# -*- coding: utf-8 -*-
#!/usr/bin/python3

import os
import logging


# log
logger = logging.getLogger('SynDiv')
formatter = logging.Formatter('[%(asctime)s] %(message)s')
handler = logging.StreamHandler()  # output to the console
handler.setFormatter(formatter)
logger.addHandler(handler)


# Determine whether the file exists, if it is empty, return the size of the file. -1->does not exist  0->size is 0  1->skip this step
def getsize(
        file_path: str
):
    """
    :param file_path:    The file path that needs to be judged
    :return: filesize    -1->does not exist  0->size is 0  1->skip this step
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
