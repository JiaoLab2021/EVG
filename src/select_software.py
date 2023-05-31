#!/usr/bin/python3
# coding=gb2312

# select_software
def main(real_depth, read_len, fasta_base):
    # 选择的软件列表
    select_software_list = []

    # 选择程序
    if read_len >= 10000:  # 读长大于10kb的时候，为三代数据，跑GraphAligner
        select_software_list.append("GraphAligner")
    else:
        select_software_list.append("Paragraph")
        select_software_list.append("BayesTyper")

        if fasta_base > 200000000:  # 如果基因组大于200Mb，用giraffe跑
            select_software_list.append("VG-Giraffe")
        else:  # 否则用vg_map跑
            select_software_list.append("VG-MAP")

        if read_len > 130 and real_depth > 5:  # 读长大于120bp，且测序深度大于5×的时候，再选择GraphTyper2
            select_software_list.append("GraphTyper2")

    return select_software_list
