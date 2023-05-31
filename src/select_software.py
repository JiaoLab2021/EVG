# -*- coding: utf-8 -*-
#!/usr/bin/python3

# select_software
def main(real_depth, read_len, fasta_base):
    # Selected software list
    select_software_list = []

    # select program
    if read_len >= 10000:  # When the read length is greater than 10kb, run GraphAligner for three generations of data
        select_software_list.append("GraphAligner")
    else:
        select_software_list.append("Paragraph")
        select_software_list.append("BayesTyper")

        if fasta_base > 200000000:  # If the genome is larger than 200Mb, use giraffe to run
            select_software_list.append("VG-Giraffe")
        else:  # Otherwise run with vg_map
            select_software_list.append("VG-MAP")

        if read_len > 130 and real_depth > 5:  # When the read length is greater than 120bp, and the sequencing depth is greater than 5×, then choose GraphTyper2
            select_software_list.append("GraphTyper2")

    return select_software_list
