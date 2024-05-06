#!/usr/bin/env python3

# -*- coding: utf-8 -*-


__data__ = "2024/05/06"
__version__ = "1.1.8"
__author__ = "Zezhen Du"
__email__ = "dzz0539@gmail.com or dzz0539@163.com"

import os
import argparse


# Solve the error when the denominator is 0
def safe_division(numerator, denominator):
    try:
        return str('%.3f' % (numerator / denominator))
    except ZeroDivisionError:
        return "0"

# Calculate all_ratio
def cal_ratio_fun(all_number, recall_number, genotype_recall_number, call_number, type_cal):
    recall_ratio = safe_division(recall_number, all_number)
    recall_precision_ratio = safe_division(recall_number, call_number)
    genotype_recall_ratio = safe_division(genotype_recall_number, all_number)

    if type_cal == "recall":
        genotype_precision_ratio = safe_division(genotype_recall_number, recall_number)
    elif type_cal == "call":
        genotype_precision_ratio = safe_division(genotype_recall_number, call_number)
    else:
        genotype_precision_ratio = "0"

    return recall_ratio, recall_precision_ratio, genotype_recall_ratio, genotype_precision_ratio

def cal_ratio(len_list, ratio_len_list, type_cal, ratio_type_list, file):
    stats = {
        'all': {'all': 0, 'genotype': 0, 'recall': 0, 'call': 0},
        'sv': {'all': 0, 'genotype': 0, 'recall': 0, 'call': 0},
        'del': {'all': 0, 'genotype': 0, 'recall': 0, 'call': 0},
        'indel': {'all': 0, 'genotype': 0, 'recall': 0, 'call': 0},
        'snp': {'all': 0, 'genotype': 0, 'recall': 0, 'call': 0},
        'ins': {'all': 0, 'genotype': 0, 'recall': 0, 'call': 0},
    }

    # Add header
    if not ratio_len_list:
        ratio_len_list = [["recall", "precision", "genotype", "genotype_precision"]]
    else:
        ratio_len_list[0] += ["recall", "precision", "genotype", "genotype_precision"]

    # Loop quantity list calculation ratio
    for i in range(1, len(len_list)):
        # The last four elements: all, genotype, recall, call
        tmp_list = len_list[i][-4:]
        # The total number under this length
        all_number = float(tmp_list[0])
        # Number of software recalls
        genotype_number = float(tmp_list[1])
        recall_number = float(tmp_list[2])
        call_number = float(tmp_list[3])

        # Calculate proportion
        recall_ratio, recall_precision_ratio, genotype_recall_ratio, genotype_precision_ratio = cal_ratio_fun(
            all_number, 
            recall_number,
            genotype_number, 
            call_number,
            type_cal
        )
        if len(ratio_len_list) - 1 < i:
            ratio_len_list.append([recall_ratio, recall_precision_ratio, genotype_recall_ratio, genotype_precision_ratio])
        else:
            ratio_len_list[i] += [recall_ratio, recall_precision_ratio, genotype_recall_ratio, genotype_precision_ratio]

        # Calculate the total
        stats['all']['all'] += all_number
        stats['all']['genotype'] += genotype_number
        stats['all']['recall'] += recall_number
        stats['all']['call'] += call_number
        
        # Find the length of variation to determine the type
        length_left = float(len_list[i][0].strip().split("|")[0])
        length_right = float(len_list[i][0].strip().split("|")[1])
        if length_right < -49:
            stats['del']['all'] += all_number
            stats['del']['genotype'] += genotype_number
            stats['del']['recall'] += recall_number
            stats['del']['call'] += call_number

            stats['sv']['all'] += all_number
            stats['sv']['genotype'] += genotype_number
            stats['sv']['recall'] += recall_number
            stats['sv']['call'] += call_number
        elif -49 <= length_left < 0 and -49 <= length_right < 0:
            stats['indel']['all'] += all_number
            stats['indel']['genotype'] += genotype_number
            stats['indel']['recall'] += recall_number
            stats['indel']['call'] += call_number
        elif length_left == length_right == 0:
            stats['snp']['all'] += all_number
            stats['snp']['genotype'] += genotype_number
            stats['snp']['recall'] += recall_number
            stats['snp']['call'] += call_number
        elif 0 < length_left <= 49 and 0 < length_right <= 49:
            stats['indel']['all'] += all_number
            stats['indel']['genotype'] += genotype_number
            stats['indel']['recall'] += recall_number
            stats['indel']['call'] += call_number
        else:
            stats['ins']['all'] += all_number
            stats['ins']['genotype'] += genotype_number
            stats['ins']['recall'] += recall_number
            stats['ins']['call'] += call_number

            stats['sv']['all'] += all_number
            stats['sv']['genotype'] += genotype_number
            stats['sv']['recall'] += recall_number
            stats['sv']['call'] += call_number

    # Calculate the proportion of different types of mutations and add to the output list
    for type in stats:
        recall_ratio, recall_precision_ratio, genotype_recall_ratio, genotype_precision_ratio = cal_ratio_fun(
            stats[type]['all'], 
            stats[type]['recall'],
            stats[type]['genotype'], 
            stats[type]['call'],
            type_cal
        )
        ratio_type_list.append(f"{type}_{file.split('/')[0]}\t{recall_ratio}\t{recall_precision_ratio}\t{genotype_recall_ratio}\t{genotype_precision_ratio}")

    return ratio_len_list, ratio_type_list

def out_merge(file, ratio_all_list, len_list, type_cal, ratio_len_list, ratio_type_list):
    with open(file, 'r') as f:
        informations = f.read().strip().split("\n")

        # recall
        if len(informations) == 42:
            # snp+indel+sv
            genotype_recall_all = float(informations[1].strip().split(":")[1])
            recall_all = float(informations[3].strip().split(":")[1])
            call_all = float(informations[5].strip().split(":")[1])
            all_all = float(informations[7].strip().split(":")[1])

            # Calculate proportion
            recall_ratio_all, recall_precision_ratio_all, genotype_recall_ratio_all, genotype_precision_ratio_all = cal_ratio_fun(
                all_all, 
                recall_all,
                genotype_recall_all, 
                call_all, 
                type_cal
            )

            # sv
            genotype_recall_sv = float(informations[35].strip().split(":")[1])
            recall_sv = float(informations[37].strip().split(":")[1])
            call_sv = float(informations[39].strip().split(":")[1])
            all_sv = float(informations[41].strip().split(":")[1])

            # Calculate proportion
            recall_ratio_sv, recall_precision_ratio_sv, genotype_recall_ratio_sv, genotype_precision_ratio_sv = cal_ratio_fun(
                all_sv, 
                recall_sv, 
                genotype_recall_sv, 
                call_sv, 
                type_cal
            )

            ratio_all_list.append('\t'.join(["all_" + file.split("/")[0], recall_ratio_all, recall_precision_ratio_all, genotype_recall_ratio_all, genotype_precision_ratio_all]))
            ratio_all_list.append('\t'.join(["sv_" + file.split("/")[0], recall_ratio_sv, recall_precision_ratio_sv, genotype_recall_ratio_sv, genotype_precision_ratio_sv]))

            # len_list.append
            line_num = 0
            for info in informations[9:33]:
                parts1 = info.split(":")[0].strip().replace("/", "|")
                parts2 = info.split(":")[1].strip().split("/")
                new_list = [parts1, parts2[6], parts2[0], parts2[2], parts2[4]]

                if len(len_list) < 24:
                    len_list.append(new_list)
                else:
                    len_list[line_num] += new_list[1:]

                line_num += 1
        else:
            print("The number of lines in the file is incorrect (42): {}".format(str(len(informations))))
            exit(1)

        ratio_len_list, ratio_type_list = cal_ratio(len_list, ratio_len_list, type_cal, ratio_type_list, file)

    return ratio_all_list, len_list, ratio_len_list, ratio_type_list


# save()
def save_result(output_file, ratio_all_list, len_list, out_file_list, ratio_len_list, ratio_type_list):
    with open(output_file, "w") as f:
        out_txt = "\t".join(out_file_list) + "\n" + "\n".join(ratio_all_list) + "\n\n"
        # quantity
        for i in len_list:
            out_txt += "\t".join(i) + "\n"
        # in proportion to length
        out_txt += "\n"
        for i in ratio_len_list:
            out_txt += "\t".join(i) + "\n"
        # proportion by type
        out_txt += "\n" + "\n".join(ratio_type_list) + "\n"

        f.write(out_txt)


def main():
    # log
    print(f"data: {__data__}")
    print(f"version: {__version__}")
    print(f"author: {__author__}")
    print(f"If you encounter any issues related to the code, please don't hesitate to contact us via email at {__email__}.\n")

    # Input parameters
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("input files list")
    required.add_argument("-i", "--input", dest="file_dir", help="Folder of recall results", required=True, type=str, nargs='+')

    parser.add_argument("-o", "--out", dest="output", help="output file name", default='recall.table')
    parser.add_argument("-t", "--type", dest="type_cal", help="Accuracy is calculated based on the number of [call/recall]", default='call', choices=['call', 'recall'])

    args = parser.parse_args()

    out_file_list = args.file_dir
    output_file = args.output

    # overall proportion
    ratio_all_list = ["file\trecall\trecall_precision\tgenotype\tgenotype_precision"]
    len_list = []

    # proportion of length
    ratio_len_list = []

    # Proportion of categories
    ratio_type_list = ["file\trecall\trecall_precision\tgenotype\tgenotype_precision"]

    for i in out_file_list:
        out_file = os.path.join(i, "vcf_evulate.out")
        ratio_all_list, len_list, ratio_len_list, ratio_type_list = out_merge(
            out_file, 
            ratio_all_list, 
            len_list,
            args.type_cal,
            ratio_len_list,
            ratio_type_list
        )

    save_result(output_file, ratio_all_list, len_list, out_file_list, ratio_len_list, ratio_type_list)


if __name__ == '__main__':
    main()
