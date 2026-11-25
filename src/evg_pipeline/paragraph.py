"""Paragraph: per-batch sample lists and multigrmpy.py genotyping."""

import os

from evg_pipeline.common import run_cmd, write_list


def write_paragraph_lists(batches, reads, bams, genotype_dir):
    """Write one Paragraph sample list (id, BAM, depth, read length) per batch; return the list paths."""
    list_files = []
    for number, batch in enumerate(batches):
        list_file = os.path.join(genotype_dir, f"bam_for_Paragraph_{number}.txt")
        text = "id\tpath\tdepth\tread length\n"
        for sample_name in batch:
            text += sample_name + "\t" + bams[sample_name] + "\t" + str(reads[sample_name]["depth"]) + "\t" \
                + str(reads[sample_name]["length"]) + "\n"
        write_list(list_file, text)
        list_files.append(list_file)
    return list_files


def genotype_paragraph(reference, vcf, list_file, bam_files, threads, out_dir):
    """Genotype one batch of samples (BAMs listed in list_file) with multigrmpy.py; return genotypes.vcf.gz."""
    cmd = '''tmpn=`mktemp -u paragraph_XXXXX` && tmpd="/tmp/${tmpn}" && mkdir ${tmpd} && echo ${tmpd} && export TMP=${tmpd} &&''' \
        + f" multigrmpy.py -i {vcf} -m {list_file} -r {reference} --threads {threads}" + ''' --scratch-dir ${tmpd}''' \
        + f" -M 1000 -o {out_dir}" + " && rm -r ${tmpd}"
    vcf_out_file = os.path.join(out_dir, "genotypes.vcf.gz")
    run_cmd(cmd, out_dir, files=[vcf_out_file] + bam_files + [bam + ".bai" for bam in bam_files])
    return vcf_out_file
