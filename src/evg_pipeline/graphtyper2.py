"""GraphTyper2: BAM list, graphtyper genotype_sv and concatenation of the region VCFs."""

import os

from evg_pipeline.common import run_cmd, write_list


def write_graphtyper_list(bams, genotype_dir):
    """Write the GraphTyper2 BAM list (one BAM per line, samples-file order); return its path."""
    list_file = os.path.join(genotype_dir, "bam_for_GraphTyper2.txt")
    text = ""
    for bam_file in bams.values():
        text += bam_file + "\n"
    write_list(list_file, text)
    return list_file


def genotype_graphtyper2(reference, vcf, list_file, bam_files, region_file, threads, out_dir):
    """Genotype all samples (BAMs listed in list_file) with graphtyper genotype_sv and concatenate the region VCFs; return graphtyper.vcf."""
    vcf_out_file = os.path.join(out_dir, "graphtyper.vcf")
    cmd = f"rm -rf {vcf_out_file} {out_dir}/*/" \
          f" && graphtyper genotype_sv {reference} {vcf} --sams={list_file} --region_file={region_file} --output={out_dir} --threads {threads}" \
          f" && for i in $(ls {out_dir}/*/ | grep 'vcf.gz' | grep -v 'tbi' | head -n1); do zcat {out_dir}/*/$i | grep '#' > {vcf_out_file}; done && zcat {out_dir}/*/*.vcf.gz >> {vcf_out_file}"
    run_cmd(cmd, out_dir, files=bam_files + [bam + ".bai" for bam in bam_files])
    return vcf_out_file
