"""Merge of one sample's genotype calls with graphvcf merge."""

import os

from evg_pipeline.common import CODE_DIR, run_cmd


def merge_sample(vcf, sample_name, calls, threads, out_dir):
    """Merge one sample's (software, VCF) calls against the converted VCF and index the results with tabix; return <sample>.vcf.gz."""
    merge_vcf_file = os.path.join(out_dir, f"{sample_name}.vcf.gz")
    high_vcf_file = f"h.{sample_name}.vcf.gz"
    cmd = f"rm -f {high_vcf_file} {high_vcf_file}.tbi && {os.path.join(CODE_DIR, 'graphvcf')} merge -t {threads} -v {vcf} "
    for software, call_file in calls:
        cmd += f"--{software} {call_file} "
    cmd += f"-n {sample_name} -o {merge_vcf_file} && tabix -f -p vcf {merge_vcf_file}" \
           f" && if [ -e {high_vcf_file} ]; then tabix -f -p vcf {high_vcf_file}; fi"
    run_cmd(cmd, out_dir, files=[merge_vcf_file + ".tbi", f"{high_vcf_file}*"])
    return merge_vcf_file
