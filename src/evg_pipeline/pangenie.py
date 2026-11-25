"""PanGenie: non-overlapping VCF, PanGenie-index and per-sample genotyping."""

import os

from evg_pipeline.common import run_cmd


def build_pangenie_index(reference, vcf, threads, index_dir):
    """Drop overlapping variants (convert.vcf) and build the PanGenie index."""
    vcf_out_file = os.path.join(index_dir, "convert.vcf")
    awk_cmd = r"""awk 'BEGIN{FS="\t"; OFS="\t"; start=0; end=0} /^#/ {print; next} $1!=chr {chr=$1; start=0; end=0} {if(start==0) {start=$2; end=$2+length($4)-1} else if($2>end) {start=$2; end=$2+length($4)-1} else{next}} {print}'"""
    if vcf.endswith((".gz", ".GZ")):
        run_cmd(f"zcat {vcf} | {awk_cmd} 1>{vcf_out_file}", index_dir)
    else:
        run_cmd(f"cat {vcf} | {awk_cmd} 1>{vcf_out_file}", index_dir)
    run_cmd(f"PanGenie-index -r {reference} -v {vcf_out_file} -t {threads} -o {os.path.join(index_dir, 'out')}", index_dir,
            files=["out*"])


def genotype_pangenie(sample_name, read1, read2, threads, index_dir, out_dir):
    """Concatenate one sample's reads and genotype them with PanGenie; return the VCF path."""
    read_file = os.path.join(out_dir, f"{sample_name}.fq")
    if read1.endswith((".gz", ".GZ")):
        cmd = f"zcat {read1} 1>{read_file}"
    else:
        cmd = f"cat {read1} 1>{read_file}"
    if read2:
        if read2.endswith((".gz", ".GZ")):
            cmd += f" && zcat {read2} 1>>{read_file}"
        else:
            cmd += f" && cat {read2} 1>>{read_file}"
    cmd += f" && PanGenie -s {sample_name} -i {read_file} -f {index_dir}/out -t {threads} -j {threads} -o {os.path.join(out_dir, sample_name)} && rm {read_file}"
    vcf_out_file = os.path.join(out_dir, f"{sample_name}_genotyping.vcf")
    try:
        run_cmd(cmd, out_dir, files=[os.path.join(index_dir, "out*"), vcf_out_file])
    finally:
        if os.path.exists(read_file):
            os.remove(read_file)
    return vcf_out_file
