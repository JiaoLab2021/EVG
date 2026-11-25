"""GraphAligner: vg construct/index, GraphAligner alignment, vg pack and vg call."""

import os

from evg_pipeline.common import run_cmd
from evg_pipeline.vg import call_vg


def build_graphaligner_index(reference, vcf, threads, index_dir):
    """Build the vg graph, XG index and snarls used by GraphAligner."""
    temp_dir = os.path.join(index_dir, "temp")
    prefix = os.path.join(index_dir, "out")
    cmd = f"vg construct -t {threads} -r {reference} -v {vcf} 1>{prefix + '.vg'} && mkdir -p {temp_dir} && vg index -t {threads} -b {temp_dir} -x {prefix + '.xg'} {prefix + '.vg'} && rm -rf {temp_dir} && vg snarls -t {threads} {prefix + '.xg'} 1>{prefix + '.snarls'}"
    run_cmd(cmd, index_dir)


def genotype_graphaligner(sample_name, read, depth, threads, index_dir, out_dir):
    """Align one sample with GraphAligner and genotype it with vg; return the VCF path."""
    gam_file = os.path.join(out_dir, f"{sample_name}.gam")
    run_cmd(f"GraphAligner -t {threads} -g {os.path.join(index_dir, 'out.vg')} -x vg -f {read} -a {gam_file}", out_dir)
    return call_vg(sample_name, min(0, int(depth / 2)), threads, index_dir, out_dir)
