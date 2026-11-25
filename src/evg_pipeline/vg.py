"""VG-MAP and VG-Giraffe: vg autoindex, vg map / vg giraffe, vg pack and vg call."""

import logging
import os

from evg_pipeline.common import get_version, parse_version, run_cmd

logger = logging.getLogger("EVG")


def get_giraffe_files(index_dir, cwd):
    """Return the minimizer and zipcodes paths written by this vg version (zipcodes is "" before vg 1.63.0; newer names if unknown)."""
    version = parse_version(get_version("vg", cwd))
    if not version:
        logger.warning("Could not read the vg version; using the index file names of vg 1.63.0 and later.")
    if not version or version >= (1, 63, 0):
        return os.path.join(index_dir, "out.shortread.withzip.min"), os.path.join(index_dir, "out.shortread.zipcodes")
    return os.path.join(index_dir, "out.min"), ""


def build_vg_index(reference, vcf, workflow, threads, index_dir):
    """Build the vg map or vg giraffe index (vg autoindex, which reuses existing out.* files, so they are removed first) and its snarls."""
    temp_dir = os.path.join(index_dir, "temp")
    prefix = os.path.join(index_dir, "out")
    cmd = f"rm -rf {temp_dir} {prefix}.* && mkdir -p {temp_dir} && vg autoindex -T {temp_dir} -t {threads} -R XG --workflow {workflow} -r {reference} -v {vcf} -p {prefix} && vg snarls -t 10 {prefix + '.xg'} > {prefix + '.snarls'} && rm -rf {temp_dir}"
    run_cmd(cmd, index_dir, files=["out.*"])


def call_vg(sample_name, min_mapq, threads, index_dir, out_dir):
    """Genotype one sample's GAM with vg pack and vg call; return the VCF path."""
    xg_file = os.path.join(index_dir, "out.xg")
    snarls_file = os.path.join(index_dir, "out.snarls")
    gam_file = os.path.join(out_dir, f"{sample_name}.gam")
    pack_file = os.path.join(out_dir, f"{sample_name}.pack")
    vcf_file = os.path.join(out_dir, f"{sample_name}.vcf")
    run_cmd(f"vg pack -t {threads} -Q {min_mapq} -x {xg_file} -g {gam_file} -o {pack_file}", out_dir)
    run_cmd(f"vg call -t {threads} -s {sample_name} {xg_file} -k {pack_file} -r {snarls_file} 1>{vcf_file}", out_dir)
    return vcf_file


def genotype_vg(software, sample_name, read1, read2, depth, threads, index_dir, out_dir):
    """Align one sample with vg map (VG-MAP) or vg giraffe (VG-Giraffe) and genotype it; return the VCF path."""
    xg_file = os.path.join(index_dir, "out.xg")
    gam_file = os.path.join(out_dir, f"{sample_name}.gam")
    if software == "VG-MAP":
        cmd = f"vg map -t {threads} -g {os.path.join(index_dir, 'out.gcsa')} -x {xg_file} -f {read1}"
        if read2:
            cmd += f" -f {read2}"
    else:
        minimizer_file, zipcodes_file = get_giraffe_files(index_dir, out_dir)
        cmd = f"vg giraffe -t {threads} -x {xg_file} -Z {os.path.join(index_dir, 'out.giraffe.gbz')} -m {minimizer_file} -d {os.path.join(index_dir, 'out.dist')} -f {read1}"
        if read2:
            cmd += f" -f {read2}"
        if zipcodes_file:
            cmd += f" -z {zipcodes_file}"
    run_cmd(cmd + f" 1>{gam_file}", out_dir, files=[os.path.join(index_dir, "out.gcsa.lcp")] if software == "VG-MAP" else [])
    return call_vg(sample_name, min(3, int(depth / 2)), threads, index_dir, out_dir)
