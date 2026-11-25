"""Input preparation: read and genome counts, downsampling, reference and VCF conversion."""

import logging
import os

from evg_pipeline.common import CODE_DIR, run_cmd, run_jobs

logger = logging.getLogger("EVG")


def count_fasta(fasta, threads, cwd):
    """Return the number of bases in a FASTA file (fastAQ count)."""
    stdout = run_cmd(f"{os.path.join(CODE_DIR, 'fastAQ')} count -t {threads} -i {fasta}", cwd)
    fasta_base = 0
    for line in stdout.split("\n"):
        if "readBase" in line:
            fasta_base = int(line.strip().split(":")[1])
    if fasta_base == 0:
        raise RuntimeError(f"The size of the reference genome is 0. Please check the file: {fasta}.")
    logger.info(f"  - {os.path.basename(fasta)}: {fasta_base} bp")
    return fasta_base


def count_fastq(read1, read2, threads, cwd):
    """Return the number of bases and the read length of one sample (fastAQ count)."""
    fastaq = os.path.join(CODE_DIR, "fastAQ")
    if read2:
        cmd = f"{fastaq} count -t {threads} -i {read1} -i {read2}"
    else:
        cmd = f"{fastaq} count -t {threads} -i {read1}"
    stdout = run_cmd(cmd, cwd)
    fastq_base = 0
    read_len = 0
    for line in stdout.split("\n"):
        if "readBase" in line:
            fastq_base = int(line.strip().split(":")[1])
        if "readLen" in line:
            read_len = int(line.strip().split(":")[1])
    if fastq_base == 0 or read_len == 0:
        logger.warning(f"The fastq file is wrong, please check the parameters: {read1} {read2}")
    return fastq_base, read_len


def sample_reads(cmd, read_out, out_dir):
    """Run one fastAQ sample command that writes read_out (the two mates of a pair run in parallel)."""
    run_cmd(cmd, out_dir, files=[read_out])


def downsample_reads(sample_name, read1, read2, ref_base, depth, threads, out_dir):
    """Downsample one sample to the target depth (fastAQ sample, BGZF output) into reads/<sample>/; return its reads, depth and read length."""
    fastaq = os.path.join(CODE_DIR, "fastAQ")
    read_base, read_len = count_fastq(read1, read2, threads, out_dir)
    if read_base == 0:
        raise RuntimeError(f"Sequencing depth is zero. Please verify the input files: {read1} {read2}")
    need_ratio = round(ref_base * depth / read_base, 3)
    read_depth = read_base / ref_base
    if need_ratio == 0:
        raise RuntimeError(f"Cannot downsample {sample_name} ({read_depth:.2f}x) to {depth:g}x: the fraction rounds to 0.")

    # ------ Downsampling ------
    if need_ratio >= 1:
        logger.info(f"  - {sample_name}: {read_depth:.2f}x, read length {read_len} bp, not downsampled (below {depth:g}x)")
        read1_out = read1
        read2_out = read2
    else:
        read1_out = f"sample.{need_ratio}.R1.{os.path.basename(read1)}"
        if not read1_out.endswith((".gz", ".GZ")):
            read1_out += ".gz"
        sample_dir = os.path.join(out_dir, "reads", sample_name)
        os.makedirs(sample_dir, exist_ok=True)
        read1_out = os.path.join(sample_dir, read1_out)
        if read2:
            read2_out = f"sample.{need_ratio}.R2.{os.path.basename(read2)}"
            if not read2_out.endswith((".gz", ".GZ")):
                read2_out += ".gz"
            read2_out = os.path.join(sample_dir, read2_out)
            cmd1 = f"{fastaq} sample -t {max(1, threads // 2)} -i {read1} -f {need_ratio} -o {read1_out}"
            cmd2 = f"{fastaq} sample -t {max(1, threads - threads // 2)} -i {read2} -f {need_ratio} -o {read2_out}"
            run_jobs([(sample_reads, cmd1, read1_out, out_dir), (sample_reads, cmd2, read2_out, out_dir)],
                     2 if threads > 1 else 1, progress=False)
        else:
            read2_out = ""
            sample_reads(f"{fastaq} sample -t {threads} -i {read1} -f {need_ratio} -o {read1_out}", read1_out, out_dir)
        logger.info(f"  - {sample_name}: {read_depth:.2f}x, read length {read_len} bp, downsampled to {depth:g}x (fraction {need_ratio})")

    # ------ Sequencing depth ------
    sample_depth = min(read_base / ref_base, depth)
    if sample_depth == 0:
        raise RuntimeError(f"Sequencing depth is zero. Please verify the input files: {read1} {read2}")
    return {"read1": read1_out, "read2": read2_out, "depth": sample_depth, "length": read_len}


def convert_reference(reference, out_dir):
    """Write convert.<reference> with non-ATGCN bases set to N, index it with samtools faidx; return its path."""
    out_name = "convert." + os.path.basename(reference)
    cmd = '''awk '{if ($1~/>/) {print $0} else {$0=toupper($0); gsub(/[^ATGCNatgcn]/,"N"); print $0}}' ''' \
        + reference + " 1>" + out_name
    run_cmd(cmd, out_dir)
    run_cmd(f"samtools faidx {out_name}", out_dir, files=[out_name + ".fai"])
    return os.path.join(out_dir, out_name)


def convert_vcf(reference, vcf, mode, threads, out_dir):
    """Convert the VCF (graphvcf convert, BGZF output), keep SVs only in fast mode, index both; return the vcf, sv and region paths."""
    vcf_name = f"convert.{os.path.basename(vcf.replace('.gz', ''))}"
    vcf_gz = os.path.join(out_dir, vcf_name + ".gz")
    region_file = os.path.join(out_dir, "CHROMOSOME.NAME")

    # ------ graphvcf convert and tabix ------
    cmd = f"{os.path.join(CODE_DIR, 'graphvcf')} convert -t {threads} -r {reference} -v {vcf} -l 350 -o {vcf_name}.gz" \
          f" && tabix -f {vcf_name}.gz"
    run_cmd(cmd, out_dir, files=[region_file, vcf_gz, vcf_gz + ".tbi"])
    if mode != "fast":
        return {"vcf": vcf_gz, "sv": vcf_gz, "region": region_file}

    # ------ SV subset (fast mode) ------
    sv_name = "sv." + vcf_name
    cmd = f"zcat {vcf_name}.gz" + ''' | awk -F '\t' 'BEGIN{FS=OFS="\t"} {if($0~/#/) {print $0} else if(length($4)<50 && length($5)<50) {pass} else {print $0}}' > ''' + sv_name \
        + f" && bgzip -@ {threads} -f {sv_name} && tabix -f {sv_name}.gz"
    run_cmd(cmd, out_dir, files=[sv_name + ".gz.tbi"])
    return {"vcf": vcf_gz, "sv": os.path.join(out_dir, sv_name + ".gz"), "region": region_file}
