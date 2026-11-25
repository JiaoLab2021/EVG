"""BWA index and read mapping (inputs of GraphTyper2, BayesTyper and Paragraph)."""

import os

from evg_pipeline.common import run_cmd

INDEX_SUFFIXES = (".amb", ".ann", ".bwt", ".pac", ".sa")


def build_bwa_index(reference, work_dir):
    """Build the BWA index next to the converted reference."""
    run_cmd(f"bwa index {reference} -p {reference}", work_dir, files=[reference + suffix for suffix in INDEX_SUFFIXES])


def map_reads(sample_name, read1, read2, reference, threads, bwa_dir):
    """Map one sample with bwa mem, sort and index the BAM; return the BAM path."""
    bam_file = os.path.join(bwa_dir, sample_name + ".bam")
    cmd = f"bwa mem -R '@RG\\tID:foo\\tSM:{sample_name}\\tLB:library1' -t {threads} {reference} {read1} {read2} | samtools view -b -S | samtools sort -@ {threads} -o {bam_file} && samtools index {bam_file}"
    run_cmd(cmd, bwa_dir, files=[bam_file + ".bai"] + [reference + suffix for suffix in INDEX_SUFFIXES])
    return bam_file
