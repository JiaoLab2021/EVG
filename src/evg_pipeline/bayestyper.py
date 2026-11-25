"""BayesTyper: per-batch sample lists, k-mer counting, clustering, genotyping and unit merge."""

import os

from evg_pipeline.common import run_cmd, write_list


def write_bayestyper_lists(batches, bams, genotype_dir):
    """Write one BayesTyper sample list (name, sex, BAM link) per batch; return the list paths."""
    list_files = []
    for number, batch in enumerate(batches):
        list_file = os.path.join(genotype_dir, f"bam_for_BayesTyper_{number}.tsv")
        batch_dir = os.path.join(genotype_dir, "BayesTyper", str(number))
        text = ""
        for sample_name in batch:
            text += sample_name + "\tF\t" + os.path.join(batch_dir, os.path.basename(bams[sample_name])) + "\n"
        write_list(list_file, text)
        list_files.append(list_file)
    return list_files


def genotype_bayestyper(reference, vcf, list_file, bams, threads, out_dir):
    """Genotype one batch of samples with kmc and BayesTyper; return the merged bayestyper.vcf.gz."""
    # ------ k-mer counting ------
    with open(list_file, "r") as f:
        lines = f.readlines()
    for line in lines:
        fields = line.strip().split()
        sample_name = fields[0]
        bam_link = fields[2]
        run_cmd(f"ln -sf {bams[sample_name]} {bam_link}", out_dir)
        run_cmd(f"kmc -k55 -ci1 -t1 -fbam {bam_link} {bam_link} {out_dir}", out_dir,
                files=[bam_link + ".kmc_pre", bam_link + ".kmc_suf"])
        run_cmd(f"bayesTyperTools makeBloom -k {bam_link} -p {threads}", out_dir,
                files=[bam_link + ".bloomData", bam_link + ".bloomMeta"])

    # ------ Clustering ------
    cluster_dir = os.path.join(out_dir, "bayestyper_cluster_data")
    cmd = f"rm -rf {cluster_dir} {os.path.join(out_dir, 'bayestyper_unit_*')}" \
          f" && bayesTyper cluster -v {vcf} -s {list_file} -g {reference} -p {threads} -o {os.path.join(out_dir, 'bayestyper')}"
    run_cmd(cmd, out_dir, files=[os.path.join(cluster_dir, name) for name in (
        "intercluster_regions.txt.gz", "multigroup_kmers.bloomData", "multigroup_kmers.bloomMeta", "parameter_kmers.fa.gz")]
        + [os.path.join(out_dir, "bayestyper_unit_*", "variant_clusters.bin")])

    # ------ Genotyping ------
    units = sorted(int(folder[len("bayestyper_unit_"):]) for folder in os.listdir(out_dir)
                   if folder.startswith("bayestyper_unit_") and folder[len("bayestyper_unit_"):].isdigit())
    unit_dirs = [os.path.join(out_dir, f"bayestyper_unit_{unit}") for unit in units]
    if not unit_dirs:
        raise RuntimeError("BayesTyper does not generate any folders starting with 'bayestyper_unit_'.")
    for unit_dir in unit_dirs:
        cmd = f"bayesTyper genotype -v {os.path.join(unit_dir, 'variant_clusters.bin')} -c {cluster_dir} -s {list_file} -g {reference} -o {os.path.join(unit_dir, 'bayestyper')} -z -p {threads} --noise-genotyping"
        run_cmd(cmd, out_dir, files=[os.path.join(unit_dir, "bayestyper.vcf.gz")])

    # ------ Merge units ------
    vcf_files = [os.path.join(unit_dir, "bayestyper.vcf.gz") for unit_dir in unit_dirs]
    vcf_out_file = os.path.join(out_dir, "bayestyper.vcf.gz")
    run_cmd(f"zcat {vcf_files[0]} | grep '^#' | gzip -c > {vcf_out_file}"
            f" && zcat {' '.join(vcf_files)} | awk '!/^#/' | sort -k 1,1 -k 2,2n -t $'\t' | gzip -c >> {vcf_out_file}", out_dir)
    return vcf_out_file
