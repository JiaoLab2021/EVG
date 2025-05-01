# EVG

<!-- [![GitHub Downloads](https://img.shields.io/github/downloads/JiaoLab2021/EVG/total.svg?style=social&logo=github&label=Download)](https://github.com/JiaoLab2021/EVG/releases) -->
<!-- [![BioConda Install](https://img.shields.io/conda/dn/duzezhen/evg.svg?style=flag&label=BioConda%20install)](https://anaconda.org/DuZeZhen/evg) -->
[![GitHub last commit](https://img.shields.io/github/last-commit/JiaoLab2021/evg.svg?label=Last%20commit&logo=github&style=flat)](https://github.com/JiaoLab2021/EVG/releases)
[![Build Status](https://github.com/JiaoLab2021/EVG/actions/workflows/ci.yaml/badge.svg)](https://github.com/JiaoLab2021/EVG/actions)

## Introduction

A comprehensive benchmark of graph-based genetic variant genotyping algorithms on plant genomes for creating an accurate ensemble pipeline

![pipeline.jpg](fig/pipeline.jpg)

## Requirements

[tabix_url]: https://github.com/samtools/htslib
[bwa_url]: https://github.com/lh3/bwa
[samtools_url]: https://github.com/samtools/samtools
[VG_url]: https://github.com/vgteam/vg
[GraphAligner_url]: https://github.com/maickrau/GraphAligner
[Paragraph_url]: https://github.com/Illumina/paragraph
[BayesTyper_url]: https://github.com/bioinformatics-centre/BayesTyper
[GraphTyper2_url]: https://github.com/DecodeGenetics/graphtyper
[PanGenie_url]: https://github.com/eblerjana/pangenie

Please note the following requirements before building and running the software:

* `Linux` operating system
* cmake version `3.12` or higher
* Python version `3.9`
* C++ compiler that supports `C++17` or higher, and the `zlib` library installed (we recommend using GCC version `"7.3.0"` or newer) for building `graphvcf` and `fastAQ`
* The following dependencies must also be installed: [tabix][tabix_url], [bwa][bwa_url], [samtools][samtools_url], [VG][VG_url], [GraphAligner][GraphAligner_url], [Paragraph][Paragraph_url], [BayesTyper][BayesTyper_url], [GraphTyper2][GraphTyper2_url], [PanGenie][PanGenie_url]

## Recent major updates:

(2025/04/30, v1.2.1)

* Updated Giraffe indexing and alignment commands for vg ≥1.63.0.
* Pinned BayesTyper to 1.5=h176a8bc_0 due to bugs in newer conda versions.

(2024/06/25, v1.2.0)

* If a sample's genotype information is missing in the VCF file, the previous version would throw a segmentation fault. In version `v1.2.0`, it will be replaced with `0|0`.

## Installation

**Install via Anaconda**

The easiest way to install EVG is through Anaconda, but please note that in this case, the Python version must be `3.9`. Conda will automatically set the Python version for you, so please ensure that your system can install Python `3.9`.

```shell
# Create a new environment named evg_env
conda create -n evg_env
# Activate the environment
conda activate evg_env
# Install EVG with all dependencies
conda install -c bioconda -c conda-forge -c kdm801 -c duzezhen evg
```

**Building on Linux**

Use the following script to build the software:

1. First, obtain the source code.

```shell
git clone https://github.com/JiaoLab2021/EVG.git
cd EVG
```

2. Next, compile the software and add the current directory to your system's `PATH` environment variable. Please make sure that `EVG`, `graphvcf`, and `fastAQ` are all in the same folder, as `EVG` will call these two programs from its own directory.

```shell
cmake ./
make
chmod +x EVG.py
ln -sf EVG.py EVG
echo 'export PATH="$PATH:'$(pwd)'"' >> ~/.bashrc
source ~/.bashrc
```

3. Assuming that you have installed all the required software dependencies, please make sure they have been added to your environment path or activated in the corresponding code environment. If you haven't installed them yet, you can use the following code to install all the dependencies:

```shell
# Create a new environment named evg_env
conda create -n evg_env
# Activate the environment
conda activate evg_env
# Install software using conda
conda install -c bioconda -c conda-forge -c kdm801 tabix bwa samtools vg graphaligner paragraph bayestyper graphtyper kmc pangenie
# "ModuleNotFoundError: No module named 'pysam.bcftools'", therefore it is recommended to upgrade pysam in this case
conda update pysam
```

**Note**

The default version of `PanGenie` installed by conda is `2.1.0`, but `EVG` requires version `3.0` or higher. If you choose `PanGenie` as your downstream tool, please remove the current `PanGenie` from your conda environment and manually install the latest version of `PanGenie`, then add it to your environment variables.

**Test**

To verify that the software has been installed correctly, perform a test run using the following steps:

```shell
EVG -h
graphvcf -h
fastAQ -h
tabix -h
bwa
samtools
vg -h
GraphAligner -h
paragraph -h
bayesTyper -h
graphtyper -h
PanGenie -h
kmc -h
jellyfish -h
# test
cd test
EVG -r test.fa -v test.vcf.gz -s sample.txt --software VG-MAP VG-Giraffe GraphAligner Paragraph BayesTyper GraphTyper2 PanGenie &>log.txt &
```

## Usage

**Input Files**

* Reference Genome
* VCF File of Population Variants
* Sample File:

```shell
# Sample File
sample1 sample1.r1.fq.gz sample1.r2.fq.gz
sample2 sample2.r1.fq.gz sample2.r2.fq.gz
...
sampleN sampleN.r1.fq.gz sampleN.r2.fq.gz
```

Please note that the Sample file must be formatted exactly as shown above, where each sample is listed with its corresponding read files.

**Running**

For convenience, let's assume the following file names for the input:

* `refgenome.fa`
* `input.vcf.gz`
* `sample.txt`

`EVG` automatically selects suitable software based on the genome, mutation and sequencing data. If desired, users can also use the `"--software"` command to specify their preferred software. The default running command is as follows:

```shell
EVG -r refgenome.fa -v input.vcf.gz -s sample.txt
```

The results are stored in the `merge/` folder, and each file is named after the corresponding sample listed in `sample.txt`: `sample1.vcf.gz`, `sample2.vcf.gz`, ..., `sampleN.vcf.gz`.

```shell
$ tree merge/
merge/
├── test1.vcf.gz
└── test2.vcf.gz

0 directories, 2 files
```

**Parameter**

* `--depth`: This parameter specifies the maximum sequencing data depth allowed for downstream analysis. If this value is exceeded, EVG will randomly downsample reads to the specified level in order to speed up the run. The default downsampling level is set at 15×, but it can be adjusted to meet specific requirements.
* `--mode`: This parameter determines the operating mode of `EVG`. In fast mode, only certain software is utilized to genotype SNPs and indels, while precise mode employs all software to genotype all variants.
* `--force`: If there are pre-existing files in the running directory of `EVG`, this parameter can be used to forcibly empty the folder. Otherwise, the software will encounter an error and exit.
* `--restart`: This parameter allows the software to resume from where it left off if it unexpectedly stops, enabling a breakpoint restart. Note that software completion is determined by file existence. It's recommended to manually check for incomplete or empty files before using this parameter and delete them.

**graphvcf**

If you already have results from different genotyping software and do not need to use EVG, you can directly use `graphvcf` to merge your results.

```shell
graphvcf merge -v merged.vcf.gz --Paragraph xx.vcf.gz --BayesTyper xx.vcf.gz --VG-Giraffe xx.vcf.gz -n sample1 -o sample.vcf.gz
```

[graphvcf_url]: https://github.com/JiaoLab2021/EVG/wiki/graphvcf-usage

Detailed instructions for using `graphvcf` can be found on the [Wiki page][graphvcf_url].

## Citation

[evg_article]: https://doi.org/10.1186/s13059-024-03239-1
[bwa_article]: https://arxiv.org/abs/1303.3997v2
[vg-map_article]: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1941-7
[vg-giraffe-article]: https://www.science.org/doi/10.1126/science.abg8871
[GraphAligner_article]: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02157-2
[Paragraph_article]: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1909-7#citeas
[BayesTyper_article]: https://www.nature.com/articles/s41588-018-0145-5
[GraphTyper2_article]: https://www.nature.com/articles/s41467-019-13341-9
[samtools_article]: https://academic.oup.com/gigascience/article/10/2/giab008/6137722?login=false
[PanGenie_article]: https://www.nature.com/articles/s41588-022-01043-w

When using the following tools, please cite the corresponding articles:

*  `EVG`:

    *  Du, ZZ., He, JB. & Jiao, WB. [A comprehensive benchmark of graph-based genetic variant genotyping algorithms on plant genomes for creating an accurate ensemble pipeline.][evg_article] Genome Biol 25, 91 (2024).

*  `vg map`: 

    *  Hickey, G., Heller, D., Monlong, J. et al. [Genotyping structural variants in pangenome graphs using the vg toolkit.][vg-map_article] Genome Biol 21, 35 (2020).

*  `vg giraffe`: 

    *  Jouni Sirén et al. [Pangenomics enables genotyping of known structural variants in 5202 diverse genomes.][vg-giraffe-article] Science 374, abg8871 (2021).

*  `GraphAligner`: 

    *  Rautiainen, M., Marschall, T. [GraphAligner: rapid and versatile sequence-to-graph alignment.][GraphAligner_article] Genome Biol 21, 253 (2020).

*  `Paragraph`: 

    *  Chen, S., Krusche, P., Dolzhenko, E. et al. [Paragraph: a graph-based structural variant genotyper for short-read sequence data.][Paragraph_article] Genome Biol 20, 291 (2019).

    *  Li, H. [Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM.][bwa_article] arXiv: Genomics. (2013).

*  `BayesTyper`: 

    *  Sibbesen, J.A., Maretty, L., [The Danish Pan-Genome Consortium. et al. Accurate genotyping across variant classes and lengths using variant graphs.][BayesTyper_article] Nat Genet 50, 1054–1059 (2018).

    *  Li, H. [Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM.][bwa_article] arXiv: Genomics. (2013).

*  `GraphTyper2`: 

    *  Eggertsson, H.P., Kristmundsdottir, S., Beyter, D. et al. [GraphTyper2 enables population-scale genotyping of structural variation using pangenome graphs.][GraphTyper2_article] Nat Commun 10, 5402 (2019).

    *  Li, H. [Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM.][bwa_article] arXiv: Genomics. (2013).

    *  Danecek, P., Bonfield, J. K., Liddle, J. et al. [Twelve years of SAMtools and BCFtools.][samtools_article] GigaScience, Volume 10, Issue 2, February 2021, giab008

*  `PanGenie`: 

    *  Ebler, J., Ebert, P., Clarke, W.E. et al. [Pangenome-based genome inference allows efficient and accurate genotyping across a wide spectrum of variant classes.][PanGenie_article] Nat Genet 54, 518–525 (2022).

## License

MIT