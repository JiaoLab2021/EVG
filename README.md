# EVG

## Introduction

A comprehensive benchmark of graph-based genetic variant genotyping algorithms on plant genomes for creating an accurate ensemble pipeline

![pipeline.jpf](./pipeline.jpf)

## Requirements

[VG_url]: https://github.com/vgteam/vg
[GraphAligner_url]: https://github.com/maickrau/GraphAligner
[Paragraph_url]: https://github.com/Illumina/paragraph
[BayesTyper_url]: https://github.com/bioinformatics-centre/BayesTyper
[GraphTyper2_url]: https://github.com/DecodeGenetics/graphtyper
[PanGenie_url]: https://github.com/eblerjana/pangenie

Please note the following requirements before building and running the software:

* `Linux` operating system
* cmake version `3.12` or higher
* Python version `3.8` or higher
* C++ compiler that supports `C++14` or higher, and the `zlib` library installed (we recommend using GCC version `"4.9"` or newer) for building `graphvcf` and `fastAQ`
* The following dependencies must also be installed: [VG][VG_url], [GraphAligner][GraphAligner_url], [Paragraph][Paragraph_url], [BayesTyper][BayesTyper_url], [GraphTyper2][GraphTyper2_url], [PanGenie][PanGenie_url]

## Installation

**Download Releases**

The simplest means of obtaining `EVG` is by downloading the static binary executable.

[Download][url]

**Building on Linux**

If you cannot or do not want to use a pre-built release of `EVG`, building the software from source code is a feasible option.

1. First, obtain the source code.

```shell
git clone https://github.com/JiaoLab2021/EVG.git
cd EVG
```

2. Next, compile the software and add the current directory to your system's `PATH` environment variable.

```shell
cmake ./
make
chmod +x EVG.py
ln -sf EVG.py EVG
echo 'export PATH="$PATH:'$(pwd)'"' >> ~/.bashrc
source ~/.bashrc
```

3. To verify that the software has been installed correctly, perform a test run using the following steps:

```shell
EVG -h
graphvcf -h
fastAQ -h
```

## Usage

**Input Files**

* Reference Genome
* VCF File of Population Variants
* Sample File:

```shell
# Sample File
sample1 path_to_sample1_read1 path_to_sample1_read2
sample2 path_to_sample2_read1 path_to_sample2_read2
...
sampleN path_to_sampleN_read1 path_to_sampleN_read2
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

The results are stored in the `genotype/` folder, and each file is named after the corresponding sample listed in `sample.txt`: `sample1.vcf.gz`, `sample2.vcf.gz`, ..., `sampleN.vcf.gz`.

**Parameter**

* `--depth`: This parameter specifies the maximum sequencing data depth allowed for downstream analysis. If this value is exceeded, EVG will randomly downsample reads to the specified level in order to speed up the run. The default downsampling level is set at 15×, but it can be adjusted to meet specific requirements.
* `--mode`: This parameter determines the operating mode of `EVG`. In fast mode, only certain software is utilized to genotype SNPs and indels, while precise mode employs all software to genotype all variants.
* `--force`: If there are pre-existing files in the running directory of `EVG`, this parameter can be used to forcibly empty the folder. Otherwise, the software will encounter an error and exit.
* `--restart`: This parameter allows the software to resume from where it left off if it unexpectedly stops, enabling a breakpoint restart. Note that software completion is determined by file existence. It's recommended to manually check for incomplete or empty files before using this parameter and delete them.
