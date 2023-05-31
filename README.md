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
