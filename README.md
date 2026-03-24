# FASTQ → BAM → Variants Pipeline (BWA / Samtools / BCFtools)

This repo contains a simple bash pipeline for processing paired-end FASTQ files through alignment, variant calling, and consensus generation.
It is designed to run locally on macOS or on Ubuntu.

## Overview

For each paired-end sample (`*_R1*.fastq` / `*_R2*.fastq`), the pipeline performs:

1. Read alignment with **BWA MEM**
2. BAM conversion, sorting, and indexing (**samtools**)
3. Alignment statistics (`samtools flagstat`)
4. Variant calling (**bcftools mpileup + call**, haploid)
5. Generation of:
   - BCF and indexed VCF
   - CSV summary of variants
   - Consensus FASTA
   - BigWig coverage file (for IGV)


## Requirements

The following tools must be installed and available in `$PATH`:

- **bwa**
- **samtools**
- **bcftools**
- **fastp**
- **htslib**
- **bedtools**
- **tabix**
- **parallel** (GNU parallel)
- **bedGraphToBigWig** (UCSC tool, not available via Homebrew, see instructions below)

## Recommended installation (Homebrew)

Clone the repository locally to the directory where you want your code to be.

```bash
git clone https://github.com/EIT-GBI/ngs-seq-local-pipeline.git
cd ngs-seq-local-pipeline
```

Make sure you have Homebrew installed on your MAC. If not, you can install it with the EIT Self Service portal. 

Then install the required tools with Homebrew:

```bash
brew install bwa samtools bcftools fastp htslib tabixpp parallel bedtools fastqc
```

To generate bigwig files, we also need the UCSC tool 'bedGraphToBigWig'. This is not available via Homebrew, but you can download the precompiled binary for Mac from their website:

```bash
mkdir -p tools/ucsc
# if on Mac do:
curl -L -o tools/ucsc/bedGraphToBigWig \
  https://hgdownload.cse.ucsc.edu/admin/exe/macOSX.arm64/bedGraphToBigWig

# if on Linux do:
curl -L -o tools/ucsc/bedGraphToBigWig \
  https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig

chmod +x tools/ucsc/bedGraphToBigWig
```

Alternatively, download the binary and place it in the `tools/ucsc` directory of the repository, then make it executable. The script assumes the path to the tool is `tools/ucsc/bedGraphToBigWig` relative to the script location.

## Installation on Cluster -> Ubuntu 22.04 LTS

We are going to use Spack to install the required tools on the cluster. A description of how to use Spack is available at the link below, but we'll try to abstract away the details and just give you the commands to run to get the tools installed - https://spack.readthedocs.io/en/latest/getting_started.html.

First, clone the Spack repository and set it up:

```bash
git clone --depth=2 https://github.com/spack/spack.git

# Source the script for your (bash) shell
. spack/share/spack/setup-env.sh
```

Spack enviroments are somewhat similar to conda environments. One can create an environment with a specific set of packages and then activate it to use those packages. We have created a Spack environment file for the pipeline with the required tools and their versions specified - `spack.lock`. You can find it in the `spack_env/` directory of the repository. To create and activate the environment, run: 

```bash
spack env create ngs-pipeline spack_env/spack.lock
spack env activate ngs-pipeline
# To check that the tools are available, you can run:
spack find
``` 
[EXTRA]. You will notice that there is also a `spack_env/spack.yaml` file. One could use this file to create the same environment, although the `spack.lock` file is more specific and ensures that the exact versions of the tools are installed. If you want to create the environment from the `spack.yaml` file, run:

```bash
spack env create ngs-pipeline spack_env/spack.yaml
spack env activate ngs-pipeline
```


Here we named our environment `ngs-pipeline`, but you can choose any name you like, as long as activate the correct envirnment before running the pipeline. Which brings us to the next point - every time you want to run the pipeline, make sure to activate the Spack environment first with `spack env activate ngs-pipeline`. This will ensure that all the required tools are available in your `$PATH` and the pipeline can run without any issues.

```bash
spack env activate ngs-pipeline
bash ./run_pipeline.sh /path/to/my_config.txt
```

## Usage

Add info about the pipeline in config.txt:

```text
# Sample name
SAMPLENAME=sample1

# Input directory for FASTQ files
INPUT_DIR=/path/to/fastq/files

# Output directory
OUTPUT_DIR=/path/to/output/directory

# Reference genome FASTA
REFERENCE=/path/to/reference/genome.fasta

# Number of threads
THREADS=4

# Minimum mapping quality
MIN_MAPQ=20

# Minimum base quality
MIN_BASEQ=20

# Coverage threshold for consensus
COVERAGE_THRESHOLD=10

# Nr of threads for various steps
THREADS_PER_SAMPLE=2

# Maximum number of samples to process in parallel
SAMPLES_PARALLEL=4

# Memory for sorting BAM files
SORT_MEM=512M

# Patterns to identify R1 FASTQ files
FASTQ_PATTERN_R1=(
  "*_R1*.fastq.gz"
  "*_R1*.fastq"
  "*_1.fastq.gz"
  "*_1.fastq"
)
# Whether to run FastQC on raw reads
RUN_FASTQC=1
```

Then run:
```bash
# Use default config.txt in current directory
bash ./run_pipeline.sh

# Or specify a config file path
bash ./run_pipeline.sh /path/to/my_config.txt
```
I recommend copying the `config.txt` file to your data directory and editing it there. Then specify the path to that config when running the pipeline like in the example above. We plan to update the pipeline from time to time to add more features, so keeping a copy of the config file with the settings you used for each run will be helpful for reproducibility. And your config will not be overwritten every time you do git clone or pull.

Pipeline is parallelized with GNU parallel or xargs (this is commented out at the moment), so multiple samples can be processed simultaneously. All logs and outputs will be saved in the specified output directory.


## TO DO:
- for consensus files add a check to see if there's any coverage at all
- check if all packages can be indeed installed with homebrew 
