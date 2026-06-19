# :dna: Long-read → BAM → Variants Pipeline (minimap2 / Samtools / BCFtools)

This repo contains a simple bash pipeline for processing **long-read** sequencing data
(Oxford Nanopore or PacBio HiFi) through alignment, variant calling, and consensus generation.
It is designed to run locally on macOS or on Ubuntu.

## :mag: Overview

Reads are single-end. Set `PLATFORM` in the config to `nanopore` or `pacbio-hifi`; this
selects the minimap2 preset (`map-ont` / `map-hifi`) and the bcftools mpileup profile
(`ont` / `pacbio-ccs`). Input can be FASTQ, or — for PacBio HiFi — an unaligned BAM (uBAM)
straight off CCS, which is streamed through `samtools fastq` into the aligner without
writing an intermediate FASTQ.

For each sample (`*.fastq[.gz]`, `*.fq[.gz]`, or `*.bam`), the pipeline performs:

1. Read length/quality filtering with **chopper**
2. Read alignment with **minimap2** (platform preset)
3. BAM sorting and indexing (**samtools**)
4. Read/alignment QC with **NanoPlot** and stats (`samtools flagstat`)
5. Variant calling (**bcftools mpileup + call**, haploid, platform-tuned via `-X`)
6. Generation of:
   - BCF and indexed VCF
   - CSV summary of variants
   - Consensus FASTA
   - BigWig coverage file (for IGV)

> **Note:** bcftools' model is short-read oriented. The `-X` profile tunes indel/error
> parameters per platform, but for production long-read calling consider a dedicated
> caller (Clair3, Medaka, or DeepVariant for HiFi).


## :hammer: Requirements

The following tools must be installed and available in `$PATH`:

- **minimap2**
- **samtools**
- **bcftools**
- **chopper** (long-read filtering; bioconda or `cargo install chopper`)
- **NanoPlot** (long-read QC; `pip install NanoPlot`)
- **htslib**
- **bedtools**
- **tabix**
- **parallel** (GNU parallel)
- **bedGraphToBigWig** (UCSC tool, not available via Homebrew, see instructions below)

## :computer: Local installation (Homebrew)

Clone the repository locally to the directory where you want your code to be.

```bash
git clone https://github.com/EIT-GBI/ngs-seq-local-pipeline.git
cd ngs-seq-local-pipeline
```

Make sure you have Homebrew installed on your MAC. If not, you can install it with the EIT Self Service portal. 

Then install the required tools with Homebrew:

```bash
brew install minimap2 samtools bcftools htslib tabixpp parallel bedtools
# chopper and NanoPlot are not in Homebrew:
pip install NanoPlot
brew install chopper   # or: conda install -c bioconda chopper / cargo install chopper
```

To generate bigwig files, we also need the UCSC tool 'bedGraphToBigWig'. This is not available via Homebrew, but you can download the precompiled binary from their website.

If you are working on a Mac do:
```bash
mkdir -p tools/ucsc
curl -L -o tools/ucsc/bedGraphToBigWig \
  https://hgdownload.cse.ucsc.edu/admin/exe/macOSX.arm64/bedGraphToBigWig
chmod +x tools/ucsc/bedGraphToBigWig
```

If you are working on Linux do:
```bash
mkdir -p tools/ucsc
curl -L -o tools/ucsc/bedGraphToBigWig \
  https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
chmod +x tools/ucsc/bedGraphToBigWig
```

Alternatively, download the binary and place it in the `tools/ucsc` directory of the repository, then make it executable. The script assumes the path to the tool is `tools/ucsc/bedGraphToBigWig` relative to the script location.

## :cloud: Installation on Cluster -> Ubuntu 22.04 LTS

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

> **Note (long-read tools):** `spack.lock` is from the previous short-read pipeline and is
> now stale — regenerate it (`spack concretize -f` in the activated env) after updating
> `spack.yaml`. Also, **chopper** and **NanoPlot** are not reliably packaged in spack;
> install them separately (`pip install NanoPlot`; `conda install -c bioconda chopper` or
> `cargo install chopper`). `filtlong` is a spack-packaged fallback filter if needed.


Here we named our environment `ngs-pipeline`, but you can choose any name you like, as long as activate the correct envirnment before running the pipeline. Which brings us to the next point - every time you want to run the pipeline, make sure to activate the Spack environment first with `spack env activate ngs-pipeline`. This will ensure that all the required tools are available in your `$PATH` and the pipeline can run without any issues.

```bash
spack env activate ngs-pipeline
bash ./run_pipeline.sh /path/to/my_config.txt
```

## :running: How to run the pipeline

Add info about the pipeline in config.txt:

```text
# Sample name
SAMPLE=sample1

# Sequencing platform: nanopore | pacbio-hifi
PLATFORM=nanopore

# Input folder containing reads (FASTQ, or uBAM for PacBio HiFi)
INPUT_DIR=/path/to/reads

# Output folder to store results
OUTPUT_DIR=/path/to/output/directory

# Reference genome FASTA
REFERENCE_GENOME=/path/to/reference/genome.fasta

# Number of threads
THREADS=8

# Minimum mapping quality
MIN_MAPQ=20

# Minimum base quality
MIN_QUAL=20

# Minimum depth (coverage threshold for consensus)
MIN_DEPTH=10

# Long-read pre-alignment filtering (chopper)
MIN_READ_LEN=1000   # drop reads shorter than this (bp)
MIN_READ_QUAL=10    # drop reads with mean quality below this (Phred)

# Nr of threads per sample
THREADS_PER_SAMPLE=4

# Maximum number of samples to process in parallel
SAMPLES_PARALLEL=4

# Memory for sorting BAM files
SORT_MEM=512M

# Patterns to identify input read files (single-end long reads + uBAM)
INPUT_PATTERN=(
  "*.fastq.gz"
  "*.fastq"
  "*.fq.gz"
  "*.fq"
  "*.bam"
)

# Whether to run NanoPlot QC
RUN_NANOPLOT=1
```

Then run:
```bash
# Use default config.txt in current directory
bash ./run_pipeline.sh

# Or specify a config file path
bash ./run_pipeline.sh /path/to/my_config.txt
```
I recommend copying the `config.txt` file to your data directory and editing it there. Then specify the path to that config when running the pipeline like in the example above. We plan to update the pipeline from time to time to add more features, so keeping a copy of the config file with the settings you used for each run will be helpful for reproducibility. And your config will not be overwritten every time you do git clone or pull.

Pipeline is parallelized across samples with GNU parallel, so multiple samples can be processed simultaneously. All logs and outputs will be saved in the specified output directory.


## :clipboard: TO DO:
- for consensus files add a check to see if there's any coverage at all
