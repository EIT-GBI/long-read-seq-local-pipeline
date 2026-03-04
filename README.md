# FASTQ → BAM → Variants Pipeline (BWA / Samtools / BCFtools)

This repo contains a simple bash pipeline for processing paired-end FASTQ files through alignment, variant calling, and consensus generation.  
It is designed to run locally on MAC OS for now

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
- **deeptools** (for `bamCoverage`)
- **tabix**
- **parallel** (GNU parallel)


### Recommended installation (Homebrew)

```bash
brew install bwa samtools bcftools fastp htslib tabix deeptools parallel
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
Pipeline is parallelized with GNU parallel or xargs (this is commented out at the moment), so multiple samples can be processed simultaneously. All logs and outputs will be saved in the specified output directory.