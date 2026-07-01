# Container for the long-read pipeline

One image with every tool `run_pipeline.sh` calls: **minimap2, samtools, bcftools,
tabix/bgzip (htslib), bedtools, bedGraphToBigWig (ucsc), chopper, NanoPlot, GNU parallel**,
plus **vsearch** and **seqkit** for downstream sequence work.
Built with micromamba + bioconda, same recipe as the `nf-dnaseq` modules.

The cluster (Sandpit / Tokyo Slurm) has **no Docker**. Pick one of two paths:

- **A. Build the `.sif` on the cluster** with Apptainer, from `longread.def` — no
  Docker and no registry needed. Simplest.
- **B. Build with Docker on your Mac, push to GHCR**, then pull/run on the cluster
  via Apptainer or Pyxis. Use this if you want a shared registry image.

## A. Build on the cluster (Apptainer, no Docker)

```bash
# on the Tokyo login node, from the repo root
srun -p cpu apptainer build \
    /mnt/lustre/containers/eit-gbi/longread-pipeline.sif \
    container/longread.def
# add --fakeroot after `apptainer build` if your site requires it.
```

Store the `.sif` under `/mnt/lustre/containers` (the doc's preferred spot for large
shared images). Then run the pipeline — bind every filesystem `config.txt` touches
(reads, reference, outdir all live under `/mnt/gbi-shared`):

```bash
srun -p cpu apptainer exec \
    --bind /mnt/gbi-shared \
    /mnt/lustre/containers/eit-gbi/longread-pipeline.sif \
    bash run_pipeline.sh config.txt
```

## B. Build on your Mac, push to GHCR

The cluster is x86_64, so the image is always built for `linux/amd64`.

```bash
container/build.sh                 # build locally (amd64, emulated on Apple Silicon)
container/build.sh --push          # build + push to GHCR
# log in first: echo $PAT | docker login ghcr.io -u <user> --password-stdin
```

Override the tag with `IMAGE=ghcr.io/eit-gbi/longread-pipeline:v1 container/build.sh --push`.
Make the GHCR package **public**, or the cluster needs registry credentials to pull.

Then on the cluster, either flavour works:

```bash
# Apptainer: pulls the OCI image and caches a .sif under scratch
export APPTAINER_CACHEDIR="$HOME/scratch/apptainer-cache"
srun -p cpu apptainer exec --bind /mnt/gbi-shared \
    docker://ghcr.io/eit-gbi/longread-pipeline:latest \
    bash run_pipeline.sh config.txt

# Pyxis/Enroot: note the host#path syntax and explicit mounts/workdir
srun -p cpu \
    --container-image="ghcr.io#eit-gbi/longread-pipeline:latest" \
    --container-mounts=/mnt/gbi-shared:/mnt/gbi-shared \
    --container-workdir="$PWD" \
    bash run_pipeline.sh config.txt
```

## One pipeline tweak needed

`run_pipeline.sh` hardcodes `bedGraphToBigWig` at `tools/ucsc/bedGraphToBigWig`.
Inside the container it's already on `$PATH`, so either:

- change that line to just `BEDGRAPHTOBIGWIG="bedGraphToBigWig"`, **or**
- `mkdir -p tools/ucsc && ln -s "$(command -v bedGraphToBigWig)" tools/ucsc/bedGraphToBigWig`
  (run inside the container / on a bound path).

Everything else runs unchanged.

## Versions

Pinned in the `Dockerfile`. `bcftools`, `bedtools`, `ucsc-bedgraphtobigwig` match the
`nf-dnaseq` modules. If a `micromamba` solve fails, bump the offending pin — bioconda
sometimes drops old patch builds.
