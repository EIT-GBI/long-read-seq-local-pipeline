# Container for the long-read pipeline

One image with every tool `run_pipeline.sh` calls: **minimap2, samtools, bcftools,
tabix/bgzip (htslib), bedtools, bedGraphToBigWig (ucsc), chopper, NanoPlot, GNU parallel**.
Built with micromamba + bioconda, same recipe as the `nf-dnaseq` modules.

## Build & push (from your Mac)

The cluster is x86_64, so the image is always built for `linux/amd64`.

```bash
# build locally (amd64, runs under emulation on Apple Silicon)
container/build.sh

# build + push to GHCR (log in first: echo $PAT | docker login ghcr.io -u <user> --password-stdin)
container/build.sh --push
```

Override the tag with `IMAGE=ghcr.io/eit-gbi/longread-pipeline:v1 container/build.sh`.

## Run on the cluster (Apptainer)

Apptainer pulls OCI images straight from GHCR and caches a `.sif`. Point the cache at
the same dir your Nextflow profile uses:

```bash
export APPTAINER_CACHEDIR=/mnt/gbi-shared/home/cristian-eitgbi/singularity_cache

# one-off pull -> .sif (optional; exec will pull on first use too)
apptainer pull "$APPTAINER_CACHEDIR/longread-pipeline.sif" \
    docker://ghcr.io/eit-gbi/longread-pipeline:latest

# run the whole pipeline inside the container.
# --bind every filesystem your config.txt touches (reads, reference, outdir).
apptainer exec \
    --bind /mnt/gbi-shared \
    "$APPTAINER_CACHEDIR/longread-pipeline.sif" \
    bash run_pipeline.sh config.txt
```

If your data and repo live under different roots, bind each one, e.g.
`--bind /mnt/gbi-shared --bind /scratch`.

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
