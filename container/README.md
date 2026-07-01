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
mkdir -p "$HOME/scratch/containers"
srun -p cpu apptainer build --ignore-fakeroot-command \
    "$HOME/scratch/containers/longread-pipeline.sif" \
    container/longread.def
```

Two Sandpit-specific gotchas are baked into that command:

- **`--ignore-fakeroot-command`** — the login node isn't in `/etc/subuid`, so
  Apptainer injects the host's `libfakeroot.so` (built against GLIBC_2.38) into the
  Debian-12 base (GLIBC 2.36) and the build dies with `GLIBC_2.38 not found`. Our
  `%post` only runs `micromamba install` into `/opt/conda`, which doesn't need
  fakeroot, so skipping it is safe.
- **Write to a path you own** (`$HOME/scratch/...`). `/mnt/lustre/containers/...` is
  admin-managed and not user-writable — building straight there fails at the final
  step with `permission denied` (the build itself succeeds; only the SIF write
  fails). To make it a *shared* image, build to scratch, then ask the cluster admins
  to place the finished `.sif` under `/mnt/lustre/containers/eit-gbi/`.

If the build still fails, use the registry-free Docker route below (D).

Store the `.sif` under `/mnt/lustre/containers` (the doc's preferred spot for large
shared images). Then run the pipeline — bind every filesystem `config.txt` touches
(reads, reference, outdir all live under `/mnt/gbi-shared`):

```bash
srun -p cpu apptainer exec \
    --bind /mnt/gbi-shared \
    "$HOME/scratch/containers/longread-pipeline.sif" \
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

> **GHCR access.** GHCR packages are **private by default**. Both cluster flavours
> below pull anonymously, so an un-pushed or private image fails with
> `403 Forbidden` at `https://ghcr.io/token`. Either make the package public
> (GitHub → Packages → the package → *Package settings* → change visibility to
> Public), or give Enroot/Apptainer a credential: put a GHCR PAT (scope
> `read:packages`) in `~/.config/enroot/.credentials` as
> `machine ghcr.io login <user> password <PAT>`. If you just want to explore, skip
> all of this and use the registry-free `.sif` from section **A**.

Then on the cluster, either flavour works (image must be pushed + accessible):

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

## D. Registry-free: Docker on Mac → .sif on cluster

If the cluster fakeroot build (A) fails and you don't want to use GHCR, ship the
image as a tarball. Building the `.sif` from a Docker archive runs **no `%post`**,
so there's no fakeroot/glibc problem.

```bash
# on your Mac (image already built by container/build.sh)
docker save ghcr.io/eit-gbi/longread-pipeline:latest | gzip > longread-pipeline.tar.gz
scp longread-pipeline.tar.gz sandpit-tokyo-login:'~/scratch/containers/'

# on the cluster
mkdir -p "$HOME/scratch/containers"
srun -p cpu apptainer build \
    "$HOME/scratch/containers/longread-pipeline.sif" \
    docker-archive:"$HOME/scratch/containers/longread-pipeline.tar.gz"
```

Then run it exactly as in section A.

## Interactive use

For poking around, testing commands by hand, or debugging the reference/mapping
before committing to a full run, drop into a shell **inside** the container. Do it
on a compute node (via `srun --pty`), not the login node.

### Apptainer (from the `.sif`)

```bash
# allocate an interactive CPU session and open a shell in the container
srun -p cpu --pty apptainer shell \
    --bind /mnt/gbi-shared \
    "$HOME/scratch/containers/longread-pipeline.sif"
```

You land at an `Apptainer>` prompt with every tool on `$PATH` and your bound
filesystems visible. Your current directory is preserved, so you can work straight
from the repo:

```text
Apptainer> cd ~/code/long-read-seq-local-pipeline
Apptainer> minimap2 --version
Apptainer> samtools view -H /mnt/gbi-shared/.../hifi_reads/....bam | head
Apptainer> seqkit stats /mnt/gbi-shared/.../reads.fastq
# run the whole pipeline, or just paste individual steps to experiment
Apptainer> bash run_pipeline.sh config.txt
Apptainer> exit
```

Ask for more resources on the interactive session as needed, e.g.
`srun -p cpu -c 8 --mem=16G --pty apptainer shell --bind /mnt/gbi-shared <sif>`.

Bind extra roots with repeated `--bind` (e.g. `--bind /mnt/gbi-shared --bind /mnt/lustre`).

### Pyxis/Enroot (from a registry image)

```bash
srun -p cpu --pty \
    --container-image="ghcr.io#eit-gbi/longread-pipeline:latest" \
    --container-mounts=/mnt/gbi-shared:/mnt/gbi-shared \
    --container-workdir="$PWD" \
    bash
```

This opens a normal `bash` prompt inside the container. First launch imports the
image (slower); reruns on the same warm worker are fast.

> Tip: an interactive shell is the quickest way to settle the reference question —
> mount the reads + candidate references and try `minimap2`/`seqkit`/`vsearch` by
> hand until reads map, then bake the winning reference into `config.txt`.

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
