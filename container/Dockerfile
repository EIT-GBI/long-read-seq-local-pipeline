# All-in-one image for the long-read (Nanopore / PacBio HiFi) bash pipeline.
# Single monolithic image (run_pipeline.sh calls every tool in one process),
# built the same way as the nf-dnaseq modules: micromamba + bioconda/conda-forge.
#
# NOTE: bioconda ships linux-64 (x86_64) packages. The GBI cluster is x86_64,
# so ALWAYS build with `--platform=linux/amd64` (see container/README.md),
# even from an Apple-Silicon Mac.
FROM mambaorg/micromamba:1.5.8

USER root

# Pinned where it's safe; htslib and ucsc-bedgraphtobigwig are left unpinned so
# the solver picks builds with a consistent libzlib (the old ucsc=445 build
# needs libzlib<1.3, which conflicts with modern htslib/samtools). htslib is
# what provides the `tabix`/`bgzip` binaries.
RUN micromamba install -y -n base -c bioconda -c conda-forge \
        minimap2=2.30 \
        samtools=1.23.1 \
        bcftools=1.23.1 \
        htslib \
        bedtools=2.31.1 \
        ucsc-bedgraphtobigwig \
        chopper=0.9.0 \
        nanoplot=1.44.1 \
        python-kaleido=0.2.1 \
        vsearch \
        seqkit \
        parallel \
        procps-ng \
    && micromamba clean --all --yes

ENV PATH=/opt/conda/bin:$PATH \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8

LABEL org.opencontainers.image.title="long-read-seq-local-pipeline" \
      org.opencontainers.image.description="minimap2, samtools, bcftools, htslib/tabix, bedtools, bedGraphToBigWig, chopper, NanoPlot, vsearch, seqkit, GNU parallel" \
      org.opencontainers.image.source="https://github.com/EIT-GBI/ngs-seq-local-pipeline"

# Sanity-check at build time: every tool on PATH, and the two that can break at
# import/runtime (NanoPlot's kaleido dep, parallel's need for `ps`) actually run.
RUN for t in minimap2 samtools bcftools tabix bgzip bedtools bedGraphToBigWig chopper NanoPlot vsearch seqkit parallel; do \
        command -v "$t" >/dev/null || { echo "MISSING: $t" >&2; exit 1; }; \
    done \
    && NanoPlot --version >/dev/null \
    && echo "ok" | parallel echo >/dev/null \
    && echo "All pipeline tools present and runnable."

CMD ["bash"]
