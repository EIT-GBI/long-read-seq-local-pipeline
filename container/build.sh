#!/usr/bin/env bash
# Build (and optionally push) the all-in-one pipeline image.
# Always targets linux/amd64 because bioconda has no linux-aarch64 builds for
# several of these tools and the GBI cluster is x86_64.
set -euo pipefail

IMAGE="${IMAGE:-ghcr.io/eit-gbi/longread-pipeline:latest}"
HERE="$(cd "$(dirname "$0")" && pwd)"

PUSH=0
[[ "${1:-}" == "--push" ]] && PUSH=1

# This script builds with Docker/buildx and is meant to run on a machine that
# HAS Docker (e.g. your Mac). The Sandpit cluster has no Docker — build the .sif
# there with Apptainer instead (no registry needed):
#   srun -p cpu apptainer build \
#     /mnt/lustre/containers/eit-gbi/longread-pipeline.sif container/longread.def
if ! command -v docker >/dev/null 2>&1; then
  echo "ERROR: 'docker' not found." >&2
  echo "  - On your Mac: install/start Docker Desktop, then re-run this script." >&2
  echo "  - On the Sandpit cluster (no Docker): build the .sif with Apptainer:" >&2
  echo "      srun -p cpu apptainer build \\" >&2
  echo "        /mnt/lustre/containers/eit-gbi/longread-pipeline.sif container/longread.def" >&2
  exit 1
fi

if [[ "$PUSH" == "1" ]]; then
  # buildx can build for a foreign platform and push in one step.
  docker buildx build --platform=linux/amd64 -t "$IMAGE" --push "$HERE"
else
  # --load brings the amd64 image into the local docker image store
  # (it will run under emulation on Apple Silicon, which is fine for testing).
  docker buildx build --platform=linux/amd64 -t "$IMAGE" --load "$HERE"
fi

echo "Built: $IMAGE"
