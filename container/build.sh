#!/usr/bin/env bash
# Build (and optionally push) the all-in-one pipeline image.
# Always targets linux/amd64 because bioconda has no linux-aarch64 builds for
# several of these tools and the GBI cluster is x86_64.
set -euo pipefail

IMAGE="${IMAGE:-ghcr.io/eit-gbi/longread-pipeline:latest}"
HERE="$(cd "$(dirname "$0")" && pwd)"

PUSH=0
[[ "${1:-}" == "--push" ]] && PUSH=1

if [[ "$PUSH" == "1" ]]; then
  # buildx can build for a foreign platform and push in one step.
  docker buildx build --platform=linux/amd64 -t "$IMAGE" --push "$HERE"
else
  # --load brings the amd64 image into the local docker image store
  # (it will run under emulation on Apple Silicon, which is fine for testing).
  docker buildx build --platform=linux/amd64 -t "$IMAGE" --load "$HERE"
fi

echo "Built: $IMAGE"
