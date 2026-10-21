#!/usr/bin/env bash
# Mirror every campaign artefact off the rented pods onto the Mini.
#
# Runpod pods are TERMINATED, not stopped, when the prepaid balance reaches zero,
# and Auto-Pay is off on this account. Pod-local disk is therefore not storage -
# it is a scratch buffer. Anything that is not here is one billing event from gone.
set -o pipefail   # deliberately NOT -e: a pod that is gone must not abort the other pull
DEST="${1:-$HOME/dev/adaptyv-c4/results/pods}"
mkdir -p "$DEST"
ts=$(date -u +%Y-%m-%dT%H:%M:%SZ)

pull () {            # pull <ssh-host> <remote-run-dir> <local-name>
  local host="$1" run="$2" name="$3"
  local out="$DEST/$name"
  mkdir -p "$out"
  # the ranked designs, the per-campaign CSVs and the worker logs; NOT the
  # trajectory dirs, which are large and regenerable
  rsync -az --timeout=120 \
        --include='3_Ranked/***' --include='*.csv' --include='*.json' \
        --include='workers/***' --include='*.log' \
        --exclude='1_Trajectories/***' --exclude='compile_cache/***' --exclude='*' \
        "$host:$run/" "$out/" 2>&1 | tail -3 || echo "  rsync failed for $name"
  local n=0
  if [ -d "$out/3_Ranked" ]; then
    n=$( { find "$out/3_Ranked" -type f \( -name '*.cif' -o -name '*.pdb' \) 2>/dev/null || true; } | wc -l | tr -d ' ')
  fi
  echo "  $name: $n ranked designs mirrored to $out"
}

echo "[$ts] pulling campaign results"
pull adaptyv-c4 /workspace/runs/c4_E1_human  c4_E1_human
pull adaptyv-c3 /workspace/runs/c3_glp1r     c3_glp1r
# the target inputs too, so a run is reproducible without the pod
rsync -az --timeout=60 --include='*.pdb' --include='*.json' --exclude='*' \
      adaptyv-c4:/workspace/runs/ "$DEST/c4_inputs/" 2>/dev/null || true
echo "[$ts] done"
