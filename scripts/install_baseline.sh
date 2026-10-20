#!/usr/bin/env bash
# Baseline install for the Challenge 4 entry: Miniforge + BindCraft (AF2 + ProteinMPNN + PyRosetta).
# Everything lands on the /workspace network volume so it survives pod stop/terminate.
set -euo pipefail
export DEBIAN_FRONTEND=noninteractive
W=/workspace
mkdir -p "$W/logs" "$W/src"

echo "=== [1/4] miniforge ==="
if [ ! -x "$W/conda/bin/conda" ]; then
  wget -q -O /tmp/miniforge.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
  bash /tmp/miniforge.sh -b -p "$W/conda"
fi
export PATH="$W/conda/bin:$PATH"
conda --version
conda config --set always_yes true --set auto_activate_base false

echo "=== [2/4] system deps ==="
apt-get update -qq
apt-get install -y -qq build-essential dssp 2>/dev/null || apt-get install -y -qq build-essential || true

echo "=== [3/4] BindCraft clone ==="
cd "$W/src"
if [ ! -d BindCraft ]; then
  git clone --depth 1 https://github.com/martinpacesa/BindCraft
fi
cd BindCraft
git log --oneline -1

echo "=== [4/4] BindCraft install (conda env + AF2 weights ~5.3 GB) ==="
bash install_bindcraft.sh --cuda '12.8' --pkg_manager 'conda'

echo "=== DONE $(date -u +%FT%TZ) ==="
