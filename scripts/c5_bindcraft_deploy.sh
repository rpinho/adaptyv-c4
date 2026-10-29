#!/bin/bash
# Deploy BindCraft multi-template C5 design run to a RunPod A40.
# Usage: bash c5_bindcraft_deploy.sh <ssh-spec>
#   e.g. bash c5_bindcraft_deploy.sh "-p 22114 root@194.68.245.6"
set -euo pipefail

SSH_TARGET="${1:?usage: c5_bindcraft_deploy.sh <ssh-spec>}"
SSH="ssh -T -i ~/.ssh/runpod_ed25519 -o StrictHostKeyChecking=no $SSH_TARGET"
SCP="scp -i ~/.ssh/runpod_ed25519 -o StrictHostKeyChecking=no"
REPO=~/dev/adaptyv-c4

echo "=== C5 BindCraft multi-template deployment ==="

# 1. Install BindCraft on the pod
echo "[1/5] Installing BindCraft..."
$SSH 'bash -s' << 'INSTALL_EOF'
set -euo pipefail
cd /root

if [ -d /root/BindCraft ] && python3 -c "import colabdesign" 2>/dev/null; then
    echo "BindCraft already installed, skipping..."
else
    # Clone repo
    [ -d /root/BindCraft ] || git clone https://github.com/martinpacesa/BindCraft.git /root/BindCraft
    cd /root/BindCraft

    # Install dependencies (no conda on RunPod — use pip directly)
    pip install --break-system-packages \
        pandas matplotlib 'numpy<2.0.0' biopython scipy pdbfixer seaborn tqdm \
        py3dmol dm-haiku 'flax<0.10.0' dm-tree joblib ml-collections immutabledict optax 2>&1 | tail -5

    # Install ColabDesign
    pip install --break-system-packages git+https://github.com/sokrypton/ColabDesign.git --no-deps 2>&1 | tail -3

    # Install PyRosetta
    pip install --break-system-packages pyrosetta --find-links https://west.rosettacommons.org/pyrosetta/quarterly/release.cxx11thread.serialization 2>&1 | tail -3

    # Check JAX GPU
    python3 -c "import jax; print('JAX devices:', jax.devices())"
    python3 -c "import colabdesign; print('ColabDesign OK')"
    python3 -c "import pyrosetta; print('PyRosetta OK')"

    # Download AF2 weights if not present
    if [ ! -d /root/BindCraft/params ]; then
        echo "Downloading AF2 weights..."
        mkdir -p /root/BindCraft/params
        wget -q -O /root/BindCraft/params/alphafold_params_2022-12-06.tar \
            "https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar"
        cd /root/BindCraft/params && tar xf alphafold_params_2022-12-06.tar && cd /root/BindCraft
    fi

    # Make dssp and DAlphaBall executable
    chmod +x /root/BindCraft/functions/dssp /root/BindCraft/functions/DAlphaBall.gcc
fi

echo "BindCraft install complete."
INSTALL_EOF

# 2. Upload target PDBs
echo "[2/5] Uploading target PDBs..."
$SSH 'mkdir -p /root/c5_design/targets /root/c5_design/settings /root/c5_design/output'

# Parse SSH_TARGET to get the scp target (handle -p port)
SCP_PORT=""
SCP_HOST=""
if echo "$SSH_TARGET" | grep -q "\-p"; then
    PORT=$(echo "$SSH_TARGET" | grep -oP '(?<=-p )\d+')
    HOST=$(echo "$SSH_TARGET" | awk '{print $NF}')
    SCP_PORT="-P $PORT"
    SCP_HOST="$HOST"
else
    SCP_HOST="$SSH_TARGET"
fi

$SCP $SCP_PORT \
    "$REPO/targets/c5_pmhc/state1_g12d_complex.pdb" \
    "$REPO/targets/c5_pmhc/state2_wt_complex.pdb" \
    "$REPO/targets/c5_pmhc/state3_peptide_free.pdb" \
    "$SCP_HOST:/root/c5_design/targets/"

# 3. Upload config files
echo "[3/5] Uploading configuration..."
$SCP $SCP_PORT \
    "$REPO/scripts/c5_bindcraft_target.json" \
    "$REPO/scripts/c5_bindcraft_advanced.json" \
    "$SCP_HOST:/root/c5_design/settings/"

# 4. Upload runner script
echo "[4/5] Uploading runner..."
$SCP $SCP_PORT \
    "$REPO/scripts/c5_bindcraft_runner.py" \
    "$SCP_HOST:/root/c5_design/"

# 5. Verify
echo "[5/5] Verifying deployment..."
$SSH 'ls -la /root/c5_design/targets/ && ls -la /root/c5_design/settings/ && python3 -c "import jax; print(\"GPU:\", jax.devices())"'

echo "=== Deployment complete ==="
echo "To run: ssh ... 'cd /root/c5_design && nohup python3 c5_bindcraft_runner.py > run.log 2>&1 &'"
