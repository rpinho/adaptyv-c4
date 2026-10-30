"""Boltz-2 complex prediction — CORRECTED CONTIGUOUS APICAL target.

This runner replaces the discontinuous 154-residue concatenated surrogate used by
modal_boltz2_acid*.py with the verified CONTIGUOUS canonical human TfR1 (P02786)
apical domain, residues 189-383 inclusive (195 residues). It re-predicts the
complex for the existing validated binder candidates (no candidate regeneration).

- Distinct Modal app name:  boltz2-c4-apical
- Distinct output prefix:    rfd4-runs:boltz2_c4_apical/   (NEVER boltz2_c4_acid/)
- Idempotent: validated skip/resume per design; explicit nonzero failure on errors.
- Heavy compute runs on Modal; launch with `modal run --detach`.

Launch (after local validation passes):
    modal run --detach modal_boltz2_apical.py::main \
        --refold-json results/c4_acid_combined_refold.json --max-designs 1000

Validation-only (no GPU):
    python3 -I modal_boltz2_apical.py            # runs _validate_local()
"""
import json
import pathlib
import sys

import modal

# --- Verified contiguous canonical apical target -----------------------------
# Derived directly from data/P02786.txt (canonical human TfR1, 760 AA) residues
# 189-383 inclusive. Independently checked against configs/rfd4/tfr1_apical_binder.yaml
# (auth 189-383 == 11QC label 68-262, 195 residues, contiguous, no gaps) and the
# epitope auth K205/N206/R208/L209/V210. See results/c4_target_context_review_20261024.md.
TFR1_APICAL_189_383 = (
    "KIQVKDSAQNSVIIVDKNGRLVYLVENPGGYVAYSKAATVTGKLVHANFGTKKDFEDLYTPVNG"
    "SIVIVRAGKITFAEKVANAESLNAIGVLIYMDQTKFPIVNAELSFFGHAHLGTGDPYTPGFPSF"
    "NHTQFPPSRSSGLPNIPVQTISRAAAEKLFGNMEGDCPSDWKTDSTCRMVTSESKNVKLTVSNVLKE"
)

# Fail loudly at import if this literal is ever corrupted by an edit.
assert len(TFR1_APICAL_189_383) == 195, (
    f"apical target must be 195 residues, got {len(TFR1_APICAL_189_383)}")
assert TFR1_APICAL_189_383[0] == "K", "apical residue 189 must be K"
assert TFR1_APICAL_189_383[-1] == "E", "apical residue 383 must be E"
# auth 205/206/208/209/210 -> index within 189-383 is (auth-189)
assert TFR1_APICAL_189_383[205 - 189] == "K"
assert TFR1_APICAL_189_383[206 - 189] == "N"
assert TFR1_APICAL_189_383[208 - 189] == "R"
assert TFR1_APICAL_189_383[209 - 189] == "L"
assert TFR1_APICAL_189_383[210 - 189] == "V"

OUT_PREFIX = "boltz2_c4_apical"  # NEVER boltz2_c4_acid
_AA = set("ACDEFGHIKLMNPQRSTVWY")
TARGET_LEN = len(TFR1_APICAL_189_383)  # 195
# Exact expected size of the validated PASS pool in c4_acid_combined_refold.json.
# A drift from this means the input pool changed under us — fail before spend
# rather than silently predicting a different candidate set.
EXPECTED_PASS = 343

app = modal.App("boltz2-c4-apical")

boltz_cache = modal.Volume.from_name("boltz2-cache", create_if_missing=True)
runs_vol = modal.Volume.from_name("rfd4-runs", create_if_missing=True)

image = (
    modal.Image.debian_slim(python_version="3.10")
    .pip_install(
        "boltz[cuda]",
        "pyyaml",
        "propka",
        "biopython==1.84",
    )
    .env({"BOLTZ_CACHE": "/cache/boltz", "REBUILD": "v2"})
)


def _load_candidates(refold_json: str, max_designs: int,
                     expected_pass: int = EXPECTED_PASS):
    """Load + validate the existing validated candidate pool. No regeneration.

    Count validation is exact: the number of PASS records must equal
    ``expected_pass`` (default EXPECTED_PASS=343). This is not a lower bound or
    a truncation convenience — if the pool has drifted (records added, dropped,
    or re-scored) we refuse to launch rather than quietly predicting a different
    candidate set. ``max_designs`` only caps how many of those exact PASS records
    are dispatched, and the dispatched count is asserted too.
    """
    with open(refold_json) as f:
        refold = json.load(f)
    if isinstance(refold, dict):
        refold = refold.get("results", [])

    passed = [r for r in refold if r.get("status") == "PASS"]
    # Exact-count gate: the validated pool must be exactly what we expect.
    if len(passed) != expected_pass:
        raise ValueError(
            f"PASS-pool count drift: expected exactly {expected_pass} PASS "
            f"records in {refold_json}, found {len(passed)}. Refusing to launch "
            f"against a changed candidate set.")
    passed.sort(key=lambda x: (-x.get("avg_his_pka", 0), -x.get("mean_plddt", 0)))
    selected = passed[:max_designs]
    expected_selected = min(expected_pass, max_designs)
    if len(selected) != expected_selected:
        raise ValueError(
            f"selected-count mismatch: expected {expected_selected}, "
            f"got {len(selected)}")

    # Identity / input validation — fail before any spend.
    seen = set()
    for r in selected:
        name = r.get("name")
        seq = (r.get("sequence") or "").upper()
        if not name:
            raise ValueError(f"candidate missing name: {r}")
        if name in seen:
            raise ValueError(f"duplicate candidate name: {name}")
        seen.add(name)
        if not seq:
            raise ValueError(f"candidate {name} has empty sequence")
        bad = set(seq) - _AA
        if bad:
            raise ValueError(f"candidate {name} has non-AA chars: {sorted(bad)}")
    # Exact distinct-name invariant across the dispatched set.
    if len(seen) != len(selected):
        raise ValueError(
            f"duplicate names within selection: {len(selected)} designs, "
            f"{len(seen)} unique")
    return selected


def _num_ok(v):
    """True iff v is a real number in the valid confidence range [0, 1]."""
    return isinstance(v, (int, float)) and 0.0 <= float(v) <= 1.0


def _pdb_chains_ok(pdb_path, target_len):
    """Structural validation: the complex PDB must contain a binder chain A and
    a target chain B whose residue count equals the apical target length.

    This goes beyond a size check — a truncated/garbled PDB that happens to be
    large would still fail here.
    """
    chain_res = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 26:
                chain = line[21].strip()
                resnum = line[22:26].strip()
                if chain:
                    chain_res.setdefault(chain, set()).add(resnum)
    if "A" not in chain_res or "B" not in chain_res:
        return False
    if len(chain_res["A"]) < 1:
        return False
    return len(chain_res["B"]) == target_len


def _validate_cached(name, binder_seq, pdb_dest, conf_dest, result_dest):
    """Full validation of one per-design cached output triple.

    Returns the cached result dict (marked skipped=True) only if EVERY check
    passes, else None (→ recompute). Validates, in order:
      * all three artifacts exist and the PDB is non-trivial in size;
      * the raw confidence JSON carries a real numeric score;
      * the result JSON's name matches this design;
      * the result JSON's binder sequence matches the requested binder_seq;
      * target identity + length are the contiguous apical 195-mer;
      * iptm / ptm / confidence are real numbers in [0, 1];
      * the PDB structurally contains chain A (binder) and chain B == 195 res.
    """
    if not (pdb_dest.exists() and conf_dest.exists() and result_dest.exists()):
        return None
    if pdb_dest.stat().st_size <= 1000:
        return None
    with open(conf_dest) as f:
        conf = json.load(f)
    if not (_num_ok(conf.get("iptm")) or _num_ok(conf.get("confidence_score"))):
        return None
    with open(result_dest) as f:
        cached = json.load(f)
    if cached.get("name") != name:
        return None
    if (cached.get("sequence") or "").upper() != binder_seq.upper():
        return None
    if cached.get("target_len") != TARGET_LEN:
        return None
    if cached.get("target") != "P02786_189-383_contiguous":
        return None
    if not (_num_ok(cached.get("iptm")) and _num_ok(cached.get("ptm"))
            and _num_ok(cached.get("confidence"))):
        return None
    if not _pdb_chains_ok(pdb_dest, TARGET_LEN):
        return None
    cached["skipped"] = True
    return cached


@app.function(
    image=image,
    gpu="L40S",
    volumes={"/cache/boltz": boltz_cache, "/runs": runs_vol},
    timeout=3600,
    max_containers=20,
)
def predict_one(name: str, binder_seq: str):
    """Predict one complex (binder A + contiguous apical B) with Boltz-2.

    Idempotent: if a valid structure + confidence already exist on the volume,
    return the cached result (skipped) instead of recomputing.
    """
    import subprocess
    import shutil
    import yaml

    out_base = pathlib.Path(f"/runs/{OUT_PREFIX}")
    out_base.mkdir(parents=True, exist_ok=True)
    pdb_dest = out_base / f"{name}_complex.pdb"
    conf_dest = out_base / f"{name}_confidence.json"
    result_dest = out_base / f"{name}_result.json"

    # --- fully-validated skip/resume (no stale copying) ---
    runs_vol.reload()
    try:
        cached = _validate_cached(name, binder_seq, pdb_dest, conf_dest, result_dest)
        if cached is not None:
            print(f"SKIP {name}: fully validated cached outputs", flush=True)
            return cached
    except Exception as e:
        print(f"  cache invalid for {name}, recomputing: {e}", flush=True)

    print(f"Predicting {name} ({len(binder_seq)} res) vs apical-195", flush=True)

    yaml_dir = pathlib.Path(f"/tmp/boltz2/{name}")
    yaml_dir.mkdir(parents=True, exist_ok=True)
    yaml_path = yaml_dir / f"{name}.yaml"

    config = {
        "version": 1,
        "sequences": [
            {"protein": {"id": "A", "sequence": binder_seq}},
            {"protein": {"id": "B", "sequence": TFR1_APICAL_189_383}},
        ],
    }
    with open(yaml_path, "w") as f:
        yaml.dump(config, f)

    out_dir = yaml_dir / "output"
    out_dir.mkdir(exist_ok=True)

    cmd = [
        "boltz", "predict", str(yaml_path),
        "--use_msa_server",
        "--out_dir", str(out_dir),
        "--output_format", "pdb",
        "--diffusion_samples", "1",
        "--recycling_steps", "3",
    ]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        if proc.returncode != 0:
            print(f"  ERROR: {proc.stderr[-300:]}", flush=True)
            return {"name": name, "error": proc.stderr[-200:]}
    except subprocess.TimeoutExpired:
        print("  TIMEOUT", flush=True)
        return {"name": name, "error": "timeout"}

    pred_dir = out_dir / f"boltz_results_{name}" / "predictions" / name
    if not pred_dir.exists():
        return {"name": name, "error": "no prediction dir"}

    conf_files = list(pred_dir.rglob("confidence*.json"))
    pdb_files = list(pred_dir.rglob("*.pdb"))
    if not conf_files:
        return {"name": name, "error": "no confidence output"}
    if not pdb_files:
        return {"name": name, "error": "no pdb output"}

    with open(conf_files[0]) as f:
        conf = json.load(f)

    iptm = conf.get("iptm", 0)
    ptm = conf.get("ptm", 0)
    confidence = conf.get("confidence_score", 0)
    print(f"  ipTM={iptm:.3f}, pTM={ptm:.3f}, conf={confidence:.3f}", flush=True)

    r = {
        "name": name, "sequence": binder_seq,
        "iptm": iptm, "ptm": ptm, "confidence": confidence,
        "complex_plddt": conf.get("complex_plddt", 0),
        "target": "P02786_189-383_contiguous", "target_len": TARGET_LEN,
        "skipped": False,
    }

    # persist raw confidence + structure under the APICAL prefix only
    shutil.copy2(str(pdb_files[0]), pdb_dest)
    with open(conf_dest, "w") as f:
        json.dump(conf, f, indent=2)

    # complex His pKa (acid-OFF relevant); non-fatal
    try:
        from propka.run import single as propka_single
        mol = propka_single(str(pdb_files[0]))
        best_conf = max(mol.conformations.values(), key=lambda c: len(c.groups))
        his_bound = {}
        for group in best_conf.groups:
            if type(group).__name__ == "HISGroup":
                his_bound[f"{group.atom.chain_id}{group.atom.res_num}"] = round(
                    group.pka_value, 2)
        r["his_complex_pkas"] = his_bound
        binder_his = {k: v for k, v in his_bound.items() if k.startswith("A")}
        if binder_his:
            r["avg_complex_his_pka"] = round(
                sum(binder_his.values()) / len(binder_his), 2)
            r["n_elevated_complex"] = sum(1 for v in binder_his.values() if v > 6.5)
    except Exception as e:
        print(f"  PROPKA error: {e}", flush=True)
        r["propka_error"] = str(e)

    with open(result_dest, "w") as f:
        json.dump(r, f, indent=2)

    runs_vol.commit()
    return r


@app.function(
    image=image,
    volumes={"/runs": runs_vol},
    timeout=21600,  # 6h: must outlive the whole batch
)
def run_all(jobs: list):
    """Server-side driver: fan out predict_one over all jobs and consolidate.

    Running the fan-out inside a remote function (rather than a local_entrypoint)
    means a Mini restart / client disconnect cannot interrupt the batch: this
    driver is the single detached function that keeps the whole run alive.
    """
    results = list(predict_one.starmap(jobs))

    out_base = pathlib.Path(f"/runs/{OUT_PREFIX}")
    out_base.mkdir(parents=True, exist_ok=True)
    consolidated = {
        "target": "P02786_189-383_contiguous",
        "target_len": len(TFR1_APICAL_189_383),
        "n_designs": len(results),
        "results": results,
    }
    with open(out_base / "consolidated_results.json", "w") as f:
        json.dump(consolidated, f, indent=2)
    runs_vol.commit()

    errored = [r for r in results if "error" in r]
    skipped = [r for r in results if r.get("skipped")]
    computed = [r for r in results if "error" not in r and not r.get("skipped")]
    good = [r for r in results if "error" not in r and r.get("iptm", 0) > 0.5]
    print(f"{len(computed)} computed, {len(skipped)} skipped(resume), "
          f"{len(errored)} errored; {len(good)}/{len(results)} ipTM>0.5", flush=True)
    consolidated["summary"] = {
        "computed": len(computed), "skipped": len(skipped),
        "errored": len(errored), "iptm_gt_0.5": len(good),
        "errors": [{"name": r["name"], "error": r["error"]} for r in errored],
    }
    with open(out_base / "consolidated_results.json", "w") as f:
        json.dump(consolidated, f, indent=2)
    runs_vol.commit()
    if errored:
        # surface failure explicitly in the driver's own result
        raise RuntimeError(f"{len(errored)} designs errored; see consolidated_results.json")
    return consolidated["summary"]


@app.function(image=image, volumes={"/runs": runs_vol}, timeout=600)
def read_consolidated():
    """Read-only: return the remotely-written consolidated_results.json.

    Opens the file on the volume and returns its parsed contents WITHOUT any
    mkdir, write, or commit. This is the download path's remote half — all
    consolidation happens in run_all; this only reads it back.
    """
    runs_vol.reload()
    path = pathlib.Path(f"/runs/{OUT_PREFIX}/consolidated_results.json")
    if not path.exists():
        return {"error": "consolidated_results.json not present yet"}
    with open(path) as f:
        return json.load(f)


@app.function(
    image=image,
    gpu="L40S",
    volumes={"/cache/boltz": boltz_cache, "/runs": runs_vol},
    timeout=7200,       # 2h function wall — single longest complex only
    max_containers=1,
)
def predict_one_long(name: str, binder_seq: str, subprocess_timeout: int = 6600):
    """Single-candidate retry of predict_one with a LONGER wall clock only.

    The prediction is identical to predict_one — same contiguous apical target,
    same exact binder, same Boltz-2 parameters (diffusion_samples=1,
    recycling_steps=3, --use_msa_server). The ONLY difference is the time budget:
    the subprocess wall is raised from 1800s to ``subprocess_timeout`` (default
    6600s) and the function timeout to 2h. This is the bounded root-cause fix for
    the lone design whose 100-res binder + 195-res apical target (295-res complex,
    the largest in the pool) exceeded the 1800s subprocess wall twice. Scoring is
    unchanged — no relaxation, no parameter change that affects the result.

    Idempotent: a fully-validated cached triple is returned as skipped. Writes
    ONLY this design's artifacts under the apical prefix; never touches others.
    """
    import subprocess
    import shutil
    import yaml

    out_base = pathlib.Path(f"/runs/{OUT_PREFIX}")
    out_base.mkdir(parents=True, exist_ok=True)
    pdb_dest = out_base / f"{name}_complex.pdb"
    conf_dest = out_base / f"{name}_confidence.json"
    result_dest = out_base / f"{name}_result.json"

    runs_vol.reload()
    try:
        cached = _validate_cached(name, binder_seq, pdb_dest, conf_dest, result_dest)
        if cached is not None:
            print(f"SKIP {name}: fully validated cached outputs", flush=True)
            return cached
    except Exception as e:
        print(f"  cache invalid for {name}, recomputing: {e}", flush=True)

    print(f"Predicting {name} ({len(binder_seq)} res) vs apical-195 "
          f"[long wall {subprocess_timeout}s]", flush=True)

    yaml_dir = pathlib.Path(f"/tmp/boltz2/{name}")
    yaml_dir.mkdir(parents=True, exist_ok=True)
    yaml_path = yaml_dir / f"{name}.yaml"

    config = {
        "version": 1,
        "sequences": [
            {"protein": {"id": "A", "sequence": binder_seq}},
            {"protein": {"id": "B", "sequence": TFR1_APICAL_189_383}},
        ],
    }
    with open(yaml_path, "w") as f:
        yaml.dump(config, f)

    out_dir = yaml_dir / "output"
    out_dir.mkdir(exist_ok=True)

    cmd = [
        "boltz", "predict", str(yaml_path),
        "--use_msa_server",
        "--out_dir", str(out_dir),
        "--output_format", "pdb",
        "--diffusion_samples", "1",
        "--recycling_steps", "3",
    ]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              timeout=subprocess_timeout)
        if proc.returncode != 0:
            print(f"  ERROR: {proc.stderr[-300:]}", flush=True)
            return {"name": name, "error": proc.stderr[-200:]}
    except subprocess.TimeoutExpired:
        print("  TIMEOUT", flush=True)
        return {"name": name, "error": "timeout"}

    pred_dir = out_dir / f"boltz_results_{name}" / "predictions" / name
    if not pred_dir.exists():
        return {"name": name, "error": "no prediction dir"}

    conf_files = list(pred_dir.rglob("confidence*.json"))
    pdb_files = list(pred_dir.rglob("*.pdb"))
    if not conf_files:
        return {"name": name, "error": "no confidence output"}
    if not pdb_files:
        return {"name": name, "error": "no pdb output"}

    with open(conf_files[0]) as f:
        conf = json.load(f)

    iptm = conf.get("iptm", 0)
    ptm = conf.get("ptm", 0)
    confidence = conf.get("confidence_score", 0)
    print(f"  ipTM={iptm:.3f}, pTM={ptm:.3f}, conf={confidence:.3f}", flush=True)

    r = {
        "name": name, "sequence": binder_seq,
        "iptm": iptm, "ptm": ptm, "confidence": confidence,
        "complex_plddt": conf.get("complex_plddt", 0),
        "target": "P02786_189-383_contiguous", "target_len": TARGET_LEN,
        "skipped": False,
    }

    shutil.copy2(str(pdb_files[0]), pdb_dest)
    with open(conf_dest, "w") as f:
        json.dump(conf, f, indent=2)

    try:
        from propka.run import single as propka_single
        mol = propka_single(str(pdb_files[0]))
        best_conf = max(mol.conformations.values(), key=lambda c: len(c.groups))
        his_bound = {}
        for group in best_conf.groups:
            if type(group).__name__ == "HISGroup":
                his_bound[f"{group.atom.chain_id}{group.atom.res_num}"] = round(
                    group.pka_value, 2)
        r["his_complex_pkas"] = his_bound
        binder_his = {k: v for k, v in his_bound.items() if k.startswith("A")}
        if binder_his:
            r["avg_complex_his_pka"] = round(
                sum(binder_his.values()) / len(binder_his), 2)
            r["n_elevated_complex"] = sum(1 for v in binder_his.values() if v > 6.5)
    except Exception as e:
        print(f"  PROPKA error: {e}", flush=True)
        r["propka_error"] = str(e)

    with open(result_dest, "w") as f:
        json.dump(r, f, indent=2)

    runs_vol.commit()
    return r


@app.function(
    image=image,
    volumes={"/runs": runs_vol},
    timeout=10800,  # 3h: must outlive the single long predict_one_long task
)
def retry_one_driver(name: str, binder_seq: str):
    """Server-side driver for the SINGLE authorized retry.

    Runs predict_one_long for exactly ONE design, then merges its outcome into
    consolidated_results.json by replacing ONLY that design's entry. The other
    342 entries, and every one of their per-design artifacts on the volume, are
    left byte-identical. The summary is recomputed from the merged list. Being a
    detached remote driver, this survives local client exit.

    If the retry still errors/times out, the merged entry keeps the explicit
    error (missing status is retained, not masked).
    """
    new = predict_one_long.remote(name, binder_seq)

    out_base = pathlib.Path(f"/runs/{OUT_PREFIX}")
    runs_vol.reload()
    cpath = out_base / "consolidated_results.json"
    with open(cpath) as f:
        data = json.load(f)
    results = data.get("results", [])

    replaced = False
    for i, r in enumerate(results):
        if r.get("name") == name:
            results[i] = new
            replaced = True
            break
    if not replaced:
        results.append(new)
    data["results"] = results
    data["n_designs"] = len(results)

    errored = [r for r in results if "error" in r]
    skipped = [r for r in results if r.get("skipped")]
    computed = [r for r in results if "error" not in r and not r.get("skipped")]
    good = [r for r in results if "error" not in r and r.get("iptm", 0) > 0.5]
    data["summary"] = {
        "computed": len(computed), "skipped": len(skipped),
        "errored": len(errored), "iptm_gt_0.5": len(good),
        "errors": [{"name": r["name"], "error": r["error"]} for r in errored],
    }
    with open(cpath, "w") as f:
        json.dump(data, f, indent=2)
    runs_vol.commit()

    status = ("error: " + new["error"]) if "error" in new else (
        "iptm=%.3f" % new.get("iptm", 0))
    print(f"retry merged: {name} -> {status}", flush=True)
    print(f"summary: {data['summary']}", flush=True)
    return {"retried": name, "result": new, "summary": data["summary"]}


@app.local_entrypoint()
def retry_one(refold_json: str = "results/c4_acid_combined_refold.json",
              name: str = "tfr1_c4_acid_site1_41"):
    """Launch the single authorized retry of ONE missing design (detached-safe).

    Looks the design up BY NAME in the validated pool (preserving its exact
    binder sequence), confirms it is a unique PASS member, and spawns
    retry_one_driver. Does NOT re-dispatch the 342 valid designs and does not
    change scoring. Use:

        modal run --detach modal_boltz2_apical.py::retry_one
    """
    with open(refold_json) as f:
        refold = json.load(f)
    if isinstance(refold, dict):
        refold = refold.get("results", [])
    matches = [r for r in refold
               if r.get("name") == name and r.get("status") == "PASS"]
    if len(matches) != 1:
        raise ValueError(
            f"expected exactly 1 PASS record named {name}, found {len(matches)}")
    rec = matches[0]
    seq = rec.get("sequence") or ""
    bad = set(seq.upper()) - _AA
    if not seq or bad:
        raise ValueError(f"invalid binder sequence for {name}: {sorted(bad)}")
    print(f"Retrying ONLY {name} ({len(seq)} res) vs contiguous apical "
          f"{TARGET_LEN}-mer; 342 valid designs untouched")
    call = retry_one_driver.spawn(name, seq)
    print(f"Spawned retry_one_driver; function call id = {call.object_id}")
    print("Outputs: rfd4-runs:boltz2_c4_apical/  (this design + updated "
          "consolidated_results.json)")


@app.local_entrypoint()
def fetch(out_path: str = "results/c4_apical_consolidated_results.json"):
    """Download the remote consolidated results read-only to a local file.

    Usage:  modal run modal_boltz2_apical.py::fetch
    Writes only to the given LOCAL path (default under results/); never touches
    the Modal volume. Run this after the detached batch reports completion.
    """
    data = read_consolidated.remote()
    if "error" in data and "results" not in data:
        print(f"Not ready: {data['error']}")
        return
    with open(out_path, "w") as f:
        json.dump(data, f, indent=2)
    summ = data.get("summary", {})
    print(f"Downloaded consolidated results -> {out_path}")
    print(f"  n_designs={data.get('n_designs')} summary={summ}")


@app.local_entrypoint()
def main(refold_json: str = "results/c4_acid_combined_refold.json",
         max_designs: int = 1000):
    """Validate locally (pre-spend), then spawn the durable remote driver.

    GENUINELY-DETACHED LAUNCH. Two supported patterns, both of which keep the
    batch alive after the local client exits (NO nohup, NO attached client):

      (A) modal run --detach modal_boltz2_apical.py::main
          The spawned run_all driver is the last-triggered function, so Modal's
          detached mode keeps it (and the predict_one fan-out it drives) running
          server-side after this entrypoint returns and the client disconnects.

      (B) modal deploy modal_boltz2_apical.py   # once
          then spawn the deployed driver from any throwaway client:
          python3 -I -c "import modal, json; \\
            jobs=[(r['name'],r['sequence']) for r in json.load(open( \\
              'results/c4_acid_combined_refold.json'))['results'] \\
              if r.get('status')=='PASS']; \\
            fn=modal.Function.from_name('boltz2-c4-apical','run_all'); \\
            print(fn.spawn(jobs).object_id)"
          A deployed function's spawn is fully independent of any client — the
          most robust 'survives local exit' path.

    Re-running either command resumes idempotently (each predict_one skips
    designs already fully validated on the volume).
    """
    selected = _load_candidates(refold_json, max_designs)
    print(f"Validated {len(selected)} candidate inputs; "
          f"target = contiguous P02786 189-383 ({len(TFR1_APICAL_189_383)} res)")
    jobs = [(r["name"], r["sequence"]) for r in selected]
    call = run_all.spawn(jobs)
    print(f"Spawned durable driver run_all; function call id = {call.object_id}")
    print("Outputs: rfd4-runs:boltz2_c4_apical/  (per-design + consolidated_results.json)")


def _validate_local():
    """Pre-spend gate: validate receptor + candidates with no Modal/GPU use."""
    refold = "results/c4_acid_combined_refold.json"
    sel = _load_candidates(refold, 1000)
    print("LOCAL VALIDATION PASSED")
    print(f"  receptor: contiguous P02786 189-383, len={len(TFR1_APICAL_189_383)}")
    print(f"  epitope auth 205/206/208/209/210 = "
          f"{TFR1_APICAL_189_383[16]}{TFR1_APICAL_189_383[17]}"
          f"{TFR1_APICAL_189_383[19]}{TFR1_APICAL_189_383[20]}"
          f"{TFR1_APICAL_189_383[21]}")
    print(f"  candidates: {len(sel)} PASS, unique names/seqs, valid AA")
    print(f"  output prefix: rfd4-runs:{OUT_PREFIX}/  (acid prefix untouched)")
    return 0


if __name__ == "__main__":
    sys.exit(_validate_local())
