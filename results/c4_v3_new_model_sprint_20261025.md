# C4-V3 NEW-MODEL SPRINT — concrete shortlist, exact inputs, canary gates (3h window)

**AUTHORIZATION REQUIRED: 1 A40 pod, ~$0.60 total (approx 1 GPU-hour at $0.59/h). Nothing launched.**
2026-10-25 17:58 ET. Owner `01a170d2`. Supersedes the inventory half of `c4_v3_alt_model_plan_20261025.md`
(`2264285f`) — that plan's *shortlist* stands; two of its *blockers* were overstated and are corrected below.

**Current spend: `$0.00` this window.** Modal is at **~$19.6 / $20 (exhausted)** — wave 2 `ap-MQxEbfl`
is `stopped`, Tasks 0. **Modal is NOT a candidate.** RunPod: 8 pods, all `EXITED`, none billing.

---

## 1. The answer in one line

**AF2-multimer v3 on one RunPod A40 is the only run worth making:** it is the only installed,
already-validated, *genuinely different engine* that can model a complex. Input is **built and
byte-verified**. Canary 3 designs ≈ **8 min, ≈$0.08**; full 20 ≈ **50 min, ≈$0.49**.

## 2. Candidate comparison — what is actually available, measured today

| candidate | installed & reachable | complex-capable | different engine from Boltz-2 | ready now | verdict |
|---|---|---|---|---|---|
| **AF2-multimer v3** (ColabFold 1.5.5) | ✅ RunPod image + localcolabfold on pod | ✅ | ✅ **yes** | ✅ input + runner exist | **DO THIS** |
| Boltz-2 (extra seed) | ✅ Modal / local | ✅ | ❌ same engine | ✅ | **reject** — measures sampling noise, not model variance |
| ESMFold v1 | ✅ Modal `hf_cache` | ❌ **monomer only** | ✅ | ⚠️ | **reject for complexes**; valid only as binder-folds-alone control |
| Chai-1 | ❌ not installed anywhere | ✅ | ✅ | ❌ new image + weights + validation | **out of scope** (>3h, unvalidated) |
| Protenix / AF3-class | ❌ not installed; AF3 weights gated | ✅ | ✅ | ❌ licensing + image | **out of scope** |
| Boltz-2 on **correct** apical | ✅ done 2026-10-24 | ✅ | — | ✅ **already the V3 basis** | baseline, not a new model |

**Why not Modal:** the $20 sprint cap is spent (~$19.6). A 20-design AF2m batch on L40S at the measured
$2–3/design would cost **$40–60** and blow the cap ~2–3×. **RunPod A40 at $0.59/h is ~3.3× cheaper per
GPU-hour and ~34× faster per design in wall-clock** (validated: 116–160 s for all 5 models).

## 3. Exact inputs — verified, not assumed

| item | value |
|---|---|
| input CSV | `results/c4_v3_af2m_input_20261025.csv` — sha **`b7e05e7fe9198b62`**, **20 designs** |
| schema | `name,binder,target,source` — **identical to** `c3_screen_input_20261024.csv` (drop-in) |
| target | **195 aa**, `TFR1_APICAL_189_383`, **byte-identical for all 20 rows** (verified against `modal_boltz2_apical.py`) |
| binder lengths | 65–100 aa |
| complex total | **260–295 aa** (vs 157–211 for every prior AF2m run: **+40 %**) |
| ranking source | `results/c4_ranked_v3_apical_20261025/` (v3 apical ranking, 20 rows) |

### Two blockers from the 17:24 plan that are now RESOLVED or overstated

1. **"Target constant is hard-coded per runner" — mostly wrong.** In `modal_af2m_c3_complex.py` the target
   is a **CSV column** (`binder`/`target` per row). Only one line changes:
   `CSV_PATH = pathlib.Path(".../results/c3_screen_input_20261024.csv")` → the C4 file. **No code rewrite.**
   My CSV has the identical schema, so this is a one-line path swap plus a new app name.
2. **"`pod_run.sh` is on the pod, not in the repo" — real, and the repo-resident substitute needs ONE edit.**
   `scripts/c3_r3_pod_run.py` has the correct settings already (`multimer_v3`, `--num-models 5`,
   `--num-recycle 3`, `--rank multimer`, `mmseqs2_uniref_env`) and takes the CSV as `argv[1]`, so it is
   drop-in for C4. **But its `CUTOFF_ET = 1792911600` (2026-10-25 03:00 ET) is in the past**, so it would
   print `CUTOFF_NO_NEW_ADMISSION` and exit cleanly having done nothing. **This must be raised first** —
   it is step 0 in §7, and it is the single easiest way to lose an hour.

**Genuine remaining risk (unchanged):** sequences are **260–295 aa vs 211 max ever run on this path**
(+40 %). MSA depth/time may grow non-linearly. **This is exactly what the canary is for.**

## 4. The 3-hour plan

| step | what | wall | $ | gate |
|---|---|---|---|---|
| 0 | *(optional, $0)* ESMFold binder-alone on 20 | 10 min CPU/Modal-free | $0 | skip if tight |
| **1** | **Start 1 A40 pod** (restart `fak1k8jnuxwtc5`) | 2 min | — | pod reaches SSH |
| 2 | Deploy C4 CSV + runner, smoke 1 design | 5 min | ~$0.03 | finite iPTM; **chain B parses as exactly 195 aa** |
| **3** | **CANARY: 3 designs** (1 clean, 1 flagged, 1 long binder) | ~8 min | **~$0.08** | ↓ see gates |
| 4 | **FULL 20** only if canary passes | ~50 min | **~$0.49** | all 20 finite |
| 5 | Pull scores to Mini, compare vs Boltz-2 | 10 min | $0 | agreement table |
| 6 | **Stop pod** | 1 min | — | must be explicit |

**Total ≈ 75–90 min wall, ≈ $0.55–0.60.** Comfortably inside 3 h, with a 2× overshoot still under $1.20.

### Stop/go gates (pre-registered — decide before data exists)

**GO to full 20 only if ALL hold:**
1. All 3 canary designs return a **finite iPTM** (no NaN, no timeout).
2. **Chain B is exactly 195 residues** in each predicted PDB — proves the target swap took.
3. Binder chain matches its input sequence (hash or exact string).
4. Per-design wall **≤ 6 min** (else 20 designs exceed the window).

**STOP immediately if any of:**
- any NaN/absent iPTM, or chain B ≠ 195 aa → **the port is wrong, do not escalate**;
- OOM / MSA timeout on the canary → sequences too long for this path (record it, do not "try harder");
- projected total cost would exceed **$2.00** → stop, this is meant to be ~$0.60.

**Also stop the pod** (not just the job) at the end. An idle A40 is $0.59/h.

## 5. What this buys, stated honestly

V3's entire ranking rests on **one predictor**. AF2m gives the single most informative cheap fact
available: **does a second engine agree that these 20 designs bind the apical construct?** Disagreement
is actionable (down-weight, or flag which rows are model-dependent). Agreement is corroboration but
**not** proof of binding — this is still in-silico.

It does **not** create new designs and does **not** change C3. It is C4 evidence only.

## 6. The one thing I need from Ricardo

> **"authorize Trial 1"** (~$0.08, 3-design canary) — or **"Trials 1+2"** (~$0.60, canary + full 20).

Plus confirmation that the **2026-10-24 assistance-scope blocker** does not apply to a C4
correct-construct run. I am treating both as *not yet granted*, so **nothing is started**.

Until then, the $0 work continues: worker A runs disjoint infrastructure checks (pod image contents,
colabfold version, MSA cache presence) that de-risk step 2 without spending anything.

## 7. Exact commands (for the moment authorization lands)

```bash
cd ~/dev/adaptyv-c4

# 0. RUNNER PREREQUISITE (must be done first — see §9)
#    scripts/c3_r3_pod_run.py refuses every design because CUTOFF_ET=1792911600
#    (2026-10-25 03:00 ET) is already in the past. Raise it before use:
#      CUTOFF_ET = int(time.time()) + 3*3600        # or pass a new cutoff
#    Without this the runner prints CUTOFF_NO_NEW_ADMISSION and exits 0 with no work.

# 1. start the least-recently-used A40 (existing pod; do NOT create a new one)
python3 scripts/runpod_control.py --tool pod-action \
  --arguments '{"id":"fak1k8jnuxwtc5","body":{"action":"start"}}'

# 2. confirm colabfold + weights + GPU present (read-only, $0)
ssh fak1k8jnuxwtc5-644110ca@ssh.runpod.io \
  'nvidia-smi --query-gpu=name --format=csv,noheader; which colabfold_batch; ls /workspace'

# 3. stage input + runner (repo-resident, reproducible)
scp results/c4_v3_af2m_input_20261025.csv \
    fak1k8jnuxwtc5-644110ca@ssh.runpod.io:/workspace/c4_v3.csv
scp scripts/c3_r3_pod_run.py \
    fak1k8jnuxwtc5-644110ca@ssh.runpod.io:/workspace/c4_pod_run.py

# 4. CANARY — 3 designs only
ssh fak1k8jnuxwtc5-644110ca@ssh.runpod.io \
  'cd /workspace && python c4_pod_run.py c4_v3.csv \
     tfr1_c4_acid_site2_120,tfr1_c4_acid_site2_177,tfr1_c4_acid_site1_194 /workspace/c4canary'

# 5. GATE CHECK before the full batch (see §4 gates):
#    - all 3 scores.json have a finite "iptm"
#    - each rank-001 PDB has exactly 195 residues in chain B
#    - binder chain matches its input sequence
#    - per-design wall <= 6 min

# 6. full 20 — ONLY after the gates pass
ssh fak1k8jnuxwtc5-644110ca@ssh.runpod.io \
  'cd /workspace && python c4_pod_run.py c4_v3.csv ALL /workspace/c4full'

# 7. ALWAYS stop the pod (state transition stops billing; an idle A40 is $0.59/h)
python3 scripts/runpod_control.py --tool pod-action \
  --arguments '{"id":"fak1k8jnuxwtc5","body":{"action":"stop"}}'

# 8. verify it actually stopped
python3 scripts/runpod_control.py --tool get-pod --arguments '{"id":"fak1k8jnuxwtc5"}'
```

## 8. Status

**Launched: nothing. Spent: $0.00.** All 8 pods `EXITED`. Modal at cap, no tasks. No CSV altered,
nothing submitted, nobody contacted. Frozen artifacts untouched.
