# C4 corrected-apical — retry + artifact comparison, 2026-10-24

Author: Claude (Opus). Scope: one bounded retry of the single missing C4 design, then
an artifact-based comparison of candidate portfolios on the corrected contiguous apical
target. No Adaptyv submission, no external contact, no candidate regeneration, no change
to scoring. Historical artifacts preserved (hashes below). The separately-running C3 full
AF2M screen (`ap-kyeVNsRjGqksuudoNSczUN`) was not touched.

## 1. State at start

- `results/c4_apical_consolidated_results.json`: 343 entries — 341 skipped + 1 computed =
  **342 valid**, 1 errored (`tfr1_c4_acid_site1_41` → `timeout`). Matches the volume:
  `rfd4-runs:boltz2_c4_apical/` held 342 `_result.json`, none for `site1_41`.
- No active C4 app. All `boltz2-c4-*` apps (`ap-tNtBKluXq1k8gupDcvSSie`,
  `ap-bRyOn8T0naTy0TFg2X5odJ`, `ap-vIswggqjslZar8Gzod4Nt1`, `ap-YrMCFoNiXtXF8XrjjHtCGb`,
  `ap-Y95jGjM5vBmGc6ihEOo59l`) were **stopped** at start. Only the C3 screen was live.

## 2. Root cause of the lone timeout + bounded fix

`tfr1_c4_acid_site1_41` has a **100-residue binder — the maximum in the 343 PASS pool**
(pool length min/median/max = 60/80/100). Paired with the 195-residue apical target it is
the **largest complex in the batch (295 res)**. It hit the `predict_one` 1800 s subprocess
wall on both the initial run and the earlier idempotent retry (`ap-bRyOn8T0naTy0TFg2X5odJ`).
The failure is a **wall-clock budget** problem, not a scoring, input, or target-geometry
problem.

Operational fix (root-cause, single longer task — **no indiscriminate whole-pool rerun**):
added `predict_one_long` + `retry_one_driver` + `retry_one` to `modal_boltz2_apical.py`.
`predict_one_long` is byte-for-byte identical to `predict_one` in target, binder, and
Boltz-2 parameters (`diffusion_samples=1`, `recycling_steps=3`, `--use_msa_server`) — the
**only** differences are the time budget (subprocess wall 1800 s → 6600 s; function timeout
3600 s → 7200 s) and `max_containers=1`. `predict_one` itself was left unmodified so the 342
existing artifacts keep their exact provenance. The driver runs the one prediction, then
merges its outcome into `consolidated_results.json` by replacing **only** that one entry and
recomputing the summary; the other 342 entries and all their per-design artifacts are left
untouched. If the retry errors/times out again, the merged entry keeps the explicit error
(missing status retained, not masked).

### Exact commands

Local pre-spend validation gate (no GPU; passed — 343 PASS, contiguous P02786 189-383
195-mer, epitope auth 205/206/208/209/210 = KNRLV):
```bash
cd .
/usr/bin/python3 modal_boltz2_apical.py
```

Single authorized detached retry (binder preserved exactly, same contiguous apical target):
```bash
cd .
modal run --detach modal_boltz2_apical.py::retry_one
```

Check / download (read-only):
```bash
modal app list | grep x8qUME
modal volume ls rfd4-runs boltz2_c4_apical/ | grep site1_41
modal run modal_boltz2_apical.py::fetch
```

### Retry run identifiers

- **App ID: `ap-x8qUMEomla0Dhcf7x8Zzb8`** (`boltz2-c4-apical`, ephemeral/detached).
- Durable driver call id: **`fc-01M5QM4HY72W1JDNVQPVZFZQ6E`** (`retry_one_driver`).
- Launched 2026-10-24 ~14:53 EDT; verified ephemeral with 1 live task after the local
  client exited (survived client disconnect, unlike the first attached attempt).
- Launch log: `results/c4_apical_retry_20261024.launch.log`.

### FINAL MISSING STATUS

**RESOLVED: FAILED AGAIN — `tfr1_c4_acid_site1_41` IS EXPLICITLY MISSING — 2026-10-24 22:44 ET.**

The single authorized retry (`ap-x8qUMEomla0Dhcf7x8Zzb8`, function call
`fc-01M5QM4HY72W1JDNVQPVZFZQ6E`) did not produce the artifact. Verified read-only:

- `modal app list` no longer lists the app (it is gone, not running).
- `modal volume ls rfd4-runs boltz2_c4_apical/` returns **1027 entries and contains no
  `tfr1_c4_acid_site1_41_*` file** (`_result.json`, `_confidence.json` and `_complex.pdb` all absent).
- `results/c4_apical_consolidated_results.json` (mtime 17:52) still records 342 valid designs.

Per the standing rule, this is the **one and only** authorized retry and it failed:
**C4 coverage is 342/343 and no further C4 compute is authorized.** The missing design is
retained as explicitly missing — it is not imputed, backfilled, or represented by a neighbour.
Everything downstream already operates on the 342 valid designs, so none of the portfolio
comparisons or the preview ranking changes.

## 3. Artifact-based portfolio comparison (on 342 valid designs)

All 343 designs are now scored against the **same corrected contiguous apical target**
(P02786 189-383, 195 res), so their Boltz-2 metrics are mutually comparable. The two
historical 20-row portfolios are both subsets of this pool, which makes the comparison
apples-to-apples on the corrected target.

### Comparability caveat (binding — do NOT overclaim)

These are **domain-only** predictions of a binder against the **monomeric** 189-383 apical
domain. The official C4 assay construct is TfR1 **121-760, homodimer**. Nothing here
validates native-homodimer binding, acid-triggered release in the real assay, or biological
improvement. The only defensible claim is **relative ranking on corrected-apical-domain
Boltz-2 metrics (complex ipTM + complex His pKa)** — a model-metric reranking against a
cleaner target than the discontinuous surrogate the historical CSVs were selected on.

### Historical portfolios score poorly on the corrected target

| portfolio | mean ipTM | mean conf | mean complex His pKa | ipTM>0.5 | His pKa>6.5 |
|---|---|---|---|---|---|
| `submission_v9.csv` (original) | 0.418 | 0.756 | 5.80 | 8/20 | 0/20 |
| `submission_v9_preview.csv` (historical) | 0.482 | 0.753 | 6.00 | 8/20 | 1/20 |
| **`c4_apical_preview_20261024.csv` (NEW)** | **0.704** | **0.800** | **6.38** | **20/20** | — |

- Pool (342 valid): ipTM min/med/max 0.119/0.507/0.942; 177/342 ipTM>0.5; 116 ipTM>0.6.
  Complex His pKa min/med/max 3.90/5.95/7.21; only 15/342 >6.5.
- v9 member ranks within the pool by ipTM run 49→342 (median ~200); preview 29→339.
  **None of the corrected-target top-20 by ipTM appear in either historical portfolio** —
  expected, because those portfolios were selected against the 154-res discontinuous
  7-segment surrogate, not this contiguous apical domain.

### NEW portfolio — selection logic (acid-OFF aware)

The C4 objective is a **pH-dependent (acid-OFF) binder**: bind at pH 7.4, release at pH 6.0.
That needs **both** a plausibly-formed complex **and** an elevated binder-His pKa (the
protonation switch in the 6.0-7.4 window). Pure ipTM ranking ignores the switch; pure pKa
ranking selects designs that barely form a complex (the 15 pKa>6.5 designs are mostly
ipTM<0.45). So the NEW preview applies:

1. **Gate:** corrected-apical complex ipTM ≥ 0.60 (confident complex; 116 designs qualify).
2. **Rank:** descending complex His pKa (acid-OFF switch signal).
3. **Diversity:** greedy 4-mer-Jaccard < 0.60 filter to avoid near-duplicate backbones.
4. Take the top 20 (8 site1 + 12 site2; His pKa 6.26-6.62; ipTM 0.60-0.78).

File: **`results/c4_apical_preview_20261024.csv`** (NEW; 20 unique names/sequences,
`molecule_class=single_chain`). It shares **zero** designs with v9 or the old preview.

Builder scripts (local, no GPU): `/tmp/c4anal/{compare,portfolio,build,write_csv}.py`
(scratch — copy into `scripts/` if a reproducible artifact is wanted).

## 4. Recommended file ranking

1. **`results/c4_apical_preview_20261024.csv`** — strongest on both comparable
   corrected-apical axes (complex ipTM and acid-OFF His pKa), diversity-filtered.
   Recommended basis for any C4 submission **pending Ricardo's review** and subject to the
   domain-only / non-native-homodimer caveat above. NOT submitted.
2. `results/c4_rebuild_20261024/submission_v9_preview.csv` — historical; mid-pool on the
   corrected target.
3. `results/submission_v9.csv` — original; weakest on the corrected target. **Preserved
   byte-for-byte; still the only file that was ever a real submission baseline.**

This ranking is a model-metric recommendation only. It does not establish experimental
performance, platform acceptance, or native-target validity. No file was submitted or
promoted; no one was contacted.

## 5. Integrity (preserved artifacts)

- `results/submission_v9.csv`: `31d98a4133c558e1a3a9f582a607383bca8d7aad1b8bd820e684cb8fa231eadd` (unchanged)
- `results/c4_rebuild_20261024/submission_v9_preview.csv`: `00ec3d599b50224f6d88942e599150d66301406e6a8a72a225be209a5a180cb6` (unchanged)
- `results/c4_apical_preview_20261024.csv` (NEW): `f31774b1c789be9b3cbb17ead78b89ede18388d979b5abc136f59176dd571ab7`
- `modal_boltz2_apical.py`: retry functions appended; `predict_one` unmodified. New SHA
  recorded in the trailing status block once the retry resolves.

---

## TRAILING STATUS — retry outcome (updated when the detached driver resolves)

RESOLVED — retry FAILED (2026-10-24 22:44 ET). `tfr1_c4_acid_site1_41` remains explicitly
missing; C4 final coverage is **342/343**. The app is no longer listed and the volume holds no
`site1_41` artifact. No further C4 compute was or will be authorized for this design.
No CSV was written, promoted or submitted; the three C4 files keep their recorded SHA-256 values.
