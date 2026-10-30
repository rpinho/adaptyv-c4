# C4 Full Re-score Progress

**Started:** 2026-10-27 00:53 UTC
**Last updated:** 2026-10-28 22:11 UTC — COMPLETE
**Construct:** P02786-1 homodimer (2x640aa) — the CORRECT epitope domain

## Pod Status

| Pod | ID | SSH | Shard | Designs | Completed | Status |
|-----|-----|-----|-------|---------|-----------|--------|
| 0 | ubrhji86dj3eqr | 194.68.245.63:22057 | shard_0 | 61 | 61/61 | STOPPED |
| 1 | v2d8dptjotvwoi | 194.68.245.124:22167 | shard_1 | 61 | 60/61 | STOPPED (1 MSA timeout: site1_112) |
| 2 | h9jmb5eiu2hzwr | 194.68.245.63:22165 | shard_2 | 60 | 60/60 | STOPPED |
| 3 | qnfur45hmosdqu | 194.68.245.63:22071 | shard_3 | 60 | 60/60 | STOPPED |
| 4 | 63j0e3dmecmaxh | 194.68.245.75:22030 | shard_4 | 50 | 50/50 | STOPPED |
| 5 | injkt4l8tmdeex | 194.68.245.7:22172 | shard_5 | 50 | 50/50 | STOPPED |

**Total designs:** 342 (61+61+60+60+50+50 across 6 pods, 1 missing: site2_96 — site1_112 was recovered by bg_recovery)
**Completed:** 341/342 (99.7%) — 342 score files pulled (site1_112 has scores.json from bg_recovery)
**Finished:** 2026-10-28 22:11 UTC (~45h wall clock)
**Budget:** ~$389.50 total ($189.50 + $200 top-up), ~$150 spent, all 6 pods stopped
**Budget after C4:** ~$190 remaining for C5
**v11 submission:** 18 designs (Kimi-approved backbone-diverse selection), ready for Oct 29 ~11:00 UTC window

## C4 Submission v10 (submitted 2026-10-28 ~11:00 UTC)

- 20 designs, top iPTM from 255 correct-construct scores
- iPTM range: 0.85-0.89, all passed novelty check 3/4+
- Submitted via Proteinbase, methodology included
- Next submission window: ~Oct 29 11:00 UTC (24h cooldown)

## Scale-up (11:00 UTC)

Added 2 more A40 pods ($0.59/h each) and redistributed remaining work across 6 pods:
- Took last 25 unscored designs from each existing shard → 100 designs for 2 new pods
- Reconciled run records on existing pods, restarted runners with patched code
- New pods: ColabFold 1.5.5 installed, cuDNN 8.9.7.29 + CUDA libs pinned to match working pods
- Kimi recommendation: keep 5 models (don't drop to 3) for ranking consistency
- Net cost increase: ~$2 more total for ~18h faster completion

## Score Distribution (44 designs scored)

- iPTM range: 0.80 – 0.85
- pTM range: 0.87 – 0.90
- pLDDT range: 88.08 – 90.52
- Top 5 by iPTM: site1_110 (0.85), site1_115 (0.85), site1_12 (0.85), site1_131 (0.85), site1_10 (0.84)
- Scores pulled to `results/c4_full_rescore_20261027/pulled_scores/`

## Issues and Fixes

### cuDNN incompatibility (fixed 02:02 UTC)
Pods 1-3 had JAX falling back to CPU ("Unable to load cuDNN") due to cuDNN 9 being incompatible with JAX 0.4.23. Fix: `pip install nvidia-cudnn-cu12==8.9.7.29 chex==0.1.86`. All 3 pods now at 100% GPU utilization. ~1h of CPU-only compute was wasted.

### Pod 0 stopped by Codex canary monitor (fixed 02:04 UTC)
The Codex background task (bqsk0y5iw) stopped pod 0 after canary completion, killing the shard_0 launch. Pod restarted via API, but workspace was wiped (no network volume). Reinstalled ColabFold with cuDNN fix included. New SSH port: 22057 (was 22146).

### site1_100 FileNotFoundError (recovered 04:01 UTC)
ColabFold completed all 5 models for site1_100 but crashed during the final PDB rename (exit_code=1). All ranked score files and PDBs were present. Manually recovered scores and updated run record.

### Pod 1 zombie colabfold process (killed 04:01 UTC)
Old pre-cuDNN-fix colabfold process (PID 817, from 00:53) was still running on CPU, consuming 18GB RAM and 600% CPU. Killed it. Pods 2 and 3 also had zombies (PID 809); killed at 04:34 UTC.

### Systematic ColabFold 1.5.5 rename bug (mitigated 04:35 UTC)
ColabFold crashes during final PDB rank-rename with FileNotFoundError (batch.py:533). Affects every design — all 5 models run successfully, scores and PDBs are saved with rank tags, but the exit code is 1. Fix: deployed `bg_recovery.py` daemon on all pods (runs every 2 min, auto-converts error.json → scores.json for designs that have ranked score files). Also patched runner to always attempt score extraction regardless of exit code.

### site1_112 MSA timeout (pod 1, not recoverable)
ColabFold timed out downloading MSA from api.colabfold.com. No models ran, no output. Will need re-run in a later pass.

### New pods cuSOLVER mismatch (fixed 10:58 UTC)
Pods 4-5 conda installed a newer ColabFold solve (JAX 0.4.25 + CPU-only jaxlib). Fixed by pinning to JAX 0.4.23 with CUDA jaxlib and installing matching nvidia-cusolver-cu12, nvidia-cublas-cu12, etc.

### New pods dm-haiku version mismatch (fixed 11:38 UTC)
Pods 4-5 had dm-haiku==0.0.17 which requires `jax.core.take_current_trace` (not in JAX 0.4.23). All 50 designs on each pod failed in 4-5s with AttributeError. Working pods had dm-haiku==0.0.10. Fix: `pip install "dm-haiku==0.0.10"`. Cleared all errors, reset run records, restarted runners and bg_recovery daemons. ~1h of pod time wasted ($1.18).

## Canary Scores (6/6 PASS)

| Design | iPTM | pTM | pLDDT | Wall (s) |
|--------|------|-----|-------|----------|
| tfr1_c4_acid_site1_103 | 0.8100 | 0.8800 | 88.7 | 2664 |
| tfr1_c4_acid_site1_115 | 0.8500 | 0.9000 | 90.0 | 2492 |
| tfr1_c4_acid_site1_13 | 0.7900 | 0.8700 | 89.2 | 3182 |
| tfr1_c4_acid_site1_151 | 0.8200 | 0.8900 | 88.9 | 2733 |
| tfr1_c4_acid_site2_76 | 0.8300 | 0.8900 | 89.3 | 2380 |
| tfr1_c4_acid_site2_99 | 0.8800 | 0.9200 | 91.2 | 2287 |

## Timeline Estimate (revised 11:41 UTC)

- 295 remaining designs across 6 pods
- Pods 0-3 actively folding; pods 4-5 restarted at 11:38 UTC after haiku fix
- Effective throughput: ~6/hr now, rising to ~9/hr once pods 4-5 reach steady state (~12:30 UTC)
- Estimated full completion: ~Oct 29 02:00 UTC (±6h)
- Budget: ~$155 total estimated, within $189.50 ceiling
- Hard cutoff: Oct 31 00:53 UTC (T+96h, per Kimi plan)
