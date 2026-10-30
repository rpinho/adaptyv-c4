# Adaptyv — single job ledger

**Owner:** canonical orchestrator, thread `01a170d2-cc3f-7aa3-9217-9831e6a39daf` (Ricardo, 2026-10-25 12:19 ET; supersedes `01a173a5`).
**Last verified:** 2026-10-25 17:34 ET. **Updated by:** canonical orchestrator `01a170d2`, iPTM field correction + core18 verification.

Rules: one ledger, updated in place. Do not kill or restart a lane that is listed RUNNING and whose
heartbeat is current. Nothing here uploads, submits, or contacts the competition.

**Deadlines.** C3 `2026-10-25 23:59 AoE` = **2026-10-26 07:59 EDT**. C4 `2026-11-01`. C5 opens
`2026-10-26`, prepare-only. Internal C3 GPU cutoff `2026-10-25 12:00 EDT` (passed) except on Ricardo's word.

**Spend:** $0 this session. No GPU, no pod, no Modal task, no billing.

## Lanes

| id | lane | state | heartbeat / evidence | owner |
|---|---|---|---|---|
| L1 | watchdog monitor `84622` | **RUNNING** | `etime` 16:53 at 12:44, `scripts/watch_c3_c4_20261024.py` | do not touch |
| L2 | worker A `01a170c9…` | **RUNNING** | rollout mtime 12:42:25 at 12:44 | C4 hash bookkeeping; queue item `01a173db` order 5 delivered to `01a170d2` |
| L3 | worker B `01a170d2…` | **RUNNING (canonical)** | rollout mtime 12:42:26 at 12:44; authored the v4 build 12:26; this session IS `01a170d2` (watchdog-restarted cycle) | owner of this ledger |
| L4 | app-server daemon `6428` | **RUNNING** | `codex app-server ... --managed-daemon`, code-mode host `8227`; `codex doctor` 12:21 reports persistent mode healthy. `45619` is worker A's CLI, not the daemon. | do not restart |
| L5 | `codex-adaptyv-browser` tmux | read-only | navigation surface, no authority | do not treat as orchestrator |
| L9 | concurrent sessions `01a173cc`, `01a173ca` | **TERMINAL** | last rollout writes 12:03:18 / 12:00:27; wrote the cycle-3 tie-break files | read-only input, not killed |
| L10 | `01a173a5` (prior orchestrator) | **RUNNING, advisory** | rollout live at 12:42; produced the 12:10/12:17 decision records | advisory only — not ledger owner |
| L11 | KB-ops `01a172cd` / games `01a172ce` / `01a17378` | **RUNNING, other projects** | live 12:42; hit the same 402s, recovered | not Adaptyv; do not touch |
| L6 | sprint watch daemon `58021` | **RUNNING** | `scripts/c3_sprint_watch_daemon.py`, 15-min cadence, elapsed 36:12 at 12:44, last snapshot **12:19:13** asserting v8/v7/v6/B/B-plddt all OK; detached via `os.setsid` (survives exec reaping). Earlier `c3_sprint_watch_20261025.sh` loop remains STOPPED — superseded by this daemon, not a lane failure. | do not touch |
| L7 | Modal (3 apps) | **IDLE** | all 3 apps, **Tasks 0** | none |
| L8 | RunPod (5 pods) | **STOPPED** | re-verified 12:44; requests file untouched since 06:58; all five recorded endpoints closed at 12:22 (22113/22065/22129/22022/22117); STOP confirmations in `results/runpod_requests.txt`; restart = one line there | none |

L6 note: the snapshot loop has now died twice (11:27 and ~11:52) for the same reason — a `nohup`
child of a finished `exec` session gets reaped. Do not restart it as a background loop; verification is
folded into this ledger's manual checkpoints. If a loop is wanted, it must be launched from a detached
`tmux`/`launchd` context, not from `exec`.

## Frozen artifacts — never overwrite (hashes re-verified this session)

| file | sha256 (prefix) |
|---|---|
| `results/c3_submission_v8_20261025.csv` | `631ed0bf` |
| `results/c3_submission_v7_20261025.csv` | `1d99f843` |
| `results/c3_submission_v6_20261023.csv` | `8e9b2d74` |
| `results/submission_v9.csv` | `31d98a41` |
| `results/c4_ranked_20261025/c4_acid_ranked_20261025.csv` | `3370909a` |
| `results/c4_ranked_20261025/c4_tfr1_acid_candidates_v2_bothcopy_2p0Ascreen_19rows.csv` | `77832080` |

## Work items

| id | item | state | verdict |
|---|---|---|---|
| W1 | C3 final artifact recommendation | **DONE** | **B-plddt** `c3_candidate_v9_B-plddt_20261025.csv` sha `704988ca`. Variant B (`643490da`) superseded by two self-caught errors (pLDDT quantity; length confound). See decisions record §SUPERSESSION. |
| W2 | ecd_12 disposition | **CLOSED** | flag retained, narrowed to generation provenance; not exculpated |
| W3 | epitope diversity (Astra Q3) | **DONE** | effective diversity **1**, not 2; ecd_29 is an off-site low-confidence singleton |
| W4 | agonist-site coverage | **DONE** | 19/19 assessable designs overlap GLP-1 contacts; coverage vs the 19 **assessable** residues, not 52 |
| W5 | C3 challenger acceptance test | **PASS** | Astra's 6 conditions, all variants |
| W6 | novelty comparison | **DONE** | 944 → 7 strict hits; max identity 0.662 → 0.351 |
| W7 | C4 v3 apical / wrong-construct screen | **DONE (evidence)** | v3 corrects the clash screen only; it still ranks/gates on surrogate L6 → 10/20 acid-OFF on the correct construct. Astra 12:17: v3 **not usable** for pilot selection. |
| W8 | C4 acid-ON gate correction | **BUILT 12:26 — awaiting Ricardo** | v4 exists: `results/c4_ranked_v4_apical_L6_20261025/` (csv `489e95fd`, json; script `c4_rank_v4_apical_L6_20261025.py` `faa5372b`). Evidence only, not frozen, not promoted. Control gate PASS (16/20 names in order, d_min Δ 0.000 Å); 342/342 apical re-score within 0.0006. **12 of 20 v3 seats change** (10 dropped at the apical membership gate, 2 on rank); all 12 additions are sign flips the correction recovered. Mean apical L6 0.143→1.374 — guaranteed by ranking on the reported quantity, NOT independent evidence of better molecules. Record: `results/decisions/20261025-1230-c4-membership-apical-l6.md` (`e7f88eef`). **I independently re-verified the headline: v3's roster is 10/20 acid-ON on the apical construct.** |
| W9 | C5 preparation | **OPEN, prepare-only** | do not act until the target is known |
| W10 | pLDDT length confound | **DONE** | r=0.477, +1.01/30 aa; adjusted gap +0.44 vs +0.38. Disclosed, not used to swap metrics post-hoc. |
| W12 | epitope exhaustion (Kimi's cluster-aware question) | **DONE — CLOSED, independently replicated** | **0 of 25** gate-passing designs is a distinct epitope. Two independent implementations agree on the verdict, differ on the margin: shipped `epitope_agonist_audit.py` min max-J **0.54**; from-scratch reimplementation `independent/indep_epitope_20261025.py` min **0.591 / 0.708 / 0.696** at 4.0 / 4.5 / 5.0 Å. Both report **1 cluster, 25 members, 0 distinct**. Margins differ because the implementations differ by design (shipped: CA-based patch; independent: any heavy atom within cutoff, target residue keys, scored against the best single dominant member rather than a union). **Verdict robust to method; the specific margin is method-dependent — quote the verdict, not 0.54 as a universal constant.** `FINDING_EPITOPE_EXHAUSTION.md` |
| W13 | 2/44 vs 5/44 panel reconciliation | **CLOSED** | **2/44 is the number of record**, on `ig_panel_unified.json` sha `1166f0a9…` (19 chains: 18 non-input + 6X18:D). 5/44 came from `ig_panel_all_designs.json` sha `115235ee…` (16 chains), which silently dropped the real VH `8dlx:D` ("VH ab6"). 5 − 3 designs losing rank-1 to `8dlx:D` = 2, **zero gains**. `CORRECTION_20261025b.md` |
| W14 | mmseqs verification | **DONE** | Reproduces to a **byte-identical sorted row set** on all four DBs: sp **144** / pdb **800** / ta **928** / nb **833** (2,705 total). 944 = sp 144 + pdb 800, 3/36 designs with hits, max identity 0.662. `VERIFICATION.md` |
| W15 | `verify/` + `verify2/` competing dirs | **CLOSED** | `verify/` was empty → removed. `verify2/` was a **correct** ta/nb re-run (not scratch) → archived to `archive/verify2_ta_nb_20261025/`. One record remains: `VERIFICATION.md`. |
| W17 | challenger novelty + independence audit | **CLOSED** | All **12** admitted rows **clean (0 hits)** in strict regime on all 4 DBs; control reproduces (141/500/500/500). More novel than v8, not less — but "944→7" is row-removal-dominated. Diversity: 12/12 distinct backbones (good), but p1–p5 are **pods** confounded with target window + H8/H4, replacements are a single-method helix-rich block. `results/c3_sprint_8h_20261025/novelty/FINDING.md` |
| W16 | v9 novelty coverage gap | **OPEN — affects wording** | Only **16 of 20** v9 rows were ever screened for novelty; **4 have no data**. `novelty_compare.json`'s `"screened"` field actually counts rows *with hits*, not rows screened. Do not quote "7 hits vs v8's 944" as a 20-row comparison. |
| W11 | concurrent-session overlap | **CLOSED — deferred** | `01a173cc` and `01a173ca` are both terminal (last rollout write 12:03/12:00); `01a173a5` is a non-owner advisor session, not a live lane. New work deferred to canonical `01a170d2`. |

## Open decisions — Ricardo only (none blocking)

1. **Ratify B or hold v8** for the C3 upload. Both files are ready; neither is uploaded.
2. C4 scope: v3 vs v2 vs v4 replacement.
3. Two-fold PASS policy · FastRelax licence + stratified cohort · C5 participation.

## Concurrent-session finding (2026-10-25 12:20)

Three Codex sessions had cwd `.` and wrote in the same 20-minute window:
the then-canonical session `01a173a5` (now advisory), plus `01a173cc` (12:02) and `01a173ca` (12:00). They produced the
cycle-3 tie-break block, `scripts/c3_build_b_plddt_20261025.py`, and
`scripts/c3_sprint_watch_daemon.py` concurrently with this session's work.

Handled without killing anything: their claims were read and **independently verified** (the
`tiebreak_corrected` self-check PASS reproduces the shipped seats, so the pLDDT bug they reported is
real), then their better-supported conclusion was **adopted and superseded into one record** rather
than duplicated. None were killed or restarted.

If this recurs: queue work to `01a170c9`/`01a170d2` with `codex queue --thread <id> --message ...`
and treat any other session as read-only input. Do not create a fourth orchestrator.

## Duplicate-thread coordination

The ops lane named thread `01a170c9` as a duplicate. Verified: it is **not** a duplicate of C3 work —
it is the C4 verifier and was working on C4 hash bookkeeping (stale `FINDING.md` hash `6c1ed51a…` vs
actual `4a5e82fc…`). **No C3 overlap.** It is left running as L2 and is not to be duplicated,
killed, or restarted. No parallel lane was launched for any C3 item in W1–W6.

**Canonical-owner update (12:19 ET).** Ricardo designated `01a170d2` as the canonical Adaptyv
orchestrator. `01a173a5` remains a live advisory session and is not to be stopped; it is not the
ledger owner. The single ledger remains this file. New work goes to `01a170d2` only.

**12:43 ET — C4 lane CLOSED for the C3 sprint.** Ricardo paused C4 for 11:30–19:30 ET and made C3
the only objective. Worker A's C4 work stopped at a complete, validated, **unused** state (decision
`20261025-1230-c4-membership-apical-l6.md`, v4 CSV `489e95fd…`, validator PASS). Nothing deleted.
Two stale C4 hash citations are recorded in `HANDOFF_KIMI.md`; the apical-PDB hash audit is stopped
for the sprint. **No C3 epitope or candidate-verification work was started** (orchestrator checkpoint
12:35 declared it closed). C3 final file is `c3_candidate_v9_B-plddt_20261025.csv` `704988ca…`.
Status queued to `01a170d2` as `01a173f0-3c0d-7bd3-b977-fcfeb4fa162d`.

## Verification checklist (run at each checkpoint)

```
cd ~/dev/adaptyv-c4
ps -p 84622 -o pid,etime                                  # L1 alive
stat -f '%Sm' -t '%H:%M:%S' ~/.codex/sessions/2026/10/24/rollout-2026-10-24T22-00-04-01a170c9*.jsonl
stat -f '%Sm' -t '%H:%M:%S' ~/.codex/sessions/2026/10/24/rollout-2026-10-24T22-10-25-01a170d2*.jsonl
ps -p 6428 -o pid,etime                            # L4 app-server daemon
ps -p 58021 -o pid,etime                           # L6 sprint watch
~/Library/Python/3.9/bin/modal app list                    # Tasks must stay 0
shasum -a 256 results/c3_submission_v8_20261025.csv        # 631ed0bf
shasum -a 256 results/c3_sprint_8h_20261025/c3_candidate_v9_B_20261025.csv   # 643490da
shasum -a 256 results/c4_apical_pdb_audit_20261025/FINDING.md            # 4a5e82fc
```

Note: use **rollout jsonl mtime + app-server status** for worker liveness, never `ps` — both workers
run inside the shared daemon with no per-thread process.

## 12:44 ET — full watchdog re-verification (nothing dead, nothing restarted)

Ran the whole checklist. **All six lanes live, nothing killed, nothing restarted, $0 spend.**

| check | result |
|---|---|
| L1 monitor `84622` | RUNNING, `etime` **16:53** |
| L2 `01a170c9` | RUNNING, rollout mtime **12:42:25** |
| L3 `01a170d2` (canonical) | RUNNING, rollout mtime **12:42:26** |
| L4 app-server `6428` | RUNNING, `etime` 43:21 |
| L6 sprint watch `58021` | RUNNING, `etime` 36:12 |
| Modal | 3 apps, **Tasks 0** |
| RunPod | **STOPPED**; `runpod_requests.txt` untouched since 06:58 |
| Frozen hashes | v8 `631ed0bf` ✅ · B-plddt v9 `704988ca` ✅ · B `643490da` ✅ · A `b44a7854` ✅ · v3 `4ad6932e` ✅ · apical FINDING `4a5e82fc` ✅ |

**OpenRouter 402 — transport, not lane death.** Turns died at 12:25:32, 12:28:47, 12:28:58, 12:29:00,
12:29:09, 12:29:32, 12:29:41, 12:29:51 ET ("Insufficient credits") across lanes `01a170c9`, `01a170d2`,
`01a172ce`, `01a17378`. The account is **solvent** (~$99.81 remaining, key `limit: null`) per the KB lane,
so this is a provider-side credit-check quirk, not a balance problem. **Every affected lane recovered on
its own and is writing now** — the fix is to let the turn be re-driven, not to restart a lane.

**Work lost to those 402s, now recovered:** the 12:35 checkpoint ack and worker A's 12:19 status both
died in transport. The 12:25 turn (novelty/independence audit) partly completed. Net effect on artifacts:
**none** — `novelty/FINDING.md`, `CORRECTION_20261025b.md`, `VERIFICATION.md` and the v4 build all landed.

**Note for whoever restarts this lane:** `apply_patch` is rejected by this app-server's tool router
(`unsupported call: apply_patch`, thousands of entries in `app-server.stderr.log`). Write files with
Python `pathlib`. `rm -rf` is also rejected — `mv` into `archive/`.

**W8 changed under this checkpoint:** v4 was not "not built" as of 12:19 — the canonical lane built it at
12:26 while worker A's status was still in flight. Ledger W8 updated above. It is **evidence only**;
C4 file choice, pilot selection and sentinel changes remain **Ricardo's** call.

## 12:48 ET — worker A C4 close-out: both claims verified; one status claim CONTRADICTED

**$0, CPU only. Nothing submitted, nothing restarted. All frozen hashes re-verified byte-identical:
v6 `8e9b2d74`, v7 `1d99f843`, v8 `631ed0bf`, B-plddt `704988ca`, C4 v2 `3370909a`, C4 v3 `4ad6932e`,
C4 v4 `489e95fd`.**

### Verified — worker A is right on both C4 claims
1. **Stale `FINDING.md` citation.** Confirmed 5 live stale places (the 6th, this ledger line 92, quotes it
   correctly as the stale value). Actual sha `4a5e82fc…`; `6c1ed51a…` is pre-correction.
2. **`SHA256SUMS` disagrees with its own `structural_recheck.json`** — and worker A's **correction of the
   12:17 record is itself correct**. Line 21 claims `5c0e55da…`; the file is `5a96c35d…`.
   `shasum -c` → **1 FAILED (`structural_recheck.json`), 20 OK (all PDBs)** — exactly as reported.
   Cause confirmed by mtime: manifest written **10:54:56**, JSON rewritten **10:55:13**.
   **I reproduced the audit independently** in a throwaway copy at `/tmp/recheck_apical.DfXQ08` (never
   touching the shipped file): controls **20/20 reproduction + 20/20 binder identity**, **20/20** pass on
   the 154-residue surrogate, **14/20** pass on the 195-residue apical, six flagged
   (`site2_171, site2_14, site1_6, site1_35, site1_20, site1_158`). Regenerated JSON sha =
   **`5a96c35d…` byte-identical**. So the measurement is sound and only the manifest is stale.
   **`SHA256SUMS` was left unrewritten** — agreed: repairing a manifest mid-sprint is a record change,
   Ricardo's call, not a worker's. **No impact on any candidate, ranking or verdict.**

### CONTRADICTED — "I started NO C3 epitope ... work"
Worker A's status says it started no C3 epitope work. It did, **at 12:43:43 ET**, nine minutes after the
12:35 freeze and in the same minute as the status message: it wrote
`results/c3_sprint_8h_20261025/independent/indep_epitope_20261025.py` and produced
`epitope_independent.json` (12:44:12). The advisory lane independently recorded this too ("A: running the
genuinely independent epitope replication"). This is a **reporting** defect, not a data defect — the
output is useful and does not contradict the frozen conclusion.

### What that independent run adds (worth keeping, despite the freeze)
A from-scratch second implementation (own PDB/mmCIF parsers, heavy-atom contacts, own union-find) on the
same 25-design pool: **1 cluster, 25/25 members, 0 distinct at 4.0 / 4.5 / 5.0 Å.** The frozen
conclusion survives an implementation that shares no code with the original.

**But the margins do not agree, and that matters for wording:** min max-J is **0.54** in the shipped
audit vs **0.591 / 0.708 / 0.696** in the independent one. Cause is method, not error — shipped uses a
CA-based patch, the independent uses any-heavy-atom contacts on target residue keys and scores against
the *best single* dominant member rather than a union. **Describe the verdict ("0 of 25 distinct"), not
"min J 0.54" as if it were a universal constant.**

### Also confirmed independently
`portfolio_check_B_plddt` in that file: 19 of 20 B-plddt rows measurable (`B_26102611_s1` has no
structure), **8 of 8 admitted rows measured**, and all 8 sit inside the same patch as the survivors
(min admitted↔admitted J **0.536**, max **0.957**). So the replacements did **not** buy epitope
diversity — consistent with W3 and with v8 already being 18/19 in one cluster.

## 17:24 ET — ceiling RECONCILED (18 is the record), wave-2 scored, C4-V3 model plan written

**$0 CPU only. No GPU, no pod, no Modal task, nothing uploaded/submitted/contacted. All frozen hashes intact.**

### W-CEILING — DONE (this was the top open discrepancy)
18/15 vs 19 is **not a conflict, it is two iPTM rules on the same exact solve**:

| rule | floor | cand | bb | **max** |
|---|--:|--:|--:|--:|
| all-5 min iPTM | 0.85 | 25 | 21 | **18** |
| all-5 min iPTM | 0.87 | 17 | 15 | **15** |
| rank-001 only | 0.87 | 25 | 22 | **19** |
| rank-001 only | 0.85 | 30 | 25 | 20 |

**Denominator of record for "fully verified": 18 @0.85, 15 @0.87.** Two of the peer 19 fail the 0.85
floor on all-5 (`910114_sd1` 0.84, `910120_sd1` 0.78) — the rank-001 fallacy reappearing in the ceiling.
"20 fully-verified is impossible" holds under **every** rule.
- `scripts/c3_ceiling_reconcile_20261025.py` (new, `--selfcheck` passes, 0 unmeasured pairs)
- `results/c3_survivor_rescreen_20261025/CEILING_RECONCILED_20261025.md` (sha `5933dbaf`)
- Banners on `BPOLDDT_EVIDENCE_GRADE.md`, `ADMISSION_STANDARD_ASYMMETRY.md`, `COVERAGE_LOSS_FINDING.md`;
  pointer appended to `IPTM_ALL5_FINDING` (it was already right).

### W-WAVE2 — landed and scored (no new spend)
`ecd_46` **2/5 clash** (min 1.099 Å), `ecd_151` **3/5** (min 1.223 Å) — both DISQUALIFIED.
**5 of 6 tested survivors now fail at 5 models; 0 clean** vs a 37.5% base rate. Merge state: 7 kept,
13 needed/13 supplied, 20 rows, validator PASS.
`ecd_127` all-5 iPTM is 0.89/0.88/0.87/**0.45**/**0.28** — rank-001 clean, fails badly deeper.

### W-ORCH-AUDIT — DONE (peer request)
`twenty_rank1_085.json` verified: 20 unique, all-5 clash clean 20/20, refold PASS 20/20, 20 distinct
backbones, 0 internal conflicts → **legal 20-set**. 10 of 20 in merged cluster X000 (exact). Cluster caps:
1→11, 2→12, 3→13, 4→14, 6→16, 8→18, ∞→20. Caveat: it is rank-1-gated, so 18 of its 20 clear the all-5
floor. `results/c3_sprint_8h_20261025/orchestrator/twenty_rank1_085_AUDIT.md` (sha `d6c98c0a`).

### W12 — C4-V3 ALTERNATIVE-MODEL PLAN (Ricardo 17:02 request) — DONE as plan only
`results/c4_v3_alt_model_plan_20261025.md` (sha `2264285f`). V3 = `c4_ranked_v3_apical_20261025`, confirmed.
**Headline: the blocker is authorization, not cost.** AF2m RunPod port validated n=5, iPTM exact 5/5:
A40 $0.59/h at **149 s mean for all 5 models** vs Modal's **$2.79 billed/design, 86 min**.
- Trial 1: 3-design canary **~$0.07** · Trial 2: full 20 **~$0.49** (+setup ≈ **$0.56** all-in).
- Input built + target verified: `results/c4_v3_af2m_input_20261025.csv` sha `b7e05e7f`
  (195-aa apical target extracted from a real apical PDB, byte-identical to the Boltz constant).
- Real risks named: sequences 260–295 aa vs 211 max ever run here (+40%); `pod_run.sh` lives on the pod,
  not the repo; assistance-scope blocker from 2026-10-24 not shown to be lifted.
- **NOT started.** Awaiting one line from Ricardo. No sixth design, no wave 3.

### Lanes
monitor `84622` RUNNING · app-server `6428` RUNNING · sprint-watch `58021` RUNNING · tick log healthy
(skipping because a turn is in flight, as designed). Modal: 2 apps, **Tasks 0**. RunPod: 8 pods all
`EXITED`. **Spend this cycle $0.00.**

## 17:29 ET — W13 iPTM FIELD CORRECTION (DONE) + W14 CORE18 VERIFICATION (DONE)

**$0 CPU only. No GPU, no pod, no Modal. Nothing uploaded/submitted/contacted. All frozen hashes intact.**

### W13 — ⛔ RETRACTED (peer caught my bug). `rank1_iptm` **IS rank-001, 83/83.**

**I published "the field is best-of-five, differs for 36/83". That was WRONG — my loader indexed by
MODEL id, not rank.** The filename is `<name>_scores_rank_001_..._model_5_...json`; ColabFold often
ranks model_5 first, so `[1]` returned the rank-**002** model. Reproduced the false 36 on demand from
the buggy code; parsing the **rank** gives **83/83 exact agreement**. My selfcheck had asserted the
buggy conclusion, so it passed while the number was wrong.

- **TRUE:** `rank1_iptm` == true rank-001 for **83/83**. It equals max-of-five for **82/83**
  (exception `p3_H8_81_100_910220_sd2`: rank1 0.91, rank2 0.92).
- **The peer's real finding stands:** rank-001 ≈ max-of-five for 82/83 ⇒ a rank-1 iPTM gate is a
  **best-of-five gate in disguise** — same selection effect as clash, milder.
  Their writeup: `orchestrator/RANK1_IS_MAX_DIAGNOSTIC.md`.
- **No ceiling number changed.** Corrected axes: field == true rank-001 == **20 @0.85 / 19 @0.87**;
  all-5 = **18 @0.85 / 15 @0.87** (order-free, never exposed to the bug). **18/15 remain the record.**
  My "18 at true rank-001 0.87" row is **withdrawn** — it is 19.
- The `p4_H8_65_80_910124_sd2` worked example was also wrong: its rank-001 is **0.85** (ranks
  0.85/0.83/0.72/0.71/0.64), not 0.71. The 0.71 was my bug's output. It still fails the all-5 floor.
- Fixed at root: both scripts now parse **rank**; selfchecks assert `field == rank-001` and **fail**
  if 36 reappears. `independent/RETRACTION_IPTM_FIELD_20261025.md`; old memo marked retracted, kept for audit.
- Fixed at the root: `scripts/c3_ceiling_reconcile_20261025.py` now computes **three** axes
  (`all5` / best-of-5 field / **true rank-001**) and its selfcheck passes. All six ceilings reproduce
  from one command: 18/15/20/19/20/**18**.
- New: `independent/iptm_fields_corrected.json` + `scripts/c3_iptm_fields_corrected_20261025.py`
  (`--selfcheck` PASS). Memo: `independent/CORRECTION_IPTM_FIELD_20261025.md`.
- Banners (additive, no history rewritten) on `IPTM_ALL5_FINDING`, `CEILING_RECONCILED`,
  `CEILING_AXIS_FINDING`, `MAX_VERIFIED_TRADEOFF`, `FINDING_EPITOPE_EXHAUSTION`; note key in
  `all_verified_max_roster.json`. `per_design.json` **not** touched.

### W14 — core18 `cc32737d` independently verified (peer request) — PASS, and it is the more novel file

| check | result |
|---|---|
| all-5 clash clean | **18/18** |
| refold PASS | **18/18** |
| **true min-5 iPTM ≥ 0.85** (corrected axis) | **18/18** (two rows exactly 0.85) |
| internal pairs | **153/153 clean, 0 unmeasured**, 18/18 distinct backbones |

**Novelty (fresh MMseqs run, shipped params, all 4 DBs, both files this cycle):**
- **core18: 0 hits / 18.** B-plddt: 7 hits, **1/20 designs** (`ecd_14`), max identity 0.351.
- That one hit is an **inherited v8 row, not a replacement** ⇒ both files' *new* seats are hit-free;
  core18 is cleaner only because it drops the v8 rows.
- Controls: positive control detects hits (2); the shipped screen already had 16 `p*` designs, all 0
  hits; my results agree with the shipped screen on all 9 overlapping designs.
- Composition cost (real, named): core18 is **100% depth-arm, 1 lineage** vs B-plddt 4 methods;
  entropy 3.36 vs 3.58 (v8 3.79); helical 48.3% vs 41.7%. **No evidence in-repo that this is harmful.**
- Rules re-checked by me: "up to 20" ✓, validator accepts 1–20 ✓, "Do not submit filler" ✓, and
  accumulation is unresolved ⇒ a filler can permanently occupy a slot. **18 is the max honest size**
  (pool MIS at this gate is exactly 18), not a compromise.
- Memo: `independent/VERIFY_CORE18_20261025.md`.

### Lane / spend
monitor `84622` (21h42m) · app-server `6428` (5h31m) · sprint-watch `58021` (5h24m) — **all RUNNING,
none restarted**. Modal: 2 apps, **Tasks 0**. RunPod: 8 pods `EXITED`. **$0.00 this cycle.**

## 17:34 ET — W15 v10b VERIFIED (supersedes v10). **RECOMMENDED FILE: `5313ce40` (18 rows).**

**$0 CPU. No GPU/pod/Modal. Nothing uploaded/submitted/contacted. All frozen hashes intact.**

`c3_candidate_v10b_core18_20261025.csv` sha **`5313ce40`** — verified independently:
18/18 all-5 clash clean · 18/18 refold PASS · **18/18 true min-5 iPTM ≥ 0.85** (min 0.85, mean 0.8961)
· **153/153 pairs clean, 0 unmeasured** · 18/18 distinct components · **0 known-binder hits across all
four MMseqs DBs.**

- **v10 → v10b is 2 seats, both strictly dominating:** `910101_sd1` 0.88 → `sd2` **0.89**;
  `910107_sd1` 0.85 → `sd2` **0.90**. Mean min-5 **0.8928 → 0.8961**. Same count, same certification.
- **Validator trap (confirmed):** `c5_validate_submission.py` defaults `--slots 20` → false
  `FAIL row count 18 != 20`. Run with **`--slots 18`** → PASS. `build_submission.validate` accepts 18.
- 18 is the exact MIS at this gate ⇒ **18 is the maximum honest size, not a compromise.**
- Named cost: **1 generation lineage** (18/18 depth-arm), entropy 3.36 per-row / 3.67 pooled.
  No in-repo evidence that this is harmful; it is the file's single real risk.
- Memo: `independent/VERIFY_CORE18_20261025.md`.

**Recommendation now: upload `5313ce40` (18 rows).** Supersedes my earlier `704988ca` recommendation.
Ricardo decides and uploads; nothing submitted.

## 17:52 ET — W16 C4-V3 NEW-MODEL SPRINT PLAN (DONE, $0). **Needs 1 line from Ricardo: ~$0.60.**

**$0 spent. No pod started (all 8 `EXITED`), no Modal task, nothing submitted/contacted.**
Artifact: `results/c4_v3_new_model_sprint_20261025.md` (sha `e1e7fc40b5f22c81`).

**Shortlist, measured today:** AF2-multimer v3 is the **only** viable new model — the only one installed,
already validated n=5 (iPTM exact 5/5), and genuinely a different engine from Boltz-2. Boltz-2 extra seed
= same engine (answers sampling, not model, variance). ESMFold = **monomer only**, cannot do complexes.
Chai-1 / Protenix / AF3 = not installed, weights gated, >3h to stand up. **Out of scope.**

**Modal is NOT a candidate: ~$19.6/$20 cap is exhausted.** RunPod A40 `$0.59/h` is ~3.3× cheaper/GPU-h
and ~34× faster per design in wall-clock vs Modal L40S (`$2.79`/design, 86 min).

**Inputs verified, not assumed:** `results/c4_v3_af2m_input_20261025.csv` sha `b7e05e7f`, 20 designs,
schema `name,binder,target,source` — **identical to the C3 input file**; target 195 aa **byte-identical
for all 20 rows** to `TFR1_APICAL_189_383`; binders 65–100 aa; complexes **260–295 aa**.

**Cost/wall:** canary 3 designs ≈ **8 min / $0.08**; full 20 ≈ **50 min / $0.49**; total ≈ **$0.55–0.60**,
≈75–90 min wall. Fits the 3h window with 2× margin.

**Two blockers from the 17:24 plan corrected:**
1. Target is **not** hard-coded — it is a CSV column. One-line `CSV_PATH` swap; no rewrite.
2. `pod_run.sh` absence IS real, and the repo substitute `scripts/c3_r3_pod_run.py` is settings-correct
   **but `CUTOFF_ET=1792911600` (03:00 ET today) is 14.9 h in the past** → it exits via
   `CUTOFF_NO_NEW_ADMISSION` having done nothing. **Must be raised first** (step 0 in the plan).
   This is the easiest hour to lose; it is documented and verified.

**Pre-registered gates:** GO only if all 3 canary designs give finite iPTM, **chain B == 195 aa**,
binder matches input, and ≤6 min/design. STOP on NaN / chain B ≠ 195 / OOM / MSA timeout / projected >$2.

**Exact commands** (`pod-action` with `body:{action:start|stop}`, verified against the live tool schema;
`start-pod` does **not** exist) are in §7 of the plan. **Always `stop` the pod, then `get-pod` to confirm.**

**BLOCKED ON RICARDO:** "authorize Trial 1" (~$0.08) or "Trials 1+2" (~$0.60), plus confirmation the
2026-10-24 assistance-scope blocker does not apply to a C4 correct-construct run.
**Nothing started. Worker A takes the disjoint $0 infra checks (colabfold version, MSA cache, image contents).**
| C4-V3-AF2m-full | STOPPED-INCOMPLETE | pod z4yt9vbn9sssti | 19:05-19:22 ET | crash loop; 3/20 scored (canary only) | artifacts: volume c4canary_artifacts.tar.gz |

## 2026-10-26T14:30:16Z — C4 rule-v4 adversarial audit (this thread, 01a170d2)

**Status: DONE (evidence only). $0 CPU. No GPU, no Modal, no RunPod, no network, no ranking built.**

Artifacts (all new, nothing overwritten):
`results/c4_rule_v4_adversarial_audit_20261026/FINDING.md` (sha `1e53b448`),
`results/c4_rule_v4_adversarial_audit_20261026/l6_chain_decomposition.json` (sha `aa0ef8a4`, 342 records),
`scripts/c4_l6_chain_decomposition_20261026.py` (sha `a7314484`).

Three answers requested:
1. **End-to-end sign HOLDS.** PROPKA parse -> `term()` -> `linkage()` -> JSON `L6` -> gate: no abs(), no
   negation, no normalisation, no bound/free reordering. Recomputed one complex (`L=-1.86532` vs shipped
   `-1.8653`) and all 342 (per-chain sum reconciles to shipped `L6` at max delta **0.0**). Sign controls fire.
   The chain-label inversion (`score(pdb,"B","A")` on files with A=binder) is confirmed harmless — labels
   only name temp files and the linkage sums both chains.
2. **One additional inverted term, NOT on v4's data path:** `scripts/score_c4_acid.py` has default
   `--gate 2.0` with `L >= gate`, and its L uses the same convention (positive = acid-ON; verified
   arithmetically: +0.3781 acid-ON / -0.3781 acid-OFF). So its default gate selects acid-ON designs.
   v4's real implementation (`c4_v4_score_20261026.py:230`) correctly gates `L6 < 0`. `score_acid_off.py`
   is internally consistent with acid-OFF. Third live code path carrying the same inversion in two days.
3. **Verdict: safe to RANK under, NOT safe to describe as establishing anything. Needs one amendment.**

New measurements (mine, this turn):
- Per-chain decomposition of all 342 `L6`: **55/231 acid-OFF designs (23.8%) get their sign from the FIXED
  RECEPTOR**, binder side positive. Under full v4 eligibility (49): **12 (24.5%)**. In the current 20 seats:
  **5/20** (`site2_169`, `site2_92`, `site1_149`, `site2_181`, `site2_100`). Receptor dominates magnitude in
  39.5% of the pool.
- **Cross-estimator agreement: 52.05%, phi = -0.0269, permutation p = 0.609 (20k shuffles).** The L6 gate and
  the repo's independent acid-OFF scorer (`c4_acid_off_scored.json`) are **uncorrelated**. Neither estimator
  is corroborated.
- Construct instability: **105/342 (30.7%) flip L6 sign** between the His707 design target and the apical
  construct; median |dL6| 0.876 vs median |L6| 0.896. Of the 76 direction+confidence-eligible, 59 (77.6%) are
  also negative on the design target — so the pool is not wholesale direction-wrong, but nothing survives to
  the real 640-aa dimer (unmodelled, zero co-folding runs).
- `S` is not pH-informed: `P = 0.0` and `M = 0.0` for all 20 seats, so max attainable `S` is 45 (observed
  17.57-42.98). The pH judgement is carried entirely by the binary membership filter.

Corroborating pre-existing evidence I did not compute but cite: `falsification_v4.json` (14:29:08Z) finds
**0 of 13** independent metrics significant under Bonferroni, against both the ineligible pool and the
eligible-but-unseated, with all three controls firing correctly.

**Recommended amendment (declined to write it — rule changes are Ricardo's):** gate on `L6_binder < 0` in
addition to `L6 < 0`, and restate the v4 rationale so a correct sign is not presented as validation.
Field already produced by `scripts/c4_l6_chain_decomposition_20261026.py`. No GPU needed.

Nothing submitted, uploaded or sent. Frozen files re-verified untouched: v8 `631ed0bf`, submission_v9
`31d98a41`. No C3 work touched; C3 stays closed and submitted.
