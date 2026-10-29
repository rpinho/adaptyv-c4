# C5 compute notes — 2026-10-27

Planning input for the KRAS G12D pMHC challenge. Constraints and their provenance are in
`../C5_CONSTRAINTS_20261027.md` (page check 2026-10-27 19:26 UTC). Nothing here is a design
decision; those go to Kimi / Fable / Opus 5.5 / Astra.

## 1. What makes C5 structurally more expensive than C4

C4 scored each design against **one** target. C5 measures each design against **three**, in the same
conditions:

1. KRAS G12D complex — HLA-A\*11:01 + b2m + `VVVGADGVGK`
2. Wild-type complex — same HLA + b2m + `VVVGAGGVGK`
3. Peptide-free — same HLA + b2m, empty groove

Two of those are negative objectives that **rank ahead of affinity**. So the in-silico readout has to
discriminate, not just rank affinity, and the discriminating signal is one side chain (Asp vs Gly at
peptide position 6) on an otherwise byte-identical ~380-residue complex. That means three predictions
per design and tighter error bars than C3 or C4 ever needed.

## 2. The MSA amortisation — the main cost lever

Established with the `adaptyv-c4-71` session, 2026-10-27: in the current C4 pipeline the **MSA cache
is per-design, fetched from `api.colabfold.com` on each run**. That is what makes C4's shard
re-splitting safe (no locality penalty), but carried into C5 unchanged it is wasteful.

**In C5 the target side is fixed across the entire campaign.** The heavy chain (281 aa) and b2m
(99 aa) are identical in all three states; the only variable is a 10-residue peptide taking exactly
two values or absent. So there are **three distinct target-side inputs for the whole challenge**, not
three per design. Compute them once, keep them, reuse them for every design.

Only the designed binder chain varies per design, and for de novo binders single-sequence mode is
standard, so the per-design MSA cost should be close to zero rather than 3x.

Sequence hashes for cache keys (computed 10-27):

| chain | len | sha256 |
|---|---|---|
| HLA-A\*11:01 heavy chain | 281 | `77f7825c0d673572048dadac819f4a73166acc3e1f1fda5ea1b441f247c9ba19` |
| beta-2-microglobulin | 99 | `19916695f0fea2713c1617e79f23e054164c645cadbd8c368d7a9f44fdd20624` |
| peptide G12D | 10 | `6909a59fdaec58efd61388dbac885325f6bc4038d8f5324418fb816791e405af` |
| peptide WT | 10 | `c04b93c3d3845086801893b7e2a7af4a67d96ede6225e7b7cee15b0093449b6c` |

Secondary benefit: at C5 funnel scale, repeatedly hitting the free public ColabFold server is a
rate-limit risk and a fair-use problem, and it is already demonstrably unreliable on the critical
path — the C4 run has one permanent MSA failure (`site1_112`). Pre-computing the three target-side
inputs removes nearly all of that exposure.

## 3. Cost estimate, and the measurement that would sharpen it

C4's **measured** cost is **$1.83-2.79 per design** (brief, 10-26), and the observed throughput on the
10-27 re-score fleet is **~1.23 designs per pod-hour** on an A40 — about 49 minutes of wall clock per
design.

Order-of-magnitude for C5: a 200-design funnel evaluated across three states lands around
**$1,100-1,700**. That is the number to agree before any launch, not after.

**The measurement that would sharpen it, and should be taken first:** instrument a single design and
split that ~49 minutes into (a) MSA fetch/wait and (b) GPU compute. If a large share is MSA wait,
then amortising the target-side MSAs per section 2 cuts the three-state penalty well below 3x and the
estimate above is too pessimistic. If it is almost all GPU, the estimate stands. Either way this is
one design's worth of spend to size a four-figure decision.

## 4. The gate before any GPU spend

Target prep and control validation are CPU-only and cost nothing:

1. Build the three states from PDB 7OW6 — chains A-C for G12D, chains A-C with D6G for wild type,
   chains A-B for peptide-free.
2. Assert the three differ **only** as intended (one residue between states 1 and 2; one chain
   between states 1 and 3). A diff that shows anything else means the prep is wrong.
3. **Validate the readout on a positive control before generating anything.** 7OW6 contains a natural
   binder to this exact complex in the TCR that the organisers removed. Check what selectivity the
   literature actually reports for it, then confirm the differential readout recovers that preference.
   A pipeline that cannot separate G12D from wild type on a real TCR will not rank designs either.

Only after step 3 discriminates does GPU spend make sense. Per the standing rule: written expected
gain, a cap, a stop condition, and review before any launch. The brief's warning applies — prior C4
GPU work burned $1.70-2.70 and returned zero usable rows.

## 5. Capacity reality as of 2026-10-27 19:50 UTC

From a live read-only `list-pods`: 6 A40s at $0.59/h, all EU-SE-1, all on the C4 re-score. Real rate
~1.23 designs/pod-hr uniformly (the two low-count pods started 9 h 45 m later, they are not slow).
Job ETA Oct 29 ~04:30 UTC, or Oct 28 ~21:00 UTC if the site2 tail is re-split across all six.

**Finishing C4 costs ~$109-116 of the ~$123 RunPod balance.** Freed C4 pods are therefore not spare
capacity for C5 — they are $0.59/h of new spend against a near-empty balance. C5 requires its own
top-up and cap from Ricardo, or it does not run.

Infra note for C5 (not worth retrofitting to C4 mid-run): a network volume mounted across the pods
would let the MSA cache and the three target-side inputs survive pod restarts and be shared, which
also removes the upload step from any re-sharding. The two existing volumes from the old C3 pods
(`i4f8nnqs0k`, `1l49pa109j`) are **not ours** and Ricardo said to ask before touching them, so this
means a new volume and his approval.
