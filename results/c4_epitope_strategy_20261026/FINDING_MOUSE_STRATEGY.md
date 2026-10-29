# FINDING: mouse cross-reactivity is not capped by apical surface conservation

> **AMENDMENT, same day, from the native-context measurement that followed this one.**
> The headline below stands, but **the conservation numbers in this document are inflated by a
> burial artifact and must not be quoted on their own.** The designs bind predominantly to apical
> surface that exists only because the 195-aa domain was cut out of the 639-aa subunit, and buried
> residues are conserved *because* they are buried. Measured on the surface that is actually
> available in the native homodimer, apical mouse identity is **56.2%**, not 72–85%:
>
> | contacted apical positions | n | mouse-identical |
> |---|---:|---:|
> | all contacted | 148 | 72.3% |
> | **exposed in the native homodimer** | 73 | **56.2%** |
> | materially occluded in the native homodimer | 31 | 90.3% |
>
> The same split holds for the apical domain as a whole (80 natively exposed positions, 56.2%
> identical; 32 occluded positions, 90.6% identical), so this is a property of the domain, not of
> the designs. **High footprint conservation is therefore not by itself evidence of a
> cross-reactive epitope.** See `results/c4_native_context_20261026/FINDING_NATIVE_CONTEXT.md`.
> The "6 of 9" premise is still the wrong denominator, and the four-position avoid-list below is
> still correct; but the honest figure for available, conserved epitope is much lower than this
> document's unqualified medians, and the 14-design intersection reported at the end **collapses to
> 0** once native accessibility is required.

Written 2026-10-26 from the real UTC clock (see `epitope_strategy.json:created_at_utc`).
CPU only, $0, read-only on all inputs. No ranking produced. No submission, no organiser contact.

## The question asked

"Is the current design strategy (a single fixed scaffold family binding one human apical patch)
even capable of mouse cross-reactivity given that the patch differs at 6 of 9 residues?"

## The answer

**Yes, it is capable, and the "6 of 9" premise does not describe where these designs actually bind.**
The surface the 342-design pool actually contacts is 72–85% identical in mouse, and 15 designs have
footprints that are **100% mouse-identical under every optimal alignment**. Mouse cross-reactivity
in this programme is **unmeasured, not disfavoured**.

What limits it is not the surface. It is (1) the inverted pH gate documented in
`results/c4_independent_check_20261026/FINDING_L6_SIGN.md`, which steered selection *toward* the
one charge-reversing hot position, and (2) the fact that no mouse complex has ever been modelled.

## Where the "6 of 9" figure comes from, and why it is the wrong denominator

It derives from the five-residue hotspot list in
`results/mouse_epitope_conservation/RESULTS.md` — K205, N206, R208, L209, V210 — of which only
N206 and L209 are identical in mouse. That list was a *prior assumption* about the epitope.

The designs do not bind there. The twelve most-contacted apical positions across the pool are:

| human | aa | mouse | identical | % of pool contacting |
|---:|:--|:--|:--|---:|
| 298 | F | F | yes | 85.1% |
| 288 | F | F | yes | 78.4% |
| 316 | F | F | yes | 77.8% |
| 300 | H | H | yes | 75.4% |
| 295 | L | L | yes | 74.9% |
| **297** | **F** | **L** | **no** | 73.7% |
| **290** | **I** | **V** | **no** | 73.1% |
| 305 | T | T | yes | 73.1% |
| 263 | T | T | yes | 72.2% |
| 314 | P | P | yes | 71.9% |
| 264 | F | F | yes | 67.8% |
| 283 | M | M | yes | 67.8% |

None of K205/N206/R208/L209/V210 appears. The prior cycle's own mouse RESULTS.md had already
flagged this — *"the old small hotspot summary does not describe every actual v4 contact surface"* —
but the strategic framing still carried the 6-of-9 number. It should be retired.

## Measured conservation of the surface actually used

Numbering control fired: 195/195 chain B residues match human P02786 author 189–383.

| quantity | value |
|---|---:|
| apical positions (189–383) | 195 |
| identical in mouse in **every** optimal alignment | 146 (74.87%) |
| positions contacted by ≥1 design | 148 |
| of those, identical in mouse | 107 (**72.3%**) |
| "hot" positions (contacted by ≥10% of pool) | 59 |
| of those, identical in mouse | 50 (**84.75%**) |

Per-design footprint identity across 342 designs:

| statistic | value |
|---|---:|
| median | **0.845** |
| p10 / p90 | 0.731 / 0.929 |
| min / max | 0.143 / **1.000** |
| designs with a **fully** mouse-identical footprint | **15** |
| designs ≥ 0.90 identical | **81** |
| designs ≥ 0.80 identical | **264** |

The two design sites behave the same (site1 n=177 median 0.850; site2 n=165 median 0.838), so this
is a property of the apical surface, not of one scaffold sub-family.

## The substituted positions are a short, mostly conservative, avoidable list

Only 9 of the 59 hot positions are substituted, and most are conservative:

| human | aa → mouse | % pool | % of v4's 20 | character |
|---:|:--|---:|---:|:--|
| 297 | F → L | 73.7% | 95% | conservative hydrophobic |
| 290 | I → V | 73.1% | 85% | conservative hydrophobic |
| 296 | S → A | 64.9% | 75% | small polar → small apolar |
| **261** | **K → E** | 53.5% | **70%** | **charge reversal (+1 → −1)** |
| 222 | Y → F | 43.3% | 45% | conservative aromatic, loses OH |
| **292** | **N → E** | 31.6% | 40% | **neutral polar → acidic** |
| **217** | **G → E** | 22.2% | 25% | **glycine → bulky acidic** |
| 286 | T → N | 12.6% | 15% | polar → polar amide |
| **325** | **R → Q** | 12.6% | 10% | **loses positive charge** |

Four are chemically hard: **K261, N292, G217, R325**. The rest are substitutions a binder interface
can usually tolerate. This is an actionable avoid-list, not a wall.

## The inverted pH gate made the mouse problem worse

v4's seats are **worse than the pool they were drawn from** on footprint conservation:

| | v4's 20 seats | whole pool |
|---|---:|---:|
| median footprint identity | **0.8121** | 0.8450 |
| min / max | 0.6087 / 0.9524 | 0.1429 / 1.0000 |
| designs ≥ 0.90 identical | **4 / 20** | 81 / 342 |
| contacting K261 (→ E263, the charge reversal) | **14 / 20 (70%)** | 53.5% |

So ranking on `L6 > 0` also enriched for contact with the single worst mouse position. The prior
cycle observed the K261→E263 reversal in 14 of 20 and correctly declined to call it a lost salt
bridge; this shows *why* 14 of 20 had it — it is a consequence of the gate, not of the surface.

## What the pool already contains that v4's rule excluded

Quantification only. **This is not a ranking and must not be used as one** — the pre-registered
rule requires eligibility and model inventory frozen by sequence hash first, and its L6 line still
encodes the inverted direction.

| filter | n |
|---|---:|
| pool | 342 |
| footprint ≥ 0.90 mouse-identical | 81 |
| + iPTM ≥ 0.60 | 21 |
| + acid-OFF (L6 < 0, direction-correct) | **17** |
| + zero charge reversal anywhere in the footprint | **14** |
| *for contrast:* v4's actual rule (iPTM ≥ 0.60 **and** L6 > 0) | 40 |

There are **14 designs on disk** that are simultaneously direction-correct on pH, confident,
mouse-conserved in footprint, and free of any charge reversal at the interface. v4's gate admitted
**none** of them. The 15 fully-conserved-footprint designs are listed in `epitope_strategy.csv`;
13 of the 15 are acid-OFF.

## What would have to change in the design approach

**Not the scaffold family, and not the apical domain.** Both are adequate. Three changes, all
CPU-only, all on molecules already on disk:

1. **Fix the pH gate direction** (see the companion finding). This is prerequisite to everything.
2. **Add footprint mouse-conservation as an explicit selection axis**, and exclude contact with
   K261, N292, G217 and R325. The challenge scores mouse cross-reactivity as its own objective and
   nothing in v4 scored it.
3. **Model mouse directly** before any mouse claim. Everything here is sequence identity over a
   human-derived pose. A mouse complex prediction is the missing experiment, and it is the one place
   GPU would buy something real — which is a proposal for `results/decision_requests/`, not an
   action.

## What is NOT claimed

- **No binding, KD, affinity or cross-reactivity is measured.** This is sequence identity mapped
  onto human-derived contact footprints. A 100% mouse-identical footprint is a necessary-not-
  sufficient condition: identical residues can still present different rotamers, loop
  conformations, backbone geometry and glycans.
- **This does not overturn the prior cycle's pose audit.** That measured rigid-body transfer
  geometry and found 13/20 with persistent backbone overlaps — of which **12 were also confounded
  by the human-transfer control**, i.e. mostly artifacts of the transfer method rather than species
  effects. The two results are consistent: conservation is high, and the transferred poses conflict
  for reasons the transfer itself largely explains.
- **The 14-design intersection is not a candidate set.** Those designs are unscreened on clash,
  novelty, foldability, diversity and pairwise redundancy, and their iPTM values are mostly modest.
- **Conservation is measured over all six optimal BLOSUM62 alignments**; a position counts as
  identical only if identical in every optimum. The human 210–212 indel remains unresolved and
  ambiguous positions are reported separately, never resolved to improve a result.

## Reproduce

```sh
cd ~/dev/adaptyv-c4
date -u '+%Y-%m-%dT%H:%M:%SZ'
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python3 scripts/c4_epitope_strategy_20261026.py --out results/c4_epitope_reproduce_20261026
# the script exits non-zero if the 195-residue numbering control does not fire
```
