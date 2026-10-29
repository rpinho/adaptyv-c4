# FINDING: the L6 sign convention survives an external test — but the zero threshold does not

Written 2026-10-26 from the real clock. CPU only, **$0**, read-only inputs. No GPU, no submission,
no organiser contact. No frozen file touched (4/4 verify). All four controls fired.

Mandated by **Fable 5**, PART 3 of the 2026-10-26 reviewer round:

> "Audit the L6 sign convention independently... Validate it against a reference pair whose
> direction is known independently of the scorer, not against the scorer's documentation. If the
> convention is the reverse of what rule v4 now assumes, the submitted 8 are drawn from the wrong
> 231 and the set must be rebuilt from the complement."

## Why the existing controls could not answer this

`scripts/c4_linkage_audit_20261025.py` carries three synthetic sign controls (null, pKa +1.0,
pKa −1.0). All three fire, and they fired again here. But they are generated from the same formula
they validate: they confirm the implementation matches its **own documented convention**. They
cannot detect the failure mode Fable 5 named — that the documented convention is mislabelled with
respect to physical reality. Only an external reference tests that.

## The external reference pair

Two complexes of the **same receptor** whose **relative** pH behaviour is established in the
iron-uptake literature, independently of anything in this repo:

- **Transferrin–TfR1** (1SUV, 3S9L). The transferrin cycle requires the Tf–TfR1 complex to survive
  the acidified endosome (~pH 5.6) and release at the cell surface at pH 7.4.
- **HFE–TfR1** (1DE4). HFE binds TfR1 at neutral pH; the interaction is disrupted at acidic pH.

A **relative** test was used deliberately. The absolute sign for a single crystal form depends on
the ionisation state captured and on PROPKA accuracy; the ordering is the robust claim.

    PREDICTION, if the documented convention (L > 0 == tighter at LOW pH) is correct:
        L6(Tf-TfR1)  >  L6(HFE-TfR1)

## Measured

| reference | L6 | titratable groups |
|---|---:|---:|
| Tf–TfR1, 3S9L | **+1.7288** | 45 |
| Tf–TfR1, 1SUV | **+1.3135** | 30 |
| HFE–TfR1, 1DE4 | **+0.2005** | 42 |

**Ordering holds.** min(Tf) = +1.31 > max(HFE) = +0.20.

## Verdict on Fable 5's flip condition: NOT met

The convention as applied in rule v4 is **consistent** with the external reference. The acid-OFF
gate is **not** selecting the wrong half. The candidate set does **not** need rebuilding from the
complement. Fable 5's stated flip condition does not fire.

## But a narrower problem appears, and it is new

**All three references score L6 > 0 — including HFE, the acid-DISRUPTED one**, which the documented
convention would place below zero. The scorer reproduces the correct **ordering** but not the
expected **absolute sign** for HFE. That is consistent with a **positive offset in L6**.

This matters because **rule v4 gates at exactly `L6 < 0`**. This audit validates the *direction of
the ordering*; it does **not** validate that zero is the correct threshold. If a positive offset
exists, the 111 acid-ON / 231 acid-OFF split is drawn in the wrong place — the direction would be
right while the eligibility boundary is wrong, changing which designs qualify without changing
which half is conceptually correct.

Calibrating that offset would need several references with known direction spanning both signs.
Three structures, one of them a single acid-OFF example, cannot do it. **Recorded as an open
question, not resolved.**

## Controls — 4/4 fired

| control | requirement | measured | fired |
|---|---|---|---|
| SYNTHETIC_SIGN | the scorer's own controls still hold | null 0.0, acid-ON +0.5918, acid-OFF −0.4539 | yes |
| SWAP_ANTISYMMETRY | L is a property of the complex, not of which chain is labelled binder | identical on relabel | yes |
| NONZERO_GROUPS | each reference scores a non-trivial number of groups | 30 / 45 / 42 | yes |
| SAME_RECEPTOR | references are all TfR1 complexes | 3 scored | yes |

**Correction made in place:** v1 of this audit reported SYNTHETIC_SIGN as failing. That was my bug —
the control read keys `null`/`acid_on`/`acid_off` while the scorer emits `null_L`/`acid_on_L`/
`acid_off_L`. The scorer's controls had in fact fired correctly all along. v1 preserved as
`results/c4_l6_sign_audit_20261026_v1_CONTROLKEY_SUPERSEDED/`. No reported number changed.

## What is NOT claimed

- **No binding, affinity, KD or pH-response is measured.** This audits a sign convention.
- **1SUV and 3S9L are holo (iron-loaded) Tf complexes, not apo-Tf.** The sharpest textbook direction
  (apo-Tf acid-ON) applies to the apo form. The relative contrast with HFE does not depend on
  resolving that, but the absolute values should not be read as apo-Tf behaviour.
- **PROPKA values are estimates**, and this repo records a prior PROPKA calibration failure
  (`726d9b5`). One reference pair is a direction check, not a calibration.
- **The generation-side defect is untouched.** Both RFD4 configs set `C_SEQ_PROB_HIS: 0.15` and aim
  every hotspot at a carboxylate; that bias is upstream of any gate and no sign audit repairs it.
