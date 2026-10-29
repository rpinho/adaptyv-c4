# C5 provenance & methodology record — TEMPLATE

Fill this in **as work happens**, not afterwards. One file per submitted batch; never edit a
filled copy after it is referenced by a manifest (write a new version instead).

- `batch_id`:
- `created_utc`:
- `orchestrator_session`:
- `competition_track` (1 / 2 / 3 — **unconfirmed until Ricardo confirms**):
- `slot_limit` (20 or 40 — follow the track, never assume):
- `decision_record` (path to the Kimi/Astra decision authorising this batch):

## 1. Target inputs (verbatim, hashed)
| item | value | sha256 |
|---|---|---|
| target sequence (as distributed) | | |
| target numbering (first..last, offsets) | | |
| provided structure file (if any) | | |
| assay construct / oligomeric state | | |
| conditions (pH, temperature, buffer) | | |
| required objective (binding / inhibition / agonism / pH-switch) | | |
| epitope or hotspot constraints | | |

Validator: `scripts/c5_validate_submission.py` must record the same sequence hash.

## 2. Method provenance, per design family
| field | value |
|---|---|
| method (e.g. RFdiffusion / BindCraft / RFD4-Proteina) | |
| version / commit | |
| weights file + sha256 | |
| seeds (exact list) | |
| number of backbones generated | |
| sequence-design method (e.g. ProteinMPNN) + version + weights sha256 | |
| sequences sampled per backbone, temperature | |
| refold / validation model + version + settings | |
| GPU type, count, wall time, cost | |
| de novo / zero-shot: how verified | |

**De-novo rule:** no starting from an existing binder and no modification of a known binder.
Record *how* this was checked (novelty scan, known-binder FASTA exclusion, template absence).

## 3. Selection and ranking
| step | rule | threshold | survivors |
|---|---|---|---|
| hard filters (length, canonical residues, class) | | | |
| structural validation (e.g. refold pLDDT / RMSD) | | | |
| interface validation (e.g. iPTM) | | | |
| diversity (cluster rule + threshold) | | | |
| final ranking key | | | |

State explicitly which metric decides order and which are only reported.

## 4. Files
| file | sha256 | rows | rank note (is this the best version, and why) |
|---|---|---|---|
| | | | |

Rule: every new candidate set is a **NEW versioned file**. Never overwrite a prior version.
One-line ranking of which file is best and why goes in `results/overnight_log.md` too.

## 5. Validation evidence
| check | result | artifact path |
|---|---|---|
| port validated against a known result (same design, two environments) | | |
| construct/sequence identity gate (byte-exact) | | |
| validator pass on final CSV | | |
| reviewer (human) sign-off | | |

## 6. Known limitations (state plainly)
- Domain-only vs native construct:
- Monomer vs oligomer:
- Predicted vs experimental template:
- Metric that is *not* biological evidence:

## 7. Sign-off
- Researcher review completed by:            (name, UTC timestamp)
- Explicitly **not** done: submission, upload, account creation, external contact.
