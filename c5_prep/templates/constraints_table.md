# C5 constraints table — FILLABLE TEMPLATE

Fill **only** from the official C5 page and the public selection prompt, after the Oct 26
announcement. Every row needs a source (URL + fetch timestamp) and a snapshot hash. An empty
row is a **blocker**, not a default. Nothing here may be inferred from C3/C4.

## 0. Snapshot provenance
| item | value |
|---|---|
| challenge page URL | |
| fetched at (UTC) | |
| page snapshot file | |
| page snapshot sha256 | |
| selection prompt URL | |
| prompt snapshot file + sha256 | |
| fetcher (read-only) | |

## 1. Target
| field | value | source |
|---|---|---|
| Target name / accession | | |
| Exact target sequence (verbatim) | | |
| Sequence sha256 | | |
| Numbering: first / last / offsets | | |
| Provided structure file (+ sha256) | | |
| Provided structure is experimental or predicted | | |
| Construct boundaries (as distributed) | | |

## 2. Assay and objective
| field | value | source |
|---|---|---|
| Assay type | | |
| Oligomeric state of the assayed construct | | |
| Conditions (pH / temperature / buffer) | | |
| Epitope or region to engage (if specified) | | |
| Objective (binding / inhibition / agonism / pH-conditional / other) | | |
| How the objective is ranked | | |
| Anything explicitly out of scope | | |

## 3. Submission schema
| field | value | source |
|---|---|---|
| Slot limit (20 vs 40) | | |
| Track (1 / 2 / 3) — **confirm with Ricardo** | | |
| Min / max length | | |
| Allowed molecule classes (exact strings) | | |
| Required columns and ordering | | |
| Fab representation rule | | |
| Deadline (with timezone) | | |
| Testing allocation for this track | | |

## 4. Rules
| field | value | source |
|---|---|---|
| De novo / zero-shot restrictions | | |
| Known-binder restrictions | | |
| External data / model restrictions | | |
| Prior-design reuse rules | | |
| Disclosure terms | | |
| Account requirements | | |

## 5. Reuse assessment (answer only after 1–4 are filled)
| question | answer | evidence |
|---|---|---|
| Is any existing C3 artifact a genuine de novo design against *this* target? | | |
| Is any existing C4 artifact usable under the rules? | | |
| Does any prior artifact violate "already submitted" / uniqueness? | | |
| Smallest compliant experiment | | |

## 6. Gate
- [ ] All of §1–§4 filled with sources, or the decision is **DECLINE**
- [ ] Ricardo confirmed track, participation, and separate budget
- [ ] If any of the above is missing at **Oct 27 12:00 ET** → record DECLINE, do no C5 work
