# C5 stop/go checklist — TEMPLATE

Run this **before** any GPU spend and again before any CSV is called ready.
Every line must be an explicit YES/NO with evidence. Any NO in a **GO** section = stop.

## A. Authorization (all required)
- [ ] Ricardo explicitly approved C5 participation — evidence:
- [ ] Track confirmed (1 / 2 / 3) — evidence:
- [ ] Separate C5 budget approved, amount $____ — evidence:
- [ ] Account/upload decision confirmed by Ricardo — evidence:
- [ ] Disclosure terms accepted by Ricardo — evidence:
- [ ] Decision deadline (Oct 27 12:00 ET) not passed without the above → otherwise **DECLINE**

## B. Target understanding (all required)
- [ ] Official C5 page fetched read-only and **snapshotted with timestamp + sha256**
- [ ] Public selection prompt fetched, snapshotted, and summarised
- [ ] Target sequence obtained and hashed (validator records the same hash)
- [ ] Numbering convention confirmed against the provided structure (not assumed)
- [ ] Assay construct, oligomeric state, and conditions recorded verbatim
- [ ] Objective is unambiguous (binding vs inhibition vs agonism vs pH-conditional)
- [ ] Slot limit taken from the track, not assumed (20 vs 40)

## C. Method readiness (all required before scaling)
- [ ] Chosen port validated against a **known result** on the same design (numbers recorded)
- [ ] Construct/sequence identity gate is byte-exact (sha256 + `cmp`)
- [ ] Fail-closed tests green (`tests/c5/`, and the relevant design/gate tests)
- [ ] De-novo/zero-shot compliance verified, with the method of verification recorded
- [ ] Diversity rule stated numerically and applied
- [ ] Every new candidate set written as a NEW versioned file; no fallback overwritten
- [ ] Confirmed the previous version's sha256 is unchanged (nothing was overwritten in place)

## D. Pre-submission (all required; submission itself needs Ricardo's word)
- [ ] `scripts/c5_validate_submission.py` PASS on the exact final file, with `--manifest`
- [ ] Manifest reviewer sign-off recorded (human, not the agent)
- [ ] Hash of final file recorded in the handoff and the provenance record
- [ ] Rank order justified in one line per file; best file named and why
- [ ] Known limitations written plainly (domain vs native, monomer vs oligomer, model metrics ≠ biology)
- [ ] **Nothing uploaded, nobody contacted** — confirmed explicitly

## E. Kill conditions (any one = stop the branch now)
- integrity/hash/construct check fails
- GPU-hour, wall-clock, or dollar ceiling reached
- no measurable progress for ~45 min (consult Kimi/Astra once, then record the dead end)
- a required gate script is missing or its fixtures fail
- a fallback artifact would have to be overwritten to proceed
- Ricardo's mission authorization window has expired (08:00 ET Oct 25) and no new instruction exists
