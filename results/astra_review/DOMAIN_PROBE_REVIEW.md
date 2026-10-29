# Native domain probe review — 2026-10-22

**Verdict:** the reviewed implementation supports one bounded native-domain generation diagnostic. No additional design-launch blocker found beyond completing the declared current-image CPU checks, including the explicit prepared-hotspot check described below. This review does not approve or launch the later protocol2 refold job, and no GPU was allocated by the reviewer.

## Verified

- DOMAIN differs from the repaired native v9 configuration only in run identity and target name/path. Forced targeting remains false, max trajectories remains1, and the existing loss settings and acceptance filters remain unchanged. The first-length monitor requires69 for this matched seed.
- I compared the actual domain PDB with the v9 local-frame receptor. All five hotspot Cα coordinates are identical; the patch center remains15.0003Å from binder origin. Thus the new file has not reintroduced the old uncentered-target initialization mistake.
- The domain contains native contiguous author residues189–383. The existing preparer checks native identities against the hash-locked11QC reference and preserves coordinates. The monitor's new CPU branch checks its195-residue CIF roundtrip and author numbering; output records distinguish this reduced context from the639-residue native monomer.
- The Modal design branch uses the monitored worker, finite runtime and single-trajectory profile. The local controller checks all three deployed functions before a fresh allocation, uses the existing release verification, disables ordinary pH follow-up, and only extracts a sequence after a completed job with verified GPU release. There is no automatic modeling successor.
- Selection is declared before the result: first upstream-ranked design, otherwise the sole saved trajectory snapshot. The pinned upstream ranked table and CIF conventions support this path. Missing/ambiguous outputs stop extraction. The extractor verifies native domain identity/author mapping and produces a sequence-only manifest, with source and sequence hashes. It does not pass source coordinates to the later refold.
- A source-pose restoration failure is recorded separately from hard-sequence diagnostic eligibility. A restoration success is explicitly staging-only and cannot qualify a binder. Reduced-domain or surrogate scores are not promoted to native binding evidence.

Independent checks:12 local-frame/profile tests and9 apical-restoration tests passed. I did not run neural inference. Current-image preparation and actual effective sequence/seed/frame evidence must come from root's CPU build; a local profile comparison alone does not establish those runtime facts.

## Small gaps sent to root

1. The initially reviewed domain preflight did not explicitly assert that the actual prepared target contains exactly the requested five HOTSPOT flags at author205,206,208,209,210. Its runtime geometry monitor counts those Cα positions, but that is later than the promised CPU feature check and does not itself verify loss flags. Add this assertion and an input/native-identity hash record before allocation. Comparing a file only to itself after parsing is not an independent native-identity check.
2. `prepare_domain_refold.py` initially builds residue dictionaries by assignment. Reject conflicting residue names, insertion-code collisions and multiple model identities rather than letting later atom rows silently overwrite earlier residue identity. Require the complete intended binder author numbering/length, not merely any65–75 dictionary entries. Actual pinned exports are expected to be unambiguous, so these are cheap extraction safeguards, not evidence the current source is wrong.
3. Protocol2 is a preparation artifact only. The existing fixed v15 worker validates exactly two other sequences and cannot run this new one-sequence manifest. A later worker/image needs its own frozen sequence identity, four-job cap and CPU preflight. Do not silently loosen the v15 whitelist or claim that follow-up is already executable.

## Scientific scope

This changes target representation, including removed surfaces, target length and resulting inference context. It is not a clean proof of one geometric mechanism even with all scalar settings held fixed. Native domain chemistry does not establish that the isolated domain has autonomous stability. A successful local pose must survive restoration and fresh full-native-target inference; a restored dimer remains distinct from dimer-aware prediction. Preserve the same disclosed monomer approximation and unchanged physical/admission requirements.

The predeclared rejected snapshot may be diagnosed because v15 showed that low-confidence search states can encode plausible hard-sequence folds. That does not relax generator acceptance or establish binding. One trial and one selected diagnostic sequence are sufficient for this next decision; no automatic seed/scaffold/weight rescue is implied.

## Reviewed file identities

- `scripts/modal_search_policy.py`: `706f8816793bc871ef005d4a9054d6fa2c960de654996e948304966a8848633a`
- `scripts/monitored_design.py`: `4af727b7c874aaa262924ed8ca03c876e97b10435d9057af6bf1f7209195dbf9`
- `scripts/modal_campaign.py`: `31ceee1bfc3351c643c049e0e54c36cbd4e30e75de18f888a45a9ed1ccade916`
- `scripts/run_domain_probe.py`: `e36f9fa82eee2cf97a9185b9dac450541595cdd4df88913c9070deb149f637c1`
- `scripts/prepare_domain_refold.py`: `e69e995a6c89a5e10c034c21b80b7d595b2184e93cac93b884daf6e9854a30e3`

Later edits require refreshed CPU evidence. Reported prior spending remains Modal approximately$9.5173, Runpod$35.6652 including history, consultation$1.1923; other costs unreconciled. No new GPU spending from this review.
