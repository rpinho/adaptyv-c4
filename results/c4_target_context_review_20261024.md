# C4 target-context review — 2026-10-24

**Verdict: context-caveated.** Boltz chain B is a concatenation of seven discontinuous TfR1 segments from an older His707 neighbourhood crop. It is not the stated contiguous apical domain 189–383, nor a shorter terminal crop of that domain. The preview’s artifact/provenance validation remains useful, but does not establish native-target readiness.

**Official scope.** Read the locally cached [official C4 page](https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/tfr1), without web access: `scratchpad/pb/tfr1.html` (SHA256 `e76de06ad64c06aba52d40878e97e7728fe8b6341232cdef913490738dbc80e6`). It specifies human P02786-1 residues **121–760, 640 residues per subunit, homodimer**, and binding at pH 7.4 with release at pH 6.0. Its displayed human sequence exactly matches that canonical interval. The apical epitope is recommended, not mandatory; **189–383 is the local apical modelling definition**, not a requirement stated numerically in this capture. Local `data/11qc.cif` contains canonical **122–760**, 639 residues in each A/B subunit: residue 121 of the assay construct is absent. `TARGET_PREP.md` calling these exactly the same construct overlooks that one-residue difference.

**Exact mapping.** All boundaries below are inclusive. Canonical numbering comes from local UniProt `data/P02786.txt` (760 residues), independently checked against every modelled receptor residue in 11QC. In both 11QC receptor chains, **author/canonical = label_seq_id + 121**.

| Boltz target chain B positions | Canonical P02786 / 11QC author positions | 11QC label positions | Count |
|---|---|---|---:|
| 1–9 | 260–268 | 139–147 | 9 |
| 10–25 | 300–315 | 179–194 | 16 |
| 26–34 | 322–330 | 201–209 | 9 |
| 35–54 | 532–551 | 411–430 | 20 |
| 55–69 | 616–630 | 495–509 | 15 |
| 70–81 | 633–644 | 512–523 | 12 |
| 82–154 | 678–750 | 557–629 | 73 |
| Local intended apical domain, positions 1–195 | 189–383 | 68–262 | 195 |

The 154-residue literal in `modal_boltz2_acid_parallel.py:14–18` exactly equals chain A of `targets/E1_his707.pdb` in file order. All five `modal_boltz2*.py` runners contain this same target; `modal_boltz2_complex.py:17` explicitly attributes it to that PDB, despite calling it apical. `scripts/trim_target.py` deliberately preserves discontinuous spatial-neighbourhood fragments with their native numbering. Boltz instead receives them as **one uninterrupted sequence** (`modal_boltz2_acid_parallel.py:55–60`), introducing six nonnative junctions and omitting the source PDB’s 55 partner-chain residues. Only **34/195** apical residues are represented; **120/154** Boltz residues lie outside that domain. This is not a simple 41-residue truncation or a global numbering offset. Crop provenance is clear; no inspected note justifies flattening it into an apical-domain surrogate. Intent beyond that reuse is unproven.

**Other pipeline contexts.** The acid-site configs use 11QC **label**, not canonical, numbers. The existing sequence exports independently confirm these intervals for all 200 designs in each arm:

| Stage/input | Actual receptor context |
|---|---|
| `configs/rfd4/tfr1_apical_binder.yaml` | A labels 68–262 → canonical 189–383, contiguous 195-residue domain |
| `tfr1_c4_acid_site1.yaml` | A555–620 → canonical A676–741 (66); B555–575 → partner B676–696 (21) |
| `tfr1_c4_acid_site2.yaml` | A530–600 → canonical A651–721 (71); same 21-residue partner segment |
| `modal_esmfold_c4_acid.py` | Binder monomer only; no receptor input |
| Acid Boltz runners / `scripts/score_acid_off.py` | Binder A plus the concatenated 154-residue receptor B above |

Thus Boltz is a fresh prediction in a different receptor context from acid-site generation. ESMFold PASS does not validate that transition. `scripts/audit_rfd4_context.py` requires an exact 195-residue target, so its apical-context results from other pools cannot qualify these acid-site/Boltz inputs. `results/apical_context/README.md` and `scripts/restore_apical_context.py` also explicitly distinguish reduced-domain staging from native-complex qualification.

**Preview versus existing v9.** Read-only checks found:

- All **343/343** canonical rebuild PDBs have exactly the configured target at B1–154. All 40 existing local PDB copies for the old 20 selected IDs share that target too. Both CSVs therefore share the same surrogate-target limitation.
- Both CSVs contain 20 unique IDs; **6 retained, 14 added, 14 removed**, with identical retained CSV rows. Rebuild manifests still match their recorded output hashes and report 122 exact plus 221 three-decimal confidence records. These checks validate data identity/completeness, not target suitability.
- Comparisons are meaningful only within that surrogate context. Confidence and contact/burial scores cannot establish binding, acid release, native epitope accessibility, or full homodimer/glycan compatibility. Artificial junctions and omitted context can alter the predicted fold and interface. No native-context pass or failure was measured here.
- The rebuild also refreshed 20 stale structures, and the preserved historical score pool does not exactly reproduce the old CSV selection. Consequently the change to 20 T1 entries versus 18 T1/2 T2 recalculated for old IDs is not a controlled demonstration of biological improvement. The stale handoff’s 13/7 split is not the current comparison.

**Use this label:** “20-design preview; artifact/provenance checks passed; target-context caveated.” Do not call it validated-ready against the challenge construct. Neither CSV is rescued by the shared mismatch. Native-target equivalence and relevant full-context evidence remain unestablished; this review does not repair, rerank, or promote candidates.

Reviewed supporting notes: `TARGET_PREP.md`, relevant `METHOD.md` domain/context sections, `HANDOFF_FROM_CLAUDE.md`, `HANDOFF_KIMI.md`, and the Oct 24 pipeline/rebuild reports. No protein sequences are reproduced here. Only this report was written; no submissions/scripts changed, jobs started, network requests made, or people contacted.

Verified CSV SHA256:

- Existing `results/submission_v9.csv`: `31d98a4133c558e1a3a9f582a607383bca8d7aad1b8bd820e684cb8fa231eadd`
- Preview `results/c4_rebuild_20261024/submission_v9_preview.csv`: `00ec3d599b50224f6d88942e599150d66301406e6a8a72a225be209a5a180cb6`
