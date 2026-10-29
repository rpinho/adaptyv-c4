# FINDING (root cause): the C4 designs were generated against the helical domain, and every evaluation since 2026-10-24 has scored them against the apical domain

Written 2026-10-26 from the real UTC clock. CPU only, $0. Read-only on all inputs.
No ranking, no submission, no organiser contact. No frozen file touched.

**This supersedes the framing of the two companion findings in this cycle and inverts the defect
narrative in `HANDOFF_TO_CLAUDE_20261026.md` §4.** It is the root cause of the anomalies both
companions measured.

## The claim

The 342 `tfr1_c4_acid_site1_*` / `site2_*` designs were generated to bind human TfR1 **author
residues ~651–741 — the helical domain — as a two-subunit groove around His707**. Every evaluation
artifact from v3 onward (the L6 recompute, the clash screen, the iPTM values, the v4 ranking, the
mouse pose audit) was computed against **author 189–383, the apical domain**, which does not
overlap the design target at a single residue.

## How the numbering was established, decisively

`configs/rfd4/tfr1_c4_acid_site1.yaml` and `tfr1_c4_acid_site2.yaml` (committed `0ff37cd`,
2026-10-23 22:15, *"Add C4-corrected RFD4 hotspots targeting TfR1 acidic residues"*) specify their
targets numerically, so the convention had to be settled before they could be read. The `C_HOT`
selectors name **specific side-chain atoms**, which makes the residue identities a decisive test:
`[CD,OE1,OE2]` are glutamate atoms, `[CG,OD1,OD2]` are aspartate atoms.

| selector | required | 11QC **author** numbering | 11QC **label** numbering (author = label+121) |
|---|---|---|---|
| site1 A573 `[CD,OE1,OE2]` | GLU | TYR — no | author 694 = **GLU** |
| site1 A593 `[CD,OE1,OE2]` | GLU | VAL — no | author 714 = **GLU** |
| site1 A607 `[CD,OE1,OE2]` | GLU | LEU — no | author 728 = **GLU** |
| site1 A561 `[CD,OE1,OE2]` | GLU | THR — no | author 682 = **GLU** |
| site2 A538 `[CG,OD1,OD2]` | ASP | ASP — match | author 659 = **ASP** |
| site2 A561 `[CD,OE1,OE2]` | GLU | THR — no | author 682 = **GLU** |
| site2 A573 `[CD,OE1,OE2]` | GLU | TYR — no | author 694 = **GLU** |
| | | **1/7** (one coincidence) | **7/7** |

The configs use **label** numbering. Translating the contigs:

| config | contig | target in author numbering |
|---|---|---|
| site1 | `60-100, /0, A555-620, /0, B555-575` | binder 60–100 aa; **chain A 676–741**; **chain B 676–696** |
| site2 | `60-100, /0, A530-600, /0, B555-575` | binder 60–100 aa; **chain A 651–721**; **chain B 676–696** |

The apical domain is author **189–383**. **Overlap: none.** Author 651–741 lies in TfR1's helical
(dimerisation) domain, and both configs include the **partner subunit**, so the intended epitope is
a composite two-chain groove — something a single-chain apical construct cannot represent at all.

Independent corroboration from the other generator arm: `configs/c4_E1_human_v3_groove.json`,
campaign `c4_tfr1_groove_human`, target `E1_his707_tight.pdb`, `"chains": "A,B"`, hotspots
**A685, A701, A702, A707, A709**, `forced_targeting: true`, `min_hotspot_contact_final: 0.5`.
Same region, both chains, and hotspot contact **enforced**. The three histidines in author 640–760
are 684, 699 and **707** — the campaign is named for His707.

## The "154-residue stitched surrogate" was the real target

`HANDOFF_TO_CLAUDE_20261026.md` §4 defect 2 condemns the Boltz chain B used by v2 as *"a
concatenation of seven discontinuous TfR1 segments from an older His707 crop, creating six
nonnative junctions and omitting 55 partner-chain residues. Only 34 of 195 apical residues are
represented; 120 of 154 modelled residues lie outside the intended domain."*

Measured:

| file | chain A | chain B |
|---|---|---|
| `targets/E1_his707.pdb` | **154 residues, 260–750** | 55 residues, 306–700 |
| `results/c4_rebuild_20261024/remote_pdbs/*_complex.pdb` | binder, 1–86 | **154 residues** |

The 154-residue construct **is** `E1_his707.pdb` chain A — the design target crop. The "55 omitted
partner-chain residues" are exactly that file's chain B, i.e. the partner subunit the design
campaign deliberately included. The phrase *"120 of 154 modelled residues lie outside the intended
domain"* measured the crop against an **assumption** that the apical domain was intended. It was
not.

So defect #2 identified something real and misdiagnosed it. **The real defect is that a
discontinuous crop was fed to Boltz-2 as one continuous chain**, creating six nonnative junctions —
a genuine modelling error. The correct remedy was to model the real target properly, as a gapped
crop or the full dimer. The remedy applied was to switch the construct to the apical domain, which
moved every downstream number onto a domain the designs were never designed against.

Consequently defects #3 and #4, and the whole v2 → v3 → v4 progression, were corrections in the
wrong direction:

| artifact | claimed status in the handoffs | status given the design target |
|---|---|---|
| v2 clash screen on the 154-residue construct (20/20 clean) | "wrong construct" | **on the design target** |
| v3/v4 clash screen on apical 195 (14/20 clean) | "the correction" | on a non-target domain |
| `linkage_apical.json`, *"L6 recomputed on the CORRECT 195-residue apical complexes"* | correct | on a non-target domain |
| v4, *"correct on all four axes it claims"* | the ranking to use | all four axes on a non-target domain |
| mouse pose audit, 20 rows | per-design mouse evidence | contacts measured on a non-target domain |

The binder chain is identical between the two constructs (the handoff states this and the sequence
hashes agree), so this is purely a question of which receptor the binders were scored against.

## This explains both companion findings

**It explains the native-context correlation.** The companion
`FINDING_NATIVE_CONTEXT.md` measured that contact frequency across the pool correlates with how
much of a position's area vanishes in the native homodimer (Spearman ρ = +0.679, n = 195). The
mechanism offered there — that truncation exposed hydrophobic intra-subunit faces and the design
method was drawn to them — is **not** the right explanation. The correct explanation is simpler:
Boltz-2 was asked to co-fold a binder designed for the helical-domain groove against an apical
construct it was never designed to touch. With no real epitope available, the binder settles onto
whatever surface that truncated construct exposes, which is predominantly its artificial
packing faces. The geometric measurement stands; its causal story is replaced.

**It explains the histidine composition.** Both RFD4 configs set `C_SEQ_PROB_HIS: 0.15` on the
binder and point every hotspot at a **carboxylate** (Glu/Asp side-chain oxygens). That is a
deliberate His–carboxylate pairing design. It is precisely the mechanism `fa98723` identified as
wrong-direction: *"Prior scoring rewarded His-carboxylate contacts that would stabilize low-pH
binding (wrong direction)."* So the **generation** was committed to acid-ON before any scoring
decision. `fa98723` corrected the scoring; it could not retroactively change what had been
generated. The companion `FINDING_L6_SIGN.md` records that v4's gate then re-selected the acid-ON
half — the inversion is therefore present in both generation and gating, from different causes.

## Correction to my own companion findings

Both companions are amended rather than rewritten. Specifically:

- `FINDING_NATIVE_CONTEXT.md`: the ρ = +0.679 measurement stands; the hydrophobic-attraction
  mechanism is **withdrawn** and replaced by the wrong-domain explanation above. Astra's review
  independently flagged that mechanism as "plausible, not demonstrated" and named imposed hotspots
  as a candidate confound; the hotspots were indeed imposed (`forced_targeting: true`), on a
  different region than the one the contacts were measured against.
- `FINDING_MOUSE_STRATEGY.md`: its footprint-conservation numbers describe contacts on a non-target
  domain and should not be read as epitope conservation for this campaign at all. The apical
  conservation measurements themselves (146/195 identical; 56.2% among natively exposed positions)
  remain valid as statements about the apical domain. **Mouse conservation of the real epitope,
  author 651–741 across both subunits, has not been measured by anyone.**
- `FINDING_L6_SIGN.md`: the gate inversion stands on the scorer definition plus its three firing
  sign controls, and does not depend on the histidine correlation. Per Astra, that correlation is
  **not** independent corroboration — L6 is computed from the same structures and ionisable groups —
  and the claim "His is the only residue that titrates between 7.4 and 6.0" is **wrong**: buried
  glutamate pKa values span roughly 5.2–9.4 experimentally, and PROPKA assigns environment-specific
  pKa values, so Asp/Glu can titrate in that window. That sentence is retracted; the conclusion it
  was offered to support does not rest on it.

## What is NOT claimed

- **I have not re-evaluated any design against its real target.** No binder has been scored against
  author 651–741 with both subunits present. That is the missing experiment, and it is the one thing
  that would say whether any of these 342 molecules is worth anything.
- **This is not a claim that the designs are bad.** It is a claim that the evidence base describes
  the wrong receptor. The designs may be good, bad, or indifferent against their actual target; no
  measurement in this repo addresses it.
- **The apical domain is not disqualified as a C4 epitope.** The challenge states it is recommended
  but not required. A future campaign could target it deliberately — but it would need designs
  generated against it, which the E2a/E2b arm (`targets/E2a_apical_202_213.pdb`,
  `E2b_apical_340_376.pdb`) was set up for and which is not what these 342 are.
- **I have not verified which generator produced each of the 342 designs row by row.** The
  identification rests on the exact name match between the designs and the two RFD4 config stanzas,
  the binder-length agreement (configs 60–100 aa; observed 65–100 aa), and the His707 target files.
  A per-design provenance join by sequence hash to generator output has not been done and should be.

## Immediate consequences

1. **No C4 file on disk is submittable**, and the reason is now more fundamental than either
   companion finding stated: the evidence base scores the wrong receptor domain.
2. **Do not build any ranking on `linkage_apical.json`, the apical clash screen, or apical iPTM.**
   All three are measurements against a non-target domain.
3. **Selection-rule v3 requires amendment on two independent counts** before any future ranking:
   its L6 direction gate, and its eligibility clause *"exact human P02786 189-383 target identity"*,
   which pins eligibility to the wrong domain.
4. The cheapest informative next step is CPU-only and does not need a new campaign: re-evaluate the
   existing 342 binder sequences against the **real** target (author 651–741, chains A and B, from
   `data/11qc.cif` or `targets/E1_his707.pdb` with the crop correctly gapped). Whether that requires
   GPU co-folding, and at what cost, is a proposal for `results/decision_requests/` — not an action
   in this cycle.

## Reproduce

```sh
cd ~/dev/adaptyv-c4
date -u '+%Y-%m-%dT%H:%M:%SZ'
cat configs/rfd4/tfr1_c4_acid_site1.yaml configs/rfd4/tfr1_c4_acid_site2.yaml
git log --format='%h %ad %s' --date=iso -- configs/rfd4/tfr1_c4_acid_site1.yaml
cat configs/c4_E1_human_v3_groove.json

# the decisive numbering test: C_HOT atoms must match residue identity
python3 - <<'PY'
import warnings; warnings.simplefilter('ignore')
from Bio.PDB.MMCIFParser import MMCIFParser
A = next(MMCIFParser(QUIET=True).get_structure('x','data/11qc.cif').get_models())['A']
res = {r.id[1]: r.get_resname() for r in A if r.id[0]==' '}
for n, exp in [(573,'GLU'),(593,'GLU'),(607,'GLU'),(561,'GLU'),(538,'ASP')]:
    print(n, 'author:', res.get(n), ' label->author', n+121, ':', res.get(n+121), ' expected', exp)
PY

# the 154-residue construct is the design target crop, not a surrogate
python3 -c "
from collections import defaultdict
for f in ['targets/E1_his707.pdb','results/c4_rebuild_20261024/remote_pdbs/tfr1_c4_acid_site1_0_complex.pdb']:
    ch=defaultdict(list)
    for l in open(f):
        if l.startswith('ATOM') and l[12:16].strip()=='CA': ch[l[21]].append(int(l[22:26]))
    print(f); [print('  chain',c,'n=',len(v),'range',min(v),'-',max(v)) for c,v in sorted(ch.items())]
"
```
