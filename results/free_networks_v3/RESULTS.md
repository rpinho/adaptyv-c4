# Isolated-binder pH-module scout

72 models: 24 variants × 3 side-chain preparation seeds; all completed. All sequence identities/checksums verified; maximum backbone displacement 0.0000 Å. The native-code crash in the first attempt is preserved separately and fixed by using three-letter residue names with MutateResidue.

**No variant produced two selected free-state His pKas >=7.** Independently repacking the isolated HH binder gives H4=6.33/H46=6.12 for all three seeds. The old bound-conformation extraction gave H46=7.25; its favorable value did not survive free-side-chain reorganization. These are fixed-backbone replicas, not an equilibrium ensemble or biological probabilities.

R49E raises free H46 to6.40 reproducibly while remaining below the same-seed HH energy proxy; V8D raises H4 to6.42. L50E reaches H46=6.74 in two replicas but6.20 in the third. None is a strong two-His module. These are mechanism diagnostics, not improved binding candidates.

**The source docking pose fails the full-receptor audit.** This scout only tests reusable isolated-binder sequence/conformation hypotheses. No receptor affinity, pH linkage, mouse cross-reactivity or experimental stability was established. Same-seed Rosetta ddE<=3 REU is an engineering filter; it is not folding free energy. Do not submit these variants on the basis of this run.

Do not expand this nearest-neighbor carboxylate search: its proposed free-pKa benefit is too small on this scaffold. Advance fresh, physically accessible full-context binders, then install geometrically explicit networks.

GPU spend from this CPU scout: $0. Separate active GPU spend and shutdown state: ../gpu_campaign/supervisor.json.
