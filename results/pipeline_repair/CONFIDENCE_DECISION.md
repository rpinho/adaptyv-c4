# One confidence-objective sensitivity test after V11

V11 completed50 gradient calls and failed screen(pLDDT.36/iPTM.13),0 accepted.
It approached the requested patch but produced CA overlaps and malformed binder
geometry, then ended42.044A away. Every observed prediction used375 surrogate
target changes; native restoration was never reached. This is failure, not an
intermediate candidate. Actual resource release verified07:45:29EDT.

Declare exactly one `c4_forced_confidence_v12` trajectory, changing only
`weights_plddt_loss`0.1→0.5 from V11. Preserve forced targeting, full target
backbone, seed/69-residue first binder, other losses, all acceptance gates and
native-restoration/admission rules. Hypothesis: increased confidence objective
pressure may change the poor-fold/patch-acquisition tradeoff. It is an objective
sensitivity test, not a proven geometry repair. Scalar component magnitudes do
not establish gradient dominance. High pLDDT alone cannot qualify a binder.

Astra verified the actual pinned schema: top-level `weights_plddt_loss` is the
correct key; nested `losses.plddt_loss.weight` is rejected. No allocation before
the installed-image CPU `read_settings`/`build_losses` check proves effective
weight.5, baseline.1 and unchanged other built objectives/parameters. Existing
prepared-input frame/sequence preflight and exact config-delta regression also
must pass. Current local suite58 tests passes.

Same30min remote /35min local caps, single-use max1/min0 GPU, no successor, actual
release before CPU follow-up. No early surrogate output may enter pH scoring.
This one-seed test cannot estimate hit rates or establish causal improvement.
Failure leads to CPU restoration/control work described in Astra's ALTERNATIVES,
not an automatic weight/length sweep or relaxation of qualification gates.

Existing continued modeling/spend authorization applies. H100GPU$3.9492/hr;
requested CPU/RAM-inclusive runtime upper$4.649328/hr,35min envelope$2.7121 plus
setup/storage. Prior Modal runtime upper$5.8023551, newly supervised Runpod$1.7952
additional to historical$33.87, consultations$1.1922999. Actual next-run state
belongs in results/modal_forced_confidence_v12/supervisor.json after allocation.
