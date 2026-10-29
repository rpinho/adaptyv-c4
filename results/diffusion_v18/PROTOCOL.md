# One different interface-generation pilot

After native sequence-gradient search and a contiguous-domain trial failed,
use official RFdiffusion complex conditioning, not another seed/weight rescue.
Exactly one new69-residue backbone, native contiguous A189-383, unchanged five
hotspots205/206/208/209/210. No scaffold/motif/known binder/partial diffusion.
Seed20261022,50 resolved reverse steps, final_step1, partial_TNone. CPU actual
sampler/checkpoint/features/masks/output numbering pass. CPU proof records the
full immutable config and measured checkpoint SHA. The model came from the
upstream official HTTPS URL and content length483619179. Its URL directory is
not an authoritative MD5 checksum; the first build's mistaken comparison was
removed and this measured SHA is pinned in the deployed image proof. It is not
an independently published checksum. Missingpip and Hydra None-replay issues
were fixed before allocation; original failed CPU logs are retained.

Official image digest fa8bbf2fd51291c861369adbf246c4909c5757763d1478cafacadcc76c4584c1,
source86507b6538f51fce57b5a72477165f03999ed7ae, old PyTorch1.12.1/CUDA11.6.
A10080GB only. Actual neural GPU forward compatibility is not proven by CPU
checks. The worker refuses CPU fallback and uses the pinned sample_init,
sample_step and writepdb interfaces. JSON config/mappings/step records replace
legacy TRB; PDB B factors are zeros and are not model confidence. Native target
A189-383, generated binderB1-69. Do not apply legacy output renumbering.

Single-use container,30min worker,35min provider,40min monitor; actual zero
runners/inputs/backlog before CPU audit. Stop on malformed peptide/backbone,
non-native/deformed domain, less than50percent of requested hotspot CAs within8A
of binder CAs, restored complete A/B backbone overlap or observed-NAG overlap.
Passing this is staging only, never binding/sequence qualification. No automatic
replacement backbone on failure. If it passes, produce exactly two ProteinMPNN
samples atT.1/seed20261022, binderC only on fully restored fixed A/B receptor,
originalv48_020weights, omitCX. CPU MPNN setup/hash proof prepared separately.
Both sequences then require their own frozen-manifest CPU image preflight and
up to8 coordinate-free AF2 jobs with native639 monomer+CPU dimer restoration.
No model initial binder coordinates, admission bypass, pH or submission.

A100 GPU2.4984/hr + requested4CPU/64GiB upper.700128 =3.198528/hr.
Generation45min planning reserve2.398896; possible validation45min3.486996;
combined5.885892 plus setup/storage, not an invoice ceiling under API failure.
Prior completed Modal9.9065633; Runpod tracked estimate35.6652; reported external
consultations1.8786524; billing/setup/storage remain unreconciled. Prices checked
Oct22 at https://modal.com/pricing . Existing continued spending authorization.
Official source: https://github.com/RosettaCommons/RFdiffusion and
https://github.com/dauparas/ProteinMPNN . Astra reviews live in ../astra_review/.
