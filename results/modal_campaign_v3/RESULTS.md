# Native-target control: completed, no accepted designs

Three independent seeds produced 0 accepted designs from 6 trajectories. Five stopped at screen; one passed screen (pLDDT0.63) and stopped at refine (pLDDT0.55). All autotuned fields were empty. No ProteinMPNN redesign or pKa-scoring candidate was produced.

The last worker finished 2026-10-21T22:19:39.214725+00:00. Each batch independently verified zero GPU runners before its CPU follow-up. A fresh Modal stats read at2026-10-21 23:58 UTC also reported zero runners, inputs and backlog. Nothing was modeling at that checkpoint.

Estimated Modal runtime upper total **$2.8760**; prior newly supervised Runpod runtime **$1.7952**. Setup, startup/teardown, storage and older campaigns are separate; these are runtime estimates, not reconciled bills.

| Batch | Trajectories | Accepted | Stop stages | Runtime upper estimate |
|---|---:|---:|---|---:|
| 1 | 2 | 0 | screen, screen | $0.8703 |
| 2 | 2 | 0 | screen, screen | $0.9058 |
| 3 | 2 | 0 | refine, screen | $1.1000 |

The control does not support a claim of improved scaffold quality. Six attempts are too few to establish that the epitope is impossible, but they do not justify repeatedly renting the same unchanged configuration. Saved trajectory structures make diagnosis possible without another GPU. Full native target, opposite-subunit restoration, clash checks and final quality gates remain required.

The next bounded hypothesis tests upstream's large-complex initialization on the same full native monomer, patch, lengths and seeds. It injects the target coordinates at the start of the structure module; the old run used zero initialization. This is a computational setup hypothesis, not a demonstrated improvement. See ../modal_campaign_v4/DECISION.md.
