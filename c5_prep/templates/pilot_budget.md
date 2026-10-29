# C5 small-pilot budget — TEMPLATE (no spend until Ricardo approves)

**Default: $0. Nothing here authorises compute.** Ricardo-only unknowns: track, participation,
separate budget, account/upload, disclosure terms. Until all are confirmed, the standing decision
at **Oct 27 12:00 ET** is DECLINE and this file stays a template.

## 1. Measured unit costs (from this repo's own runs — do not re-derive)

| lane | unit cost | source |
|---|---|---|
| RunPod A40, 1 pod | $0.59/h | pod list, 2026-10-24 |
| AF2M screen, one 111-aa complex on an **exclusive** A40 | 105 s wall → **$0.017/design** | `results/c3_af2m_port_20261025/validation_20261025.md` |
| RFD4 backbone generation, one backbone | 45–77 s → **$0.007–0.013/backbone** | round-2 bins A/B/C, 36/36 |
| ProteinMPNN, one backbone (2 seqs) | ~5 s → **~$0.001/backbone** | p4 `mpnn.manifest` |
| ESMFold refold, one monomer | ~30–60 s → **~$0.005–0.010/design** | p3 `esmout` |
| Modal L40S | $1.95/h — **~3.3× the A40**; only for work the pods cannot do | Modal billing snapshot |

Contended-GPU runs cost 2–3× the exclusive numbers (285–350 s vs 105 s observed). Budget the
exclusive figure and treat contention as a risk, not a plan.

## 2. Pilot scope (fill when the unknown is answered)
| field | value |
|---|---|
| smallest compliant experiment | |
| designs to generate | |
| designs to screen | |
| pods needed (≤5) | |
| wall-clock ceiling | |
| GPU-hour ceiling | |
| dollar ceiling | |
| hard stop (date/time ET) | |
| what result would justify scaling | |
| what result means stop | |

## 3. Arithmetic check (fill; must be internally consistent)
```
backbones      N_b  × $0.010  = $___
MPNN           N_b  × $0.001  = $___
refold         N_r  × $0.008  = $___
screen         N_s  × $0.017  = $___
pod idle/setup overhead (assume 20%) = $___
                                    TOTAL = $___
```
Ceiling must be ≥ TOTAL and must be approved **before** the first pod is started.

## 4. Stop rules (non-negotiable)
- Stop at the GPU-hour ceiling, the wall-clock ceiling, or the dollar ceiling — whichever hits first.
- A stopped pod is removed, not left idle; idle A40 = $0.59/h of pure waste.
- RunPod only; **no new Modal spend** without a separate explicit approval.
- Any failed integrity check (target hash, construct identity, validator) stops the branch.
- No measurable progress for ~45 min → stop debugging, consult Kimi/Astra, record the dead end.

## 5. Sign-off
| item | who | when (UTC) |
|---|---|---|
| participation approved | Ricardo | |
| track confirmed | Ricardo | |
| separate budget approved (amount $___) | Ricardo | |
| pilot plan reviewed | Kimi/Astra | |
