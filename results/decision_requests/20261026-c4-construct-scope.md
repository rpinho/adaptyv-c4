# DECISION REQUEST — C4 construct scope (Ricardo only)

Raised 2026-10-26. **No GPU has been spent on this and none will be without your explicit yes.**
Nothing is blocked except which construct any future C4 work is computed on.

## The question, in one line

Do we keep evaluating designs on the **195-aa apical surrogate** (cheap, everything already on
disk), or pay for the **official 640-aa construct** (correct, expensive)?

## Why it is now cheaper than the last time I asked

The binary framing above is stale. The native-context measurement shows the surrogate's dominant
error is **intra-subunit — 13.0% of footprint area lost to the design's own subunit** — not the
partner chain (4.9%) and not glycans (1.5%). Membrane-proximal geometry is 36.62 A away and is not
a factor.

That matters because it means a **640-aa monomer** may recover most of the correction. The partner
subunit contributes ~4.9%, so if the goal is "stop scoring a surface that does not exist", a monomer
gets most of the way there at a fraction of the homodimer cost.

Three options, cheapest first:

| option | what it costs | what it fixes | what it leaves broken |
|---|---|---|---|
| **A. Stay on 195-aa surrogate** | $0 (already paid) | nothing new | 13% intra-subunit error, partner absent, 44.5% median footprint retention |
| **B. 640-aa monomer** | new GPU, ~1-2x per-design cost | intra-subunit error (13%) | partner interface (4.9%) |
| **C. Official 640-aa homodimer + glycans** | new GPU, full cost | everything measured here | nothing known |

Measured per-design GPU cost is **$1.83-2.79/design**. A 20-design ranking under B is therefore
roughly **$37-56**; under C, roughly double that. I have not run either and will not.

## What I recommend

**B**, and only after the current $0 work lands. Reason: the evidence says the error we can actually
fix cheaply is intra-subunit, and B fixes it while C's extra cost buys a 4.9% term we have not shown
to change any ranking. But B is not free and the per-design numbers above are estimates from history,
not quotes.

## What I need from you

One of: **A**, **B**, **C**, or **"none yet — keep doing $0 work"**. Reply with the letter.

If you pick B or C, also confirm a spend ceiling. Current cumulative C3 sprint spend is $8.96/$20
and C4 has spent $2 on a GPU run that produced zero usable rows, so I am treating any GPU request as
requiring an explicit number from you.

## What is true regardless of your answer

No C4 file on disk is submittable today. The pH gate was inverted (now corrected in rule v4), the
pool is aimed mostly at a different surface, only 19 of 343 designs touch the anchor residue His707,
and zero designs pass all three objectives together. Changing the construct does not fix any of
that — it only makes future measurements honest. The generation defect (`C_SEQ_PROB_HIS: 0.15`,
carboxylate-aimed hotspots) is upstream of the construct and no construct choice repairs it.

---

# AMENDMENT 2026-10-26T14:56Z — the evidence behind this request moved. The question is unchanged.

A later $0 measurement this same day
(`results/c4_real_epitope_mouse_20261026/FINDING_REAL_EPITOPE.md`) measured the **real** design
epitope — author 651-741 — for mouse conservation and native accessibility. Two numbers change how
this request should be read. **Your options are still A / B / C / "none yet", and I have not acted.**

## 1. B over C is now supported by a measured zero, not by inference

I recommended B partly because the partner subunit looked like a small term. At the real epitope it
is small where it counts and **exactly zero where it counts most**:

| term | band-wide (A^2) | at the ten named generator hotspots |
|---|---:|---|
| lost to the design's own subunit | 2644.8 (34.2%) | 50-160 A^2 at five of them |
| lost to the partner subunit | 994.7 (12.9%) | **0.0 A^2 at all ten** |
| lost to observed glycans | 101.7 (1.3%) | negligible |

B (640-aa monomer) therefore captures the **entire** hotspot-level correction. C's extra cost buys
a term measured at zero at every hotspot. If you pick a paid option, B is the one the evidence
supports.

## 2. But the ceiling on what either option can discover is now known, and it is low

The same measurement found that the pH-switch geometry this campaign was built on does not exist in
the native receptor:

- **All three histidines in the band are natively inaccessible.** H699 is completely buried
  (0.0 A^2). H684 retains 21.8 A^2. **H707 — the residue the campaign is named for — presents only
  28.5 A^2 (rel. SASA 0.127), and that is true in the band alone as well as in the full homodimer.
  No construct changes it.**
- The only four natively exposed hotspots are **all carboxylates** (D659, E694, E714, E728) — the
  acid-ON anchors `fa98723` identified as the wrong pH direction.
- **0 of the 20 direction-correct v4 seats contact H707 at all**; 19 of 343 designs do.
- The median design makes only **5** contacts on natively-exposed real-band surface.

Good news in the same measurement, for completeness: mouse conservation is **not** a blocker at the
real epitope — 83.5% band-wide, **80.0%** among natively exposed positions (apical was 56.2%), and
all ten named hotspots plus all three histidines are mouse-identical.

## 3. What that means for your choice

Paying $37-56 under B would measure *this pool* against *this epitope* more honestly. It would not
change that the epitope offers no accessible receptor histidine, that the available anchors push
the pH response the wrong way, and that no seat touches the anchor residue. **"None yet — keep
doing $0 work" is a more defensible answer than it was when I filed this request**, and I am
flagging that rather than quietly leaving a stale recommendation standing.

I am not changing my recommendation on my own authority, because the construct choice is yours. If
you still want a paid run, B. If you want my honest read now: **none yet.**

Nothing has been spent. No GPU has been started. Cumulative C4 GPU spend remains the historical $2.
