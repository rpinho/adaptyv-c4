#!/usr/bin/env python
"""How much pH selectivity can one histidine actually buy?

Contribution of a single titratable group to L = log10 Ka(pH_low) - log10 Ka(pH_high),
as a function of its free-state pKa and the depression it suffers on binding. Used to set
the per-residue design target instead of guessing it.
"""
import math
PH_LOW, PH_HIGH = 6.0, 7.4
def term(pb, pf, ph):
    return math.log10(1 + 10 ** (pb - ph)) - math.log10(1 + 10 ** (pf - ph))
def contrib(pf, d):
    pb = pf - d
    return term(pb, pf, PH_LOW) - term(pb, pf, PH_HIGH)
print("per-His contribution to L (pH %.1f vs %.1f); more negative = more pH-selective\n" % (PH_LOW, PH_HIGH))
deltas = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]
print("pKa_free " + "".join("  d=%-5.1f" % d for d in deltas))
best = (0, None, None)
for pf in [5.5, 6.0, 6.3, 6.6, 7.0, 7.4, 7.8, 8.2]:
    row = [contrib(pf, d) for d in deltas]
    for d, v in zip(deltas, row):
        if v < best[0]:
            best = (v, pf, d)
    print("%-9.1f" % pf + "".join("  %+7.2f" % v for v in row))
print("\nbest single-group contribution in this grid: %+0.2f at pKa_free %.1f, depression %.1f"
      % best)
for n in (2, 3, 4):
    print("  %d such histidines -> L = %+0.2f %s" % (n, n * best[0],
          "(clears the -2 gate)" if n * best[0] <= -2 else "(does NOT clear -2)"))
print("\nwith a modest 1.0-unit depression at the optimum: %+0.2f each -> need %d of them"
      % (contrib(6.6, 1.0), math.ceil(2 / abs(contrib(6.6, 1.0)))))
