# C4 linkage audit — bound vs unbound (CPU only)

**Status: evidence only. Nothing promoted, no candidate replaced, nothing submitted.**

## Controls
```
{"null_L": 0.0, "acid_on_L": 0.5918, "acid_off_L": -0.4539, "null_ok": true, "acid_on_ok": true, "acid_off_ok": true, "all_ok": true}
```

Positive L = tighter binding at the lower pH (C4 conditional-binding direction).

## Result
| metric | value |
|---|---|
| preview rows | 20 |
| scored | 20 |
| missing complex PDB | 0 |
| L(6.0/7.4) > 0 | **11** |
| L(6.0/7.4) < 0 | **9** |
| L(5.5/7.4) > 0 | 10 |
| L(5.5/7.4) < 0 | 10 |

## Ranked by L(6.0/7.4), most acid-ON first

| design | H | L 5.5/7.4 | L 6.0/7.4 | direction |
|---|---|---|---|---|
| tfr1_c4_acid_site1_156 | 5 | 2.083 | 1.68 | ACID_ON |
| tfr1_c4_acid_site1_137 | 3 | 2.188 | 1.58 | ACID_ON |
| tfr1_c4_acid_site2_74 | 3 | 1.464 | 0.901 | ACID_ON |
| tfr1_c4_acid_site2_171 | 3 | 1.312 | 0.858 | ACID_ON |
| tfr1_c4_acid_site1_56 | 6 | 0.981 | 0.796 | ACID_ON |
| tfr1_c4_acid_site2_155 | 4 | 1.002 | 0.713 | ACID_ON |
| tfr1_c4_acid_site1_1 | 3 | 1.432 | 0.685 | ACID_ON |
| tfr1_c4_acid_site1_144 | 10 | 0.113 | 0.538 | ACID_ON |
| tfr1_c4_acid_site2_137 | 5 | 0.555 | 0.318 | ACID_ON |
| tfr1_c4_acid_site2_62 | 8 | -0.124 | 0.036 | ACID_ON |
| tfr1_c4_acid_site2_1 | 4 | -0.148 | 0.008 | ACID_ON |
| tfr1_c4_acid_site2_92 | 7 | 0.138 | -0.003 | ACID_OFF |
| tfr1_c4_acid_site1_174 | 6 | -0.174 | -0.181 | ACID_OFF |
| tfr1_c4_acid_site2_197 | 4 | -0.537 | -0.429 | ACID_OFF |
| tfr1_c4_acid_site2_52 | 7 | -1.2 | -0.603 | ACID_OFF |
| tfr1_c4_acid_site2_4 | 5 | -1.615 | -0.918 | ACID_OFF |
| tfr1_c4_acid_site1_177 | 7 | -1.693 | -1.063 | ACID_OFF |
| tfr1_c4_acid_site2_169 | 8 | -2.261 | -1.253 | ACID_OFF |
| tfr1_c4_acid_site1_143 | 8 | -2.933 | -1.78 | ACID_OFF |
| tfr1_c4_acid_site2_73 | 3 | -2.921 | -1.911 | ACID_OFF |

## Interpretation

PROPKA bound-vs-free linkage on model structures; detects sign contradictions, does not establish experimental pH-conditional binding or release. Elevated bound pKa alone is not sufficient evidence of switching.

Sign controls pass, so a negative L here means the modelled complex binds *more* tightly at the higher pH — the C3/endosomal-release direction, not C4's objective. This audit cannot prove release; it can only show whether the modelled electrostatics point the wrong way.
