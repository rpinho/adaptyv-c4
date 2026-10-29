# NOTE: `results/c4_ranked_v4_apical_L6_20261025/` is superseded and carries the INVERTED gate

Written 2026-10-26 from the real clock. The directory was deliberately **not renamed**: five
scripts reference that exact path and the integrity manifest records per-file hashes for it.
A rename was attempted and fully reverted. All five reviewers chose "leave in place and mark".

## What it is

The first "v4" ranking. Its rule is `iPTM>=0.60 AND apical L6>0`, which selects **tighter binding
at pH 6.0**, the opposite of challenge objective 2. It seats 20 acid-ON designs and has a field
named `n_eligible_acid_on_apical`.

| file | sha256 |
|---|---|
| `c4_acid_ranked_v4_apical_L6_20261025.csv` | `489e95fdca8e4d032acec444a42dc2387c122a168c674b0b32145010c100068b` |
| `c4_acid_ranked_v4_apical_L6_20261025.json` | `1a251c820fb1b6cad79c7ba60f6cfdde8583029b1b49b79dfd02396d2f65cdb7` |

## Do not confuse it with the other "v4"

`results/c4_ranking_v4_acid_off_20261026/v4_ranking_acid_off.json`
sha256 `75236da7cc1ae106c5649cc310a347e8564c4fbd525cfc7f4779abba6d613b3f`
is the acid-OFF gate. **It is also withdrawn as a ranking** by selection rule v5
(`results/decisions/20261026-c4-selection-rule-v5.md`, commit `dcf5f57`), which removes L6 from
eligibility and ranking altogether.

**Neither "v4" is a ranking to use.** Neither is submittable.
