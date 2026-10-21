> **Coverage diagnostic, superseded.** Two terminal-His variants contain partial linkage sums. Use the audited replacement in `../campaign/`.

# CPU switch-installation results

Campaign status: **completed**. Completed backbones: **25/25**.

- Configured cap: 3 ranked pairs per backbone from 8 shortlisted sites.
- HH pairs attempted: 75; retained by neutral-energy filters and scored: 38.
- Retained HH with a retained QQ twin: 27.
- Retained HH with L <= -2: 0; with a retained QQ twin: 0.
- HH pairs with both selected His free-state pKa >= 7: 0.
- Runpod GPU spend this session: **$0.00**. Last recorded historical pod spend: **$33.87**. Both pods left stopped. Persistent-volume storage not freshly reconciled.

References and variants use target AB and binder C. Historical ordering used A only. The whole complex is first repacked and side chains minimised with the backbone fixed; all variants then use the same local protocol as their reference. Retention requires negative reference and variant dG and ddG <= 3 REU. ddG is Rosetta REU; L is PROPKA Wyman linkage at pH 6.0 versus 7.4. Neither establishes experimental affinity or the assayed pH ratio.

## Best retained HH variants

| Backbone | Sites | L | ΔL vs relaxed AB reference | ddG (REU) | QQ retained | Selected free pKa |
|---|---|---:|---:|---:|---|---|
| l63_59b0f47ecf134650 | CD4;CT46 | -1.530 | -1.179 | +1.74 | True | C46:7.25;C4:6.33 |
| l76_06cf7c724b17b60d | CQ53;CA56 | -1.359 | -0.539 | -6.17 | False | C53:6.48;C56:5.96 |
| l76_06cf7c724b17b60d | CQ53;CS43 | -1.316 | -0.496 | -9.48 | True | C53:6.50;C43:6.43 |
| l63_59b0f47ecf134650 | CD4;CD42 | -1.252 | -0.901 | -1.93 | True | C42:6.29;C4:6.33 |
| l56_8c4eb2498f38b9a8 | CE15;CQ27 | -1.197 | -0.551 | +2.00 | True | C27:7.03;C15:5.74 |
| l70_54a5b53444de2b49 | CA12;CR26 | -1.189 | -0.643 | +2.66 | False | C12:6.52;C26:6.32 |
| l78_cc1057879d35a8ec | CQ7;CE4 | -1.189 | -0.507 | +2.64 | False | C4:6.30;C7:6.20 |
| l70_54a5b53444de2b49 | CA12;CR27 | -1.187 | -0.641 | +2.90 | False | C12:6.52;C27:6.66 |
| l98_ac2cbaba03c60f1c | CE95;CR80 | -1.128 | -0.596 | +2.23 | True | C95:6.37;C80:6.43 |
| l80_24ee63503d07a260 | CE57;CE66 | -1.006 | -0.607 | +2.83 | False | C57:6.49;C66:6.34 |

## Limits and stop point

The configured ranked-pair cap is a screen, not an exhaustive combinatorial search. The coordinate ranking is a heuristic; selected acidic residues can themselves be removed by mutation and neighbouring sites can interact. Measured mutant pKa and the full acid budget, not pair-score sums, must decide whether the mechanism improved.

Ranking used source coordinates before preparation. Rotamer changes can materially change PROPKA pKa; single-conformer sensitivity limits absolute pKa and energy claims.

No APBS estimator, orthogonal refolding, mouse validation, account creation or submission was run. This CPU campaign stops after its configured backbone/pair count. No rented pod is needed.

Files: `campaign.json` (full state), per-design manifests and PDBs, `variants.csv` (including rejected twins), `residues.csv` (both tails of every scored linkage sum).
