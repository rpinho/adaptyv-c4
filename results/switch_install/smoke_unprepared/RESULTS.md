> **INVALIDATED diagnostic run. Do not select these variants.** The binding reference used inconsistent preparation. See [INVALIDATED.md](INVALIDATED.md).

# CPU switch-installation results

Campaign status: **completed**. Completed backbones: **1/1**.

- Configured cap: 1 ranked pairs per backbone from 8 shortlisted sites.
- HH pairs attempted: 1; retained by ddG filter and scored: 1.
- Retained HH with L <= -2: 0; with a retained QQ twin: 0.
- HH pairs with both selected His free-state pKa >= 7: 0.
- Runpod GPU spend this session: **$0.00**. Last recorded historical pod spend: **$33.87**. Both pods left stopped. Persistent-volume storage not freshly reconciled.

References and variants use target AB and binder C. Historical ordering used A only. All variants use the same local fixed-backbone repack/minimise protocol as their reference. ddG is Rosetta REU; L is PROPKA Wyman linkage at pH 6.0 versus 7.4. Neither establishes experimental affinity or the assayed pH ratio.

## Best retained HH variants

| Backbone | Sites | L | ΔL vs relaxed AB reference | ddG (REU) | QQ retained | Selected free pKa |
|---|---|---:|---:|---:|---|---|
| l58_9f94b53776050aad | CH24;CY32 | -0.704 | -0.215 | -0.35 | True | C24:6.76;C32:6.30 |

## Limits and stop point

The configured ranked-pair cap is a screen, not an exhaustive combinatorial search. The coordinate ranking is a heuristic; selected acidic residues can themselves be removed by mutation and neighbouring sites can interact. Measured mutant pKa and the full acid budget, not pair-score sums, must decide whether the mechanism improved.

No APBS estimator, orthogonal refolding, mouse validation, account creation or submission was run. This CPU campaign stops after its configured backbone/pair count. No rented pod is needed.

Files: `campaign.json` (full state), per-design manifests and PDBs, `variants.csv` (including rejected twins), `residues.csv` (both tails of every scored linkage sum).
