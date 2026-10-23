# Applicable scoring mode is unestablished

The present predictor implements the general protein novelty scale. Adaptyv
explicitly uses a separate scale when its structural typer identifies an
antibody fragment. That scale depends on IMGT-numbered CDRH3 identity plus
global sequence similarity; a conserved framework alone does not establish a
low official novelty level. [Primary source](https://www.adaptyvbio.com/blog/novelty).

ANTIBODY_SCOPE_EXAMPLE.json records an actual, full-query-coverage match from
this batch to a PLAbDab heavy-sequence record. It is evidence that the gap is
relevant to these inputs, not proof that the organizer would type that query
as an antibody. A sequence-database match is not the structural typer.

Accordingly every numeric output from this run is a **provisional general-scale
calculation**. Neither that example nor a general local level 1/2 establishes
official rejection for a query whose antibody status is unknown. General local
3/4 likewise does not establish acceptance. No design is qualified, discarded,
ranked or submitted by this analysis. Antibody typing and the CDR-specific
pipeline remain unmeasured; no CDR novelty score has been computed.

The full batch and fixed general-scale thresholds are retained unchanged, so
this warning does not select a more favorable policy or alter measured scores.

Subsequent local CPU sequence-annotation evidence is recorded separately in
../antibody_scope_audit/RESULTS.md. It does not establish the organizer's type.
