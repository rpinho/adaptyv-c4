# Adaptive review policy, before inspecting higher-prior outputs

The first400-design geometry shortlist has20 native-context staging passes,
but only4 have a His ring within4.5A of resolved receptor atoms and only3 have
a His nitrogen within4.5A of a native Lys/Arg group. Fifteen share exactly
the same primary geometry tuple. Geometry alone leaves many ties and does
not address the requested pH mechanism. Preserve that first ranking unchanged.

After the higher-prior arms finish and release, use the canonical geometry
ranker to take up to20 rows per arm into an80-row analysis pool, with its
existing family cap3 within each arm. This prevents file-name ties from
silently excluding an entire chemistry arm. The final hypothesis ranking
applies a global family cap3 again.
Audit that pool using the unchanged native-context gate and staging helper.
Attach the descriptive His contacts in annotate_shortlist.py. No new pH,
affinity or novelty model, energetic score, protonation calculation or assay
claim is introduced.

Use an explicit second mode in the SAME canonical ranker for the final
20-row **hypothesis review**: interpretable native contexts with His-N/Lys-Arg
proximity first, then any His-ring/receptor proximity, then other hypotheses.
Within each tier, use the measured hotspot/clash/bond geometry ordering and
the declared diversity cap. Require interpretable context for contact-based
ranking, retaining every unavailable case in provenance. Do not impose a
minimum global His count or fixed affinity-ratio requirement.

These tiers are experimental selection hypotheses. Contacts do not establish
protonation, affinity changes, their sign, or pH selectivity. Acidic-neighbor
contacts remain visible for review; no additive energy is inferred. Preserve
the geometry-only pool and both ranker provenance files so this change is
explicit. Report performance only as measured composition/geometry/context,
not as improved binding or a competition pass rate.

The unchanged canonical submission emitter may produce a format-valid review
CSV from the ordered annotations. Organizer novelty/type checks, Twin-Strep
interpretation, replacement-vs-accumulation and final submission remain with
Ricardo. This policy was adopted after the initial pair, before examining
the higher-prior sequence/metric outcomes. The two extra GPU arms already
have explicit scoped approval; no further paid successor is authorized.
