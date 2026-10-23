# Adapter correction before running the candidate batch

The first control attempt stopped at our negative-control assertion. Our adapter
expected a null chain type for an unassigned query. ANARCII2.0.8 actually returns
the explicit sentinel `F`, an error `Score less than cut off.`, and no numbering.
That behavior is present in the publisher's inference/model_runner.py. Both
positive controls had the expected H/K assignment; the native negative was
unassigned, not a false positive. Raw outputs, failing assertion, complete old
adapter and source/model receipts are retained in attempt_01_null_sentinel/.

The adapter now requires the publisher's actual unassigned representation for
that negative control. It treats only error-free H/K/L outputs as assigned
variable domains and requires their numbered residues to match the source.
No score cutoff, model weight, model setting or scientific pass condition changed.
The candidate batch had not run when the representation error was found.
