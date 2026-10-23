# Finite local run and unattended reporting

Ricardo requested overnight operation while the existing 400-design CPU run was
active. No GPU authorization or generator-lane ownership is inferred from this.

The existing CPU controller runs in its own process session without a controlling
terminal. Its process IDs and the detached observer's IDs are in PERSISTENCE.json.
The observer was verified reparented to launchd (PPID 1), in its own session,
without a controlling terminal. Closing the Remote Control connection alone does
not terminate these processes. caffeinate prevents idle sleep while the existing
controller lives, with a maximum 12-hour assertion; it does not change saved
system power settings.

watch_reporting.py observes this one existing workflow for at most 12 hours.
On completion it runs finish_reporting.py, which checks input order, sequence
hashes and lengths, JSON/CSV agreement, actual stored identity arithmetic,
recomputed local structure/level rules, all five sequence-search raw hashes,
and protected-file fingerprints. It then writes EVIDENCE.csv, VALIDATION.json
and the measured section of RESULTS.md. REPORTING_STATUS.json records progress
or the actual error. No batch is repeated and no new analysis, generator,
cloud worker, GPU, submission or account action is launched by the observer.

This is a finite queued computation, not an indefinitely improving model or a
claim that the conversational agent stays active after a final response. It
cannot survive a Mac reboot, forced process termination, exhausted disk or an
unhandled scientific/tool failure. Failed work is recorded, not called complete.
An explicit persistent Goal is a separate mechanism for agent continuation
between turns; activation was asked about without blocking this CPU batch.
Official reference: https://developers.openai.com/cookbook/examples/codex/using_goals_in_codex

Recorded tracked spend: Runpod $56.589988 + own Modal estimate $10.718947 =
~$67.31; new GPU/cloud compute $0.00. Other-lane/account billing is not measured.
