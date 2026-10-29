# Archive receipt — codex-adaptyv lane stopped

Performed 2026-10-26, on Ricardo's explicit instruction ("yes, kill all three and the idle tmux
sessions"). Reversible: every stopped Codex session has a rollout on disk and can be resumed by id.

## Correction made BEFORE acting

My earlier report named PID 2184 as "tmux codex-adaptyv-c4 -> codex resume 01a15f2f". **That was
wrong and would have been destructive.** PID 2184 is the **tmux server itself** — the parent of
every session on this machine, including this Claude session (PID 4091), the VCC lane (332), the
games lane (51034), the general lane (31221), and the KB/games Codex workers (39623, 40552). Its
command line is merely the `tmux new-session` that first started the server. The
`codex-adaptyv-c4` session no longer exists and Codex session `01a15f2f` is **not running**.
PID 2184 was NOT killed.

## Stopped

| PID | codex session | tmux session | model/route | notes |
|---|---|---|---|---|
| 45619 | `01a170c9-52f7-79b1-a5cb-8cb0f59d5e8b` | `codex-adaptyv-worker` | deepseek-v4.1-flash via OpenRouter | worker A; ran `--dangerously-bypass-approvals-and-sandbox` |
| 48408 | `01a173a5-524b-7a03-b0be-a6f851bd3f2d` | `codex-adaptyv` | deepseek-v4.1-flash via OpenRouter | non-owner advisor lane (ledger W11); same bypass flags |
| 86282 | (codex `agents`) | `codex-adaptyv-browser` | deepseek-v4.1-flash via OpenRouter | agents listener |

tmux sessions killed: `codex-adaptyv`, `codex-adaptyv-worker`, `codex-adaptyv-browser`,
`kimi-consult` (held only a `sleep`).

## To restore

```sh
codex resume 01a170c9-52f7-79b1-a5cb-8cb0f59d5e8b -C ~/dev/adaptyv-c4
codex resume 01a173a5-524b-7a03-b0be-a6f851bd3f2d -C ~/dev/adaptyv-c4
# previously-ended, listed for completeness:
codex resume 01a15f2f-94e9-77c3-a384-0a24d7871728 -C ~/dev/adaptyv-c4
codex resume 01a170d2-cc3f-7aa3-9217-9831e6a39daf -C ~/dev/adaptyv-c4   # canonical orchestrator, self-archived
```

## NOT touched

- PID 2184 (tmux server) and every unrelated lane: VCC, family-learning-games, KB, catalog, general.
- Codex workers 39623 / 40552 — these carry KB_OPS and GAMES_OPS prompts, not Adaptyv.
- The three read-only guards: `com.ricardo.adaptyv-c4-cpu-review` (6h),
  `com.ricardo.adaptyv-c3-overnight-hash` (900s), `com.ricardo.adaptyv-c3-gpu-freeze` (120s).
- `com.ricardo.codex-adaptyv-tick.plist` — present on disk, deliberately NOT loaded and NOT
  reloaded: ADAPTYV_THREAD points at the archived thread and reloading would queue cycles into a
  dead lane.
- No frozen submission file was read, written or staged.
