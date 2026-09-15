# Known run-to-run failure variation (machine-local; copy to known_failure_variation.md)

Tests whose FAILURE varies from one run to the next on this installation, even with no code change. The roster comparison flags them in `needs_classification` on every comparison; a classification may cite an entry here instead of re-reading the same error texts, PROVIDED the current text matches the entry's signature. A test whose text does not match its signature is read in full and the entry is revised. An entry is evidence about the test's own behaviour, never a reason to dismiss a change: a change can still turn a varying failure into a different one.

Each entry: the test, what stays the same in every failure (the signature), what varies, the evidence that established it, and the date.

## Entries

- **`phenix/regression/tst_chat_window_mcp.py`** — signature: `AssertionError: mcp_stub not in any tool_result block` at `exercise_end_to_end_real_subprocess`, return code -6 after 3 of 3 attempts. Varies: whether a second exception (`FileNotFoundError` from `shutil.rmtree` on a `meta.json.*.tmp` file) appears during cleanup. Evidence: four server logs read 2026-09-14 (Sept 12 and Sept 13 pairs; the extra exception appeared on the baseline side in one pair and the candidate side in the other). Established 2026-09-14.
- **`phenix/regression/tst_mcp_server.py`** — signature: `subprocess.TimeoutExpired` launching `python -m phenix.mcp`, return code 1 after 3 of 3 attempts. Varies: which sub-exercise hits the time-out first (`max_concurrent_jobs` at 60 s; `cli_version` or `cli_help` at 30 s). Evidence: the same four logs, 2026-09-14. Established 2026-09-14.
