# smoke corpus — MANIFEST

Chosen by SET COVER over a 47-feature coverage vector, not by hand.
Every channel, every screen rule, every refusal rule, every source mode,
agent/GUI/error shape, identified/abstained, and the specific whole-log
cases behind past defects. A Tier-3 test asserts the cover is still
complete; changing this set should be reviewable, so say why.

| log | md5 | why it is here |
|---|---|---|
| `agent/autobuild_4_bromodomain.log` | `0b59e25b` | 999.90 sentinel cycle with Built=0, a control skip, and the rebuild_in_place decision |
| `agent/autosol_2_p9.log` | `d20786a2` | largest agent-shape log; multiple completion events; embedded refine table |
| `agent/find_reference_14.log` | `a30dafb9` | the user complaint the consuming project exists to answer: 18 exclusions in 2 reason groups, 3 of them anonymous |
| `agent/molprobity_4_esterase.log` | `216db14a` | validation-style output; measurement labels with no stage table |
| `agent/predict_and_build_2_7mjs.log` | `f9c0295a` | wizard log: the cryo-EM branch phase; says phenix.refine 38 times and never names itself |
| `agent/refine_5_beta_blip.log` | `db7d3f6b` | the refine stage table: 28 stages, the `end` summary-row trap, 3 stages worsening r_free |
| `agent/xtriage_1_p9.log` | `40ceb4ea` | the agent FINAL QUALITY METRICS block; xtriage measurement labels |
| `agent_gui/autobuild_4_bromodomain.log` | `5ed1962d` | GUI-shape twin of `agent/autobuild_4_bromodomain.log` — the paired-shape invariant needs both members |
| `agent_gui/autosol_2_p9.log` | `0ea579cf` | GUI-shape twin of `agent/autosol_2_p9.log` — the paired-shape invariant needs both members |
| `agent_gui/find_reference_14.log` | `a30dafb9` | GUI-shape twin of `agent/find_reference_14.log` — the paired-shape invariant needs both members |
| `agent_gui/molprobity_4_esterase.log` | `f4ef2b3b` | GUI-shape twin of `agent/molprobity_4_esterase.log` — the paired-shape invariant needs both members |
| `agent_gui/predict_and_build_2_7mjs.log` | `d654fc40` | GUI-shape twin of `agent/predict_and_build_2_7mjs.log` — the paired-shape invariant needs both members |
| `agent_gui/refine_5_beta_blip.log` | `7fe427d7` | GUI-shape twin of `agent/refine_5_beta_blip.log` — the paired-shape invariant needs both members |
| `agent_gui/xtriage_1_p9.log` | `8548a257` | GUI-shape twin of `agent/xtriage_1_p9.log` — the paired-shape invariant needs both members |
| `work/err/xtriage_1_err.log` | `8137da71` | a 5-line failed run: no completion record, no error keyword, nothing to find |
| `work/ok/auto_sharpen_398.log` | `de243d68` | small long-tail program: headerless, abstains on identification |
| `work/ok/fem_407.log` | `4639ed17` | 100 BARE CR characters -- splitlines() gives 539 lines where grep -n gives 439 (D13) |
| `work/ok/map_to_model_412.log` | `22b42066` | long-tail cryo-EM tool; numeric-row screen rule at volume |
| `work/ok/sec17-sad__ai_agent_2.log` | `4219dc43` | the agent's own log: quotes children's wrapper blocks mid-file, so source=unknown (D41) |
