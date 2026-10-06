# Project method template

Guide instructions: use only relevant sections, replace placeholders with
observed or explicitly proposed values, and omit irrelevant example rows.
Reuse an existing method's structure. This template is not authority and
its example capabilities are not requirements. Follow `../SETUP.md` for
recovery, verification, authorization, saving and migration.

## Project and instructions

- Project/checkout: <actual target and account/machine scope>
- Current authority: <actually loaded project instruction path>
- Canonical project method: <this method's actual path>
- Related repositories: <paths and permitted scope, if relevant>
- Shared defaults: <existing source and project overrides, if any>
- Setup record: <record containing decisions, prior bytes and verification>

## Applicable work

| Activity | Applicable now, deferred, or not applicable | Reason and requirements still missing |
| --- | --- | --- |
| <e.g. local tests> | <state> | <reason or remaining check> |

## Working locations

| Purpose | Path and host/account | Source and verification state |
| --- | --- | --- |
| Source checkout | <path> | <source/check> |
| Isolated worktrees / temporary work | <authorized or proposed path> | <source/check; cleanup responsibility> |
| Durable records / results | <path retained after cleanup> | <source/check> |
| Search roots for discovery | <named roots; never the home directory or other applications' data> | <source/check> |

## Environment and commands

| Purpose | Exact command or definition source | Host, shell, working directory and environment | Status and authorization reference |
| --- | --- | --- | --- |
| <applicable build/test/check> | <known command; do not guess> | <actual context> | <verification separately from permission> |

Record expected runtime when known. Record required executable locations,
PATH/activation and hook requirements. An alias needs its definition, not
just its name. If a value is unavailable, say what was searched and what
will establish it. Keep a recovered value while marking it unverified.

## Remote or CI work — only if relevant

For each relevant host/service: account and access method (no secrets),
installation, workspace, commands and shell, log retrieval, unattended or
user-run steps, approved resource cap, lock/owner/load rules, rebuild and
cleanup responsibility. Distinguish hardware facts from resource approval.
Name the fresh checks needed before use. Do not borrow another host's
paths or authorization. Omit this section for a local-only project.

## Coordination and action boundaries

Record relevant branch/worktree and sharing rules; who may update a working
installation; current restrictions on tests, integration and publication;
and pointers to applicable decisions/grants. Recording a command does not
authorize running it. Preserve no-push and stage-specific restrictions.

Record the destination's tag and release convention with its source, and
the accepted push shape: explicit `<commit>:refs/heads/<branch>`, tag
following disabled, no development or recovery tags. Name each
installation kept in step after a publication: canonical root, access
host (hostnames sharing one filesystem are one installation),
repositories, permitted operations, verification, reservation and load
rule, recovery boundary (expected HEAD, recorded operation state, no
intervening work lost) and stop conditions. Publication alone authorizes
no installation update. Record how suite and outside-reading waivers are
given: in the Developer's words, per batch and scope; "commit and
publish" waives neither.

## Sources and unresolved items

| Setting or rule | Value / proposed default and scope | Source | State / needed check or decision |
| --- | --- | --- | --- |
| <recovered item> | <retain value even if stale> | <path, section, date> | <recorded / verified with check / unknown / conflict / deferred / not applicable> |

Keep both sides of an unresolved conflict. For a verified fact, identify
the host/account/shell, date and evidence in the setup record. For a user
choice, identify its actual decision rather than labeling it a measured
fact. Name which operation each unresolved item affects and the next step.
