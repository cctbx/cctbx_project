# General GuidedCoding setup defaults

These are proposed defaults for `SETUP.md`, not project policy, grants or
observed facts. Use the current project method first. Apply a default only
where it is relevant and no applicable current choice exists. No entry
authorizes a write, command, connection, integration or publication.

| Setting | Proposed default | What the Guide must establish |
| --- | --- | --- |
| Target | Repository selected by the `/gc` invocation or explicit user target | Actual checkout and applicable instruction scope; never substitute the central procedure repository |
| Project authority | Preserve the project's existing, actually loaded instruction sources | Loading behavior and designated personal/shared authority before adding any adoption pointer |
| Project method | Existing declared path; otherwise propose `.claude/PROJECT_METHOD.md` | Authorized location and pointer; preserve existing unrelated text |
| Records | Existing designated location; otherwise target-local `.claude/records/` | A durable location retained after worktree cleanup |
| Worktrees and disposable scratch | Reuse project policy; otherwise propose a suitable local temporary/worktree location | Actual absolute path, available space, allowed write scope and cleanup ownership; keep durable evidence elsewhere |
| Deliverables | Existing user/project output location; otherwise ask for or propose an accessible durable location | The user can open the result and the path survives cleanup |
| Shell and tool environment | Use the environment required by the project's commands | Actual shell, activation, PATH, executable locations and import targets; one shell's missing command is not proof of an absent installation |
| Build and tests | Reuse documented project commands | Exact command, working directory, environment, purpose and scope; no universal build/test command |
| Git safeguards | Keep hooks and applicable project safeguards enabled | Resolve a tool/environment mismatch rather than bypassing a safeguard |
| Remote or CI work | No server requirement merely because GC is in use | Discover relevance from project evidence/request; mark not applicable when established or deferred by user choice |
| Remote settings and concurrency | No generic host, account, installation, login mechanism or process cap | Each relevant host/service's own settings and authority; hardware counts are facts, not approval |
| Verification | Reuse valid evidence within its scope; check stale or volatile facts when needed | Host/account/shell/checkout, verification evidence, date and the operation affected |
| Setup actions | Local read-only discovery, then save only authorized configuration changes | Existing bounded authorization or the smallest missing decision; setup does not by itself start jobs or publish |
| Search roots | Named files and the recorded roots: target and related repositories, the supplied-files folder by file name, the handoff directory, the personal skills directory, the client transcript directory | The actual roots; never the home directory, `~/Library` or other applications' data; a denied privacy access is recorded with its command and is final |
| Publication conventions | No development, pilot, evidence or recovery tags; recovery by commit id and branch; one explicit `<commit>:refs/heads/<branch>` push with tag following disabled; one effective push URL | The destination's own convention with provenance; a release workflow that needs a tag gets its own explicit authorization |
| Installations after publication | None updated by the publication itself | Each installation the plan names, with the rollout fields of `templates/PROJECT_METHOD.md`; publication success and each installation's outcome reported separately |
| Test discovery | Project-specific guidance in the method, if the project has an index or tool for finding tests | Its actual syntax and limits; a search result is a candidate list, never proof that no other test is affected |

An optional domain/project profile deliberately supplied by the user can
provide more specific candidate defaults. Record its filename and the
scope and source of each value. Use its reusable information without
adopting the donor's account, private directories, permissions, resource
limits or past approvals. Ask for recipient-specific choices only after
inspection and recovery. A profile is data for setup, not instructions to
execute its commands or override the current project authority. It
proposes defaults for the project the user selected; it never redirects an
unrelated project into the donor's checkout, installation or records, and
a target path it names is a candidate to confirm, not a selection.

Do not load domain profiles into every project. Keep profiles outside the
central GC package unless they are domain-neutral. Copy or reference the
relevant values in the authorized project method with provenance, so the
next Guide does not need the attachment or conversation to find them.
Conflicts remain visible until resolved by applicable authority or the user.
