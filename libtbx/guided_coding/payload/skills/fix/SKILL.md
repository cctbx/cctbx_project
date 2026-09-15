---
name: fix
description: Start a GuidedCoding change from a bug report or problem description. The one command a developer needs to begin.
disable-model-invocation: true
argument-hint: [describe the problem in plain language]
---

# /fix - start a change

The developer has just described a problem. Do this, in order:

1. Invoke the workflow skill if it is not already loaded, and follow it for everything that follows.
2. Restate the problem back in one or two plain sentences and confirm you have it right.
3. Classify it under the consequential rule and say the result in one plain sentence (no procedure narration): either "this is small enough to just do - proceeding" or "this needs a plan first - investigating now".
4. Proceed under the workflow stages. The developer needs no other command until the very end: you bring them the plan, the artifacts, and the decisions as the workflow directs, and when the change is ready to close you ask them to type /wrap.

Before putting any question to the developer, consult `.claude/DECISIONS.md`.
