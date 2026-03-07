# PHENIX AI Agent — User Guide

---

## 1. What Is the AI Agent?

The PHENIX AI Agent automates macromolecular structure
determination. You give it your experimental data and
(optionally) some guidance, and it figures out which
PHENIX programs to run, in what order, with what
settings — and then runs them for you automatically.

Think of it as an experienced crystallographer sitting
next to you at the computer. It looks at your data,
decides what to do first, checks the results, and
decides what to do next. If something goes wrong, it
tries a different approach. When it's done, it tells
you what it found and where to look.

### What it can do

The AI Agent can solve structures using three major
approaches:

**X-ray crystallography — molecular replacement (MR).**
If you have diffraction data and a similar structure
(a "search model"), the agent can place your model in
the unit cell using Phaser, refine it, rebuild it, and
produce a nearly complete structure. If you don't have
a search model but you have the protein sequence, the
agent can predict a model using AlphaFold and then
proceed with molecular replacement automatically.

**X-ray crystallography — experimental phasing
(SAD/MAD).** If your protein contains selenium,
sulfur, or another anomalous scatterer, the agent can
use the anomalous signal in your data to determine
phases from scratch using AutoSol, then build and
refine the model.

**Cryo-EM.** If you have a cryo-EM map (with or
without an initial model), the agent can sharpen the
map, build a model into it, dock an existing model,
and refine against the map.

In all cases, the agent handles the details:
selecting the right programs, choosing the right
input files, setting appropriate parameters for your
resolution and data quality, and tracking progress
toward a good model.

### What it cannot do

The AI Agent is a powerful assistant, but it has
limitations you should understand:

- **It does not collect data.** You must provide
  your own diffraction data or cryo-EM map.
- **It does not replace your scientific judgment.**
  The agent can build and refine a model, but you
  must inspect the result, check the electron density,
  and decide whether the model is correct. Automated
  model building can sometimes place sidechains
  incorrectly or miss features in the density.
- **It will not always succeed.** Some structures are
  genuinely difficult — poor data quality, low
  resolution, high disorder, no good search model.
  When the agent cannot make progress, it will stop
  and tell you why. This is useful information, not
  a failure.
- **It does not know your biology.** The agent works
  with structural data, not biological context. It
  won't know that a particular ligand should be in the
  active site, or that two chains should form a
  specific interface — unless you tell it in your
  advice.

### How it works (briefly)

Each time the agent runs a PHENIX program, it looks
at the results and decides what to do next. This
cycle of "run a program, check the results, decide
the next step" repeats until the structure is solved
or the agent determines it cannot make further
progress.

At a higher level (when you choose "expert" analysis
depth), the agent creates a multi-stage strategy plan
at the beginning of the session. For example, for a
molecular replacement problem it might plan: (1) check
data quality, (2) place the search model, (3) initial
refinement, (4) model rebuilding, (5) final
refinement. It then tracks progress against this plan,
advancing through stages as goals are met, or trying
a different approach if things aren't working.

You can see the agent's reasoning, its plan, and its
progress in the Agent Progress panel as it runs.

---

## 2. What You Need to Provide

The agent needs your experimental data and, depending
on your experiment, a few additional files. It also
accepts optional advice text where you can tell it
what you know about your project.

### For X-ray molecular replacement

This is the most common scenario. You have diffraction
data and either a search model or a sequence.

**Required — one of these combinations:**

- Diffraction data + search model:
  - A reflection file (`.mtz`, `.sca`, or `.hkl`)
  - A coordinate file for the search model (`.pdb`)

- Diffraction data + sequence (no search model):
  - A reflection file (`.mtz`, `.sca`, or `.hkl`)
  - A sequence file (`.fa`, `.fasta`, `.seq`, or `.dat`)
  - The agent will predict a model using AlphaFold
    and proceed with MR automatically

**Optional:**

- A ligand restraints file (`.cif`) if your structure
  contains a non-standard ligand (e.g., a drug, ATP,
  a cofactor). The agent will fit the ligand into the
  density after refinement.
- Advice text describing what you know: the resolution
  you want to use, whether the data is twinned, which
  ligand to fit, etc.

**Example:**
> Files: `lysozyme_data.mtz`, `search_model.pdb`
> Advice: "Solve by molecular replacement and refine
> to 1.8 Angstrom resolution"

### For X-ray experimental phasing (SAD/MAD)

If your protein was expressed with selenomethionine
or your crystal contains heavy atoms, the anomalous
signal in your diffraction data can be used to
determine phases without a search model.

**Required:**

- A reflection file with anomalous data (`.mtz`,
  `.sca`, or `.hkl`)
- A sequence file (`.fa`, `.fasta`, `.seq`, or `.dat`)
- Advice telling the agent what anomalous scatterer
  to look for — this is important because the agent
  needs to know the atom type

**Optional:**

- A ligand restraints file (`.cif`) if applicable
- Resolution limit or other guidance in the advice

**Example:**
> Files: `p9.sca`, `sequence.dat`
> Advice: "Solve by SAD using selenium as the
> anomalous atom at 2.5 Angstrom resolution"

**Important:** For experimental phasing, you must tell
the agent which heavy atom to use. Unlike molecular
replacement (where the agent can often figure out the
approach from the files alone), SAD/MAD phasing
requires you to specify the anomalous scatterer in
your advice (e.g., "use selenium", "Se-Met protein",
"use S-SAD with sulfur").

### For cryo-EM

If you have a cryo-EM density map and want to build
or refine an atomic model.

**Map input — any of these work:**

The most common input is a pair of **half-maps**
(the two independent reconstructions from your
processing software). From half-maps, the agent can
estimate resolution and validate the final model.
You can also provide a single full map, or both
half-maps and a full map — the agent will use
whatever you give it.

- Two half-maps (`.mrc` or `.ccp4`) — most common
  and recommended
- A single full map (`.mrc` or `.ccp4`)
- Both half-maps and a full map

**Required — one of these in addition to maps:**

- An existing model (`.pdb`) — the agent will refine
  it against the map
- A sequence file (`.fa`, `.fasta`, `.seq`, or
  `.dat`) — the agent will build a model from
  scratch or predict one using AlphaFold

**Optional:**

- A ligand restraints file (`.cif`) if applicable
- Advice about expected contents, symmetry, or
  resolution

**Examples:**
> Files: `half_map_1.mrc`, `half_map_2.mrc`,
> `model.pdb`
> Advice: "Refine the model against the 3.2 Angstrom
> cryo-EM map"

> Files: `half_map_1.mrc`, `half_map_2.mrc`,
> `sequence.fa`
> Advice: "Build a model into the 2.8 Angstrom map"

### Advice: what to write and when

The "Advice for the Agent" text box is optional but
can significantly improve results. The agent reads
your advice and uses it to make better decisions
about which programs to run and how to configure them.

**When advice is most helpful:**

- When the agent can't tell what you want from the
  files alone (e.g., which heavy atom for SAD phasing)
- When you know something special about your data
  (twinned, low completeness, unusual symmetry)
- When you want to limit what the agent does
  ("just run xtriage and stop")
- When you want to fit a specific ligand ("fit ATP
  in the active site")

**When advice is less necessary:**

- For straightforward MR cases where you provide
  data + search model — the agent will do the right
  thing without advice
- For simple cryo-EM refinement with data + model

We'll cover how to write effective advice in detail
in Section 7.

### File formats the agent accepts

The agent recognizes files by their extension:

| File type | Extensions | What it's used for |
|-----------|------------|-------------------|
| Reflection data | `.mtz`, `.sca`, `.hkl` | X-ray diffraction data |
| Map | `.mrc`, `.ccp4`, `.map` | Cryo-EM density maps (full maps or half-maps) |
| Coordinates | `.pdb` | Search models, existing models |
| Sequence | `.fa`, `.fasta`, `.seq`, `.dat` | Protein/nucleic acid sequence |
| Ligand restraints | `.cif` | Non-standard ligand definitions |

You can provide multiple files of the same type —
for example, a protein model and a ligand model, or
both a data file and a pre-calculated map. The agent
will figure out what each file is and how to use it.

**Tip:** If you have a directory with your files and
a `README` text file describing the project, you can
point the agent at the directory. It will read the
README as advice and pick up all the data files
automatically.

[SCREENSHOT: The GUI input panel showing the file
list with several files added, the advice text box
with example text, and the settings panel below]

---

## 3. Setting Up a Run

This section walks you through setting up and
launching an AI Agent run in the PHENIX GUI.

### Step 1: Open the AI Agent

In the PHENIX GUI, select **AI Agent** from the
program list (under "Automation" or by searching).
The AI Agent setup panel will appear.

[SCREENSHOT: PHENIX GUI program list with AI Agent
highlighted]

### Step 2: Add your files

Click the **Add files** button and select your input
files. You can also drag and drop files into the file
list. The agent accepts any combination of the file
types listed in Section 2.

You can remove a file by selecting it in the list and
clicking **Remove selected**.

**Tip:** You don't need to tell the agent which file
is which — it will figure out that `.mtz` is your
reflection data, `.pdb` is your model, and so on.
Just add everything and let the agent sort it out.

**Auto-discovery:** If you set `input_directory` (or
the GUI detects it from a README file), the agent
will automatically discover all crystallographic
files in that directory. This is especially useful
in rules-only mode, where the agent can't use an
LLM to extract file names from a README. You'll see
a log message like:
```
[FILES] Auto-discovered 3 file(s) from /path/to/dir:
  data.mtz
  model.pdb
  sequence.fa
```

**Ligand detection:** If your input model already
contains ligands (HETATM records that aren't water),
the agent will detect this automatically and note it:
```
[FILES] Input model 1aba.pdb contains ligand(s): ATP, MG
```
This enables polder omit map calculations without
requiring a ligandfit step in the session.

[SCREENSHOT: File list with data.mtz, model.pdb, and
seq.fa added]

### Step 3: Write advice (optional)

The **Advice for the Agent** text box is where you
can give the agent guidance in plain English. Good
advice is short and specific:

- "Solve by molecular replacement and refine"
- "Solve at 2.5 Angstrom, use Se as anomalous atom"
- "Refine and fit the ATP ligand"
- "Just run xtriage to check data quality"

You can leave this blank for straightforward cases
(e.g., you provided data + search model and just want
MR + refinement). See Section 7 for detailed advice
on writing effective advice.

**README shortcut:** If your input directory contains
a file named `README`, `README.txt`, or `README.md`,
the agent will read it automatically and use its
contents as advice. Many PHENIX tutorials include a
README — just point the agent at the tutorial
directory.

### Step 4: Choose settings

The settings panel has several controls. For most
runs, the defaults are fine — you only need to change
them if you have a specific reason.

**Analysis Depth** (`thinking_level`)

This is the most important setting. It controls how
deeply the agent analyzes your data at each step:

| Setting | What it does | When to use it |
|---------|-------------|---------------|
| None | Fastest. Minimal analysis, no expert reasoning. | Quick tests, simple workflows |
| Basic | Adds AI reasoning about each step. | Light analysis |
| Advanced | Full structural validation, expert knowledge base, detailed reasoning. | When you don't need strategic planning |
| Expert | (Default) Everything in Advanced, plus a multi-stage strategy plan with progress tracking and model placement detection. | Most runs — recommended |

For your first run, **Expert** (the default) is the
best choice. It gives you the full planning system
with goal tracking, and it includes the model
placement gate that prevents the agent from running
unnecessary molecular replacement on models that
already fit the data.

**Max cycles**

The maximum number of programs the agent will run
(default: 20). Each "cycle" is one PHENIX program
execution. A typical structure determination takes
5-15 cycles. You can increase this for complex cases
or decrease it to limit the agent's work.

**Restart mode**

- **Fresh** (default for new runs): Start a new
  analysis from scratch.
- **Resume**: Continue from where the previous run
  left off. Use this when the agent stopped and you
  want to try again with different advice or settings.

Note: After using "Display session and stop" or
"Remove last N cycles", the restart mode
automatically switches to **Resume** for the next
run. This is the expected behavior — you just
inspected or edited a session, so the natural next
action is to continue it.

**Use LLM**

When checked (default), the agent uses an AI language
model to make intelligent decisions. When unchecked,
the agent uses simple rules to choose programs — this
is faster and works offline, but makes less
sophisticated decisions.

**Provider**

Which AI provider to use (Google, OpenAI, or Ollama).
Google is the default and works well for most cases.
Ollama runs locally on your machine if you have it
installed.

**Verbosity**

How much detail to show in the output:
- **Quiet**: Only errors and the final result.
- **Normal** (default): Key decisions, metrics, and
  reasoning at each step.
- **Verbose**: Full detail including file selection
  and internal state. Useful for debugging.

**Open sub-job windows**

When checked (default), the agent opens a window
for each PHENIX program it runs so you can watch
the progress. Uncheck this if you find the extra
windows distracting.

### Step 5: Job control (optional)

These controls are for managing previous runs:

**Display session and stop**

View the results of a previous run without running
the agent again:
- **None** (default): Run the agent normally.
- **Basic**: Show a compact summary of the previous
  run (one line per cycle).
- **Detailed**: Show full detail including files,
  metrics, and reasoning for each cycle.

This is useful when you want to review what the
agent did in a previous session.

**Remove last N cycles**

Roll back the last N cycles from the previous
session. Use this when the agent made a bad decision
in its last few steps and you want to resume from
an earlier point. Set this to the number of cycles
to remove, then run with **Restart mode = Resume**
and updated advice.

### Step 6: Click Run

Click **Run** to start the agent. The Agent Progress
panel will appear, showing the agent's work in real
time. See Section 4 for how to read this display.

[SCREENSHOT: The complete setup panel with all
sections visible, ready to run]

---

## 4. What You See While It Runs

Once you click Run, the PHENIX GUI switches to the
**Agent Progress** tab. This is a live display that
shows you what the agent is doing, what it's
thinking, and how the structure determination is
progressing.

### The plan header (Expert mode)

Since **Expert** is the default analysis depth, the
first thing you'll see is the agent's strategy plan:

```
[CONFIG] thinking_level=expert — goal-directed
  planning active

==================================================
 STRATEGY PLAN: MR + refinement (standard X-ray)
==================================================
 ○ Stage 1: Analyze data quality  [phenix.xtriage]
          Goal: Data quality analysis complete
 ○ Stage 2: Find MR solution  [phenix.phaser,
              phenix.predict_and_build]
          Goal: TFZ >8, LLG >100
 ○ Stage 3: Initial refinement  [phenix.refine]
          Goal: R-free <0.30
 ○ Stage 4: Rebuild problem regions
              [phenix.autobuild, phenix.refine]
          Goal: R-free <0.28
 ○ Stage 5: Final refinement  [phenix.refine]
          Goal: R-free <0.25
==================================================
```

Each line is a stage of the plan. The symbols show
the status of each stage:

| Symbol | Meaning |
|--------|---------|
| ○ | Pending — hasn't started yet |
| ● | Active — currently being worked on |
| ✓ | Complete — goal was met |
| ⊘ | Skipped — not needed (e.g., model already placed) |
| ✗ | Failed — goal wasn't met |

The **Goal** line tells you what the agent is aiming
for in each stage. For example, "R-free <0.30" means
the agent will keep refining until R-free drops below
0.30, or until it runs out of cycles for that stage.

**Model placement detection:** If the agent discovers
that your model already fits the data (via
`model_vs_data` or a successful refinement), it will
automatically skip any pending MR or phasing phases:

```
[PLACEMENT] Model is placed — phenix.model_vs_data
  CC=0.52 (cycle 2). Destructive programs (phaser,
  autosol) will be suppressed.
[PLACEMENT] Skipped plan stage 'molecular_replacement'
  — model already placed
```

This prevents the agent from destroying a good model
by running unnecessary molecular replacement.

If you chose **Advanced** or lower, you won't see the
plan header — the agent still works correctly, it
just doesn't show the strategic overview.

### Cycle display

Each program the agent runs appears as a numbered
cycle:

```
──── Stage: Initial refinement (cycle 2, up to 3) ────
  Goal: R-free <0.30 (current: 0.312)

Cycle 5: phenix.refine

  Reasoning:
    R-free has dropped from 0.35 to 0.31 over the
    last two cycles. Continuing refinement with
    ordered solvent to improve the model further.

  $ phenix.refine model.pdb data.mtz
  Running...
```

Here's what each part tells you:

**Stage context** (Expert mode): Shows which stage
of the plan this cycle belongs to, how many cycles
have been used in this stage, and how close the
current metric is to the goal.

**Cycle number and program**: "Cycle 5: phenix.refine"
means this is the 5th program the agent has run, and
it chose phenix.refine.

**Reasoning**: The agent's explanation of why it chose
this program. This helps you understand the agent's
decision-making.

**Command**: The actual PHENIX command being run.
File paths are shortened to just the file name for
readability.

**Running...**: The program is executing. Depending
on the program and your data, this may take seconds
(xtriage) to hours (autobuild at high resolution).

### Expert Assessment

When using **Advanced** or **Expert** analysis depth,
the agent shows a detailed assessment after each
cycle:

```
  [Expert Assessment (guide_step, high)]
  Structure:
    R-free: 0.312 (improved from 0.350)
    R-work: 0.275
    Clashscore: 8.2
    Ramachandran favored: 96.5%
  Problems: 3 residues in disallowed regions

  The model is improving steadily. R-free dropped
  significantly this cycle. Consider running one
  more round of refinement before attempting
  model rebuilding, as the current map quality
  should support autobuild.

  [Expert Guidance]
  Continue refinement. The R-free trajectory
  suggests convergence within 1-2 more cycles.
```

**How to read the metrics:**

| Metric | What it means | Good | Concerning |
|--------|-------------|------|-----------|
| R-free | Agreement between model and data (lower is better) | < 0.25 | > 0.35 |
| R-work | Like R-free but on working data (always lower than R-free) | < 0.20 | > 0.30 |
| Clashscore | Number of atomic clashes per 1000 atoms (lower is better) | < 5 | > 20 |
| Ramachandran favored | Percentage of residues in ideal backbone geometry | > 97% | < 90% |
| Map CC | Map-model correlation for cryo-EM (higher is better) | > 0.7 | < 0.5 |

Don't worry if you don't understand all the metrics —
the Expert Assessment text explains what they mean for
your specific case. Focus on R-free: if it's going
down, things are working.

**The confidence and action indicators** (shown in
brackets after "Expert Assessment") tell you the
agent's overall assessment:

- **Confidence**: `high`, `medium`, `low`, or
  `hopeless` — how confident the agent is about the
  current state.
- **Action**: `guide_step` (continue normally),
  `let_run` (agent is just monitoring), `stop`
  (agent thinks it's done), `pivot` (agent wants
  to try a different approach).

### Stage transitions (Expert mode)

When the agent completes a stage of its plan, you'll
see a transition block:

```
==================================================
 ✓ STAGE COMPLETE: initial_refinement
   R-free reached target (0.28 < 0.30)
 → ADVANCING TO: model_rebuilding
==================================================
```

Or if the agent decides to try a different approach:

```
==================================================
 ⚠ RETREAT: initial_refinement → molecular_replacement
   R-free stuck above 0.45 after 3 cycles
   Strategy blacklisted: phaser_default
==================================================
```

A retreat is not a failure — it means the agent
recognized that its current approach isn't working
and is backtracking to try something different. This
is exactly what an experienced crystallographer would
do.

### Cycle results

After each program finishes, you'll see the result:

```
  >> Sub-job: [OK] phenix.refine
Cycle 5: phenix.refine -> COMPLETE
  Result: OK
```

Or if something went wrong:

```
  >> Sub-job: [FAILED] phenix.phaser
Cycle 3: phenix.phaser -> FAILED
  Result: FAILED: No MR solution found
```

A single failed cycle doesn't necessarily mean the
whole run has failed — the agent may recover
automatically by trying different parameters or a
different approach.

### When the agent finishes

When the agent completes, you'll see either a summary
of what was accomplished, or an explanation of why it
stopped. The **Results** tab will contain a structured
summary of the entire run, including all cycles and
their outcomes.

[SCREENSHOT: The Agent Progress panel during a run,
with annotations pointing to the plan header, a cycle
display, an Expert Assessment block, and a stage
transition block]

---

## 5. When Things Go Well

A successful run follows a recognizable pattern:
metrics improve over cycles, the agent advances
through its plan stages, and the final model has good
geometry and reasonable R-free (for X-ray) or map
correlation (for cryo-EM).

### What success looks like

You'll see R-free dropping cycle by cycle:

```
Cycle 3: phenix.refine   R-free: 0.350
Cycle 5: phenix.refine   R-free: 0.295
Cycle 7: phenix.refine   R-free: 0.262
Cycle 8: phenix.refine   R-free: 0.248
```

In Expert mode, the agent will advance through plan
stages with green checkmarks:

```
 ✓ Stage 1: data_assessment
 ✓ Stage 2: molecular_replacement
 ✓ Stage 3: initial_refinement
 ● Stage 4: model_rebuilding  ← currently here
 ○ Stage 5: final_refinement
```

When the agent finishes, it produces a summary of the
final model quality and where to find your files.

### Where to find your output

All output files are in the agent's working directory
(usually `ai_agent_directory/` inside your project
folder, or the directory shown in the log at the end
of the run). Key files:

- **Refined model**: The most recent `.pdb` file
  produced by the last refinement cycle. Look for the
  file with the highest cycle number.
- **Map coefficients**: The `.mtz` file from the last
  refinement, containing 2Fo-Fc and Fo-Fc map
  coefficients for viewing in Coot.
- **Session log**: `session.json` contains the full
  history of every cycle, decision, and metric.
- **Structure report** (Expert mode):
  `structure_determination_report.txt` — a text
  summary of the entire determination.
- **Session summary** (Expert mode):
  `session_summary.json` — a machine-readable summary
  with final metrics, stage outcomes, and hypothesis
  results.

### What to do next: the Top 5 Coot checklist

The agent produces a model, but you are the scientist.
Before depositing or publishing, always inspect the
model manually. Here are the five most important
things to check in Coot:

**1. Inspect the difference density map.** Display the
Fo-Fc map (usually the green/red map). Look for large
positive peaks (green blobs) near your model — these
may indicate missing atoms, alternate conformations,
unmodeled ligands, or ions. Large negative peaks (red)
may indicate atoms that shouldn't be there.

**2. Check Ramachandran outliers.** The agent reports
these in the Expert Assessment. In Coot, go to
Validate → Ramachandran Plot and navigate to each
outlier. Decide whether the local electron density
supports the unusual backbone conformation, or
whether the residue should be refit.

**3. Verify ligand placement.** If the agent fit a
ligand, examine it in the density. The ligand should
sit cleanly in a blob of difference density. If the
real-space correlation coefficient (RSCC) reported by
the agent is below 0.7, the ligand placement may be
questionable — inspect it closely.

**4. Look at B-factor distribution.** Regions with
very high B-factors (coloring the model by B-factor
in Coot makes this easy to spot) may be poorly ordered
or incorrectly built. Consider whether these regions
should be truncated or rebuilt.

**5. Check crystal contacts and special positions.**
Use Coot's symmetry display to check that the model
makes sensible crystal contacts. Verify that any
molecules on special positions are handled correctly
(correct occupancy, no clashes with symmetry mates).

---

## 6. When Things Go Wrong

Structure determination doesn't always succeed, and
that's normal. The AI Agent is designed to recognize
problems and either try a different approach or stop
with a clear explanation. This section covers the most
common issues and what to do about them.

### "No matching plan template"

**What it means:** In Expert mode, the agent couldn't
find a strategy plan that matches your combination of
files and advice. This typically means you have an
unusual setup that the agent doesn't have a pre-built
strategy for.

**What to do:**
- Check that you've provided the right files. For
  example, SAD phasing needs advice specifying the
  anomalous atom type.
- The agent will still run in "reactive mode" (choosing
  programs one at a time without a plan), so this isn't
  fatal — you just won't see the phase-tracking
  features.

### "R-free stuck above 0.40"

**What it means:** The model isn't fitting the data
well. After molecular replacement, R-free values above
0.40 often indicate a problem with the MR solution
itself.

**Possible causes:**
- Wrong search model (too different from your protein)
- Wrong space group
- Data quality issues (low completeness, twinning)
- Wrong number of copies in the asymmetric unit

**What to do:**
- Check the xtriage output for warnings about space
  group, twinning, or data quality
- Try a different search model (a closer homolog, or
  use AlphaFold prediction by providing only the
  sequence)
- If the space group is uncertain, try alternative
  space groups
- Consult with a more experienced crystallographer

### "Trying a different approach" (RETREAT)

**What it means:** The agent decided that its current
strategy isn't working and is going back to an earlier
phase to try something different. For example, if
refinement isn't converging, the agent might go back
and try a different molecular replacement solution.

**Is this bad?** No — this is actually one of the
agent's most useful behaviors. It's doing exactly what
an experienced crystallographer would do: recognizing
that the current approach is stuck and trying
something different. You'll see a message like:

```
 ⚠ RETREAT: initial_refinement → molecular_replacement
   R-free stuck above 0.45 after 3 cycles
```

The agent will then try a different strategy for
molecular replacement (e.g., different search model
processing) before attempting refinement again.

### "Safety Stop"

**What it means:** The agent detected something
fundamentally wrong with the workflow — something that
can't be fixed by trying a different program. This is
a protective measure to prevent wasting time on a
situation that needs human intervention.

**Common causes:**
- The experiment type appeared to change mid-run
  (e.g., the data looked like X-ray at first but
  then cryo-EM files appeared)
- A critical input file is missing or corrupted
- The workflow reached an impossible state

**What to do:** Read the error message carefully —
it will tell you what went wrong. Fix the underlying
issue (usually a file problem) and start a fresh run.

### "Strategy Stop" / "Agent stopped"

**What it means:** The agent finished its plan or
decided that no further improvement is possible. In
Expert mode, you'll see the final phase status:

```
 ■ GATE STOP
   All phases complete. Final R-free: 0.231
```

This is usually a good outcome — the agent achieved
its goals. Check the final metrics and proceed to
manual inspection in Coot.

If the agent stopped because it *couldn't* make
progress (rather than because it finished), the
message will explain why:

```
 ■ GATE STOP
   R-free not improving after 3 retreat attempts
```

### "Loop Detection Stop"

**What it means:** The agent detected that it was
running the same program repeatedly without making
progress. This safety mechanism prevents the agent
from wasting cycles (and your time) on a stuck
workflow.

**What to do:** Look at the last few cycles' metrics.
If R-free is flat, the model may have converged at a
local minimum. Try providing different advice or a
different search model.

### Troubleshooting flowchart

When the agent stops and you're not sure why, follow
this decision tree:

```
The agent stopped. What do I do?
    |
    +-- Is there a "Safety Stop" message?
    |   +-- Yes: Read the error message.
    |   |   Fix the input problem, then re-run.
    |   +-- No:
    |       |
    +-- Did the agent produce a model?
    |   +-- Yes: Check the R-free value.
    |   |   +-- R-free < 0.30: Probably a good
    |   |   |   result. Inspect in Coot (Section 5).
    |   |   +-- R-free 0.30-0.40: Usable but may
    |   |   |   need manual improvement.
    |   |   +-- R-free > 0.40: MR solution may be
    |   |       wrong. Try different search model.
    |   +-- No:
    |       |
    +-- Did the agent retreat multiple times?
    |   +-- Yes: Data may have issues. Check
    |   |   xtriage for twinning, anisotropy,
    |   |   or low completeness.
    |   +-- No:
    |       |
    +-- Is this a SAD/MAD experiment?
    |   +-- Yes: Check that you specified the
    |   |   anomalous atom type in your advice.
    |   |   Check the anomalous signal strength
    |   |   in the xtriage output.
    |   +-- No:
    |       |
    +-- Try re-running with:
        - More specific advice
        - A different search model
        - Ask a colleague for help
```

---

## 7. Giving Better Advice

The advice you give the agent can make a significant
difference in the quality and speed of the results.
Here's how to write advice that helps.

### Before and after: advice examples

| What you want | Vague (less effective) | Specific (more effective) |
|---------------|----------------------|--------------------------|
| Solve by MR | "Solve my structure" | "Solve by MR using the AlphaFold prediction. Target R-free below 0.25." |
| SAD phasing | "Use anomalous data" | "Solve by SAD using selenium as the anomalous atom at 2.5 A resolution." |
| Fit a ligand | "There's a ligand" | "Fit ATP in the active site after refinement." |
| Just check data | "Check my data" | "Run xtriage and stop." |
| Twinned data | "Something is weird with the data" | "Data appears to be twinned with twin law h,-k,-l." |
| Cryo-EM rebuild | "Build a better model" | "Rebuild the model into the 3.0 A cryo-EM map, focus on chain A." |

### What makes advice effective

- **Name the method.** "Molecular replacement", "SAD
  phasing", "cryo-EM refinement" — this helps the
  agent choose the right overall strategy.
- **Name specific programs if you know them.**
  "Use phaser for MR" or "run autobuild after
  refinement" — this gives the agent a strong hint.
- **Mention ligands by their three-letter code.**
  "Fit ATP" or "include NAD as a ligand" — the agent
  needs to know the ligand identity to fit it.
- **Set the resolution.** "Refine to 2.5 Angstrom"
  — this helps the agent choose appropriate parameters.
- **Mention special circumstances.** "The data is
  twinned", "the protein is a Se-Met derivative",
  "the crystal has P3₁21 symmetry" — this prevents
  the agent from having to guess.
- **Tell it when to stop.** "Run xtriage and stop"
  or "stop after one round of refinement" — the agent
  respects explicit stop conditions.

### What doesn't help

- **Too vague:** "Make it work" or "solve the
  structure" — these give the agent almost nothing to
  work with.
- **Raw PHIL parameters:** "refinement.main.
  number_of_macro_cycles=10" — the agent doesn't parse
  PHIL strings from advice. Use plain language instead:
  "use 10 macro-cycles of refinement".
- **Impossible requests:** Asking for things that
  PHENIX programs can't do.
- **Very long text:** Keep advice to 3-4 sentences.
  The agent extracts structured directives from your
  advice — long rambling text is harder to parse
  accurately.

---

## 8. Resuming and Modifying Runs

### Resuming a stopped run

If the agent stopped and you want to try again with
different advice:

1. In the setup panel, set **Restart mode** to
   **Resume**
2. Update your advice text (e.g., "try a different
   space group" or "use the predicted model instead")
3. Click **Run**

The agent will pick up where it left off, using the
new advice for its next decisions. All previous cycles
are preserved.

### Rolling back bad cycles

If the agent made a bad decision in its last few
cycles (e.g., it fit a ligand into the wrong density
and then refined with it), you can remove those cycles
before resuming:

1. Set **Remove last N cycles** to the number of
   cycles you want to undo (e.g., 3)
2. Set **Restart mode** to **Resume**
3. Update your advice if needed
4. Click **Run**

The agent will remove the last N cycles from the
session history and then continue from the earlier
state with your new advice.

### Viewing a previous run

To review what the agent did without running it again:

1. Set **Display session and stop** to **basic**
   (compact summary) or **detailed** (full detail)
2. Click **Run**

The agent will display the session history and then
stop. No new programs are run. This is useful when
you want to check the status of a run that finished
while you were away, or to review the agent's
decisions before deciding whether to resume.

---

## 9. Understanding the Output Files

After the agent finishes, all output files are in the
agent's working directory. Here's what you'll find:

### Key files to look for

**Your refined model:** The `.pdb` file from the last
successful refinement cycle. This is in the
sub-directory for that cycle's job (e.g.,
`phenix_refine_01/`). The Results panel in the GUI
will also show you where each output file is.

**Map coefficients:** The `.mtz` file from the last
refinement cycle, containing calculated maps. Open
this in Coot alongside the refined model to inspect
the electron density.

**Session data:** `ai_agent_directory/session.json`
contains the complete history of every decision the
agent made, every metric it measured, and every file
it produced. This is primarily for the agent's own
use (for resuming), but it can be useful for
understanding what happened during a run.

### Expert mode additional outputs

When running with Expert analysis depth (the default),
the agent produces additional output files:

**HTML structure report:**
`ai_agent_directory/structure_report.html`

A self-contained HTML report with:
- Outcome status (determined / stopped / incomplete)
- Final model metrics (R-free, R-work, geometry)
- Inline SVG trajectory chart showing R-free or
  Map CC progression with retreat markers
- Phase timeline table
- Output file locations

Click the **Open Structure Report** button in the
Results tab to view this report in your browser.

**Text structure report:**
`ai_agent_directory/structure_determination_report.txt`

A human-readable text summary with:
- Data characteristics (resolution, space group)
- Final model quality (R-free, R-work, geometry)
- Phase timeline (which stages completed, how many
  cycles each took)
- Hypotheses tested (if any)

**Session summary JSON:**
`ai_agent_directory/session_summary.json`

A machine-readable summary containing:
- Final metrics (R-free, R-work, clashscore, etc.)
- Data characteristics
- Phase outcomes (complete, skipped, failed)
- Metric trajectory (R-free at each cycle)
- Outcome classification (complete, stopped,
  incomplete, reactive)

This file is useful if you want to process results
programmatically or compare multiple runs.

### Sub-job directories

Each PHENIX program the agent runs creates its own
sub-directory within the working directory. These
contain the full output of each program — log files,
output models, map files, etc. The agent's Results
panel shows which sub-directory corresponds to each
cycle.

---

## 10. Glossary

**B-factor (temperature factor).** A number assigned
to each atom that describes how much it vibrates or
is disordered. Higher B-factors mean more disorder.
Measured in Å².

**Clashscore.** The number of serious atomic clashes
(atoms too close together) per 1000 atoms in the
model. Lower is better. A good model has a clashscore
below 5.

**Cryo-EM (cryo-electron microscopy).** A structural
biology technique that images frozen protein samples
using an electron microscope, producing 3D density
maps.

**Electron density map.** A 3D map showing where
electrons (and therefore atoms) are located in the
crystal or cryo-EM sample. The atomic model is built
to fit this map.

**Ligand.** A small molecule bound to a protein —
for example, a drug, a cofactor (like ATP or NAD),
or a metal ion. Ligands require special restraint
files (`.cif`) to be modeled correctly.

**Map CC (map correlation coefficient).** A measure
of how well the atomic model matches the cryo-EM
density map. Ranges from 0 (no agreement) to 1
(perfect agreement). Values above 0.7 are good.

**Map coefficients.** Numbers stored in an `.mtz`
file that can be used to calculate electron density
maps. The "2Fo-Fc" map shows the overall density;
the "Fo-Fc" (difference) map shows where the model
disagrees with the data.

**Molecular replacement (MR).** A method for
determining the arrangement of molecules in a crystal
by using a known similar structure as a starting
point. The similar structure is called the "search
model."

**Phase (of a plan).** When the agent runs in Expert
mode, it divides the structure determination into
phases — for example, "data assessment", "molecular
replacement", "refinement". Each phase has a specific
goal.

**R-free.** The most important quality metric for an
X-ray crystal structure. It measures how well the
model predicts diffraction data that was NOT used
during refinement. Ranges from 0 (perfect) to about
0.6 (random). A good structure has R-free below 0.25;
below 0.30 is acceptable; above 0.40 suggests a
serious problem.

**R-work.** Like R-free, but measured on the data
that WAS used during refinement. R-work is always
lower than R-free. The gap between them (R-free minus
R-work) should be small (< 0.05); a large gap
suggests overfitting.

**Ramachandran plot.** A plot of protein backbone
angles (phi and psi) for each residue. Most residues
should fall in "favored" regions of the plot. Outliers
may indicate errors in the model.

**Resolution.** How much detail is visible in the
experimental data, measured in Ångströms (Å). Lower
numbers mean more detail: 1.5 Å is high resolution
(individual atoms visible), 3.0 Å is moderate, 4.0 Å
is low.

**Retreat.** When the AI Agent decides its current
approach isn't working and goes back to an earlier
phase to try a different strategy.

**SAD (Single-wavelength Anomalous Dispersion).** A
method for determining crystal structure phases using
the anomalous (wavelength-dependent) scattering of
heavy atoms in the crystal, such as selenium in
selenomethionine-labeled proteins.

**Safety Stop.** When the AI Agent detects a
fundamental problem that requires human intervention,
such as a missing input file or a contradictory
workflow state.

**Search model.** A known protein structure (usually
a homolog of your protein) used as the starting point
for molecular replacement.

**Space group.** The symmetry of the crystal lattice.
Determines how many copies of the molecule are in the
unit cell and how they are related.

**Strategy Stop.** When the AI Agent completes its
plan or determines that no further improvement is
possible.

---

*This guide covers PHENIX AI Agent v114. For the
latest information, see the PHENIX documentation
at https://phenix-online.org.*

*Screenshots will be added in a future revision
(Phase B).*
