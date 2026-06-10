"""
Intent Classifier for PHENIX AI Agent (v115).

Classifies user advice into one of four intent categories
that drive the agent's stopping behavior:

  solve             → full workflow to convergence
  solve_constrained → full workflow, method locked
  task              → single program, stop after
  tutorial          → follow described steps, stop when done

Usage:
    from agent.intent_classifier import classify_intent

    result = classify_intent(user_advice)
    # result = {
    #   "intent": "task",
    #   "confidence": "high",
    #   "task_program": "phenix.xtriage",
    #   "method_constraint": None,
    #   "reason": "Explicit single-program command",
    # }
"""

from __future__ import absolute_import, division, print_function

import re


# ================================================================
# Intent categories
# ================================================================

SOLVE = "solve"
SOLVE_CONSTRAINED = "solve_constrained"
TASK = "task"
TUTORIAL = "tutorial"

# ================================================================
# Program name patterns (for task detection)
# ================================================================

# Canonical program names and their common aliases
_PROGRAM_ALIASES = {
    "phenix.xtriage": [
        r"xtriage", r"data\s+quality",
        r"twinning", r"twinning\s+analysis",
    ],
    "phenix.phaser": [
        r"phaser", r"molecular\s+replacement\b",
        r"\bMR\b",
    ],
    "phenix.refine": [
        r"\brefine\b", r"refinement",
    ],
    "phenix.real_space_refine": [
        r"real.?space.?refin",
    ],
    "phenix.autobuild": [
        r"autobuild", r"auto.?build",
        r"model\s+building",
    ],
    "phenix.autosol": [
        r"autosol", r"auto.?sol",
    ],
    "phenix.ligandfit": [
        r"ligand\s*fit", r"fit.*ligand",
        r"place.*ligand", r"ligand\s+fitting",
    ],
    "phenix.polder": [
        r"polder", r"omit\s+map",
    ],
    "phenix.molprobity": [
        r"molprobity", r"validat(?:e|ion)",
    ],
    "phenix.mtriage": [
        r"mtriage",
    ],
    "phenix.map_to_model": [
        r"map.?to.?model",
    ],
    "phenix.map_sharpening": [
        r"map.?sharpen", r"sharpen.*map",
        r"auto.?sharpen",
    ],
    "phenix.dock_in_map": [
        r"dock.*(?:in|into)\s+map",
    ],
    "phenix.resolve_cryo_em": [
        r"resolve.?cryo",
        r"resolve_cryo_em",
        r"density\s+modif",
    ],
    "phenix.map_symmetry": [
        r"map.?symmetry",
    ],
    "phenix.predict_and_build": [
        r"predict.?and.?build",
    ],
    "phenix.process_predicted_model": [
        r"process.?predict",
    ],
    "phenix.ensemble_refinement": [
        r"ensemble.?refin",
    ],
}

# ================================================================
# Task trigger patterns
# ================================================================

# These patterns indicate a single-program command.
# The user wants ONE thing done, then stop.
_TASK_TRIGGERS = [
    # Imperative: "run X", "just run X", "only run X"
    r"(?:just\s+|only\s+)?run\s+(?:the\s+)?",
    # Imperative: "do X", "just do X"
    r"(?:just\s+|only\s+)?do\s+(?:the\s+)?",
    # Imperative: "check X", "check for X"
    r"check\s+(?:for\s+)?(?:the\s+)?",
    # Imperative: "calculate X"
    r"calculate\s+(?:the\s+)?(?:a\s+)?",
    # Imperative: "analyze X"
    r"analyze\s+(?:the\s+)?",
]

# ================================================================
# Method constraint patterns (solve_constrained)
# ================================================================

# Order matters: more specific patterns first.
# MR-SAD must come before MR to avoid false match.
_METHOD_CONSTRAINTS = [
    ("mr_sad", [
        r"(?:solve|determine).*(?:with|using|by|via)"
        r"\s+MR[- ]?SAD",
        r"(?:use|try)\s+MR[- ]?SAD",
    ]),
    ("predict_and_build", [
        r"(?:solve|determine).*(?:with|using|by|via)"
        r"\s+predict.?and.?build",
        r"(?:use|try)\s+predict.?and.?build",
    ]),
    ("molecular_replacement", [
        r"(?:solve|determine).*(?:with|using|by|via)"
        r"\s+(?:molecular\s+replacement|phaser)",
        r"(?:use|try)\s+molecular\s+replacement",
        # MR alone — but not MR-SAD (already matched above)
        r"(?:solve|determine).*(?:with|using|by|via)"
        r"\s+MR(?![- ]?SAD)\b",
        r"(?:use|try)\s+MR(?![- ]?SAD)\b",
    ]),
    ("sad_phasing", [
        r"(?:solve|determine).*(?:with|using|by|via)"
        r"\s+(?:SAD|Se-?SAD|S-?SAD)",
        r"(?:use|try)\s+SAD\s+phasing",
        r"(?:solve|determine).*(?:with|using|by)"
        r"\s+anomalous",
    ]),
    ("mad_phasing", [
        r"(?:solve|determine).*(?:with|using|by|via)"
        r"\s+MAD",
        r"(?:use|try)\s+MAD\s+phasing",
    ]),
]

# ================================================================
# Tutorial detection patterns
# ================================================================

# Preprocessor signature headers (strong signal)
_PREPROCESSOR_SIGS = [
    r"input\s+files?\s+found",
    r"experiment\s+type",
    r"key\s+parameters?",
    r"special\s+instructions?",
]

# Tutorial-style language (moderate signal)
_TUTORIAL_PHRASES = [
    r"in\s+this\s+tutorial",
    r"this\s+tutorial\s+(?:will|shows|demonstrates)",
    r"the\s+goal\s+(?:of\s+this|is\s+to)",
    r"you\s+will\s+learn",
    r"step\s+\d+\s*[:.]\s",
    r"(?:first|next|then|finally)[,:]?\s+run\s",
    r"the\s+purpose\s+(?:of\s+this|is)",
    r"this\s+exercise",
    r"this\s+example\s+(?:will|shows|demonstrates)",
]

# ================================================================
# Solve detection patterns
# ================================================================

_SOLVE_PHRASES = [
    r"solve\s+(?:the\s+)?structure",
    r"determine\s+(?:the\s+)?structure",
    r"process\s+(?:these\s+|the\s+)?files",
    r"figure\s+(?:out|it\s+out)",
    r"what\s+(?:can|should)\s+(?:you|I|we)\s+do",
    r"do\s+(?:your|the)\s+(?:best|thing)",
]

# Phrases that signal a multi-step pipeline rather than a
# single-program task.  When a task is detected AND one of
# these patterns is present, the advice describes a sequence
# of programs and should be classified as solve, not task.
_PIPELINE_INDICATORS = [
    r"\bfollowed\s+by\b",
    r"\bafter\s+(?:that|this|which)\b",
    r"\bsubsequently\b",
    # Note: r"\bnext[,\s]" excluded -- fires on "the next step
    # is xtriage" (next as prefix, not sequencing connector).
    # "then <action-verb>\" — covers \"then build\", \"then refine\",
    # "then run", "then dock", etc. without matching "then stop".
    r"\bthen\s+(?:build|refine|dock|fit|run|do|apply|use)\b",
    # Fix 2 (v116): "and build/autobuild/refine" — OpenAI frequently writes
    # "Run AutoSol ... and build an initial model" which fires the autosol
    # task trigger but is clearly a multi-step workflow.  Without this,
    # intent=task stops the agent after autosol, skipping autobuild/refinement.
    r"\band\s+(?:build|autobuild|refine)\b",
    r"\bwith\s+autobuild\b",
    # "model building" in the advice always implies a full pipeline
    # (e.g. "perform substructure determination, phasing, and initial model building")
    r"\bmodel\s+building\b",
]

def classify_intent(user_advice):
    """Classify user advice into an intent category.

    Args:
        user_advice: The user's advice/instructions string.
            Can be raw user input, preprocessed README, or
            empty/None.

    Returns:
        dict with keys:
            intent: One of "solve", "solve_constrained",
                "task", "tutorial".
            confidence: "high" | "medium" | "low"
            task_program: str or None — for task intent,
                the specific program requested.
            method_constraint: str or None — for
                solve_constrained, the method to use.
            reason: str — human-readable explanation.

    Never raises.
    """
    # --- No advice → solve ---
    if not user_advice:
        return _result(
            SOLVE, "high", None, None,
            "No advice provided — solve the structure")

    # Type guard: non-string input treated as no advice
    if not isinstance(user_advice, str):
        return _result(
            SOLVE, "low", None, None,
            "Non-string advice — defaulting to solve")

    if not user_advice.strip():
        return _result(
            SOLVE, "high", None, None,
            "No advice provided — solve the structure")

    advice = user_advice.strip()
    advice_lower = advice.lower()

    # --- 0. Preprocessor output is ALWAYS tutorial ---
    # If we see preprocessor signature headers, the text is
    # a structured README summary, not a user command.
    # Check this FIRST before any pattern matching.
    sig_count = sum(
        1 for sig in _PREPROCESSOR_SIGS
        if re.search(
            r"(?i)^[ \t]*%s\s*:" % sig,
            advice, re.MULTILINE))
    if sig_count >= 2:
        return _result(
            TUTORIAL, "high", None, None,
            "Preprocessor output detected "
            "(%d signature headers)" % sig_count)

    # --- 0b. Strong tutorial opening phrases ---
    # "In this tutorial..." is always tutorial, regardless of
    # any method mentions in the text.  These are descriptive
    # texts, not user commands.
    _STRONG_TUTORIAL_OPENERS = [
        r"\bin\s+this\s+tutorial\b",
        r"\bthis\s+tutorial\s+(?:will|shows|demonstrate)",
        r"\bthe\s+goal\s+of\s+this\b",
        r"\bthis\s+exercise\b",
        r"\bthis\s+example\s+(?:will|shows)",
    ]
    for opener in _STRONG_TUTORIAL_OPENERS:
        if re.search(opener, advice_lower):
            return _result(
                TUTORIAL, "high", None, None,
                "Strong tutorial opener detected")

    # --- 0c. Step-by-step instructions ---
    # Multiple "Step N:" lines indicate a tutorial procedure,
    # even though individual steps may look like task commands.
    step_count = len(re.findall(
        r"(?:^|\n)\s*step\s+\d+\s*[:.]\s",
        advice_lower))
    if step_count >= 2:
        return _result(
            TUTORIAL, "high", None, None,
            "Step-by-step instructions detected "
            "(%d steps)" % step_count)

    # --- 1. Check for single-program task (highest priority) ---
    task_prog = _detect_task(advice_lower)
    if task_prog:
        # Before declaring this a single-step task, check for
        # pipeline indicators.  "Sharpen the map, then build a
        # model" names a task program but is a multi-step workflow
        # — the agent should not stop after sharpening.
        is_pipeline = any(
            re.search(p, advice_lower)
            for p in _PIPELINE_INDICATORS)
        if not is_pipeline:
            return _result(
                TASK, "high", task_prog, None,
                "Explicit single-program command: %s"
                % task_prog)

    # --- 2. Check for method-constrained solve ---
    method, constraint_reason = _detect_method_constraint(
        advice_lower)
    if method:
        return _result(
            SOLVE_CONSTRAINED, "high", None, method,
            constraint_reason)

    # --- 3. Check for tutorial-style text ---
    is_tutorial, tut_reason = _detect_tutorial(
        advice, advice_lower)
    if is_tutorial:
        return _result(
            TUTORIAL, "medium", None, None,
            tut_reason)

    # --- 4. Check for explicit solve phrases ---
    for pattern in _SOLVE_PHRASES:
        if re.search(pattern, advice_lower):
            return _result(
                SOLVE, "high", None, None,
                "Explicit solve request")

    # --- 5. Default: if there's substantial text,
    #         treat as tutorial; if short, treat as solve ---
    # Short advice (< 50 chars) without clear intent
    # is likely a brief instruction like "use nproc=4"
    # which should not change the workflow.
    # Longer text likely describes a procedure.
    line_count = len([
        l for l in advice.split("\n")
        if l.strip()])
    if line_count >= 5:
        return _result(
            TUTORIAL, "low", None, None,
            "Multi-line text — likely procedural "
            "description (%d lines)" % line_count)

    # Default to solve — if the user gave brief advice
    # but no clear intent, let the agent decide freely.
    return _result(
        SOLVE, "low", None, None,
        "No clear intent pattern — defaulting to solve")


# ================================================================
# Detection helpers
# ================================================================

def _detect_task(advice_lower):
    """Detect single-program task intent.

    Returns the program name if a task is detected,
    None otherwise.

    Three detection paths:
    1. Trigger + alias: "run xtriage", "check for twinning"
    2. Trigger + phenix.X: "run phenix.map_sharpening"
    3. Direct alias: "sharpen the map", "auto_sharpen"
       (short imperative commands without a trigger prefix)
    """
    # Path 1: trigger prefix + program alias
    for trigger in _TASK_TRIGGERS:
        m = re.search(trigger, advice_lower)
        if not m:
            continue
        rest = advice_lower[m.end():]
        # 1a: Check program aliases
        for program, aliases in _PROGRAM_ALIASES.items():
            for alias in aliases:
                if re.match(alias, rest, re.IGNORECASE):
                    return program
        # 1b: Check "phenix.X" directly after trigger
        pm = re.match(r"phenix\.(\w+)", rest)
        if pm:
            from agent.directive_extractor import (
                _fix_program_name)
            fixed = _fix_program_name(pm.group(0))
            if fixed:
                return fixed

    # Path 2: direct alias match (no trigger prefix)
    # Only for short advice (under 80 chars) that
    # doesn't contain solve/structure language.
    if (len(advice_lower) < 80
            and "structure" not in advice_lower
            and "solve" not in advice_lower):
        for program, aliases in _PROGRAM_ALIASES.items():
            for alias in aliases:
                if re.search(alias, advice_lower,
                             re.IGNORECASE):
                    return program

    return None


def _detect_method_constraint(advice_lower):
    """Detect a method constraint in solve-style advice.

    Returns (method_name, reason) or (None, None).
    """
    for method, patterns in _METHOD_CONSTRAINTS:
        for pattern in patterns:
            if re.search(pattern, advice_lower,
                         re.IGNORECASE):
                return method, (
                    "Solve-constrained: method=%s"
                    % method)
    return None, None


def _detect_tutorial(advice, advice_lower):
    """Detect tutorial-style text.

    Returns (is_tutorial, reason) tuple.
    """
    # Strong signal: preprocessor signature headers
    sig_count = sum(
        1 for sig in _PREPROCESSOR_SIGS
        if re.search(
            r"(?i)^[ \t]*%s\s*:" % sig,
            advice, re.MULTILINE))
    if sig_count >= 2:
        return True, (
            "Preprocessor output detected "
            "(%d signature headers)" % sig_count)

    # Moderate signal: tutorial-style phrases
    phrase_count = sum(
        1 for phrase in _TUTORIAL_PHRASES
        if re.search(phrase, advice_lower))
    if phrase_count >= 2:
        return True, (
            "Tutorial-style language detected "
            "(%d phrases)" % phrase_count)

    # Single strong tutorial phrase is enough
    for phrase in _TUTORIAL_PHRASES[:3]:
        if re.search(phrase, advice_lower):
            return True, (
                "Tutorial phrase detected: "
                + phrase[:40])

    return False, ""


# ================================================================
# Result builder
# ================================================================

def _result(intent, confidence, task_program,
            method_constraint, reason):
    """Build an intent classification result dict."""
    return {
        "intent": intent,
        "confidence": confidence,
        "task_program": task_program,
        "method_constraint": method_constraint,
        "reason": reason,
    }
