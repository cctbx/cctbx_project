"""
Natural-language → PHIL parameter converter.

Reads ``natural_language`` entries from the ``strategy_flags`` section of
programs.yaml and converts plain-English phrases in user guidelines into
valid PHIL ``key=value`` strings that can be appended to PHENIX commands.

Because the flag template (e.g. ``main.number_of_macro_cycles={value}``)
comes directly from programs.yaml — the same source as the command builder
uses — the generated PHIL is always structurally correct for that program.

Usage::

    from agent.nl_to_phil import extract_nl_params

    params = extract_nl_params(
        guidelines="Run only one macro-cycle of refinement.",
        program_name="phenix.refine",
    )
    # → ["main.number_of_macro_cycles=1"]

Design
------
Each strategy_flag entry in programs.yaml may contain a ``natural_language``
list of matchers:

    cycles:
      flag: "main.number_of_macro_cycles={value}"
      type: int
      natural_language:
        - pattern: "(?:run\\s+only\\s+)?({count})\\s+macro[- ]cycles?"
          value_group: 1      # regex group that captures the count
          word_to_int: true   # convert "one"→1, "two"→2, …

Special token ``{count}`` in a pattern is expanded to a combined
word-number alternation: ``(?:one|two|…|\\d+)``.

If ``value`` is specified instead of ``value_group``, that literal is
used directly (useful for boolean flags like ``simulated_annealing=True``).

Deduplication
-------------
When multiple strategy_flags map to the same flag template (e.g. both
``cycles`` and ``macro_cycles`` expand to ``main.number_of_macro_cycles``),
only the first match is returned to avoid duplicate injection.
"""

from __future__ import absolute_import, division, print_function

import re

# ---------------------------------------------------------------------------
# Word-to-integer lookup used when word_to_int: true
# ---------------------------------------------------------------------------
_WORD_TO_INT = {
    'one': 1, 'two': 2, 'three': 3, 'four': 4, 'five': 5,
    'six': 6, 'seven': 7, 'eight': 8, 'nine': 9, 'ten': 10,
    'eleven': 11, 'twelve': 12,
}

# Combined alternation for the {count} placeholder
_COUNT_ALT = '(?:' + '|'.join(list(_WORD_TO_INT.keys()) + [r'\d+']) + ')'


def _expand_pattern(pattern):
    """Replace the ``{count}`` placeholder with the word+digit alternation."""
    return pattern.replace('{count}', _COUNT_ALT)


def _coerce_value(raw, flag_type, word_to_int=False):
    """
    Coerce a matched string to the correct PHIL type.

    Args:
        raw:         The raw matched string (e.g. "one", "5", "True").
        flag_type:   "int", "float", "boolean", "string", or None.
        word_to_int: If True, convert English words to integers first.

    Returns:
        str representation ready for PHIL (e.g. "1", "3.5", "True", "my_str"),
        or None if coercion fails.
    """
    if word_to_int:
        lower = raw.lower().strip()
        if lower in _WORD_TO_INT:
            raw = str(_WORD_TO_INT[lower])
        elif raw.isdigit():
            raw = raw  # already numeric
        else:
            return None  # unrecognised word

    t = (flag_type or '').lower()
    if t == 'int':
        try:
            return str(int(raw))
        except (ValueError, TypeError):
            return None
    elif t == 'float':
        try:
            return str(float(raw))
        except (ValueError, TypeError):
            return None
    elif t == 'boolean':
        if raw.lower() in ('true', '1', 'yes'):
            return 'True'
        if raw.lower() in ('false', '0', 'no'):
            return 'False'
        return None
    else:
        # string — return as-is, no quoting needed for simple identifiers
        return raw


def _load_strategy_flags(program_name):
    """
    Return the strategy_flags dict for *program_name* from programs.yaml,
    or an empty dict if not found.
    """
    try:
        try:
            from libtbx.langchain.knowledge.yaml_loader import get_program
        except ImportError:
            from knowledge.yaml_loader import get_program
        prog = get_program(program_name)
        if prog and isinstance(prog, dict):
            return prog.get('strategy_flags') or {}
    except Exception:
        pass
    return {}


def extract_nl_params(guidelines, program_name):
    """
    Scan *guidelines* for natural-language phrases that match entries in the
    program's ``strategy_flags.natural_language`` rules and return a list of
    valid PHIL ``key=value`` strings.

    Args:
        guidelines:    User advice / guidelines string.
        program_name:  Full program name (e.g. "phenix.refine").

    Returns:
        list[str]: PHIL params (e.g. ["main.number_of_macro_cycles=1"]).
                   Empty list if nothing matched.
    """
    strategy_flags = _load_strategy_flags(program_name)
    if not strategy_flags:
        return []

    results = []
    seen_flag_templates = set()  # avoid duplicate injection for aliased flags

    for flag_name, flag_def in strategy_flags.items():
        if not isinstance(flag_def, dict):
            continue
        nl_rules = flag_def.get('natural_language')
        if not nl_rules:
            continue

        flag_template = flag_def.get('flag', '')
        flag_type = flag_def.get('type', 'string')

        # Skip if we already emitted a param for this flag template
        if flag_template in seen_flag_templates:
            continue

        for rule in nl_rules:
            if not isinstance(rule, dict):
                continue

            # Expand {count} placeholder in pattern
            raw_pattern = rule.get('pattern', '')
            if not raw_pattern:
                continue
            expanded = _expand_pattern(raw_pattern)

            m = re.search(expanded, guidelines, re.IGNORECASE)
            if not m:
                continue

            # Determine value
            if 'value' in rule:
                # Literal value specified in YAML
                phil_val = rule['value']
            elif 'value_group' in rule:
                grp = rule['value_group']
                try:
                    raw_val = m.group(grp)
                except IndexError:
                    continue
                phil_val = _coerce_value(
                    raw_val,
                    flag_type,
                    word_to_int=bool(rule.get('word_to_int', False)),
                )
                if phil_val is None:
                    continue
            else:
                continue

            # Render the PHIL param from the flag template
            phil_param = flag_template.replace('{value}', phil_val)

            results.append(phil_param)
            seen_flag_templates.add(flag_template)
            break  # first matching rule wins for this flag

    return results
