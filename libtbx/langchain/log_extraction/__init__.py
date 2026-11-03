"""Deterministic reader for Phenix log files.

    from log_extraction import scan
    structure = scan(open(path, "rb").read().decode("utf-8", "replace"))

Text in, structure out, every item carrying the line it came from. It reads
what the program says about itself; it does not interpret, rank, diagnose or
call a model.

See docs/EXTRACTOR_ARCHITECTURE.md for the whole picture, docs/DESIGN.md to
extend it, and docs/EXTRACTOR_REQUIREMENTS_v2.md for the spec.
"""

from __future__ import absolute_import, division, print_function

from log_extraction.log_structure_extractor import (      # noqa: F401
  scan,
  LogStructure,
  Section, Phase, Stage, Cycle, Decision, ControlSkip, Exclusion,
  Measurement, LabeledValue, CompletionRecord, Unparsed,
  Candidate, Identification, Region,
  SOURCE_PROGRAM, SOURCE_AGENT, SOURCE_UNKNOWN,
  UNCLAIMED, GENERIC_ONLY, RULE_EXCLUDED,
  ITEM_KINDS, DIAGNOSTIC_KINDS, BASIC_KINDS,
  split_lines, strip_agent_prefix, screen_line,
)

__version__ = "0.1.0"
