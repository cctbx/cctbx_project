"""Where the test corpora are, without anyone having to export anything.

RESOLUTION ORDER, first hit wins:

  1. the environment variable, if set        -- an override, not a requirement
  2. `corpus/` beside this package           -- the shipped default
  3. a few conventional siblings             -- if someone moved it one level up

Nothing here is needed by `scan()`. This exists so that
`phenix.python tests/run_all_tests.py` runs the corpus-level invariants out
of the box. Those invariants are most of the value -- fixture-only tests
encode the same assumptions as the code -- and a suite that needs two exports
before it checks anything real is a suite that will be run without them.
"""

from __future__ import absolute_import, division, print_function

import os

HERE = os.path.dirname(os.path.abspath(__file__))

# corpus/
#   work/{ok,err}      253 production logs, GUI shape        (the working set)
#   agent/             the frozen wave-1 twenty, agent shape
#   agent_gui/         the same twenty, stripped
DEFAULT_ROOT = os.path.join(HERE, "corpus")

_CANDIDATE_ROOTS = (
  DEFAULT_ROOT,
  os.path.join(os.path.dirname(HERE), "corpus"),          # langchain/corpus
  os.path.join(os.path.dirname(os.path.dirname(HERE)), "log_corpus"),
)


def _first_existing(relative, variable):
  """-> an absolute path, or None. The variable wins if it is set."""
  override = os.environ.get(variable)
  if override:
    return override if os.path.isdir(override) else None
  for root in _CANDIDATE_ROOTS:
    candidate = os.path.join(root, relative)
    if os.path.isdir(candidate):
      return candidate
  return None


def working_corpus():
  """253 GUI-shape logs: 230 that ran, 23 that failed."""
  return _first_existing("work", "PHENIX_LOG_CORPUS")


def agent_corpus():
  """The frozen wave-1 twenty, agent-shape, with program_truth.json."""
  return _first_existing("agent", "PHENIX_LOG_CORPUS_AGENT")


def agent_gui_corpus():
  """The same twenty with the agent wrapper stripped."""
  return _first_existing("agent_gui", "PHENIX_LOG_CORPUS_GUI")


def describe():
  """One line per corpus, for a runner to print when something is missing."""
  lines = []
  for name, path in (("working", working_corpus()),
                     ("agent-shape", agent_corpus()),
                     ("agent GUI-shape", agent_gui_corpus())):
    lines.append("  %-16s %s" % (name, path or "NOT FOUND"))
  return "\n".join(lines)


def missing_corpus_message():
  """What to tell someone whose corpus is not where it should be."""
  return (
    "The corpus was not found. It ships in this package at\n"
    "    %s\n"
    "If it is somewhere else, either move it there or set\n"
    "PHENIX_LOG_CORPUS / PHENIX_LOG_CORPUS_AGENT.\n\n"
    "Looked in:\n%s" % (DEFAULT_ROOT, describe()))
