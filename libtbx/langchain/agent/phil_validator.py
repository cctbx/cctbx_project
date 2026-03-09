"""
PHIL Strategy Validator (Fix 4, v115).

Validates LLM-generated strategy parameters against a program's
known strategy_flags from programs.yaml.  Strips unrecognized
parameters before they reach the command builder, preventing
'unrecognized PHIL parameter' errors at runtime.

Usage:
    from agent.phil_validator import validate_phil_strategy

    cleaned, stripped = validate_phil_strategy(
        "phenix.autobuild",
        {"resolution": 2.5, "obs_labels": "I(+)"}
    )
    # cleaned = {"resolution": 2.5}
    # stripped = [("obs_labels", "I(+)")]
"""

from __future__ import absolute_import, division, print_function
import logging

logger = logging.getLogger(__name__)

# Keys that appear in LLM strategy dicts but are NOT PHIL
# parameters — they are consumed by the build pipeline or
# represent file codes.  These must never be stripped.
_BUILD_PIPELINE_KEYS = frozenset({
  "ligand",           # 3-letter ligand code (file, not PHIL)
  "output_prefix",    # consumed by CommandBuilder
  "output_name",      # consumed by CommandBuilder
  "output_root",      # consumed by CommandBuilder
})

# Cache for loaded strategy_flags to avoid repeated YAML reads
_STRATEGY_FLAGS_CACHE = {}


def _load_strategy_flags(program):
  """Load strategy_flags for a program from programs.yaml.

  Returns:
      set of allowed flag names, or None if program not found
      or YAML loader unavailable.
  """
  if program in _STRATEGY_FLAGS_CACHE:
    return _STRATEGY_FLAGS_CACHE[program]

  # Try the libtbx yaml_loader first (production path)
  try:
    try:
      from libtbx.langchain.knowledge.yaml_loader \
        import get_program
    except ImportError:
      from knowledge.yaml_loader import get_program
    pdef = get_program(program)
    if pdef:
      flags = pdef.get("strategy_flags", {})
      allowed = set(flags.keys()) if flags else set()
      _STRATEGY_FLAGS_CACHE[program] = allowed
      return allowed
  except Exception:
    pass

  # Fallback: load programs.yaml directly (standalone/test)
  try:
    import os
    import yaml
    _dir = os.path.dirname(os.path.abspath(__file__))
    _yaml_path = os.path.join(
      _dir, "..", "knowledge", "programs.yaml")
    _yaml_path = os.path.normpath(_yaml_path)
    if os.path.isfile(_yaml_path):
      with open(_yaml_path) as f:
        data = yaml.safe_load(f)
      pdef = data.get(program)
      if pdef and isinstance(pdef, dict):
        flags = pdef.get("strategy_flags", {})
        allowed = set(flags.keys()) if flags else set()
        _STRATEGY_FLAGS_CACHE[program] = allowed
        return allowed
  except Exception:
    pass

  _STRATEGY_FLAGS_CACHE[program] = None
  return None


def validate_phil_strategy(program, strategy):
  """Validate LLM strategy dict against program's strategy_flags.

  Strips keys that are not recognized by the target program.
  This prevents unrecognized-PHIL-parameter errors at runtime.

  Args:
      program: Program name (e.g., "phenix.autobuild")
      strategy: Dict of PHIL params from LLM intent

  Returns:
      (cleaned_strategy, stripped_list) tuple.
      cleaned_strategy: dict with unrecognized keys removed.
      stripped_list: list of (key, value) tuples that were
          removed.  Empty list if nothing stripped.

  Never raises — returns strategy unchanged if validation
  cannot be performed (YAML not available, unknown program,
  empty strategy, or no strategy_flags defined).
  """
  if not strategy or not program:
    return strategy, []

  allowed = _load_strategy_flags(program)
  if allowed is None:
    # Program not found or YAML loader unavailable
    return strategy, []

  if not allowed:
    # Program has no strategy_flags defined — can't validate.
    # Pass through to avoid blocking programs that just don't
    # have flags declared yet.
    return strategy, []

  # Also allow build-pipeline keys
  allowed = allowed | _BUILD_PIPELINE_KEYS

  cleaned = {}
  stripped = []
  for key, value in strategy.items():
    if key in allowed:
      cleaned[key] = value
    else:
      stripped.append((key, value))

  if stripped:
    names = ", ".join(
      "%s=%s" % (k, v) for k, v in stripped)
    logger.info(
      "PHIL VALIDATION: stripped %d unrecognized "
      "param(s) for %s: %s",
      len(stripped), program, names)

  return cleaned, stripped
