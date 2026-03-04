"""
Strategy Memory for the Thinking Agent (v113).

Accumulates scientific understanding across cycles.
Pure Python — no LLM or external dependencies.
Persisted via session_info round-trip (JSON-safe dict).
"""

from __future__ import absolute_import, division, print_function


class StrategyMemory(object):
  """Persistent memory of crystallographic strategy across cycles.

  Tracks data quality observations, phasing strategy decisions,
  R-free history, and expert concerns. Updated each cycle by the
  think node and persisted via session_info.
  """

  def __init__(self):
    self.data_quality = ""       # e.g. "good", "twinned", "anisotropic"
    self.phasing_strategy = ""   # e.g. "MR", "SAD", "direct_methods"
    self.concerns = []           # List of strings, e.g. ["weak anomalous"]
    self.decisions = []          # List of (cycle, decision_text) tuples
    self.r_free_history = []     # R-free values in order
    self.programs_run = []       # Program names in order

  def to_dict(self):
    """Serialize to JSON-safe dict."""
    return {
      "data_quality": self.data_quality,
      "phasing_strategy": self.phasing_strategy,
      "concerns": list(self.concerns),
      "decisions": list(self.decisions),
      "r_free_history": list(self.r_free_history),
      "programs_run": list(self.programs_run),
    }

  @classmethod
  def from_dict(cls, d):
    """Reconstruct from dict. Tolerant of missing/extra keys."""
    mem = cls()
    if not isinstance(d, dict):
      return mem
    mem.data_quality = d.get("data_quality", "")
    mem.phasing_strategy = d.get("phasing_strategy", "")
    mem.concerns = list(d.get("concerns", []))
    mem.decisions = list(d.get("decisions", []))
    mem.r_free_history = list(d.get("r_free_history", []))
    mem.programs_run = list(d.get("programs_run", []))
    return mem

  def update(self, assessment):
    """Update memory from a think-node assessment dict.

    Args:
      assessment: Dict with keys like 'analysis', 'guidance',
        'data_quality', 'phasing_strategy', 'concerns'.
        All keys optional.
    """
    if not isinstance(assessment, dict):
      return
    if assessment.get("data_quality"):
      self.data_quality = assessment["data_quality"]
    if assessment.get("phasing_strategy"):
      self.phasing_strategy = assessment["phasing_strategy"]
    new_concerns = assessment.get("concerns", [])
    if isinstance(new_concerns, list):
      for c in new_concerns:
        if c and c not in self.concerns:
          self.concerns.append(c)
    # Cap concerns to avoid unbounded growth
    if len(self.concerns) > 10:
      self.concerns = self.concerns[-10:]

  def record_outcome(self, program, result, metrics=None):
    """Record the outcome of a cycle.

    Args:
      program: Program name, e.g. "phenix.refine"
      result: "SUCCESS" or "FAILED"
      metrics: Optional dict with extracted metrics
    """
    self.programs_run.append(program)
    # Cap to last 30 programs
    if len(self.programs_run) > 30:
      self.programs_run = self.programs_run[-30:]

    if metrics and isinstance(metrics, dict):
      r_free = metrics.get("r_free")
      if r_free is not None:
        try:
          r_free = float(r_free)
          self.r_free_history.append(r_free)
          # Cap to last 20 values
          if len(self.r_free_history) > 20:
            self.r_free_history = self.r_free_history[-20:]
        except (ValueError, TypeError):
          pass

  def record_decision(self, cycle, text):
    """Record an expert decision.

    Args:
      cycle: Cycle number
      text: Short decision text
    """
    self.decisions.append([cycle, text[:200]])
    # Cap to last 10 decisions
    if len(self.decisions) > 10:
      self.decisions = self.decisions[-10:]

  def metrics_stalled(self):
    """Check if R-free has stalled (no improvement for 3+ values).

    Returns:
      True if the last 3 R-free values show no improvement
      (all within 0.001 of the first of the three).
    """
    rfree = self.r_free_history
    if len(rfree) < 3:
      return False
    recent = rfree[-3:]
    # Stalled if none improved by more than 0.001
    return all(r >= recent[0] - 0.001 for r in recent[1:])

  def summary(self):
    """One-line summary for debug logging."""
    parts = []
    if self.data_quality:
      parts.append("quality=%s" % self.data_quality)
    if self.phasing_strategy:
      parts.append("strategy=%s" % self.phasing_strategy)
    if self.r_free_history:
      parts.append("r_free=%.3f" % self.r_free_history[-1])
    if self.concerns:
      parts.append("concerns=%d" % len(self.concerns))
    return ", ".join(parts) if parts else "(empty)"
