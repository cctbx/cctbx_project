"""
Expert Knowledge Base — Runtime Loader

Usage pattern:
  Static YAML loaded once at startup into a Python dict.
  At each decision point, filter by category + tags to get
  the 3-5 most relevant entries, then inject into the
  thinking agent's prompt alongside the validation report.

  No vector DB, no hashing, no build step.
  Edit the YAML, restart, done.
"""

import yaml
from collections import defaultdict


class ExpertKnowledgeBase:
  """Simple tag-based lookup for the expert knowledge base."""

  def __init__(self, yaml_path):
    with open(yaml_path) as f:
      self.entries = yaml.safe_load(f)
    # Index by category
    self.by_category = defaultdict(list)
    # Index by tag → entries (for IDF weighting)
    self.by_tag = defaultdict(list)
    for entry in self.entries:
      cat = entry.get('category', '')
      self.by_category[cat].append(entry)
      for tag in entry.get('tags', []):
        self.by_tag[tag].append(entry)
    # Precompute IDF weights: tags on fewer entries
    # get higher weight. A tag on 1 entry has weight
    # 1.0; a tag on all entries has weight ~0.
    n = max(len(self.entries), 1)
    self._tag_idf = {}
    for tag, elist in self.by_tag.items():
      # IDF = log(N / count), but we use a simpler
      # inverse-frequency: N / count, then normalize
      # so the max weight is 1.0.
      self._tag_idf[tag] = n / max(len(elist), 1)
    # Normalize to [0, 1]
    max_idf = max(self._tag_idf.values()) if (
      self._tag_idf
    ) else 1.0
    for tag in self._tag_idf:
      self._tag_idf[tag] /= max_idf

  def query(self, category=None, tags=None,
            max_results=5):
    """Find relevant entries by category and/or tags.

    Scores entries using IDF-weighted tag matching:
    tags that appear on fewer entries contribute more
    to the relevance score. This prevents entries with
    only broad tags (like 'refinement') from outranking
    entries that match on specific diagnostic tags
    (like 'plateau' or 'twinning').

    Args:
      category: str or None. Filter to this category.
      tags: list of str or None. Entries matching more
            tags rank higher.
      max_results: int. Maximum entries to return.

    Returns:
      List of entries sorted by relevance score.
    """
    if category:
      candidates = self.by_category.get(
        category, []
      )
    else:
      candidates = self.entries

    if not tags:
      return candidates[:max_results]

    # Score by IDF-weighted tag overlap
    scored = []
    tag_set = set(tags)
    for entry in candidates:
      entry_tags = set(entry.get('tags', []))
      overlap = tag_set & entry_tags
      if not overlap:
        continue
      score = sum(
        self._tag_idf.get(t, 0.5)
        for t in overlap
      )
      scored.append((score, entry))

    scored.sort(key=lambda x: -x[0])
    return [
      entry for _, entry in scored[:max_results]
    ]

  def format_for_prompt(self, entries, max_chars=2000):
    """Format entries as compact text for the thinking agent prompt.

    Stays within a character budget so we don't blow up the
    context window.
    """
    lines = []
    chars = 0
    for entry in entries:
      block = self._format_entry(entry)
      if chars + len(block) > max_chars:
        break
      lines.append(block)
      chars += len(block)
    return '\n'.join(lines)

  def _format_entry(self, entry):
    """Format a single entry compactly.

    Handles both v1 schema (assessment/action) and
    v2 schema (llm_error/correct).
    """
    eid = entry.get('id', '?')
    question = entry.get('question', '').strip()

    # Detect schema version
    if 'llm_error' in entry:
      # v2: focused on what LLMs get wrong
      correct = str(
        entry.get('correct', '')
      ).strip()
      llm_err = str(
        entry.get('llm_error', '')
      ).strip()
      return (
        f"[{eid}] {question}\n"
        f"LLM pitfall: {llm_err}\n"
        f"Correct: {correct}\n"
      )

    # v1: assessment/action format
    assessment = entry.get('assessment', '')
    if isinstance(assessment, dict):
      parts = []
      for k, v in assessment.items():
        parts.append(
          f"  {k}: {str(v).strip()}"
        )
      assessment_text = '\n'.join(parts)
    else:
      assessment_text = str(assessment).strip()
    action = entry.get('action', '')
    if isinstance(action, dict):
      parts = []
      for k, v in action.items():
        parts.append(
          f"  {k}: {str(v).strip()}"
        )
      action_text = '\n'.join(parts)
    else:
      action_text = str(action).strip()
    return (
      f"[{eid}] {question}\n"
      f"Assessment: {assessment_text}\n"
      f"Action: {action_text}\n"
    )

  @property
  def stats(self):
    """Summary statistics."""
    cats = defaultdict(int)
    confs = defaultdict(int)
    for e in self.entries:
      cats[e.get('category', '?')] += 1
      confs[e.get('confidence', '?')] += 1
    return {
      'total': len(self.entries),
      'by_category': dict(cats),
      'by_confidence': dict(confs),
    }


# ---------------------------------------------------------------
# Example: how the agent uses it at a decision point
# ---------------------------------------------------------------

def example_usage():
  kb = ExpertKnowledgeBase(
    'expert_knowledge_base_v2.yaml'
  )
  print(f"Loaded {kb.stats['total']} entries")
  print(
    f"Categories: {kb.stats['by_category']}"
  )
  print()

  # Scenario: refinement stuck, R-free high
  context_tags = [
    'r_free', 'r_free_stuck', 'plateau',
    'refinement', 'convergence',
  ]
  entries = kb.query(
    category='stopping',
    tags=context_tags,
    max_results=3,
  )
  prompt_text = kb.format_for_prompt(
    entries, max_chars=2000
  )
  print("=== Stopping decision ===")
  print(prompt_text)
  print()

  # Scenario: considering adding waters
  entries = kb.query(
    category=None,
    tags=[
      'waters', 'ordered_solvent',
      'r_free', 'refinement',
    ],
    max_results=3,
  )
  prompt_text = kb.format_for_prompt(
    entries, max_chars=2000
  )
  print("=== Water addition decision ===")
  print(prompt_text)


if __name__ == '__main__':
  example_usage()
