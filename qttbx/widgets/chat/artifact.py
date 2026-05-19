"""``Artifact`` dataclass and renderer registry.

Renderers are callables ``(artifact, storage) -> QWidget``; the
``ArtifactPanel`` looks them up by ``artifact.kind``. The default image
renderer is registered lazily by ``ArtifactPanel`` so importing this
module stays Qt-free.
"""

from dataclasses import dataclass, field


@dataclass
class Artifact:
  kind: str                          # 'image' is the only kind currently rendered
  payload: dict = field(default_factory=dict)
  caption: str = ""
  source: str = ""                   # e.g., 'user', 'assistant', 'mcp:phenix'


_renderers = {}


def register_renderer(kind, factory):
  """Register a renderer ``factory(artifact, storage) -> QWidget`` for
  ``kind``."""
  _renderers[kind] = factory


def renderer_for(kind):
  return _renderers.get(kind)
