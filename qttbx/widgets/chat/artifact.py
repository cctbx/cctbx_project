"""``Artifact`` dataclass and renderer registry.

Renderers are callables ``(artifact, storage) -> QWidget``; the
``ArtifactPanel`` looks them up by ``artifact.kind``. The default image
renderer is registered lazily by ``ArtifactPanel`` so importing this
module stays Qt-free.
"""

from dataclasses import dataclass, field


@dataclass
class Artifact:
  """A renderable artifact looked up by ``kind`` in the renderer registry.

  Parameters
  ----------
  kind : str
      Renderer key; ``'image'`` is the only kind currently rendered.
  payload : dict, optional
      Renderer-specific data. Defaults to an empty dict.
  caption : str, optional
      Human-readable caption.
  source : str, optional
      Origin of the artifact, e.g. ``'user'``, ``'assistant'``,
      ``'mcp:phenix'``.
  """
  kind: str                          # 'image' is the only kind currently rendered
  payload: dict = field(default_factory=dict)
  caption: str = ""
  source: str = ""                   # e.g., 'user', 'assistant', 'mcp:phenix'


_renderers = {}


def register_renderer(kind, factory):
  """Register a renderer factory for an artifact kind.

  Parameters
  ----------
  kind : str
      Artifact kind the factory renders.
  factory : callable
      Renderer of signature ``factory(artifact, storage) -> QWidget``.
  """
  _renderers[kind] = factory


def renderer_for(kind):
  """Return the renderer factory registered for ``kind``, or ``None``.

  Parameters
  ----------
  kind : str
      Artifact kind to look up.

  Returns
  -------
  callable or None
      The registered renderer factory, or ``None`` if none is
      registered for ``kind``.
  """
  return _renderers.get(kind)
