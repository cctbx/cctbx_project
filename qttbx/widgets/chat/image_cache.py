"""LRU image cache for the chat UI.

Keyed by ``(conv_id, sha256)`` for full-size images. Capacity is in
entry count, not bytes — ~200 full-size ``QImage`` instances is a
reasonable proxy for ~200 MB at typical chat image sizes.

``get_image`` is a module-level helper backed by a single default
``ImageCache`` instance set via ``set_default_cache``, which lets
multiple widgets share a single cache without each holding a
reference."""

from collections import OrderedDict

from qttbx.qt import QtGui


_DEFAULT_CAPACITY = 200


class ImageCache:
  """LRU cache of ``QImage`` instances keyed by an opaque tuple.

  Parameters
  ----------
  capacity : int, optional
      Maximum number of entries to retain. The least-recently-used
      entry is evicted once this count is exceeded.
  """

  def __init__(self, capacity=_DEFAULT_CAPACITY):
    self._capacity = capacity
    self._entries = OrderedDict()        # key -> QImage

  def __len__(self):
    return len(self._entries)

  def get(self, key):
    """Return the cached image for ``key`` (marking it MRU), or ``None``."""
    if key not in self._entries:
      return None
    self._entries.move_to_end(key)       # mark most-recently-used
    return self._entries[key]

  def put(self, key, image):
    """Insert or refresh ``key`` as MRU, evicting the LRU over capacity."""
    if key in self._entries:
      self._entries.move_to_end(key)
      self._entries[key] = image
      return
    self._entries[key] = image
    while len(self._entries) > self._capacity:
      self._entries.popitem(last=False)


_default_cache = ImageCache()


def set_default_cache(cache):
  """Override the module-level cache.

  Tests call this; production code should leave the default in place.

  Parameters
  ----------
  cache : ImageCache
      The cache instance to install as the module-level default.
  """
  global _default_cache
  _default_cache = cache


def get_image(storage, conv_id, sha256):
  """Load and decode a full-size image, using the default cache.

  Never raises: a placeholder gray ``QImage`` is returned on missing or
  corrupt attachments.

  Parameters
  ----------
  storage : object
      Attachment store providing ``load_attachment(conv_id, sha256)``.
  conv_id : str
      Conversation the attachment belongs to.
  sha256 : str
      Content hash identifying the attachment.

  Returns
  -------
  QtGui.QImage
      The decoded image, or a gray placeholder on failure.
  """
  key = (conv_id, sha256)
  cached = _default_cache.get(key)
  if cached is not None:
    return cached
  img = QtGui.QImage()
  try:
    data = storage.load_attachment(conv_id, sha256)
  except Exception:
    # Missing or unreadable; surface a placeholder rather than crashing
    # the widget.
    return _placeholder()
  if not img.loadFromData(data):
    return _placeholder()
  _default_cache.put(key, img)
  return img


def _placeholder(width=240, height=180):
  img = QtGui.QImage(width, height, QtGui.QImage.Format_ARGB32)
  img.fill(QtGui.QColor(220, 220, 220))
  return img
