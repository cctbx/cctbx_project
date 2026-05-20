"""LRU image cache for the chat UI (Section 12.5 of the design spec).

Keyed by (conv_id, sha256) for full-size, (conv_id, sha256, width) for
thumbnails. Capacity is in entry count, not bytes — the design caps at
~200 full-size QImages which is a reasonable proxy for the ~200MB target
at typical chat image sizes. Disk-backed thumbnail cache is a follow-up.

`get_image` / `get_thumbnail` are module-level helpers backed by a
single default ImageCache instance set via `set_default_cache` (this lets
multiple widgets share a single cache without each holding a reference)."""

from collections import OrderedDict

from qttbx.qt import QtCore, QtGui


_DEFAULT_CAPACITY = 200


class ImageCache:

  def __init__(self, capacity=_DEFAULT_CAPACITY):
    self._capacity = capacity
    self._entries = OrderedDict()        # key -> QImage

  def __len__(self):
    return len(self._entries)

  def get(self, key):
    if key not in self._entries:
      return None
    self._entries.move_to_end(key)       # mark most-recently-used
    return self._entries[key]

  def put(self, key, image):
    if key in self._entries:
      self._entries.move_to_end(key)
      self._entries[key] = image
      return
    self._entries[key] = image
    while len(self._entries) > self._capacity:
      self._entries.popitem(last=False)


_default_cache = ImageCache()


def set_default_cache(cache):
  """Override the module-level cache. Tests call this; production code
  should leave the default in place."""
  global _default_cache
  _default_cache = cache


def get_image(storage, conv_id, sha256):
  """Load + decode a full-size image. Returns a placeholder gray QImage
  on missing/corrupt attachments — never raises."""
  key = ("img", conv_id, sha256)
  cached = _default_cache.get(key)
  if cached is not None:
    return cached
  img = QtGui.QImage()
  try:
    data = storage.load_attachment(conv_id, sha256)
  except Exception:
    # Missing or unreadable; surface a placeholder rather than crashing the
    # widget — Section 12 design favors graceful degradation.
    return _placeholder()
  if not img.loadFromData(data):
    return _placeholder()
  _default_cache.put(key, img)
  return img


def get_thumbnail(storage, conv_id, sha256, width=240):
  """Width-scaled thumbnail. Aspect preserved. Cached separately from
  full-size."""
  key = ("thumb", conv_id, sha256, int(width))
  cached = _default_cache.get(key)
  if cached is not None:
    return cached
  full = get_image(storage, conv_id, sha256)
  if full.width() <= 0:
    return _placeholder()
  thumb = full.scaledToWidth(
    int(width), QtCore.Qt.SmoothTransformation)
  _default_cache.put(key, thumb)
  return thumb


def _placeholder(width=240, height=180):
  img = QtGui.QImage(width, height, QtGui.QImage.Format_ARGB32)
  img.fill(QtGui.QColor(220, 220, 220))
  return img
