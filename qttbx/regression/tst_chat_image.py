"""Image-handling tests: the LRU image cache + storage round-trip
helpers (image_cache.py) and the click-to-open modal viewer
(image_lightbox.py). Both are about turning a stored attachment into
something visible on screen."""

import os
import shutil
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from libtbx.utils import format_cpu_times, null_out

try:
  from qttbx.qt import QtCore, QtGui, QtWidgets
except ImportError:
  print("PySide2/PySide6 not available; skipping")
  print("OK")
  sys.exit(0)


from qttbx.widgets.chat.agent.storage import ConversationStorage
from qttbx.widgets.chat.image_cache import (
  ImageCache, get_image, set_default_cache)


def _qapp():
  from qttbx.widgets.font_init import init_default_app_font
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  return app


def _png_bytes(width=10, height=10, color=(255, 0, 0)):
  _qapp()
  img = QtGui.QImage(width, height, QtGui.QImage.Format_ARGB32)
  img.fill(QtGui.QColor(*color))
  buf = QtCore.QBuffer()
  buf.open(QtCore.QBuffer.WriteOnly)
  img.save(buf, "PNG")
  return bytes(buf.data())


def _solid_pixmap(w, h):
  pm = QtGui.QPixmap(w, h)
  pm.fill(QtGui.QColor(QtCore.Qt.blue))
  return pm


# ---- image cache ---------------------------------------------------------


def exercise_get_image_round_trips():
  tmp = tempfile.mkdtemp()
  try:
    storage = ConversationStorage(project_dir=Path(tmp), log=null_out())
    set_default_cache(ImageCache(capacity=8))
    att = storage.store_attachment("c1", _png_bytes(), "image/png")
    img = get_image(storage, "c1", att.sha256)
    assert isinstance(img, QtGui.QImage)
    assert not img.isNull()
    assert img.width() == 10 and img.height() == 10
  finally:
    shutil.rmtree(tmp)


def exercise_missing_attachment_returns_placeholder():
  tmp = tempfile.mkdtemp()
  try:
    storage = ConversationStorage(project_dir=Path(tmp), log=null_out())
    set_default_cache(ImageCache(capacity=8))
    img = get_image(storage, "no-such-conv", "0" * 64)
    assert isinstance(img, QtGui.QImage)
    assert not img.isNull()
  finally:
    shutil.rmtree(tmp)


def exercise_lru_evicts_oldest_at_capacity():
  tmp = tempfile.mkdtemp()
  try:
    storage = ConversationStorage(project_dir=Path(tmp), log=null_out())
    cache = ImageCache(capacity=3)
    set_default_cache(cache)
    shas = []
    for i in range(4):
      att = storage.store_attachment(
        "c1", _png_bytes(color=(i*40, 0, 0)), "image/png")
      shas.append(att.sha256)
      get_image(storage, "c1", att.sha256)
    # The OLDEST entry (the victim) must be the one evicted, and the 3
    # most-recently-used must remain -- asserting only the size would pass even
    # if the wrong entry were dropped. The cache key is (conv_id, sha256);
    # cache.get is a peek (None when absent).
    assert len(cache) == 3, len(cache)
    assert cache.get(("c1", shas[0])) is None, "oldest entry should be evicted"
    for sha in shas[1:]:
      assert cache.get(("c1", sha)) is not None, sha
  finally:
    shutil.rmtree(tmp)


# ---- image lightbox ------------------------------------------------------


def exercise_lightbox_holds_pixmap_and_closes_on_click():
  _qapp()
  from qttbx.widgets.chat.image_lightbox import ImageLightbox
  pm = _solid_pixmap(400, 300)
  dlg = ImageLightbox(pixmap=pm)
  assert dlg.pixmap is pm, "pixmap reference must be retained"
  # Simulate click anywhere -- the dialog should accept (close).
  dlg.show()
  dlg._on_clicked()
  assert dlg.isHidden()


def exercise_lightbox_is_modal():
  _qapp()
  from qttbx.widgets.chat.image_lightbox import ImageLightbox
  pm = _solid_pixmap(10, 10)
  dlg = ImageLightbox(pixmap=pm)
  assert dlg.isModal() or dlg.windowModality() != QtCore.Qt.NonModal


def exercise():
  exercise_get_image_round_trips()
  exercise_missing_attachment_returns_placeholder()
  exercise_lru_evicts_oldest_at_capacity()
  exercise_lightbox_holds_pixmap_and_closes_on_click()
  exercise_lightbox_is_modal()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
