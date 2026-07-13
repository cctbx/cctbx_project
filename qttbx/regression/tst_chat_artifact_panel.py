"""ArtifactPanel tests with the real (Plan 5) implementation."""

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
from qttbx.widgets.font_init import init_default_app_font
from qttbx.widgets.chat.artifact import Artifact
from qttbx.widgets.chat.artifact_panel import ArtifactPanel


def _png_bytes(color=(0, 0, 255)):
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  img = QtGui.QImage(40, 30, QtGui.QImage.Format_ARGB32)
  img.fill(QtGui.QColor(*color))
  buf = QtCore.QBuffer()
  buf.open(QtCore.QBuffer.WriteOnly)
  img.save(buf, "PNG")
  return bytes(buf.data())


def exercise_empty_state_has_no_current_artifact():
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  p = ArtifactPanel()
  assert p.current_artifact() is None


def exercise_add_image_artifact_auto_advances():
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  tmp = tempfile.mkdtemp()
  try:
    storage = ConversationStorage(project_dir=Path(tmp), log=null_out())
    att = storage.store_attachment("c", _png_bytes(), "image/png")
    p = ArtifactPanel(storage=storage)
    p.add_artifact(Artifact(
      kind="image", payload={"conv_id": "c", "sha256": att.sha256,
                             "mime": "image/png"}, caption="x.png"))
    assert p.current_artifact() is not None
    # The caption is the contract. The old `or "image" in ..." branch matched
    # the kind fallback, so it still passed if the caption was dropped entirely.
    assert "x.png" in p.status_text(), p.status_text()
  finally:
    shutil.rmtree(tmp)


def exercise_user_navigation_pins_index():
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  tmp = tempfile.mkdtemp()
  try:
    storage = ConversationStorage(project_dir=Path(tmp), log=null_out())
    a1 = storage.store_attachment("c", _png_bytes(color=(255, 0, 0)),
                                  "image/png")
    a2 = storage.store_attachment("c", _png_bytes(color=(0, 255, 0)),
                                  "image/png")
    p = ArtifactPanel(storage=storage)
    p.add_artifact(Artifact(kind="image", payload={
      "conv_id": "c", "sha256": a1.sha256, "mime": "image/png"}))
    p.add_artifact(Artifact(kind="image", payload={
      "conv_id": "c", "sha256": a2.sha256, "mime": "image/png"}))
    assert p.current_artifact().payload["sha256"] == a2.sha256
    p.go_prev()
    assert p.current_artifact().payload["sha256"] == a1.sha256
    # New arrivals while navigating away should NOT auto-advance.
    a3 = storage.store_attachment("c", _png_bytes(color=(0, 0, 255)),
                                  "image/png")
    p.add_artifact(Artifact(kind="image", payload={
      "conv_id": "c", "sha256": a3.sha256, "mime": "image/png"}))
    assert p.current_artifact().payload["sha256"] == a1.sha256
  finally:
    shutil.rmtree(tmp)


def exercise_show_image_focuses_existing_or_inserts():
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  tmp = tempfile.mkdtemp()
  try:
    storage = ConversationStorage(project_dir=Path(tmp), log=null_out())
    a = storage.store_attachment("c", _png_bytes(), "image/png")
    p = ArtifactPanel(storage=storage)
    # Not yet added — show_image should insert.
    p.show_image("c", a.sha256)
    assert p.current_artifact().payload["sha256"] == a.sha256
    # Calling again should focus the existing entry.
    p.show_image("c", a.sha256)
    assert len([x for x in p._artifacts
                if x.payload.get("sha256") == a.sha256]) == 1
  finally:
    shutil.rmtree(tmp)


def exercise_reload_skips_unchanged_width():
  """A resize to an unchanged width must not re-run the expensive reload.

  Every resizeEvent used to re-decode the image (get_image) and re-scale
  the full-res pixmap, making drag-resize janky. A width guard should make
  a repeated resize to the same width a no-op for the expensive path.
  """
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  tmp = tempfile.mkdtemp()
  try:
    storage = ConversationStorage(project_dir=Path(tmp), log=null_out())
    att = storage.store_attachment("c", _png_bytes(), "image/png")
    import qttbx.widgets.chat.artifact_panel as ap_mod
    calls = {"n": 0}
    real_get_image = ap_mod.get_image
    def _spy(*args, **kwargs):
      calls["n"] += 1
      return real_get_image(*args, **kwargs)
    ap_mod.get_image = _spy
    try:
      label = ap_mod._ScaledImageLabel(storage, "c", att.sha256)
      label.resize(200, 100)
      ev = QtGui.QResizeEvent(QtCore.QSize(200, 100), QtCore.QSize(120, 80))
      label.resizeEvent(ev)
      first = calls["n"]
      assert first >= 1, first
      # A second resize to the SAME width must not touch get_image again.
      label.resizeEvent(ev)
      assert calls["n"] == first, (first, calls["n"])
    finally:
      ap_mod.get_image = real_get_image
  finally:
    shutil.rmtree(tmp)


def exercise_save_name_extension_follows_mime():
  """The suggested save filename extension follows the payload mime.

  save_current used to hardcode a ``.png`` default, so a JPEG/webp/gif
  artifact was offered a wrong ``.png`` name. The default name should
  derive its extension from the payload mime, falling back to ``.png``.
  """
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  tmp = tempfile.mkdtemp()
  try:
    storage = ConversationStorage(project_dir=Path(tmp), log=null_out())
    att = storage.store_attachment("c", _png_bytes(), "image/jpeg")
    p = ArtifactPanel(storage=storage)
    art = Artifact(kind="image", payload={
      "conv_id": "c", "sha256": att.sha256, "mime": "image/jpeg"})
    name = p._suggested_save_name(art)
    assert name.lower().endswith((".jpg", ".jpeg")), name
    # A missing / unknown mime still falls back to ``.png``.
    art_nomime = Artifact(kind="image", payload={
      "conv_id": "c", "sha256": att.sha256})
    assert p._suggested_save_name(art_nomime).endswith(".png"), \
      p._suggested_save_name(art_nomime)
    # A caption still wins verbatim when present.
    art_cap = Artifact(kind="image", caption="figure.tiff", payload={
      "conv_id": "c", "sha256": att.sha256, "mime": "image/jpeg"})
    assert p._suggested_save_name(art_cap) == "figure.tiff", \
      p._suggested_save_name(art_cap)
  finally:
    shutil.rmtree(tmp)


def exercise():
  exercise_empty_state_has_no_current_artifact()
  exercise_add_image_artifact_auto_advances()
  exercise_user_navigation_pins_index()
  exercise_show_image_focuses_existing_or_inserts()
  exercise_reload_skips_unchanged_width()
  exercise_save_name_extension_follows_mime()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
