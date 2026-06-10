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
    assert "x.png" in p.status_text() or "image" in p.status_text()
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


def exercise():
  exercise_empty_state_has_no_current_artifact()
  exercise_add_image_artifact_auto_advances()
  exercise_user_navigation_pins_index()
  exercise_show_image_focuses_existing_or_inserts()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
