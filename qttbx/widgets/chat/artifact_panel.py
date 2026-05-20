"""Side artifact panel.

Holds a history of ``Artifact`` instances and renders the current one
via the renderer registry in ``artifact.py``. Ships with an image
renderer; other kinds can be added via ``register_renderer``."""

from qttbx.qt import QtCore, QtGui, QtWidgets

from qttbx.widgets.chat.artifact import (
  Artifact, register_renderer, renderer_for)
from qttbx.widgets.chat.image_cache import get_image


class ArtifactPanel(QtWidgets.QWidget):

  def __init__(self, parent=None, storage=None):
    super().__init__(parent)
    self._storage = storage
    self._artifacts = []                              # list[Artifact]
    self._index = -1
    self._user_navigated = False
    self._build_ui()
    _ensure_default_renderers_registered()

  # ---- UI ------------------------------------------------------------------

  def _build_ui(self):
    layout = QtWidgets.QVBoxLayout(self)
    layout.setContentsMargins(8, 8, 8, 8)
    self._caption_label = QtWidgets.QLabel("No artifacts yet", self)
    self._caption_label.setAlignment(QtCore.Qt.AlignCenter)
    self._caption_label.setWordWrap(True)
    layout.addWidget(self._caption_label)
    self._stage = QtWidgets.QStackedLayout()
    self._stage_widget = QtWidgets.QWidget(self)
    self._stage_widget.setLayout(self._stage)
    layout.addWidget(self._stage_widget, stretch=1)
    nav = QtWidgets.QHBoxLayout()
    self._prev_btn = QtWidgets.QPushButton("<", self)
    self._next_btn = QtWidgets.QPushButton(">", self)
    self._save_btn = QtWidgets.QPushButton("Save", self)
    self._copy_btn = QtWidgets.QPushButton("Copy", self)
    self._counter_label = QtWidgets.QLabel("0 / 0", self)
    self._prev_btn.clicked.connect(self.go_prev)
    self._next_btn.clicked.connect(self.go_next)
    self._save_btn.clicked.connect(self.save_current)
    self._copy_btn.clicked.connect(self.copy_current)
    nav.addWidget(self._prev_btn)
    nav.addWidget(self._counter_label)
    nav.addWidget(self._next_btn)
    nav.addStretch(1)
    nav.addWidget(self._save_btn)
    nav.addWidget(self._copy_btn)
    layout.addLayout(nav)
    self._update_controls()

  # ---- data ----------------------------------------------------------------

  def set_storage(self, storage):
    self._storage = storage

  def add_artifact(self, artifact):
    """Append; advance to it if the user hasn't navigated away."""
    self._artifacts.append(artifact)
    if not self._user_navigated:
      self._index = len(self._artifacts) - 1
      self._refresh()
    else:
      self._update_controls()

  def show_image(self, conv_id, sha256):
    """Focus an image (already in history if previously seen, else add
    a fresh Artifact)."""
    for i, art in enumerate(self._artifacts):
      payload = art.payload or {}
      if art.kind == "image" and payload.get("sha256") == sha256 \
         and payload.get("conv_id") == conv_id:
        self._index = i
        self._user_navigated = True
        self._refresh()
        return
    art = Artifact(kind="image", payload={
      "conv_id": conv_id, "sha256": sha256, "mime": "image/png"})
    self._artifacts.append(art)
    self._index = len(self._artifacts) - 1
    self._user_navigated = True
    self._refresh()

  def clear(self):
    self._artifacts = []
    self._index = -1
    self._user_navigated = False
    self._refresh()

  # ---- nav -----------------------------------------------------------------

  def go_prev(self):
    if self._index > 0:
      self._index -= 1
      self._user_navigated = True
      self._refresh()

  def go_next(self):
    if self._index < len(self._artifacts) - 1:
      self._index += 1
      # If user is at the latest after pressing >, treat as caught up.
      if self._index == len(self._artifacts) - 1:
        self._user_navigated = False
      self._refresh()

  def current_artifact(self):
    if 0 <= self._index < len(self._artifacts):
      return self._artifacts[self._index]
    return None

  def status_text(self):
    """Return the current caption label text (used by tests)."""
    return self._caption_label.text()

  # ---- refresh -------------------------------------------------------------

  def _refresh(self):
    while self._stage.count():
      w = self._stage.widget(0)
      self._stage.removeWidget(w)
      w.setParent(None)
      w.deleteLater()
    art = self.current_artifact()
    if art is None:
      self._caption_label.setText("No artifacts yet")
      self._update_controls()
      return
    factory = renderer_for(art.kind)
    if factory is None:
      widget = QtWidgets.QLabel(
        "No renderer for artifact kind: %s" % art.kind, self)
    else:
      widget = factory(art, self._storage)
    self._stage.addWidget(widget)
    cap = art.caption or art.payload.get("filename") or ""
    src = art.source or ""
    label = " - ".join(s for s in (cap, src) if s) or art.kind
    self._caption_label.setText(label)
    self._update_controls()

  def _update_controls(self):
    n = len(self._artifacts)
    cur = (self._index + 1) if self._index >= 0 else 0
    self._counter_label.setText("%d / %d" % (cur, n))
    self._prev_btn.setEnabled(self._index > 0)
    self._next_btn.setEnabled(0 <= self._index < n - 1)
    self._save_btn.setEnabled(self._index >= 0)
    self._copy_btn.setEnabled(self._index >= 0)

  # ---- buttons -------------------------------------------------------------

  def save_current(self):
    art = self.current_artifact()
    if art is None or self._storage is None:
      return
    if art.kind != "image":
      return
    payload = art.payload or {}
    sha = payload.get("sha256", "")
    conv_id = payload.get("conv_id", "")
    suggested = (art.caption or "%s.png" % sha[:8])
    path, _ = QtWidgets.QFileDialog.getSaveFileName(
      self, "Save image", suggested,
      "Images (*.png *.jpg *.jpeg *.webp *.gif)")
    if not path:
      return
    try:
      data = self._storage.load_attachment(conv_id, sha)
      with open(path, "wb") as fh:
        fh.write(data)
    except Exception as exc:
      QtWidgets.QMessageBox.warning(
        self, "Save failed", str(exc))

  def copy_current(self):
    art = self.current_artifact()
    if art is None or self._storage is None or art.kind != "image":
      return
    payload = art.payload or {}
    img = get_image(self._storage,
                    payload.get("conv_id", ""),
                    payload.get("sha256", ""))
    QtWidgets.QApplication.clipboard().setImage(img)


# ---- default image renderer -------------------------------------------------


def _image_renderer(artifact, storage):
  """Default renderer for kind='image'. Loads via image cache; scales to
  panel width on resize via the wrapper widget."""
  payload = artifact.payload or {}
  conv_id = payload.get("conv_id", "")
  sha = payload.get("sha256", "")
  return _ScaledImageLabel(storage, conv_id, sha)


class _ScaledImageLabel(QtWidgets.QLabel):
  def __init__(self, storage, conv_id, sha256, parent=None):
    super().__init__(parent)
    self._storage = storage
    self._conv_id = conv_id
    self._sha256 = sha256
    self.setAlignment(QtCore.Qt.AlignCenter)
    self.setMinimumSize(120, 80)
    self._reload()

  def _reload(self):
    if self._storage is None or not self._sha256:
      self.setText("(no image)")
      return
    img = get_image(self._storage, self._conv_id, self._sha256)
    pix = QtGui.QPixmap.fromImage(img)
    if self.width() > 0:
      pix = pix.scaledToWidth(
        self.width(), QtCore.Qt.SmoothTransformation)
    self.setPixmap(pix)

  def resizeEvent(self, event):
    super().resizeEvent(event)
    self._reload()


_default_renderers_registered = False


def _ensure_default_renderers_registered():
  global _default_renderers_registered
  if _default_renderers_registered:
    return
  register_renderer("image", _image_renderer)
  _default_renderers_registered = True
