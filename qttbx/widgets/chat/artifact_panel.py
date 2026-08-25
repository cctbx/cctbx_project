"""Side artifact panel.

Holds a history of ``Artifact`` instances and renders the current one
via the renderer registry in ``artifact.py``. Ships with an image
renderer; other kinds can be added via ``register_renderer``.
"""

from qttbx.qt import QtCore, QtGui, QtWidgets

from qttbx.widgets.chat.artifact import (
  Artifact, register_renderer, renderer_for)
from qttbx.widgets.chat.image_cache import get_image


class ArtifactPanel(QtWidgets.QWidget):
  """Side panel showing a navigable history of artifacts.

  Parameters
  ----------
  parent : QtWidgets.QWidget, optional
      Parent widget.
  storage : object, optional
      Attachment store used to load and save artifact payloads.
  """

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

  def add_artifact(self, artifact):
    """Append an artifact, advancing to it unless the user navigated away.

    An image already in history (same ``conv_id`` + ``sha256``, the identity
    rule ``show_image`` uses) is skipped without disturbing the current
    selection -- the live emit path and the turn-done / tool-results rescans
    may push the same image more than once.

    Parameters
    ----------
    artifact : Artifact
        The artifact to append to the history.
    """
    payload = artifact.payload or {}
    if artifact.kind == "image" and payload.get("sha256") and \
       self._find_image(payload.get("conv_id"), payload.get("sha256")) >= 0:
      return
    self._artifacts.append(artifact)
    if not self._user_navigated:
      self._index = len(self._artifacts) - 1
      self._refresh()
    else:
      self._update_controls()

  def show_image(self, conv_id, sha256):
    """Focus an image, adding a fresh ``Artifact`` if not already shown.

    If the image was previously seen it is reused from history;
    otherwise a fresh ``Artifact`` is appended.

    Parameters
    ----------
    conv_id : str
        Conversation id the image belongs to.
    sha256 : str
        Content hash identifying the image.
    """
    i = self._find_image(conv_id, sha256)
    if i < 0:
      self._artifacts.append(Artifact(kind="image", payload={
        "conv_id": conv_id, "sha256": sha256, "mime": "image/png"}))
      i = len(self._artifacts) - 1
    self._index = i
    self._user_navigated = True
    self._refresh()

  def _find_image(self, conv_id, sha256):
    """Return the history index of the image with this identity, or -1.

    Parameters
    ----------
    conv_id : str
        Conversation id the image belongs to.
    sha256 : str
        Content hash identifying the image.
    """
    for i, art in enumerate(self._artifacts):
      payload = art.payload or {}
      if art.kind == "image" and payload.get("sha256") == sha256 \
         and payload.get("conv_id") == conv_id:
        return i
    return -1

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
    """Return the currently selected ``Artifact``, or ``None`` if empty."""
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

  def _suggested_save_name(self, art):
    """Return the default filename offered by the save dialog.

    Prefers the artifact caption verbatim; otherwise builds
    ``<sha8><ext>`` where ``ext`` comes from the payload ``mime`` (via
    ``mimetypes.guess_extension``), mirroring how ``store_attachment``
    derives extensions and falling back to ``.png`` when the mime is
    missing or unknown.

    Parameters
    ----------
    art : Artifact
        The artifact whose default save name is being computed.

    Returns
    -------
    str
        The suggested default filename.
    """
    import mimetypes
    if art.caption:
      return art.caption
    payload = art.payload or {}
    sha = payload.get("sha256", "")
    mime = payload.get("mime")
    ext = (mimetypes.guess_extension(mime) if mime else None) or ".png"
    return "%s%s" % (sha[:8], ext)

  def save_current(self):
    """Prompt for a path and save the current image artifact to disk."""
    art = self.current_artifact()
    if art is None or self._storage is None:
      return
    if art.kind != "image":
      return
    payload = art.payload or {}
    sha = payload.get("sha256", "")
    conv_id = payload.get("conv_id", "")
    suggested = self._suggested_save_name(art)
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
    """Copy the current image artifact to the system clipboard."""
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
  """Render an ``image`` artifact as a width-scaling label.

  Loads the image via the image cache and scales it to the panel width
  on resize through the wrapper widget.

  Parameters
  ----------
  artifact : Artifact
      Artifact whose ``payload`` carries ``conv_id`` and ``sha256``.
  storage : object
      Attachment store the image is loaded from.

  Returns
  -------
  QtWidgets.QWidget
      A label widget that rescales the image on resize.
  """
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
    self._last_width = None
    self._reload()

  def _reload(self):
    if self._storage is None or not self._sha256:
      self.setText("(no image)")
      return
    width = self.width()
    # resizeEvent fires continuously during a drag-resize; skip the expensive
    # decode + smooth rescale (and the storage re-read) when the target width
    # is unchanged since the last reload.
    if width == self._last_width:
      return
    self._last_width = width
    img = get_image(self._storage, self._conv_id, self._sha256)
    pix = QtGui.QPixmap.fromImage(img)
    if width > 0:
      pix = pix.scaledToWidth(
        width, QtCore.Qt.SmoothTransformation)
    self.setPixmap(pix)

  def resizeEvent(self, event):
    super().resizeEvent(event)
    self._reload()


def _ensure_default_renderers_registered():
  # Idempotent one-time registration of the default image renderer. Querying
  # the registry directly is the natural idempotency check -- no separate
  # module-level "already ran" flag to keep in sync.
  if renderer_for("image") is None:
    register_renderer("image", _image_renderer)
