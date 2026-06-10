"""Read-only Markdown view.

QTextBrowser handles inline markdown natively via setMarkdown(); we add a
small stylesheet so code blocks are monospace. Streaming append is
markdown-stitched (re-set the whole document with the accumulated text)
because Qt's incremental markdown support is limited — the cost is
negligible at typical chat message sizes."""

from qttbx.qt import QtCore, QtGui, QtWidgets

# Stylesheet template; %s is filled with the platform's fixed-pitch
# family at instance-construction time (when a QApplication exists). The
# bare CSS keyword 'monospace' is avoided -- Qt has no font literally
# named 'Monospace' on macOS / Windows, so the keyword forces a one-time
# ~50 ms alias scan and prints a qt.qpa.fonts warning.
_CSS_TEMPLATE = """
pre, code { font-family: '%s'; font-size: 11pt; }
pre { background: #f4f4f4; padding: 6px; border-radius: 4px; }
table { border-collapse: collapse; }
th, td { border: 1px solid #ccc; padding: 3px 6px; }
"""

# Markdown features used for every render. The streamed text is
# assistant/tool-controlled, so we disable raw HTML: Qt's default GitHub
# dialect parses embedded <img>/<a> HTML into live rich text, and
# QTextBrowser then loads the referenced local file (e.g.
# file:///etc/passwd) at paint time. MarkdownNoHTML keeps any raw HTML as
# literal text. getattr-guarded so an older Qt without the flag still
# imports (loadResource below is the second, independent line of defense).
_MD_FEATURES = QtGui.QTextDocument.MarkdownDialectGitHub
_MD_NO_HTML = getattr(QtGui.QTextDocument, "MarkdownNoHTML", None)
if _MD_NO_HTML is not None:
  _MD_FEATURES = _MD_FEATURES | _MD_NO_HTML


class MarkdownView(QtWidgets.QTextBrowser):
  """Read-only QTextBrowser that renders accumulated markdown text."""

  def __init__(self, parent=None):
    super().__init__(parent)
    # The rendered markdown is assistant-controlled, so we must NOT let Qt
    # auto-open whatever URL a link points at (a model could emit
    # file:///... or another local-resource scheme). Disable auto-open and
    # route clicks through _on_anchor_clicked, which only forwards safe
    # web/mail schemes to the OS handler.
    self.setOpenLinks(False)
    self.setOpenExternalLinks(False)
    self.anchorClicked.connect(self._on_anchor_clicked)
    self.setReadOnly(True)
    self.setFrameStyle(QtWidgets.QFrame.NoFrame)
    mono_family = QtGui.QFontDatabase.systemFont(
      QtGui.QFontDatabase.FixedFont).family()
    self.document().setDefaultStyleSheet(_CSS_TEMPLATE % mono_family)
    self._raw = ""
    # Grow with content; the outer ConversationView is the sole
    # scroller (per the chat UI redesign — nested scrollbars are out).
    from qttbx.widgets.chat.auto_height import set_auto_height
    set_auto_height(self)

  def set_markdown(self, text):
    """Replace the document with ``text`` rendered as markdown."""
    self._raw = text or ""
    # Render via the document so the MarkdownNoHTML feature flag takes
    # effect -- the QTextBrowser-level setMarkdown() binding drops the
    # features argument, but QTextDocument.setMarkdown() honors it.
    self.document().setMarkdown(self._raw, _MD_FEATURES)

  def append_markdown(self, text):
    """Append ``text`` to the accumulated markdown and re-render."""
    self._raw = (self._raw or "") + (text or "")
    self.document().setMarkdown(self._raw, _MD_FEATURES)

  def _on_anchor_clicked(self, url):
    """Open a clicked link only if it uses a safe scheme.

    The rendered markdown is assistant-controlled, so link targets are
    untrusted. Only ``http`` / ``https`` / ``mailto`` are forwarded to the
    OS handler; ``file://``, ``smb://`` and custom app schemes are ignored
    so a click cannot open a local file or trigger an external handler.
    In-document anchors (empty scheme with a fragment) scroll within the
    view rather than opening anything.

    Parameters
    ----------
    url : QtCore.QUrl
        The anchor URL emitted by ``QTextBrowser.anchorClicked``.
    """
    scheme = (url.scheme() or "").lower()
    if scheme in ("http", "https", "mailto"):
      QtGui.QDesktopServices.openUrl(url)
    elif not scheme and url.hasFragment():
      # In-document anchor (e.g. a markdown TOC link): scroll within the
      # view instead of opening anything externally.
      self.scrollToAnchor(url.fragment())

  def clear(self):                                       # noqa: A003
    """Reset the view to empty."""
    self._raw = ""
    self.document().setMarkdown("", _MD_FEATURES)

  def loadResource(self, resource_type, name):
    """Refuse to load any external resource.

    The document is assistant/tool-controlled, so a markdown image
    (``![](file:///...)``) or any other embedded reference must never
    read a local file or fetch a URL at render time. Inline chat images
    are rendered separately by ``ImageCell`` (not through this view), so
    blocking every resource load here is safe.

    Parameters
    ----------
    resource_type : int
        The ``QTextDocument.ResourceType`` Qt is requesting.
    name : QtCore.QUrl
        The resource URL embedded in the document.

    Returns
    -------
    QtCore.QByteArray
        Always empty, so nothing is loaded.
    """
    return QtCore.QByteArray()
