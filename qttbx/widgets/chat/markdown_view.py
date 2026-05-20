"""Read-only Markdown view.

QTextBrowser handles inline markdown natively via setMarkdown(); we add a
small stylesheet so code blocks are monospace. Streaming append is
markdown-stitched (re-set the whole document with the accumulated text)
because Qt's incremental markdown support is limited — the cost is
negligible at typical chat message sizes."""

from qttbx.qt import QtGui, QtWidgets

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


class MarkdownView(QtWidgets.QTextBrowser):
  def __init__(self, parent=None):
    super().__init__(parent)
    self.setOpenExternalLinks(True)
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
    self._raw = text or ""
    self.setMarkdown(self._raw)

  def append_markdown(self, text):
    self._raw = (self._raw or "") + (text or "")
    self.setMarkdown(self._raw)

  def clear(self):                                       # noqa: A003
    self._raw = ""
    self.setMarkdown("")
