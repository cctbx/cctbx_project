"""Read-only Markdown view.

QTextBrowser handles inline markdown natively via setMarkdown(); we add a
small stylesheet so code blocks are monospace. Streaming append is
markdown-stitched (re-set the whole document with the accumulated text)
because Qt's incremental markdown support is limited — the cost is
negligible at typical chat message sizes."""

import re

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

# A lone '~' immediately before a digit is the "~2.5" ("approximately")
# shorthand the assistant emits for approximate numbers -- it is content,
# not markup. But the GitHub dialect above treats '~' as a GFM strikethrough
# marker, so two such tildes in one block (e.g. '~0.02 ... ~0.001') get
# paired and everything between them renders struck through. Escaping the
# tilde to a literal '~' before rendering breaks that pairing while keeping
# the character visible.
_APPROX_TILDE_RE = re.compile(r"(?<!~)~(?=\d)")

# A markdown inline code span: a run of backticks, the shortest run of
# characters, then a matching run of the same length. Tildes inside such a
# span (and inside fenced code blocks, bracketed line-by-line below) are
# literal code, never strikethrough, so the approximate-tilde escape skips them.
_INLINE_CODE_RE = re.compile(r"(`+)(?:.+?)\1")

# A fenced-code delimiter line: up to three leading spaces then a run of 3+
# backticks or 3+ tildes (CommonMark). Used to bracket fenced blocks.
_CODE_FENCE_RE = re.compile(r"^ {0,3}(`{3,}|~{3,})")


def _escape_outside_code_spans(line):
  """Apply the approximate-tilde escape to ``line`` but leave inline ``code``
  spans untouched -- a ``~`` inside a span is literal code, not strikethrough."""
  out = []
  pos = 0
  for m in _INLINE_CODE_RE.finditer(line):
    out.append(_APPROX_TILDE_RE.sub(r"\\~", line[pos:m.start()]))
    out.append(m.group(0))                 # code span: left verbatim
    pos = m.end()
  out.append(_APPROX_TILDE_RE.sub(r"\\~", line[pos:]))
  return "".join(out)


def _escape_approx_tildes(text):
  """Escape a lone ``~`` that precedes a digit so GFM cannot pair it into an
  unwanted strikethrough; the escaped ``\\~`` renders as a literal ``~``.

  The filter is deliberately NARROW -- it matches ONLY a ``~`` directly
  followed by a digit (the ``~2.5`` "approximately" shorthand). It is not a
  general strikethrough switch: a genuine ``~~word~~`` or ``~word~`` is left
  untouched and still renders struck through, and ordinary tildes
  (``~/path``, a bare ``~``) and ``~~`` runs are ignored. That narrowness is
  intentional -- the assistant only ever writes ``~`` before a number, so
  this neutralizes exactly the accidental pattern and nothing else, which
  lets the view stay on the GitHub dialect (tables etc. keep rendering)
  instead of disabling strikethrough wholesale or re-parsing the markdown.

  Fenced code blocks (``` / ~~~) and inline ``code`` spans are skipped: GFM
  strikethrough cannot occur inside code, so there a ``~`` is literal (a
  version ``~1.2.0``, a size ``~5GB``) and a backslash escape would render /
  copy as a spurious character.

  Parameters
  ----------
  text : str
      Raw markdown about to be handed to ``setMarkdown``.

  Returns
  -------
  str
      ``text`` with each lone ``~``-before-a-digit that lies outside code
      escaped to ``\\~``; text without that pattern is returned unchanged.
  """
  if not text:
    return text
  out = []
  fence = None             # (char, length) while inside a fenced block, else None
  for line in text.splitlines(keepends=True):
    body = line.rstrip("\r\n")
    m = _CODE_FENCE_RE.match(body)
    if fence is None:
      if m:
        delim = m.group(1)
        fence = (delim[0], len(delim))     # opening fence; its info string is code
        out.append(line)
      else:
        out.append(_escape_outside_code_spans(line))
    else:
      out.append(line)                     # inside a fenced block: never escape
      # A closing fence: same char, at least as long, nothing but whitespace
      # after the run.
      if m:
        delim = m.group(1)
        if (delim[0] == fence[0] and len(delim) >= fence[1]
            and body[m.end():].strip() == ""):
          fence = None
  return "".join(out)


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

  def append_markdown(self, text):
    """Append ``text`` to the accumulated markdown and re-render."""
    self._raw = (self._raw or "") + (text or "")
    # Render via the document so the MarkdownNoHTML feature flag takes
    # effect -- the QTextBrowser-level setMarkdown() binding drops the
    # features argument, but QTextDocument.setMarkdown() honors it.
    # _escape_approx_tildes neutralizes only the '~<digit>' "approximately"
    # shorthand (see its docstring) so it cannot form an accidental GFM
    # strikethrough; self._raw stays the true original so the chat export
    # stays faithful.
    self.document().setMarkdown(
      _escape_approx_tildes(self._raw), _MD_FEATURES)

  def searchable_cells(self):
    """This cell's searchable text: the view itself, kind ``"text"``."""
    return [("text", self)]

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
