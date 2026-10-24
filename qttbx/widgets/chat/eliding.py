"""Text widgets that shrink instead of flooring the layout.

A widget that reports its full text width as its ``minimumSizeHint`` pins the
minimum width of every container it sits in. In the chat UI that chain ends at
the ConversationView -- a ``setWidgetResizable`` QScrollArea whose whole job is
to let the bubbles track the window -- and at the ChatWindow itself. So one
long string anywhere (a title, a model id, a question header, an option
description, an MCP tool name, a tool status) stops the bubbles tracking the
window, raises a horizontal scrollbar, and leaves the window refusing to be
resized back down. A force-killed Coot bridge produced exactly that: its
~230-char MCP connection error reached a tool disclosure header.

None of that text is bounded in practice -- it comes from the model, from the
user's ``--model``, or from a tool server -- so a character budget cannot fix
it: 60 characters of a proportional font is still ~880 px. These widgets bound
the *pixels* instead:

- ``minimumSizeHint`` width 0, so the widget never floors its container;
- a horizontal size policy that can shrink, since a QSizePolicy of ``Minimum``
  (the QCheckBox / QRadioButton default) makes the layout treat ``sizeHint``
  as the floor and ignore ``minimumSizeHint`` entirely;
- ``sizeHint`` reporting the FULL text's width, so a layout with room still
  hands the widget its natural width -- the widget holds an *elided* copy of
  the text, and Qt would otherwise derive the hint from that copy and never
  grow it back;
- the text elided to whatever width the widget actually gets, re-elided from
  the widget's own ``resizeEvent``, which fires whenever the layout changes its
  width -- not only when the window is resized. Hooking a container's
  ``resizeEvent`` instead misses a sibling appearing and taking the width.

Use ``WrappingLabel`` for prose, which should wrap rather than elide. Word wrap
alone does not unfloor a label -- a wrapped QLabel's ``minimumSizeHint`` width
is its widest *unbreakable token*, so a lone path or identifier floors it just
as the whole string would -- so it needs the same zeroed minimumSizeHint.
"""

import html

from qttbx.qt import QtCore, QtWidgets

# A few pixels of breathing room added to a hugging widget's natural width, so
# the text is never painted hard against the edge and never sits exactly on the
# rounding boundary where the last glyph would elide.
_HUG_SLACK = 3


class _ElideMixin(object):
  """Elide-to-fit behaviour shared by the label and button variants.

  Mixed in ahead of exactly one Qt text-widget base (QLabel, QToolButton,
  QCheckBox, QRadioButton), so ``super()`` in these methods resolves to that
  base. Subclasses must call ``_init_elide`` from their ``__init__``.
  """

  def _init_elide(self, mode=QtCore.Qt.ElideRight, keep=0):
    self._full_text = ""
    self._elide_mode = mode
    self._keep = keep
    self._tooltip_override = None
    self._shown = ""
    # Make the widget shrinkable and change nothing else. Without ShrinkFlag a
    # layout treats sizeHint as the floor and never consults minimumSizeHint,
    # so zeroing that hint would do nothing at all (QCheckBox / QRadioButton
    # default to Minimum, QToolButton to Fixed -- none of them can shrink).
    #
    # Add the flag to whatever the base chose rather than imposing one policy
    # on all of them: each default carries intent. QToolButton is Fixed so it
    # hugs its text, and its autoRaise hover highlight covers the text and no
    # more; blanket-setting Preferred hands it a GrowFlag it never had and
    # stretches the header across the whole bubble. Fixed with the flag is
    # Maximum (hug, but yield when squeezed), Minimum with it is Preferred,
    # MinimumExpanding with it is Expanding, and every other policy has it
    # already: a QLabel is Preferred and comes through untouched. Named,
    # not computed: the bindings disagree on whether a policy is an int.
    policy = self.sizePolicy()
    horizontal = policy.horizontalPolicy()
    if horizontal == QtWidgets.QSizePolicy.Fixed:
      policy.setHorizontalPolicy(QtWidgets.QSizePolicy.Maximum)
    elif horizontal == QtWidgets.QSizePolicy.Minimum:
      policy.setHorizontalPolicy(QtWidgets.QSizePolicy.Preferred)
    elif horizontal == QtWidgets.QSizePolicy.MinimumExpanding:
      policy.setHorizontalPolicy(QtWidgets.QSizePolicy.Expanding)
    self.setSizePolicy(policy)

  # ---- public API ---------------------------------------------------------

  def set_full_text(self, text, tooltip=None):
    """Set the text to render, elided to whatever width the widget gets.

    Parameters
    ----------
    text : str
        Text to display. Coerced with ``str``, so a caller may pass an
        exception or any other object without crashing the render, and
        stored as Qt will hold it (see ``_as_qt_holds``).
    tooltip : str, optional
        Tooltip to show in place of the default rule, which carries the full
        text only when the elide actually dropped something. Pass this when
        ``text`` is itself already an abbreviation of something longer, so
        the default rule cannot see what was lost: a tool disclosure header
        shortens a long status before it ever gets here, and puts the whole
        of it on the tooltip.
    """
    full = _as_one_line(_as_qt_holds("" if text is None else str(text)))
    override = None if tooltip is None else _as_qt_holds(str(tooltip))
    if full == self._full_text and override == self._tooltip_override:
      return                 # a header refreshed with the status it shows
    self._full_text = full
    self._tooltip_override = override
    # A screen reader gets the whole text as the name: with none, Qt
    # falls back to the rendered text, which may be "Alway…". And as
    # the description, only what Qt would otherwise read wrong. Qt
    # reads the tooltip as the description when there is none, which
    # for plain text is the whole text while a cut copy is shown and
    # nothing while it fits: right as it is, where a description of
    # the text as well had a fitting widget's name read out twice. A
    # tooltip with a tag or an ampersand in it is escaped, and would
    # be read out as the markup it was escaped into: that text, or the
    # caller's tooltip, is the description. Written here, from the
    # text, and not from the elide: a widget is elided once per layout
    # pass, and the passes a new text or font sets off had it say its
    # whole text at one width and nothing at the next, two
    # DescriptionChanged events to a live screen reader for one change.
    # Only when it changes: a new status under the same tooltip has
    # nothing to add.
    self.setAccessibleName(self._full_text)
    if override is not None:
      described = override
    elif _could_be_markup(self._full_text):
      described = self._full_text
    else:
      described = ""
    described = _one_break_per_line(described)
    if described != self.accessibleDescription():
      self.setAccessibleDescription(described)
    # Shown whole first: the chrome is measured around the text shown,
    # and a widget showing none has a hint of its own -- ten px of
    # nothing, for a label -- that is no chrome. The layout then hands
    # it the width the text earns, and the resize elides against that.
    self._show(self._full_text)
    self.updateGeometry()
    self._apply_elide()

  def full_text(self):
    """Return the untruncated text, as this one-line widget holds it."""
    return self._full_text

  # ---- Qt overrides -------------------------------------------------------

  def sizeHint(self):
    """Report the full rendered text's width as the natural width.

    The widget holds an elided copy of the text, so the inherited hint would
    track the elision and the layout would never give the width back. A
    hugging widget (a QToolButton, a checkbox) is handed exactly this width,
    so it has to be enough to render the whole text: measured with the same
    metric the elide uses, plus a few px of slack so the last glyph never
    lands on the rounding boundary and gets dropped.
    """
    return QtCore.QSize(
      self._text_width(self._full_text) + self._chrome_width() + _HUG_SLACK,
      super().sizeHint().height())

  def minimumSizeHint(self):
    """Never floor the container -- the whole point of this class -- or
    by the width of its kept characters, when it was given some."""
    return QtCore.QSize(self._floor_width(),
                        super().minimumSizeHint().height())

  def resizeEvent(self, event):
    super().resizeEvent(event)
    self._apply_elide()

  def showEvent(self, event):
    super().showEvent(event)
    self._apply_elide()

  def changeEvent(self, event):
    super().changeEvent(event)
    # A font change leaves the width alone and the copy cut for the old
    # metrics: grown, it is painted clipped, with the ellipsis the part
    # that falls off; shrunk, it shows less than fits.
    if event.type() == QtCore.QEvent.FontChange:
      self._apply_elide()

  # ---- internals ----------------------------------------------------------

  def _text_width(self, text):
    """Width one line of ``text`` occupies, as the elide and painter see it.

    ``elidedText`` keeps a string whole exactly while the budget is at least
    this, which ``QFontMetrics.horizontalAdvance`` (rounded down) can undercut
    by a pixel -- enough to drop the last glyph at the hug boundary. Measuring
    both the hint and the elide with this one metric keeps them consistent.
    """
    return self.fontMetrics().size(QtCore.Qt.TextSingleLine, text).width()

  def _chrome_width(self):
    """Width this widget adds around its text: margins, frame, indicator."""
    return max(0, super().sizeHint().width()
               - self._text_width(self._shown))

  def _floor_width(self):
    """The width of the ``keep`` characters this widget never drops.

    Zero without any. Otherwise what the elide leaves at that count --
    the head of the text and an ellipsis under ElideRight, the tail
    under ElideLeft, both ends under ElideMiddle -- measured as
    ``sizeHint`` measures the whole, and never more than the whole.
    """
    keep, full = self._keep, self._full_text
    if keep <= 0 or not full:
      return 0
    if len(full) <= keep:
      return self.sizeHint().width()
    mode = self._elide_mode
    if mode == QtCore.Qt.ElideLeft:
      kept = "…" + full[len(full) - keep:]
    elif mode == QtCore.Qt.ElideMiddle:
      # Qt keeps as many either side, so an odd count is kept only with
      # room for the next even one.
      half = (keep + 1) // 2
      kept = full[:half] + "…" + full[len(full) - half:]
    elif mode == QtCore.Qt.ElideNone:
      return self.sizeHint().width()
    else:
      kept = full[:keep] + "…"
    # Never past the whole text: for narrow glyphs the ellipsis alone is
    # wider than the characters it stands in for.
    return min(self._text_width(kept) + self._chrome_width() + _HUG_SLACK,
               self.sizeHint().width())

  def _apply_elide(self):
    self._show(self._elided_text())
    # The tooltip carries the caller's text, or the whole text while a
    # cut copy is shown: an ordinary short label sprouts no redundant
    # hover label.
    if self._tooltip_override is not None:
      carried = self._tooltip_override
    elif self._shown != self._full_text:
      carried = self._full_text
    else:
      carried = ""
    self.setToolTip(self._as_tooltip(carried) if carried else "")

  def _show(self, text):
    """Render ``text``, kept as the copy the chrome is measured around.

    A button paints an ``&`` as a mnemonic -- nothing, with the next
    character underlined -- so one is doubled for the painter.
    """
    self._shown = text
    if isinstance(self, QtWidgets.QAbstractButton):
      text = text.replace("&", "&&")
    self.setText(text)

  def _elided_text(self):
    """The full text cut down to the width the widget actually has."""
    # Until it is laid out, the widget's width is a Qt default (100 px for a
    # fresh QLabel), not an assigned one -- eliding against that would cut
    # text that is about to have room. Hold the full text instead: the zeroed
    # minimumSizeHint means it still cannot floor anything meanwhile, and the
    # resizeEvent / showEvent from the first layout pass re-elides it against
    # the real width.
    if not self.isVisible():
      return self._full_text
    available = self.width() - self._chrome_width()
    # Show the whole text once it fits, rather than leaning on elidedText at
    # the boundary: at equal width elidedText can still drop the last glyph to
    # its ellipsis, which is the '(finishe…' this guards against.
    if available >= self._text_width(self._full_text):
      return self._full_text
    # Never the empty string for a text that has something to show.
    # _chrome_width measures around the CURRENT text, and a QCheckBox with
    # no text drops its label margin (and six px of height): the chrome
    # shrinks, the next budget fits an ellipsis, the chrome grows back, the
    # ellipsis no longer fits -- and the text flips between '' and '…' on
    # every pass, a resize loop at 100% CPU in a band of widths just above
    # the indicator's. An ellipsis the widget cannot fit is clipped, and
    # still carries the full text on the tooltip. An empty full text stays
    # empty: an icon-only button must not grow an ellipsis.
    #
    # elidedText keeps the string whole from the budget of its advance,
    # which an italic or bold glyph overhangs, and the whole does not
    # fit: a budget just short of the advance keeps the copy a cut one,
    # rather than the whole text shown clipped, with no tooltip.
    budget = min(available,
                 self.fontMetrics().horizontalAdvance(self._full_text) - 1)
    return (self.fontMetrics().elidedText(
      self._full_text, self._elide_mode, budget)
      or ("…" if self._full_text else ""))

  def _as_tooltip(self, text):
    """``text`` as a tooltip that reads the way this widget renders it.

    A tooltip auto-detects rich text, so markup in a string the widget
    shows literally -- every string here: a button's text is always
    plain, and a label shows markup as it is -- would render on hover: a
    tool name that the approval card displays verbatim. Convert such a
    string
    to escaped HTML, its whitespace kept: a bare ``<p>`` folded the
    indentation of a traceback a disclosure puts here, and its newlines.
    Ordinary text is left alone but for its breaks, so it stays
    comparable.
    """
    text = _one_break_per_line(text)
    if not _could_be_markup(text):
      return text
    return ('<p style="white-space:pre-wrap">%s</p>'
            % html.escape(text, quote=False))


def _as_qt_holds(text):
  """``text`` as Qt will hold it once the binding has handed it over, so
  that the widget measures what it shows.

  A lone surrogate -- what a path with a byte that is not UTF-8 decodes
  to -- has no UTF-16 form: Qt drew it as a glyph box, gave it back as a
  "?" nine pixels narrower, and a disclosure header quoting such a path
  measured its chrome from the one and its text from the other, flipping
  between whole and elided on every layout pass at 100% CPU. A NUL cut
  the string short, "phenix" standing for "phenix<NUL>start_job"; the
  replacement character marks both; a pair of surrogates is the one
  character Qt makes of it. A leading byte order mark is
  dropped on the way in -- one per pass, so a doubled one would be
  held with one -- and a reversed one, U+FFFE, has the whole string
  read with its bytes swapped: the marks go, and the reversed one is
  marked. The string terminator U+009C cuts the string short as the
  NUL does, and is marked the same way.
  """
  text = text.lstrip("\ufeff")
  if text.startswith("\ufffe"):
    text = "\ufffd" + text[1:]
  text = text.replace("\0", "\ufffd").replace("\u009c", "\ufffd")
  return text.encode("utf-16-le", "surrogatepass").decode(
    "utf-16-le", "replace")


def _as_one_line(text):
  """``text`` as a widget of one line paints it, so that the one-line
  measure counts what is painted.

  A tab is painted at a stop the measure does not know, and a line
  break -- LF, CR, or a Unicode line or paragraph separator -- starts
  a line the painter shows and the measure does not see, or cuts the
  measure short: a disclosure header with a tab in its status measured
  22 px of chrome around nothing. Each becomes a space, a CR LF pair
  one space.
  """
  text = text.replace("\r\n", " ")
  for breaker in "\t\n\r\u2028\u2029":
    text = text.replace(breaker, " ")
  return text


def _one_break_per_line(text):
  """``text`` with each kind of line break made the one LF it is, and
  its lone CRs made spaces.

  The rich-text engine makes a break of every CR as well as every LF,
  so a traceback from Windows came out double-spaced; a lone CR -- a
  progress bar's, or trailing an LF -- is painted as a space in plain
  text, and dropped it glued the words either side of it. A Unicode
  line or paragraph separator is a break to that engine and nothing to
  plain text: as an LF, the tooltip and its description agree on it.
  """
  return (text.replace("\r\n", "\n").replace("\r", " ")
          .replace("\u2028", "\n").replace("\u2029", "\n"))


def _could_be_markup(text):
  """Whether a tooltip might take ``text`` for rich text.

  A superset of Qt's test (``QtGui.Qt.mightBeRichText``, on PySide2
  and PySide6 alike), which looks no further than the first line for a
  tag it knows or a literal ``&lt;``: escaping a string Qt would have
  shown literally anyway costs at most a tooltip that word-wraps,
  while a leaked tag renders.
  """
  return "<" in text or "&" in text


class ElidingLabel(_ElideMixin, QtWidgets.QLabel):
  """QLabel that elides its text to fit and never floors its container.

  The text is shown as it is. Every label here shows plain or untrusted
  text -- a question card's header is a chip of a few words the model
  writes -- so a tag in it stays a tag, where QLabel left to decide
  rendered it, images from local paths included, and the elide could
  cut inside it.

  Parameters
  ----------
  parent : QtWidgets.QWidget, optional
      Parent widget.
  mode : QtCore.Qt.TextElideMode, optional
      Where to drop characters. ``ElideLeft`` keeps the tail visible, which
      is what a path wants.
  keep : int, optional
      Characters the label always has room for, so that it floors its
      container by their width rather than by nothing: for a label whose
      first few characters are the signal -- the kind of a tool's source
      -- and must outlast any squeeze. Under ElideRight they are the
      head of the text; under ElideLeft the tail; under ElideMiddle both
      ends. Zero, the default, floors nothing.
  """

  def __init__(self, parent=None, mode=QtCore.Qt.ElideRight, keep=0):
    super().__init__(parent)
    self.setTextFormat(QtCore.Qt.PlainText)
    self._init_elide(mode, keep)


class ElidingToolButton(_ElideMixin, QtWidgets.QToolButton):
  """QToolButton that elides its text to fit and never floors its container."""

  def __init__(self, parent=None, mode=QtCore.Qt.ElideRight):
    super().__init__(parent)
    self._init_elide(mode)


class ElidingCheckBox(_ElideMixin, QtWidgets.QCheckBox):
  """QCheckBox that elides its text to fit and never floors its container."""

  def __init__(self, parent=None, mode=QtCore.Qt.ElideRight):
    super().__init__(parent)
    self._init_elide(mode)


class ElidingRadioButton(_ElideMixin, QtWidgets.QRadioButton):
  """QRadioButton that elides its text and never floors its container."""

  def __init__(self, parent=None, mode=QtCore.Qt.ElideRight):
    super().__init__(parent)
    self._init_elide(mode)


class WrappingLabel(QtWidgets.QLabel):
  """QLabel that wraps prose and still never floors its container.

  Word wrap on its own is not enough: a wrapped QLabel reports its widest
  unbreakable token as its ``minimumSizeHint`` width, so a single long path,
  URL, or snake_case identifier floors the container exactly as the whole
  string would. Prose belongs here -- wrapped and fully readable. A tag, an
  id, or anything that must stay on one line belongs in ``ElidingLabel``.

  The text is shown as it is, like an ``ElidingLabel``'s: the prose is
  the model's, and QLabel left to decide rendered a tag it found in it,
  images from local paths included.

  Parameters
  ----------
  text : str, optional
      Initial text.
  parent : QtWidgets.QWidget, optional
      Parent widget.
  """

  def __init__(self, text="", parent=None):
    super().__init__(parent)
    self.setTextFormat(QtCore.Qt.PlainText)
    self.setWordWrap(True)
    self.setText(text)

  def minimumSizeHint(self):
    """Zero width: the wrapped text reflows to whatever width it is given."""
    return QtCore.QSize(0, super().minimumSizeHint().height())
