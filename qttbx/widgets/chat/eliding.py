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

  def _init_elide(self, mode=QtCore.Qt.ElideRight):
    self._full_text = ""
    self._elide_mode = mode
    self._tooltip_override = None
    # Make the widget shrinkable and change nothing else. Without ShrinkFlag a
    # layout treats sizeHint as the floor and never consults minimumSizeHint,
    # so zeroing that hint would do nothing at all (QCheckBox / QRadioButton
    # default to Minimum, QToolButton to Fixed -- none of them can shrink).
    #
    # OR the flag into whatever the base chose rather than imposing one policy
    # on all of them: each default carries intent. QToolButton is Fixed so it
    # hugs its text, and its autoRaise hover highlight covers the text and no
    # more; blanket-setting Preferred hands it a GrowFlag it never had and
    # stretches the header across the whole bubble. Fixed|Shrink is Maximum
    # (hug, but yield when squeezed), Minimum|Shrink is Preferred, and a
    # QLabel is already Preferred and comes through untouched.
    policy = self.sizePolicy()
    policy.setHorizontalPolicy(QtWidgets.QSizePolicy.Policy(
      int(policy.horizontalPolicy()) | int(QtWidgets.QSizePolicy.ShrinkFlag)))
    self.setSizePolicy(policy)

  # ---- public API ---------------------------------------------------------

  def set_full_text(self, text, tooltip=None):
    """Set the text to render, elided to whatever width the widget gets.

    Parameters
    ----------
    text : str
        Text to display. Coerced with ``str``, so a caller may pass an
        exception or any other object without crashing the render.
    tooltip : str, optional
        Tooltip to show in place of the default rule, which carries the full
        text only when the elide actually dropped something. Pass this when
        ``text`` is itself already an abbreviation of something longer, so
        the default rule cannot see what was lost: a tool disclosure header
        shortens a long status before it ever gets here, and puts the whole
        of it on the tooltip.
    """
    self._full_text = "" if text is None else str(text)
    self._tooltip_override = None if tooltip is None else str(tooltip)
    self.updateGeometry()
    self._apply_elide()

  def full_text(self):
    """Return the untruncated text, regardless of what is rendered."""
    return self._full_text

  # ---- Qt overrides -------------------------------------------------------

  def sizeHint(self):
    """Report the FULL text's width as the natural width.

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
    """Never floor the container -- the whole point of this class."""
    return QtCore.QSize(0, super().minimumSizeHint().height())

  def resizeEvent(self, event):
    super().resizeEvent(event)
    self._apply_elide()

  def showEvent(self, event):
    super().showEvent(event)
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
               - self._text_width(self.text()))

  def _apply_elide(self):
    self.setText(self._elided_text())
    self._apply_tooltip(self.text())

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
    if available <= 0:
      return self._full_text
    # Show the whole text once it fits, rather than leaning on elidedText at
    # the boundary: at equal width elidedText can still drop the last glyph to
    # its ellipsis, which is the '(finishe…' this guards against.
    if available >= self._text_width(self._full_text):
      return self._full_text
    return self.fontMetrics().elidedText(
      self._full_text, self._elide_mode, available)

  def _apply_tooltip(self, elided):
    if self._tooltip_override is not None:
      self.setToolTip(self._tooltip_override)
      return
    # Only carry a tooltip when something was actually dropped, so an
    # ordinary short label doesn't sprout a redundant hover label.
    self.setToolTip(self._full_text if elided != self._full_text else "")


class ElidingLabel(_ElideMixin, QtWidgets.QLabel):
  """QLabel that elides its text to fit and never floors its container.

  Parameters
  ----------
  parent : QtWidgets.QWidget, optional
      Parent widget.
  mode : QtCore.Qt.TextElideMode, optional
      Where to drop characters. ``ElideLeft`` keeps the tail visible, which
      is what a path wants.
  """

  def __init__(self, parent=None, mode=QtCore.Qt.ElideRight):
    super().__init__(parent)
    self._init_elide(mode)


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

  Parameters
  ----------
  text : str, optional
      Initial text.
  parent : QtWidgets.QWidget, optional
      Parent widget.
  """

  def __init__(self, text="", parent=None):
    super().__init__(text, parent)
    self.setWordWrap(True)

  def minimumSizeHint(self):
    """Zero width: the wrapped text reflows to whatever width it is given."""
    return QtCore.QSize(0, super().minimumSizeHint().height())
