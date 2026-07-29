"""Eliding / wrapping text widgets used across the chat UI.

These exist so that unbounded text -- a title, a model id, a question header,
an option description, an MCP tool name, a tool status -- cannot floor the
minimum width of the ConversationView or the ChatWindow. Every test here is
about that floor and about staying readable once the floor is gone.
"""

import os
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from libtbx.utils import format_cpu_times

try:
  from qttbx.qt import QtCore, QtGui, QtWidgets
except ImportError:
  print("PySide2/PySide6 not available; skipping")
  print("OK")
  sys.exit(0)


LONG = ("MCP error -32000: Connection closed. Failed to connect to the Coot "
        "RPC bridge at 127.0.0.1:44100: [Errno 61] Connection refused.")


def _qapp():
  from qttbx.widgets.font_init import init_default_app_font
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  return app


def _in_layout(widget, width):
  """Put ``widget`` in a shown container of ``width`` and lay it out."""
  app = _qapp()
  host = QtWidgets.QWidget()
  layout = QtWidgets.QHBoxLayout(host)
  layout.setContentsMargins(0, 0, 0, 0)
  widget.setParent(host)
  layout.addWidget(widget, stretch=1)
  host.resize(width, 40)
  host.show()
  app.processEvents()
  return host


def _room_for(widget):
  """Width a host must have for ``widget`` to render its full text.

  Never hardcode a "wide enough" pixel count. Qt on the Windows CI runner
  finds no font files at all -- it logs ``QFontDatabase: Cannot find font
  directory .../conda_base/Library/lib/fonts`` -- and then measures every
  glyph with the fixed-pitch fallback engine, one font pixel size each:
  17 px for the default 12 pt at the offscreen plugin's 100 dpi. LONG wants
  128 * 17 = 2176 px under that against ~1000 px in a desktop proportional
  font, so a literal that looks roomy on a developer's machine leaves the
  text elided on CI, and only there.

  Ask the widget instead: sizeHint reports the FULL text's width by
  contract, which exercise_eliding_label_size_hint_reports_the_full_text
  pins down.
  """
  return widget.sizeHint().width()


def _too_narrow_for(widget):
  """Width that forces ``widget`` to elide, on any font.

  Half of what the text asks for: past the elide threshold under every
  font, and still leaving enough glyphs that an ElideLeft tail is
  recognisable. A fixed 160 px is not -- under a fallback box font that is
  nine glyphs, which keeps no useful end of anything.
  """
  return _room_for(widget) // 2


def exercise_eliding_label_never_floors_its_container():
  """The floor is the whole point: minimumSizeHint width must be zero."""
  from qttbx.widgets.chat.eliding import ElidingLabel
  _qapp()
  label = ElidingLabel()
  label.set_full_text(LONG)
  assert label.minimumSizeHint().width() == 0, \
    label.minimumSizeHint().width()


def exercise_eliding_label_size_hint_reports_the_full_text():
  """sizeHint must describe the FULL text, not the elided copy on screen.

  The widget renders an elision of its text, so an inherited sizeHint would
  shrink along with it -- and a layout that reads that hint would never hand
  the width back when the window grows again. This is what makes the elide
  reversible.
  """
  from qttbx.widgets.chat.eliding import ElidingLabel
  _qapp()
  label = ElidingLabel()
  label.set_full_text(LONG)
  wide = label.sizeHint().width()
  host = _in_layout(label, _too_narrow_for(label))
  assert label.text() != LONG, "precondition: the text should be elided"
  # Tolerance, not equality: the chrome around the text is measured by
  # subtracting the text's advance from the inherited hint, and QLabel sizes
  # its own text a pixel differently than QFontMetrics.horizontalAdvance. The
  # property under test is that the hint still describes the whole string --
  # a collapse would drop it to the half-width the elision occupies.
  assert label.sizeHint().width() >= wide - 2, (
    "sizeHint collapsed from %d to %d once the text was elided; the layout "
    "can never grow it back" % (wide, label.sizeHint().width()))
  del host


def exercise_eliding_label_elides_to_fit_and_grows_back():
  """Narrow: elided, full text on the tooltip. Wide again: whole text back."""
  from qttbx.widgets.chat.eliding import ElidingLabel
  app = _qapp()
  label = ElidingLabel()
  label.set_full_text(LONG)
  room = _room_for(label)
  host = _in_layout(label, _too_narrow_for(label))
  assert "…" in label.text(), label.text()
  assert label.toolTip() == LONG, label.toolTip()
  assert label.full_text() == LONG, "full_text must survive the elide"
  host.resize(room, 40)
  app.processEvents()
  assert label.text() == LONG, label.text()
  assert label.toolTip() == "", \
    "a label with room to render in full needs no tooltip"


def exercise_eliding_label_re_elides_when_the_layout_resizes_it():
  """The re-elide hangs off the widget's own resizeEvent.

  A sibling appearing takes width from the label without resizing the
  container, so hooking the container's resizeEvent leaves the text stale --
  painting hard-clipped, with no ellipsis and no tooltip to recover it.
  """
  from qttbx.widgets.chat.eliding import ElidingLabel
  app = _qapp()
  label = ElidingLabel()
  label.set_full_text(LONG)
  host = _in_layout(label, _room_for(label))
  assert label.text() == LONG, "precondition: the text fits in a full row"
  # A sibling claims half the row -- a plain QLabel floors at its own text, so
  # the whole shortfall lands on the label under test. Half of LONG rather
  # than a fixed run of characters, so the split stays a half on any font: a
  # literal wide enough to squeeze a proportional font squeezes a fallback box
  # font to nothing, and a zero-width label is excused from eliding at all.
  # The host never changes size.
  sibling = QtWidgets.QLabel(LONG[:len(LONG) // 2], host)
  host.layout().addWidget(sibling)
  app.processEvents()
  assert label.width() > 0, "the sibling took the whole row; nothing to elide"
  metrics = label.fontMetrics()
  assert metrics.horizontalAdvance(label.text()) <= label.width(), (
    "stale elide: %d px of text in a %d px label"
    % (metrics.horizontalAdvance(label.text()), label.width()))
  assert label.toolTip() == LONG, label.toolTip()


def exercise_elide_left_keeps_the_tail_visible():
  """ElideLeft is what a path wants: the basename is the informative end."""
  from qttbx.widgets.chat.eliding import ElidingLabel
  _qapp()
  label = ElidingLabel(mode=QtCore.Qt.ElideLeft)
  label.set_full_text("/Users/someone/very/long/path/to/a/chat_debug.log")
  host = _in_layout(label, _too_narrow_for(label))
  assert "…" in label.text(), "precondition: the path should be elided"
  assert label.text().endswith("chat_debug.log"), label.text()
  del host


def exercise_tooltip_override_survives_a_resize():
  """An explicit tooltip must not be clobbered by the next re-elide.

  The debug-log slot renders a deliberately truncated path but wants the
  whole one on hover, so its tooltip cannot follow the default
  'only when the elide dropped something' rule.
  """
  from qttbx.widgets.chat.eliding import ElidingLabel
  app = _qapp()
  label = ElidingLabel()
  label.set_full_text("debug: …/chat_debug.log", tooltip="/full/path.log")
  room = _room_for(label)
  host = _in_layout(label, room)
  assert label.toolTip() == "/full/path.log", label.toolTip()
  host.resize(room // 2, 40)
  app.processEvents()
  assert label.toolTip() == "/full/path.log", (
    "the re-elide clobbered the caller's tooltip: %r" % label.toolTip())


def exercise_set_full_text_coerces_non_str():
  """Callers may pass an exception; rendering it must not raise."""
  from qttbx.widgets.chat.eliding import ElidingLabel
  _qapp()
  label = ElidingLabel()
  label.set_full_text(RuntimeError("bridge is gone"))
  assert "bridge is gone" in label.full_text(), label.full_text()
  label.set_full_text(None)
  assert label.full_text() == "", label.full_text()


def exercise_eliding_check_box_can_shrink_below_its_size_hint():
  """QCheckBox defaults to a Minimum horizontal policy.

  Under that policy a layout treats sizeHint as the floor and never consults
  minimumSizeHint at all -- so zeroing minimumSizeHint alone would look
  correct in isolation and still floor the card. Guards the policy fix.
  """
  from qttbx.widgets.chat.eliding import ElidingCheckBox, ElidingRadioButton
  _qapp()
  for cls in (ElidingCheckBox, ElidingRadioButton):
    btn = cls()
    btn.set_full_text(LONG)
    assert btn.minimumSizeHint().width() == 0, \
      "%s floors at %d" % (cls.__name__, btn.minimumSizeHint().width())
    policy = btn.sizePolicy().horizontalPolicy()
    assert policy != QtWidgets.QSizePolicy.Minimum, (
      "%s keeps the Minimum policy, so the layout ignores minimumSizeHint "
      "and floors at sizeHint anyway" % cls.__name__)
    narrow = _too_narrow_for(btn)
    host = _in_layout(btn, narrow)
    assert btn.width() <= narrow, "%s refused to shrink: %d px" % (
      cls.__name__, btn.width())
    del host


def exercise_becoming_shrinkable_does_not_make_a_widget_growable():
  """Adding shrink must not smuggle in grow.

  Each base picks its horizontal policy deliberately: a QToolButton is Fixed
  so it hugs its text (its autoRaise hover highlight should cover the text and
  no more), while a QCheckBox is Minimum so it fills its row. Blanket-setting
  one policy for all of them gives the button GrowFlag it never had and
  stretches it across the whole bubble. Only ShrinkFlag may be added.
  """
  from qttbx.widgets.chat.eliding import ElidingCheckBox, ElidingLabel
  from qttbx.widgets.chat.eliding import ElidingRadioButton, ElidingToolButton
  _qapp()
  grow = QtWidgets.QSizePolicy.GrowFlag
  for cls, base in ((ElidingToolButton, QtWidgets.QToolButton),
                    (ElidingCheckBox, QtWidgets.QCheckBox),
                    (ElidingRadioButton, QtWidgets.QRadioButton),
                    (ElidingLabel, QtWidgets.QLabel)):
    ours = cls().sizePolicy().horizontalPolicy()
    theirs = base().sizePolicy().horizontalPolicy()
    assert int(ours) & int(QtWidgets.QSizePolicy.ShrinkFlag), \
      "%s cannot shrink, so its zeroed minimumSizeHint is ignored" % cls.__name__
    assert bool(int(ours) & int(grow)) == bool(int(theirs) & int(grow)), (
      "%s changed its grow behaviour: base=%d ours=%d"
      % (cls.__name__, int(theirs), int(ours)))


def exercise_eliding_tool_button_hugs_its_text():
  """A Fixed-policy base must not be stretched by its layout once wrapped."""
  from qttbx.widgets.chat.eliding import ElidingToolButton
  _qapp()
  btn = ElidingToolButton()
  btn.set_full_text("▸ coot_ping (finished, 0.1s)")
  row = _room_for(btn) * 3
  host = _in_layout(btn, row)
  assert btn.width() <= btn.sizeHint().width() + 2, (
    "button grew to %d px in a %d px row; it should hug its %d px of text"
    % (btn.width(), row, btn.sizeHint().width()))
  del host


def exercise_hugging_widget_shows_its_whole_text():
  """At the width the layout hands a hugging widget, no text may be dropped.

  A hugging widget is given exactly its sizeHint width, so if sizeHint asks
  for even a pixel less than the text truly needs, the widget elides text it
  should have shown -- '(finished)' comes out '(finishe…'. horizontalAdvance
  rounds down and elidedText measures the true layout width, so the two must
  not be the pixel budget on either side of the same string. Swept over many
  lengths because whether the rounding bites depends on where each string's
  sub-pixel width lands.
  """
  from qttbx.widgets.chat.eliding import ElidingToolButton, ElidingLabel
  _qapp()
  for cls in (ElidingToolButton, ElidingLabel):
    for name in ("coot_ping", "run_python", "coot_close_unresponsive",
                 "mcp__coot__run_python_multiline",
                 "mcp__structure_tools__compute_real_space_correlation",
                 "phenix_get_program_help", "reduce", "elbow_builder_x"):
      text = "▸ %s (finished)" % name
      w = cls()
      w.set_full_text(text)
      # Ample room; the widget hugs its text.
      host = _in_layout(w, _room_for(w) * 2)
      assert "…" not in w.text(), (
        "%s dropped text with a full row to render in: %r (needs %d px, "
        "sizeHint asked for %d)"
        % (cls.__name__, w.text(),
           w.fontMetrics().size(QtCore.Qt.TextSingleLine, text).width(),
           w.sizeHint().width()))
      assert w.toolTip() == "", (
        "%s hung a tooltip on text that fits: %r" % (cls.__name__, w.text()))
      del host


def exercise_wrapping_label_wraps_and_never_floors():
  """Prose wraps and stays whole; an unbroken token still cannot floor.

  Word wrap alone leaves a wrapped QLabel reporting its widest unbreakable
  token as its minimumSizeHint width, which is why this needs the override
  even though it never elides.
  """
  from qttbx.widgets.chat.eliding import WrappingLabel
  _qapp()
  label = WrappingLabel("some ordinary prose that should wrap and stay whole")
  assert label.wordWrap(), "prose must wrap"
  assert label.minimumSizeHint().width() == 0, \
    label.minimumSizeHint().width()
  unbroken = WrappingLabel("a_single_" + "very_long_" * 12 + "token")
  assert unbroken.minimumSizeHint().width() == 0, (
    "an unbroken token floors the container at %d px"
    % unbroken.minimumSizeHint().width())
  # The text itself is never truncated -- that is the difference from
  # ElidingLabel.
  host = _in_layout(label, 120)
  assert label.text() == "some ordinary prose that should wrap and stay whole"
  del host


_EXERCISES = (
  exercise_eliding_label_never_floors_its_container,
  exercise_eliding_label_size_hint_reports_the_full_text,
  exercise_eliding_label_elides_to_fit_and_grows_back,
  exercise_eliding_label_re_elides_when_the_layout_resizes_it,
  exercise_elide_left_keeps_the_tail_visible,
  exercise_tooltip_override_survives_a_resize,
  exercise_set_full_text_coerces_non_str,
  exercise_eliding_check_box_can_shrink_below_its_size_hint,
  exercise_becoming_shrinkable_does_not_make_a_widget_growable,
  exercise_eliding_tool_button_hugs_its_text,
  exercise_hugging_widget_shows_its_whole_text,
  exercise_wrapping_label_wraps_and_never_floors,
)


def exercise_all_of_it_again_under_a_font_this_machine_lacks():
  """Every check again with metrics nothing like the developer's.

  This whole file is about how much room text takes, so it can be green on
  a machine with fonts and red on one without: the Windows CI runner ships
  no font files, falls back to a fixed-pitch engine at 17 px a glyph, and
  needs more than twice the width macOS does for the same string. That is
  invisible to anyone running the suite locally, so run it twice -- once on
  the real font, once on one wide enough that any width still written as a
  pixel literal cannot be enough.
  """
  app = _qapp()
  original = QtGui.QFont(app.font())
  wide = QtGui.QFont(original)
  # max(): a font set by pixel size rather than points reports -1 here.
  wide.setPointSize(max(24, original.pointSize() * 2))
  app.setFont(wide)
  try:
    for fn in _EXERCISES:
      fn()
  finally:
    app.setFont(original)


def exercise():
  for fn in _EXERCISES:
    fn()
  exercise_all_of_it_again_under_a_font_this_machine_lacks()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
