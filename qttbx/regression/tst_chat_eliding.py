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
  # font to nothing, and a zero-width label shows an ellipsis whatever it
  # was told, which proves nothing about the re-elide. The host never
  # changes size.
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
  # The policies by the flags they carry, as the bindings spell them;
  # computing on the enum raises on PySide6.
  P = QtWidgets.QSizePolicy
  shrinks = (P.Maximum, P.Preferred, P.Expanding, P.Ignored)
  grows = (P.Minimum, P.Preferred, P.MinimumExpanding, P.Expanding,
           P.Ignored)
  for cls, base in ((ElidingToolButton, QtWidgets.QToolButton),
                    (ElidingCheckBox, QtWidgets.QCheckBox),
                    (ElidingRadioButton, QtWidgets.QRadioButton),
                    (ElidingLabel, QtWidgets.QLabel)):
    ours = cls().sizePolicy().horizontalPolicy()
    theirs = base().sizePolicy().horizontalPolicy()
    assert ours in shrinks, (
      "%s cannot shrink, so its zeroed minimumSizeHint is ignored: %s"
      % (cls.__name__, ours))
    assert (ours in grows) == (theirs in grows), (
      "%s changed its grow behavior: base=%s ours=%s"
      % (cls.__name__, theirs, ours))
  # A base that grows and expands but cannot shrink gets to shrink too.
  label = ElidingLabel()
  label.setSizePolicy(P.MinimumExpanding, P.Preferred)
  label._init_elide()
  assert label.sizePolicy().horizontalPolicy() == P.Expanding, (
    label.sizePolicy().horizontalPolicy())


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


def exercise_wrapping_label_shows_markup_literally():
  """A wrapping label shows markup as the text it is.

  A question card's body and option descriptions are the model's
  words, and QLabel left to decide rendered a tag it found in them --
  an image from a local path included -- where the text's own angle
  brackets and Markdown were left alone.
  """
  from qttbx.widgets.chat.eliding import WrappingLabel
  _qapp()
  for text in ("Delete residues A 12-20 (<b>irreversible</b>)?",
               '<img src="file:///etc/passwd"> Proceed?'):
    label = WrappingLabel(text)
    assert label.textFormat() == QtCore.Qt.PlainText, label.textFormat()
    assert label.text() == text, label.text()
    label = WrappingLabel()
    label.setText(text)
    assert label.textFormat() == QtCore.Qt.PlainText, label.textFormat()


def exercise_eliding_button_never_renders_empty_and_settles():
  """A squeezed checkbox must settle, at every width, on a non-empty text.

  _chrome_width measures around the CURRENT text, and a QCheckBox with no
  text drops its label margin: the chrome shrinks, the next budget fits an
  ellipsis, the chrome grows back, the ellipsis no longer fits -- and the
  text flips between '' and '…' on every pass, a resize loop at 100% CPU
  in a band of widths just above the indicator's. Seen with the approval
  card's "Always allow this tool" in a 335 to 347 px view.
  """
  from qttbx.widgets.chat.eliding import ElidingCheckBox, ElidingRadioButton
  app = _qapp()
  for cls in (ElidingCheckBox, ElidingRadioButton):
    btn = cls()
    btn.set_full_text(LONG)
    room = _room_for(btn)
    host = _in_layout(btn, room)
    for width in list(range(1, 64)) + list(range(64, room, 16)):
      host.resize(width, 40)
      seen = []
      for _ in range(6):
        app.processEvents()
        seen.append(btn.text())
      assert "" not in seen and len(set(seen)) == 1, (
        "%s at %d px: text over six passes %r" % (cls.__name__, width, seen))
    del host
  # An empty full text stays empty: an icon-only button must not grow an
  # ellipsis, and the ellipsis it grew looped the same way (no text, no
  # margin; an ellipsis, a margin) in a band just above its icon's width.
  from qttbx.widgets.chat.eliding import ElidingToolButton
  for cls in (ElidingCheckBox, ElidingToolButton):
    btn = cls()
    btn.set_full_text("")
    host = _in_layout(btn, 60)
    for width in range(1, 40):
      host.resize(width, 40)
      seen = []
      for _ in range(6):
        app.processEvents()
        seen.append(btn.text())
      assert seen == [""] * 6, (
        "%s with no text at %d px: %r" % (cls.__name__, width, seen))
    del host


def exercise_first_text_into_a_shown_empty_widget_shows_whole():
  """A widget filled after it is shown shows a text that fits whole.

  Qt gives a label with no text a hint of its own, ten px of nothing,
  which the chrome measure took for chrome: the first text set into a
  shown empty label whose width was up to those ten px more than the
  text needs showed "phenix.refi…" with a tooltip until the next resize.
  ChatTopBar's title and model labels start empty and are filled while
  the bar is shown. A hugging widget is laid out again for its new
  hint, so it is read after that pass; a stretched label gets no
  resize at all, and showed what the first elide left.
  """
  from qttbx.widgets.chat.eliding import (
    ElidingCheckBox, ElidingLabel, ElidingRadioButton, ElidingToolButton)
  app = _qapp()
  text = "phenix.refine"
  for cls in (ElidingLabel, ElidingToolButton, ElidingCheckBox,
              ElidingRadioButton):
    probe = cls()
    probe.set_full_text(text)
    room = _room_for(probe)
    for width in range(room, room + 12):
      w = cls()
      host = _in_layout(w, width)
      w.set_full_text(text)
      app.processEvents()
      assert w.text() == text, (cls.__name__, width, w.text())
      assert w.toolTip() == "", (cls.__name__, width, w.toolTip())
      del host


def exercise_font_change_re_elides():
  """A font change re-elides the text for the new metrics.

  The widgets re-elided on a resize, a show and a new text, and a font
  change is none of those when the layout leaves the widget's size
  alone: the copy cut for the old font stayed. Grown, it was wider than
  the widget and painted clipped, the ellipsis the part that fell off,
  so the text read as cut mid-word; shrunk, it showed less than fits.
  """
  from qttbx.widgets.chat.eliding import ElidingLabel, ElidingToolButton
  app = _qapp()
  saved = QtGui.QFont(app.font())
  try:
    for cls in (ElidingLabel, ElidingToolButton):
      w = cls()
      w.set_full_text(LONG)
      host = _in_layout(w, _too_narrow_for(w))
      # Pinned, so that neither font's height resizes the widget: the
      # hook is then the only way to a re-elide.
      w.setFixedHeight(w.height())
      app.processEvents()
      before = w.text()
      for factor in (1.6, 0.6):
        font = QtGui.QFont(saved)
        if saved.pointSizeF() > 0:
          font.setPointSizeF(saved.pointSizeF() * factor)
        else:                              # a pixel-sized application font
          font.setPixelSize(int(saved.pixelSize() * factor))
        app.setFont(font)
        app.processEvents()
        assert w.text() == w._elided_text(), (
          cls.__name__, factor, w.text(), w._elided_text())
        assert w.text() != before, (cls.__name__, factor, w.text())
        app.setFont(saved)
        app.processEvents()
      del host
  finally:
    app.setFont(saved)


def exercise_full_text_is_what_qt_holds():
  """A string the binding rewrites on its way into Qt is stored as Qt
  holds it, so every measure compares like with like.

  A lone surrogate -- what a path with a byte that is not UTF-8 decodes
  to -- reached Qt as a glyph box 9 px wider than the "?" it came back
  as, and a disclosure header quoting such a path flipped between whole
  and elided on every layout pass, at 100% CPU, in a band of host widths
  a few px wide. A leading byte order mark is dropped on the way in, so
  a label showing the text whole still measured it as not fitting at
  one width and sprouted a tooltip, and a reversed one has the whole
  string read with its bytes swapped; a NUL cut the string short, and
  "phenix" stood for "phenix<NUL>start_job", as a string terminator
  does. A tab is painted at a stop the one-line measure does not know,
  and a line break starts a line it does not see: a header with a tab
  in its status measured 22 px of chrome around nothing.
  """
  from qttbx.widgets.chat.eliding import ElidingLabel, ElidingToolButton
  app = _qapp()
  held_as = (
    ("wrote /data/caf\udce9.pdb", "wrote /data/caf\ufffd.pdb"),
    ("\ufeffphenix_start_job", "phenix_start_job"),
    ("\ufeff\ufeffphenix_start_job", "phenix_start_job"),
    ("\ufffephenix_start_job", "\ufffdphenix_start_job"),
    ("\ud83d\ude00 smile", "\U0001f600 smile"),
    ("phenix\x00start_job", "phenix\ufffdstart_job"),
    ("grep \u009c(finished)", "grep \ufffd(finished)"))
  for cls in (ElidingLabel, ElidingToolButton):
    for text, held in held_as:
      w = cls()
      w.set_full_text(text)
      host = _in_layout(w, _room_for(w))
      assert w.full_text() == held, (cls.__name__, repr(w.full_text()))
      assert w.text() == held, (cls.__name__, repr(w.text()))
      assert w.toolTip() == "" == w.accessibleDescription(), (
        cls.__name__, repr(text), w.toolTip(), w.accessibleDescription())
      # A caller's tooltip the same way, so a clamped status with a NUL
      # in it keeps its tooltip whole.
      w.set_full_text("short", tooltip=text)
      assert w.toolTip() == held, (cls.__name__, repr(w.toolTip()))
      del host
  # A tab and every kind of line break become a space on the one line
  # the widget paints; a caller's tooltip, which has lines, keeps a tab
  # and shows each kind of break as the one break it is.
  one_line = (("grep\t(finished)", "grep\t(finished)"),
              ("grep\n(finished)", "grep\n(finished)"),
              ("grep\r\n(finished)", "grep\n(finished)"),
              ("grep\u2028(finished)", "grep\n(finished)"),
              ("grep\u2029(finished)", "grep\n(finished)"))
  for cls in (ElidingLabel, ElidingToolButton):
    for text, tip in one_line:
      w = cls()
      w.set_full_text(text)
      host = _in_layout(w, _room_for(w))
      assert w.full_text() == "grep (finished)", (cls.__name__, repr(text))
      assert w.text() == "grep (finished)", (cls.__name__, repr(w.text()))
      if cls is ElidingLabel:
        assert w._chrome_width() == 0, (repr(text), w._chrome_width())
      w.set_full_text("short", tooltip=text)
      assert w.toolTip() == tip, (cls.__name__, repr(w.toolTip()))
      assert w.accessibleDescription() == tip, (
        cls.__name__, repr(w.accessibleDescription()))
      del host
  # And a header built at any width around the fit boundary settles.
  # The band is where the whole text's hint and the elided copy's meet,
  # and a widget built in it never settled: each width gets a fresh one.
  text = "\u25b8 grep (wrote /data/caf\udce9.pdb)"
  probe = ElidingToolButton()
  probe.set_full_text(text)
  room = _room_for(probe)
  for width in range(max(1, room - 40), room + 8):
    btn = ElidingToolButton()
    btn.set_full_text(text)
    host = _in_layout(btn, width)
    seen = []
    for _ in range(6):
      app.processEvents()
      seen.append(btn.text())
    assert len(set(seen)) == 1, (
      "at %d px: text over six passes %r" % (width, seen))
    del host


def exercise_plain_text_tooltip_shows_markup_literally():
  """A tooltip auto-detects rich text; a plain-text widget's must not.

  The approval card's title is PlainText because the tool name is
  externally supplied, so markup in it has to stay literal on hover too;
  a button's text is always plain. Ordinary text is left alone, so the
  tooltip of a path still compares equal to the path. The caller's own
  tooltip goes the same way as the default: the card
  puts a registered tool name there, and a tool disclosure the whole of
  a status it clamped -- a traceback, whose indentation and line breaks
  have to come through the escape with it.
  """
  from qttbx.widgets.chat.eliding import ElidingLabel, ElidingToolButton
  _qapp()
  marked = "<b>phenix_start_job</b> " + LONG
  # What Qt would render, left to decide: a <module> is no tag it
  # knows, but an &lt; is as good as one.
  traceback = ('Traceback (most recent call last):\n'
               '  File "run.py", line 1, in <module>\n'
               '    raise ValueError("bad")\nValueError: bad')
  # (text, what the tooltip shows): a CR LF pair from Windows is one
  # break, not two, and a lone CR -- after an LF, or a progress bar's --
  # is the space a plain tooltip paints, where the rich-text engine made
  # a break of each.
  tracebacks = (
    (traceback, traceback),
    ("bad: &lt;unset&gt;", "bad: &lt;unset&gt;"),
    (traceback.replace("\n", "\r\n"), traceback),
    (traceback.replace("\n", "\n\r"), traceback.replace("\n", "\n ")),
    ("done 50%\rdone 100% <b>", "done 50% done 100% <b>"),
    ("done\u2029again <b>", "done\nagain <b>"),
    ("done\u2028again <b>", "done\nagain <b>"))

  def literal(w, tip, what, text=marked, shows=marked):
    assert tip != text, (what, tip)
    assert QtGui.Qt.mightBeRichText(tip), (what, tip)
    shown = QtGui.QTextDocumentFragment.fromHtml(tip).toPlainText()
    assert shown == shows, (what, shown)
    # A screen reader gets the text, not the markup it was escaped into:
    # Qt describes a widget by its tooltip when it has no description.
    assert w.accessibleDescription() == shows, (
      what, w.accessibleDescription())

  for cls in (ElidingLabel, ElidingToolButton):
    w = cls()
    w.set_full_text(marked)
    host = _in_layout(w, _too_narrow_for(w))
    literal(w, w.toolTip(), cls.__name__)
    w.set_full_text("short", tooltip=marked)
    literal(w, w.toolTip(), cls.__name__ + " override")
    for text, shows in tracebacks:
      w.set_full_text("short", tooltip=text)
      literal(w, w.toolTip(), cls.__name__ + " override", text, shows)
    del host
  # A CR between the "&" and the "lt;" is painted as a space: described
  # with the CR dropped, it read as "x < y".
  label = ElidingLabel()
  label.set_full_text("short", tooltip="x &\rlt; y")
  host = _in_layout(label, _too_narrow_for(label))
  literal(label, label.toolTip(), "CR", "x &\rlt; y", "x & lt; y")
  del host


def exercise_accessible_name_is_the_full_text():
  """A screen reader gets the whole string, not the elided copy.

  Qt falls back to the rendered text for a widget with no accessible name,
  which for the approval card's checkbox in a 400 px view is "Alway…".
  """
  from qttbx.widgets.chat.eliding import ElidingCheckBox, ElidingLabel
  _qapp()
  for cls in (ElidingCheckBox, ElidingLabel):
    w = cls()
    w.set_full_text(LONG)
    host = _in_layout(w, _too_narrow_for(w))
    assert w.text() != LONG, "%s did not elide" % cls.__name__
    assert w.accessibleName() == LONG, (cls.__name__, w.accessibleName())
    # The description follows the tooltip, as plain text: Qt falls back
    # to the tooltip for it, and the tooltip of an "R&D ..." is escaped.
    w.set_full_text("R&D " + LONG)
    assert w.accessibleDescription() == "R&D " + LONG, (
      cls.__name__, w.accessibleDescription())
    assert "&amp;" not in w.accessibleDescription(), (
      cls.__name__, w.accessibleDescription())
    w.set_full_text("short", tooltip="R&D override")
    assert w.accessibleDescription() == "R&D override", (
      cls.__name__, w.accessibleDescription())
    w.set_full_text("short")
    assert w.toolTip() == w.accessibleDescription() == "", (
      cls.__name__, w.toolTip(), w.accessibleDescription())
    del host
  # None for a plain text, at any width: Qt reads the tooltip for it,
  # which is the text while a cut copy is shown, and each write would
  # be a DescriptionChanged event to a live screen reader.
  label = ElidingLabel()
  calls = []
  described = label.setAccessibleDescription
  label.setAccessibleDescription = lambda v: (calls.append(v), described(v))
  label.set_full_text(LONG)
  host = _in_layout(label, _too_narrow_for(label))
  for step in range(1, 21):
    host.resize(_too_narrow_for(label) - 4 * step, 40)
    _qapp().processEvents()
  assert label.text() != LONG, label.text()
  assert calls == [], calls
  del host
  # What the label shows as it is, the reader gets as it is, markup
  # and all: read as rendered, "line one <b>bold</b>" was "line one
  # bold" and a "<module>" nothing.
  label = ElidingLabel()
  for literal in ("[<b>Plan</b>] step", "line one <b>bold</b>", "x <y",
                  "in <module>", "a &amp; b"):
    label.set_full_text(literal)
    assert label.accessibleName() == literal, label.accessibleName()


def exercise_kept_characters_floor_the_label():
  """A label given characters to keep floors its container by their width.

  The approval card's "<source> requests:" lead keeps its kind that way:
  with nothing kept, a long tool name beside it took it to a lone
  ellipsis. The floor is the kept head (ElideRight) or tail (ElideLeft)
  and an ellipsis, and never more than the whole text takes.
  """
  from qttbx.widgets.chat.eliding import ElidingLabel
  _qapp()
  text = "mcp:phenix requests:"
  for mode, kept in ((QtCore.Qt.ElideRight, "mcp:phe"),
                     (QtCore.Qt.ElideLeft, "quests:")):
    label = ElidingLabel(mode=mode, keep=7)
    label.set_full_text(text)
    floor = label.minimumSizeHint().width()
    assert 0 < floor < label.sizeHint().width(), (floor, label.sizeHint())
    host = _in_layout(label, 1)
    assert label.width() >= floor, (label.width(), floor)
    shown = label.text()
    assert shown != text and kept in shown and "\u2026" in shown, shown
    assert label.toolTip() == text, label.toolTip()
    # Grown back, it shows the whole text as any eliding label does.
    host.resize(_room_for(label), 40)
    _qapp().processEvents()
    assert label.text() == text, label.text()
    del host
  # ElideMiddle keeps as many characters either side of the ellipsis,
  # so an odd count is kept only at the width of the next even one.
  for keep in (1, 3, 7):
    label = ElidingLabel(mode=QtCore.Qt.ElideMiddle, keep=keep)
    label.set_full_text(text)
    host = _in_layout(label, 1)
    shown = label.text()
    assert "\u2026" in shown and len(shown) - 1 >= keep, (keep, shown)
    del host
  short = ElidingLabel(keep=7)
  short.set_full_text("mcp:a")
  assert short.minimumSizeHint().width() == short.sizeHint().width()
  # Nor past the whole text when the ellipsis is wider than what it
  # stands in for: in a proportional font seven i's and an ellipsis
  # outmeasure eight i's. In a fixed-pitch one they measure the same,
  # so there the clamp has nothing to do and this cannot see it.
  narrow = ElidingLabel(keep=7)
  narrow.set_full_text("i" * 8)
  fm = narrow.fontMetrics()
  kept = fm.size(QtCore.Qt.TextSingleLine, "i" * 7 + "\u2026").width()
  whole = fm.size(QtCore.Qt.TextSingleLine, "i" * 8).width()
  assert kept >= whole, (kept, whole)
  if kept > whole:
    assert narrow.minimumSizeHint().width() == narrow.sizeHint().width(), (
      narrow.minimumSizeHint(), narrow.sizeHint())
  none = ElidingLabel()
  none.set_full_text(text)
  assert none.minimumSizeHint().width() == 0, none.minimumSizeHint()


def exercise_description_is_held_as_written():
  """A description is stored as Qt holds it, so a resize finds it
  unchanged and writes nothing.

  A leading byte order mark Qt drops: held with it, every resize wrote
  the description again."""
  from qttbx.widgets.chat.eliding import ElidingLabel
  app = _qapp()
  label = ElidingLabel()
  label.set_full_text("short", tooltip="\ufeffphenix")
  host = _in_layout(label, 300)
  calls = []
  described = label.setAccessibleDescription
  label.setAccessibleDescription = lambda v: (calls.append(v), described(v))
  for width in range(290, 270, -1):
    host.resize(width, 40)
    app.processEvents()
  assert calls == [], calls
  assert label.accessibleDescription() == "phenix", (
    repr(label.accessibleDescription()))
  del host


def exercise_a_copy_is_cut_while_the_whole_text_does_not_fit():
  """Short of the room its whole text takes, a widget shows a cut copy.

  ``elidedText`` keeps a string whole from the budget of its advance,
  and an italic or bold glyph overhangs that: at those widths the whole
  text was shown with its last glyph clipped, and no tooltip.
  """
  from qttbx.widgets.chat.eliding import ElidingToolButton
  app = _qapp()
  text = "grep (finished: wrote /data/out.pdb and /data/out.mtz)"
  font = QtGui.QFont(app.font())
  font.setItalic(True)
  w = ElidingToolButton()
  w.setFont(font)
  w.set_full_text(text)
  fm = w.fontMetrics()
  whole = fm.size(QtCore.Qt.TextSingleLine, text).width()
  advance = fm.horizontalAdvance(text)
  host = _in_layout(w, 1)
  chrome = w._chrome_width()
  for width in range(advance + chrome - 3, whole + chrome + 4):
    host.resize(width, 40)
    app.processEvents()
    available = w.width() - w._chrome_width()
    assert (w.text() == text) == (available >= whole), (
      width, available, whole, advance, w.text())
    assert (w.toolTip() == "") == (w.text() == text), (width, w.toolTip())
  del host


def exercise_label_shows_markup_literally():
  """A label shows markup as the text it is, elided like any text.

  Every label here shows plain or untrusted text -- a question card's
  header is a chip of a few words the model writes -- and QLabel's
  default decided per string to render a tag it found, images from
  local paths included, through machinery whose two copies, the markup
  rendered and a copy cut, never measured alike.
  """
  from qttbx.widgets.chat.eliding import ElidingLabel
  app = _qapp()
  marked = "<b>phenix_start_job</b> " + LONG
  label = ElidingLabel()
  label.set_full_text(marked)
  assert label.textFormat() == QtCore.Qt.PlainText, label.textFormat()
  assert label.accessibleName() == marked, label.accessibleName()
  host = _in_layout(label, _room_for(label))
  assert label.text() == marked, label.text()
  host.resize(_too_narrow_for(label), 40)
  app.processEvents()
  assert label.text().startswith("<b>phenix_start_job</b> "), label.text()
  assert label.text().endswith("\u2026"), label.text()
  shown = QtGui.QTextDocumentFragment.fromHtml(label.toolTip()).toPlainText()
  assert shown == marked, (label.toolTip(), shown)
  assert label.accessibleDescription() == marked, (
    label.accessibleDescription())
  del host


def exercise_the_cards_checkbox_describes_itself_once_after_a_font_change():
  """The approval card's checkbox describes itself once for a font change.

  The application's font reaches the widgets in no particular order, so
  the checkbox elided at its old width, said its whole text, and then
  said nothing once the layout had grown it for the buttons' new font:
  two descriptions to a screen reader, at the widths around the fit,
  for a change of half again as much on Qt 5 and of two and a half on
  Qt 6. A probe at twice the size, in a nest of its own, missed both.
  """
  from qttbx.widgets.chat.agent.tools import ToolApprovalRequest
  from qttbx.widgets.chat.eliding import ElidingCheckBox
  from qttbx.widgets.chat.tool_approval import ToolApprovalCard
  app = _qapp()
  request = ToolApprovalRequest(
    request_id="r1", tool_name="mcp__coot__close_model",
    tool_source="mcp:coot", input={"model": 1}, risk="write")
  saved = QtGui.QFont(app.font())
  twice = []
  for factor in (1.5, 2.0, 2.5):
    bigger = QtGui.QFont(saved)
    bigger.setPointSizeF(factor * saved.pointSizeF())
    for width in range(80, 400, 4):
      win = QtWidgets.QWidget()
      outer = QtWidgets.QVBoxLayout(win)
      card = ToolApprovalCard(win)
      card.set_requests([request])
      outer.addWidget(card)
      outer.addStretch(1)
      w = card.findChildren(ElidingCheckBox)[0]
      win.resize(width, 400)
      win.show()
      app.processEvents()
      calls = []
      described = w.setAccessibleDescription
      w.setAccessibleDescription = lambda v, c=calls, d=described: (
        c.append(v), d(v))
      try:
        app.setFont(bigger)
        for _ in range(5):
          app.processEvents()
        writes = list(calls)
      finally:
        app.setFont(saved)
        app.processEvents()
      if len(writes) > 1:
        twice.append((factor, width, writes))
      win.close()
  assert twice == [], twice[:3]


def exercise_a_header_describes_itself_once_for_a_status_change():
  """A tool header describes itself once when its status changes.

  Its text and its font change together -- the running style is the
  italic -- and the layouts above it grow it in turn, each resize an
  elide: a header that said its whole text at one width and nothing
  at the next sent a screen reader two descriptions.
  """
  from qttbx.widgets.chat.tool_call_disclosure import ToolCallDisclosure
  app = _qapp()
  twice = []
  for width in range(120, 700, 4):
    win = QtWidgets.QWidget()
    outer = QtWidgets.QVBoxLayout(win)
    bubble = QtWidgets.QFrame(win)
    inner = QtWidgets.QWidget(bubble)
    QtWidgets.QVBoxLayout(bubble).addWidget(inner)
    outer.addWidget(bubble)
    outer.addStretch(1)
    row = ToolCallDisclosure("mcp__phenix__phenix_get_log_summary",
                             "running", parent=inner)
    QtWidgets.QVBoxLayout(inner).addWidget(row)
    row.set_args({"job_id": "abc123", "lines": 40})
    win.resize(width, 300)
    win.show()
    app.processEvents()
    w = row.header_button
    calls = []
    described = w.setAccessibleDescription
    w.setAccessibleDescription = lambda v, c=calls, d=described: (
      c.append(v), d(v))
    row.set_status("finished (2.3 s)", color="default")
    for _ in range(5):
      app.processEvents()
    if len(calls) > 1:
      twice.append((width, calls))
    win.close()
  assert twice == [], twice[:3]


def exercise_the_same_text_again_changes_nothing():
  """Text and tooltip set again as they are leave the widget alone.

  A header refreshed with the status it already shows was elided and
  described over again for nothing.
  """
  from qttbx.widgets.chat.eliding import ElidingToolButton
  app = _qapp()
  w = ElidingToolButton()
  w.set_full_text("grep (running)", tooltip="grep over /data")
  host = _in_layout(w, _too_narrow_for(w))
  writes = []
  shown = w.setText
  w.setText = lambda v: (writes.append(v), shown(v))
  described = w.setAccessibleDescription
  w.setAccessibleDescription = lambda v: (writes.append(v), described(v))
  w.set_full_text("grep (running)", tooltip="grep over /data")
  app.processEvents()
  assert writes == [], writes
  del host


def exercise_a_font_change_on_a_hugging_widget_describes_itself_once():
  """A hugging widget whose font changes shows its whole text once the
  layout has grown it, and describes itself not at all for it.

  A tool header's stylesheet goes from italic to roman as it finishes;
  under a font whose italic is the narrower, the re-elide at the italic
  width wrote the whole text as the description, then nothing once the
  layout had grown the header. The application's font the same way.
  """
  from qttbx.widgets.chat.eliding import ElidingToolButton
  app = _qapp()
  text = "grep (finished: wrote /data/out.pdb and /data/out.mtz)"
  w = ElidingToolButton()
  w.set_full_text(text)
  host = _in_layout(w, 1)
  # Twice the font in force, so the change is one under the second pass
  # of this file as well, where the application font is already large.
  big = "font-size: %dpt" % (2 * app.font().pointSizeF())
  w.setStyleSheet(big)
  app.processEvents()
  room = _room_for(w)
  w.setStyleSheet("")
  app.processEvents()
  host.resize(room + 40, 40)
  app.processEvents()
  assert w.text() == text, w.text()
  calls = []
  described = w.setAccessibleDescription
  w.setAccessibleDescription = lambda v: (calls.append(v), described(v))
  w.setStyleSheet(big)
  for _ in range(3):
    app.processEvents()
  assert w.text() == text, (w.text(), w.width(), w.sizeHint().width())
  assert calls == [], calls
  w.setStyleSheet("")
  app.processEvents()
  assert w.text() == text and calls == [], (w.text(), calls)
  saved = QtGui.QFont(app.font())
  bigger = QtGui.QFont(saved)
  bigger.setPointSizeF(2 * saved.pointSizeF())
  try:
    app.setFont(bigger)
    for _ in range(3):
      app.processEvents()
    assert w.text() == text, (w.text(), w.width(), w.sizeHint().width())
    assert calls == [], calls
  finally:
    app.setFont(saved)
  del host


def exercise_ampersand_in_a_button_is_painted_as_itself():
  """A button paints an ``&`` as a mnemonic -- nothing, and the next
  character underlined -- so one in its text is doubled for the painter,
  and the button is sized for it; a label paints it as it is.

  A status quoting a URL's query string lost each ``&`` and underlined
  the letter after it.
  """
  from qttbx.widgets.chat.eliding import (
    ElidingCheckBox, ElidingLabel, ElidingToolButton)
  app = _qapp()
  text = "GET /x?a=1&b=2 failed"
  for cls in (ElidingToolButton, ElidingCheckBox):
    w = cls()
    w.set_full_text(text)
    plain = w.__class__.__mro__[2]()        # the Qt base
    plain.setText(text.replace("&", "&&"))
    assert w.sizeHint().width() >= plain.sizeHint().width(), (
      cls.__name__, w.sizeHint().width(), plain.sizeHint().width())
    host = _in_layout(w, _room_for(w))
    assert w.text() == text.replace("&", "&&"), (cls.__name__, w.text())
    assert w.toolTip() == "", (cls.__name__, w.toolTip())
    assert w.accessibleName() == text, (cls.__name__, w.accessibleName())
    host.resize(_too_narrow_for(w), 40)
    app.processEvents()
    assert "\u2026" in w.text(), (cls.__name__, w.text())
    assert "&" not in w.text().replace("&&", ""), (cls.__name__, w.text())
    # The tooltip is escaped for its "&", and shows the text as it is.
    shown = QtGui.QTextDocumentFragment.fromHtml(w.toolTip()).toPlainText()
    assert shown == text == w.accessibleDescription(), (
      cls.__name__, shown, w.accessibleDescription())
    del host
  label = ElidingLabel()
  label.set_full_text(text)
  host = _in_layout(label, _room_for(label))
  assert label.text() == text, label.text()
  del host


def exercise_the_description_shadows_an_escaped_tooltip_only():
  """A widget's description is the tooltip it was given, or its text
  when a tooltip of that text would be escaped, and otherwise nothing:
  written when the text is set, whatever width it then gets, and
  never for a resize or a font change.

  Qt reads the tooltip as the description when there is none, which
  for plain text is the whole text while a cut copy is shown and
  nothing while it fits: right as it is, and a description of the
  text as well had a fitting widget's name read out twice. A tooltip
  with a tag or an ampersand in it is escaped, and read out as the
  markup it was escaped into, so that text is the description.
  Written from the elide instead, the description followed the width,
  and the card's checkbox and a finishing tool header sent a screen
  reader two for one change.
  """
  from qttbx.widgets.chat.eliding import ElidingLabel, ElidingToolButton
  app = _qapp()
  plain = "grep (finished: wrote /data/out.pdb and /data/out.mtz)"
  marked = "grep (finished: wrote <b>/data/out.pdb</b> & /data/out.mtz)"
  for cls in (ElidingLabel, ElidingToolButton):
    for text, carried in ((plain, []), (marked, [marked])):
      w = cls()
      calls = []
      described = w.setAccessibleDescription
      w.setAccessibleDescription = lambda v, c=calls, d=described: (
        c.append(v), d(v))
      w.set_full_text(text)
      assert calls == carried, (cls.__name__, text, calls)
      host = _in_layout(w, _room_for(w))
      # A button doubles an ampersand for its painter.
      shown = w.text().replace("&&", "&")
      assert shown == text and w.toolTip() == "", (
        cls.__name__, shown, w.toolTip())
      host.resize(_too_narrow_for(w), 40)
      app.processEvents()
      shown = w.text().replace("&&", "&")
      assert shown != text and w.toolTip() != "", (
        cls.__name__, shown, w.toolTip())
      saved = QtGui.QFont(app.font())
      bigger = QtGui.QFont(saved)
      bigger.setPointSizeF(2 * saved.pointSizeF())
      try:
        app.setFont(bigger)
        for _ in range(3):
          app.processEvents()
      finally:
        app.setFont(saved)
        app.processEvents()
      assert calls == carried, (cls.__name__, text, calls)
      w.set_full_text("short", tooltip="the whole story")
      assert calls == carried + ["the whole story"], (cls.__name__, calls)
      # A new text under the same tooltip has nothing to add.
      w.set_full_text("shorter", tooltip="the whole story")
      assert calls == carried + ["the whole story"], (cls.__name__, calls)
      del host


def exercise_text_growth_on_a_hugging_widget_writes_no_description():
  """A hugging widget handed longer plain text shows it whole once the
  layout has grown it, and writes no description for it.

  Elided against the width it still had, a tool header that finished
  wrote its whole text as its description, then nothing once the layout
  had grown it: two changes to a screen reader for one.
  """
  from qttbx.widgets.chat.eliding import ElidingToolButton
  app = _qapp()
  text = "grep (finished: wrote /data/out.pdb and /data/out.mtz)"
  w = ElidingToolButton()
  w.set_full_text(text)
  room = _room_for(w)
  w.set_full_text("grep (running)")
  host = _in_layout(w, room + 40)
  calls = []
  described = w.setAccessibleDescription
  w.setAccessibleDescription = lambda v: (calls.append(v), described(v))
  w.set_full_text(text)
  for _ in range(3):
    app.processEvents()
  assert w.text() == text, (w.text(), w.width(), w.sizeHint().width(),
                            w.font().pointSizeF(), app.font().pointSizeF())
  assert calls == [], calls
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
  exercise_wrapping_label_shows_markup_literally,
  exercise_eliding_button_never_renders_empty_and_settles,
  exercise_first_text_into_a_shown_empty_widget_shows_whole,
  exercise_font_change_re_elides,
  exercise_full_text_is_what_qt_holds,
  exercise_description_is_held_as_written,
  exercise_the_description_shadows_an_escaped_tooltip_only,
  exercise_ampersand_in_a_button_is_painted_as_itself,
  exercise_text_growth_on_a_hugging_widget_writes_no_description,
  exercise_the_same_text_again_changes_nothing,
  exercise_the_cards_checkbox_describes_itself_once_after_a_font_change,
  exercise_a_header_describes_itself_once_for_a_status_change,
  exercise_label_shows_markup_literally,
  exercise_a_copy_is_cut_while_the_whole_text_does_not_fit,
  exercise_a_font_change_on_a_hugging_widget_describes_itself_once,
  exercise_plain_text_tooltip_shows_markup_literally,
  exercise_accessible_name_is_the_full_text,
  exercise_kept_characters_floor_the_label,
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
