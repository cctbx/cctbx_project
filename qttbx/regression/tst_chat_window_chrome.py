"""Window chrome widgets: ChatTopBar (the thin title / model / debug
strip above the conversation) and EdgeRail (the thin toggle columns on
either side of the splitter). Both are chrome, both are tiny."""

import os

from libtbx.utils import format_cpu_times

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


def _qapp():
  from qttbx.qt import QtWidgets
  from qttbx.widgets.font_init import init_default_app_font
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
  init_default_app_font(app)
  return app


# ---- ChatTopBar ----------------------------------------------------------


def exercise_top_bar_shows_title_and_model_and_hides_debug_by_default():
  _qapp()
  from qttbx.widgets.chat.chat_top_bar import ChatTopBar
  bar = ChatTopBar()
  bar.set_title("test_profile / first user message")
  bar.set_model("claude-opus-4-7")
  assert "test_profile" in bar.title_label.text()
  assert "claude-opus-4-7" in bar.model_label.text()
  # No debug-log set => the slot is hidden so the bar doesn't show empty space.
  # (The old `or text() == ""` branch was always true on a fresh bar, so a
  # regression that showed the empty slot would have passed unnoticed.)
  assert bar.debug_label.isHidden()


# The title is '<profile> / <first 60 chars of first user text>', so an
# ordinary conversation produces one this long.
LONG_TITLE = ("phenix_expert / Use the sequence to get a model to phase the "
              "data. Split the")


def _grow_until_title_fits(app, bar, cap=4000):
  """Resize ``bar`` to the first width at which the title renders in full.

  No pixel constant works here: the title's rendered width is
  font-dependent (Linux CI fonts run past 900 px where macOS stays under
  700), and the bar's QHBoxLayout hands each slot its stretch share of
  the width, size hints notwithstanding, so the title only renders in
  full once its share -- not the bar -- clears the text width. Search in
  small steps so the found width barely fits, which keeps the bar inside
  the band where a later sibling can still squeeze the title.
  """
  w = bar.sizeHint().width()
  while w < cap:
    bar.resize(w, 28)
    app.processEvents()
    if bar.title_label.text() == bar.title_label.full_text():
      return w
    w += 25
  raise AssertionError("no width under %d px renders the title in full" % cap)


def exercise_top_bar_long_title_does_not_floor_the_bar_width():
  """A long title must not pin the top bar's minimum width.

  title_label is a plain QLabel: with word wrap off, its minimumSizeHint is
  the FULL text width. The bar sits in the centre column, so that floors the
  column -- and hence the whole ChatWindow -- at the natural width of the
  title, and the window silently refuses to be resized any narrower.
  """
  _qapp()
  from qttbx.widgets.chat.chat_top_bar import ChatTopBar
  bar = ChatTopBar()
  bar.set_title(LONG_TITLE)
  bar.set_model("fable")
  assert bar.minimumSizeHint().width() < 400, (
    "top bar minimum width %d is floored by the title"
    % bar.minimumSizeHint().width())


def exercise_top_bar_elides_long_title_and_keeps_the_full_text():
  """Once the title can no longer floor the width it must elide rather than
  clip, and the untruncated title stays reachable on the tooltip."""
  app = _qapp()
  from qttbx.widgets.chat.chat_top_bar import ChatTopBar
  bar = ChatTopBar()
  bar.set_title(LONG_TITLE)
  bar.set_model("fable")
  bar.resize(320, 28)
  bar.show()
  app.processEvents()
  assert bar.title_label.width() <= 320, bar.title_label.width()
  assert "…" in bar.title_label.text(), bar.title_label.text()
  assert bar.title_label.toolTip() == LONG_TITLE, bar.title_label.toolTip()
  # A title with room to render in full is left alone.
  _grow_until_title_fits(app, bar)
  assert bar.title_label.text() == LONG_TITLE, bar.title_label.text()


# A proxy / Bedrock / Fireworks model id, passed straight through from
# --model. Nothing bounds its length.
LONG_MODEL = "accounts/fireworks/models/llama-v3p1-405b-instruct"
LONG_DEBUG_PATH = ("/Users/someone/Library/Application Support/phenix/chat/"
                   "sessions/2026-07-17T09-14-22Z-a4f19c/chat_debug.log")


def exercise_top_bar_model_and_debug_labels_do_not_floor_the_bar_width():
  """The title is not the only floor in the bar.

  model_label carries the raw --model value and debug_label carries the
  debug log path; both are plain QLabels whose minimumSizeHint is their full
  text width, so either one pins the bar -- and hence the whole ChatWindow --
  exactly the way the title did. ChatWindow calls set_debug_log_path at
  construction whenever debug logging is on, so the debug floor is not an
  edge case.
  """
  _qapp()
  from qttbx.widgets.chat.chat_top_bar import ChatTopBar
  bar = ChatTopBar()
  bar.set_title(LONG_TITLE)
  bar.set_model(LONG_MODEL)
  bar.set_debug_log_path(LONG_DEBUG_PATH)
  assert bar.minimumSizeHint().width() < 400, (
    "top bar minimum width %d is floored by the model / debug labels"
    % bar.minimumSizeHint().width())


def exercise_top_bar_gives_the_title_priority_over_the_side_slots():
  """The title must not be the first thing to vanish.

  Once every slot can shrink, QHBoxLayout makes the stretch item absorb the
  whole deficit -- so without a share of its own the title collapses to zero
  width at any narrow bar while a debug path happily keeps 500 px of it. The
  title is the bar's primary content; the model id and the debug path are
  secondary chrome and should elide first.
  """
  app = _qapp()
  from qttbx.widgets.chat.chat_top_bar import ChatTopBar
  bar = ChatTopBar()
  bar.set_title(LONG_TITLE)
  bar.set_model(LONG_MODEL)
  bar.set_debug_log_path(LONG_DEBUG_PATH)
  bar.show()
  for width in (900, 700, 520):
    bar.resize(width, 28)
    app.processEvents()
    assert bar.title_label.width() > bar.model_label.width(), (
      "at a %d px bar the title has %d px and the model id has %d px"
      % (width, bar.title_label.width(), bar.model_label.width()))
    assert bar.title_label.width() >= width // 3, (
      "at a %d px bar the title is squeezed down to %d px"
      % (width, bar.title_label.width()))


def exercise_top_bar_pins_side_slot_text_to_the_right_edge():
  """A short model id / debug path must hug the bar's right edge.

  The side slots take a stretch share of surplus width (that share is what
  spreads the SHORTFALL on a narrow bar instead of collapsing the title), so
  on a wide bar each slot is wider than its text. A QLabel paints where its
  alignment says: with the default AlignLeft the model id floats mid-bar
  with the blank space on its right -- it used to sit flush against the
  right edge when the slots took no stretch (stretch=0 slots hug their
  text). The slots must right-align so the surplus lands on their LEFT.
  If the side slots ever go back to taking no stretch, slot width == text
  width and this alignment requirement is vacuous, not wrong.
  """
  app = _qapp()
  from qttbx.qt import QtCore
  from qttbx.widgets.chat.chat_top_bar import ChatTopBar
  bar = ChatTopBar()
  bar.set_title("profile / hi")
  bar.set_model("fable")
  bar.show()
  bar.resize(900, 28)
  app.processEvents()
  # Preconditions for the regression: the model slot is the rightmost
  # visible slot, ends at the bar's right content edge, and is wider than
  # its text (i.e. there IS surplus for alignment to place).
  assert bar.debug_label.isHidden()
  gap = bar.width() - 1 - bar.model_label.geometry().right()
  assert gap <= 10, (
    "model slot ends %d px short of the bar's right edge" % gap)
  fm = bar.model_label.fontMetrics()
  assert bar.model_label.width() > fm.size(
    QtCore.Qt.TextSingleLine, "fable").width(), \
    "precondition: the model slot should be wider than its text"
  assert bar.model_label.alignment() & QtCore.Qt.AlignRight, (
    "model id paints at the left of its %d px slot, leaving the blank "
    "space against the bar's right edge" % bar.model_label.width())
  # With the debug slot shown it becomes the rightmost slot and takes the
  # same surplus, so a short path has the same problem.
  bar.set_debug_log_path("/tmp/chat.log")
  app.processEvents()
  gap = bar.width() - 1 - bar.debug_label.geometry().right()
  assert gap <= 10, (
    "debug slot ends %d px short of the bar's right edge" % gap)
  assert bar.debug_label.alignment() & QtCore.Qt.AlignRight, (
    "debug path paints at the left of its %d px slot, leaving the blank "
    "space against the bar's right edge" % bar.debug_label.width())


def exercise_top_bar_re_elides_the_title_when_a_sibling_takes_its_width():
  """The title must re-elide when the layout -- not the window -- resizes it.

  title_label's width also changes on internal re-layout: set_debug_log_path
  showing the debug slot, or set_model with a longer name, shrinks the title
  with no resize of the bar itself. Hooking the elide on the bar's resizeEvent
  alone leaves a stale elide -- the title paints hard-clipped mid-character,
  with no ellipsis, and the tooltip is empty because the previous pass found
  the text fit, so the full title is unreachable even by hover.
  """
  app = _qapp()
  from qttbx.widgets.chat.chat_top_bar import ChatTopBar
  bar = ChatTopBar()
  bar.set_title(LONG_TITLE)
  bar.set_model("fable")
  bar.show()
  app.processEvents()
  # Wide enough that the title renders in full to begin with -- searched
  # for, since the needed width is font-dependent and a hardcoded 900 px
  # failed on Linux CI fonts -- and, by stopping at the first width that
  # fits, still narrow enough that the debug slot then pushes the title
  # under its natural width.
  _grow_until_title_fits(app, bar)
  assert bar.title_label.text() == LONG_TITLE, "precondition: title fits"
  # Show the debug slot: it takes width from the title, but the bar's own size
  # never changes, so no resizeEvent ever reaches the bar.
  bar.set_debug_log_path(LONG_DEBUG_PATH)
  app.processEvents()
  assert bar.title_label.text() != LONG_TITLE, \
    "precondition: the debug slot should have squeezed the title"
  label = bar.title_label
  fm = label.fontMetrics()
  assert fm.horizontalAdvance(label.text()) <= label.width(), (
    "title text is %d px wide in a %d px label -- stale elide, painting "
    "clipped mid-character"
    % (fm.horizontalAdvance(label.text()), label.width()))
  assert label.toolTip() == LONG_TITLE, (
    "full title unreachable after the re-layout; tooltip is %r"
    % label.toolTip())


def exercise_top_bar_shows_and_hides_debug_log_path():
  _qapp()
  from qttbx.widgets.chat.chat_top_bar import ChatTopBar
  bar = ChatTopBar()
  bar.set_debug_log_path("/Users/x/proj/.phenix_chat/logs/debug-20260518T103000.log")
  assert not bar.debug_label.isHidden()
  assert "debug-20260518T103000.log" in bar.debug_label.text()
  bar.set_debug_log_path(None)
  assert bar.debug_label.isHidden()


# ---- EdgeRail ------------------------------------------------------------


def exercise_left_rail_starts_collapsed_and_toggles():
  """A left-edge rail starts at 'collapsed' (right-arrow = 'show panel').
  A click emits toggled(True) and flips the arrow to left-arrow.
  Another click goes back to right-arrow and emits toggled(False)."""
  _qapp()
  from qttbx.widgets.chat.edge_rail import EdgeRail
  rail = EdgeRail(side="left", tooltip_show="Show conversations",
                  tooltip_hide="Hide conversations")
  states = []
  rail.toggled.connect(states.append)
  assert rail.button.text() == "▶"
  rail.button.click()
  assert states == [True]
  assert rail.button.text() == "◀"
  rail.button.click()
  assert states == [True, False]
  assert rail.button.text() == "▶"


def exercise_right_rail_arrow_orientation_inverted():
  """Right-edge rail uses the inverse arrows: left-arrow collapsed
  (point INTO the window to show), right-arrow expanded (point OUT
  to hide)."""
  _qapp()
  from qttbx.widgets.chat.edge_rail import EdgeRail
  rail = EdgeRail(side="right", tooltip_show="Show artifacts",
                  tooltip_hide="Hide artifacts")
  assert rail.button.text() == "◀"
  rail.button.click()
  assert rail.button.text() == "▶"


def exercise_rail_width_is_fixed_thin():
  """The rail width should NOT grow with content -- it's chrome, not
  data. 18 px matches the spec; +/- 4 px slack."""
  _qapp()
  from qttbx.widgets.chat.edge_rail import EdgeRail
  rail = EdgeRail(side="left", tooltip_show="x", tooltip_hide="y")
  assert 14 <= rail.maximumWidth() <= 22, rail.maximumWidth()


def exercise_set_expanded_is_source_of_truth_for_arrow():
  """set_expanded(bool) is the public mechanism ChatWindow uses to
  reconcile the arrow with whatever the splitter actually shows --
  e.g., after a programmatic setSizes() call, or after the user
  drags the splitter handle. Calling it must NOT emit toggled()."""
  _qapp()
  from qttbx.widgets.chat.edge_rail import EdgeRail
  rail = EdgeRail(side="left", tooltip_show="show", tooltip_hide="hide")
  states = []
  rail.toggled.connect(states.append)
  rail.set_expanded(True)
  assert rail.button.text() == "◀"
  rail.set_expanded(False)
  assert rail.button.text() == "▶"
  # set_expanded() must not re-emit -- otherwise ChatWindow's "call back
  # to sync the rail" pattern would feedback-loop.
  assert states == []


def exercise():
  exercise_top_bar_shows_title_and_model_and_hides_debug_by_default()
  exercise_top_bar_long_title_does_not_floor_the_bar_width()
  exercise_top_bar_elides_long_title_and_keeps_the_full_text()
  exercise_top_bar_model_and_debug_labels_do_not_floor_the_bar_width()
  exercise_top_bar_gives_the_title_priority_over_the_side_slots()
  exercise_top_bar_pins_side_slot_text_to_the_right_edge()
  exercise_top_bar_re_elides_the_title_when_a_sibling_takes_its_width()
  exercise_top_bar_shows_and_hides_debug_log_path()
  exercise_left_rail_starts_collapsed_and_toggles()
  exercise_right_rail_arrow_orientation_inverted()
  exercise_rail_width_is_fixed_thin()
  exercise_set_expanded_is_source_of_truth_for_arrow()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
