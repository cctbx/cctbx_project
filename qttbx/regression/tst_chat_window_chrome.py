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
  assert bar.debug_label.isHidden() or bar.debug_label.text() == ""


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
  exercise_top_bar_shows_and_hides_debug_log_path()
  exercise_left_rail_starts_collapsed_and_toggles()
  exercise_right_rail_arrow_orientation_inverted()
  exercise_rail_width_is_fixed_thin()
  exercise_set_expanded_is_source_of_truth_for_arrow()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
