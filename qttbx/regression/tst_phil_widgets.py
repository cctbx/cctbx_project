import os
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import sys
try:
  from PySide2.QtWidgets import QApplication
except ImportError:
  print("PySide2 not available; skipping")
  print("OK")
  sys.exit(0)

from libtbx.utils import format_cpu_times
from libtbx.test_utils import approx_equal

_app = None

def _get_app():
  global _app
  if _app is None:
    _app = QApplication.instance() or QApplication([])
  return _app

import libtbx.phil
from qttbx.phil import PhilModel

def _example_master():
  return libtbx.phil.parse("""
refinement {
  macro_cycles = 3
    .type = int(value_min=1, value_max=20)
  notes = ""
    .type = str
}
""")

def exercise_phil_model_get_definition():
  m = PhilModel()
  m.initialize_model(_example_master())
  d = m.get_definition("refinement.macro_cycles")
  assert d.is_definition
  assert d.name == "macro_cycles"
  assert d.type.phil_type == "int"
  assert d.type.value_min == 1
  assert d.type.value_max == 20
  print("exercise_phil_model_get_definition OK")

def exercise_phil_model_thread_assertion_passes_in_gui_thread():
  # Direct setData in the GUI thread should succeed (no AssertionError).
  m = PhilModel()
  m.initialize_model(_example_master())
  idx = m.index(0, 1)  # any valid index will do; setData asserts before doing work
  # We only care that the call doesn't AssertionError; result may be False.
  m.setData(idx, "ignored", role=2)  # Qt.EditRole == 2; avoids importing Qt here
  print("exercise_phil_model_thread_assertion_passes_in_gui_thread OK")

def exercise_phil_model_index_for_path():
  m = PhilModel()
  m.initialize_model(_example_master())
  idx = m.index_for_path("refinement.macro_cycles")
  assert idx.isValid()
  assert idx.column() == 1, "value column expected, got %d" % idx.column()
  # The PhilItem at this index should be the macro_cycles definition.
  item = idx.internalPointer()
  assert item.definition.name == "macro_cycles"

  # Bad path raises.
  raised = False
  try:
    m.index_for_path("refinement.no_such_thing")
  except ValueError:
    raised = True
  assert raised
  print("exercise_phil_model_index_for_path OK")

def exercise_phil_model_populate_after_factory_removal():
  m = PhilModel()
  m.initialize_model(_example_master())
  # Walk to refinement.macro_cycles and verify it exists with correct value.
  d = m.get_definition("refinement.macro_cycles")
  assert d.is_definition
  v = m.get_phil_extract_value("refinement.macro_cycles")
  assert v == 3
  print("exercise_phil_model_populate_after_factory_removal OK")

def exercise_phil_item_extract_path_simple():
  """PhilItem.extract_path returns dotted segments for non-multiple paths."""
  from qttbx.phil import PhilModel
  m = PhilModel()
  m.initialize_model(libtbx.phil.parse('''
refine {
  cycles = 3
    .type = int
}
'''.strip()))
  # Walk to the cycles definition's PhilItem.
  master = m._root.child(0)
  refine = master.child(0)
  cycles = refine.child(0)
  assert cycles.extract_path() == ["refine", "cycles"], cycles.extract_path()
  print("exercise_phil_item_extract_path_simple OK")


def exercise_phil_item_extract_path_with_instance_index():
  """PhilItem.extract_path emits (name, index) tuples for multi-scope segments."""
  from qttbx.phil import PhilItem
  # Build a synthetic PhilItem chain to verify extract_path's tuple
  # encoding. The actual multi-instance tree is built in Task 2.
  master = PhilItem()
  outer_scope = PhilItem(parent=master)
  class _Def(object):
    def __init__(self, name):
      self.name = name
  outer_scope.definition = _Def("ncs_group")
  outer_scope.instance_index = 2
  inner_def = PhilItem(parent=outer_scope)
  inner_def.definition = _Def("selection")
  master.appendChild(outer_scope)
  outer_scope.appendChild(inner_def)
  assert inner_def.extract_path() == [("ncs_group", 2), "selection"], \
    inner_def.extract_path()
  print("exercise_phil_item_extract_path_with_instance_index OK")


def exercise_phil_model_multi_instance_tree():
  """PhilModel builds N child items for a multi-scope with N extract instances."""
  from qttbx.phil import PhilModel
  m = PhilModel()
  master = libtbx.phil.parse('''
ncs_group {
  selection = "all"
    .type = str
}
'''.strip())
  master.objects[0].multiple = True
  m.initialize_model(master)
  # By default, an empty multi-scope has one default instance.
  master_item = m._root.child(0)
  assert master_item.childCount() == 1, master_item.childCount()
  ncs0 = master_item.child(0)
  assert ncs0.instance_index == 0, ncs0.instance_index
  print("exercise_phil_model_multi_instance_tree OK")


def exercise_phil_model_indexed_get_set():
  """get/set_phil_extract_value accept indexed paths and reach the right instance."""
  from qttbx.phil import PhilModel
  m = PhilModel()
  master = libtbx.phil.parse('''
ncs_group {
  selection = "all"
    .type = str
}
'''.strip())
  master.objects[0].multiple = True
  m.initialize_model(master)
  # Manually grow the extract to two instances.
  m._extract.ncs_group.append(master.objects[0].extract())
  # Indexed path access should reach instance 1.
  m.set_phil_extract_value([("ncs_group", 1), "selection"], "chain A")
  got = m.get_phil_extract_value([("ncs_group", 1), "selection"])
  assert got == "chain A", got
  # Instance 0 should still hold the default.
  got0 = m.get_phil_extract_value([("ncs_group", 0), "selection"])
  assert got0 == "all", got0
  print("exercise_phil_model_indexed_get_set OK")


def exercise_phil_model_persistent_index_simple():
  """persistent_index_for_path on a non-multi path returns a stable QPersistentModelIndex."""
  from qttbx.phil import PhilModel
  m = PhilModel()
  m.initialize_model(libtbx.phil.parse('''
refine {
  cycles = 3
    .type = int
}
'''.strip()))
  qpi = m.persistent_index_for_path("refine.cycles")
  # Duck-type via class name (no isinstance).
  assert type(qpi).__name__ == "QPersistentModelIndex"
  assert qpi.isValid()
  assert qpi.column() == 1
  print("exercise_phil_model_persistent_index_simple OK")


def exercise_phil_model_persistent_index_multi_scope():
  """persistent_index_for_path resolves a multi-scope segment via scope_indices."""
  from qttbx.phil import PhilModel
  m = PhilModel()
  master = libtbx.phil.parse('''
ncs_group {
  selection = "all"
    .type = str
}
'''.strip())
  master.objects[0].multiple = True
  m.initialize_model(master)
  # Grow to 3 instances.
  for _ in range(2):
    m._extract.ncs_group.append(master.objects[0].extract())
  # Re-populate so the model tree reflects the 3 instances.
  m.beginResetModel()
  m._root = type(m._root)()
  m._populate_tree(m._scope, m._root)
  m.endResetModel()
  qpi = m.persistent_index_for_path("ncs_group.selection",
                                     scope_indices=[2])
  assert qpi.isValid()
  print("exercise_phil_model_persistent_index_multi_scope OK")


def exercise_phil_model_persistent_index_requires_scope_indices():
  """persistent_index_for_path raises when scope_indices is missing on a multi-path."""
  from qttbx.phil import PhilModel
  m = PhilModel()
  master = libtbx.phil.parse('''
ncs_group {
  selection = "all"
    .type = str
}
'''.strip())
  master.objects[0].multiple = True
  m.initialize_model(master)
  try:
    m.persistent_index_for_path("ncs_group.selection")
  except ValueError as e:
    assert "scope_indices" in str(e), str(e)
    print("exercise_phil_model_persistent_index_requires_scope_indices OK")
    return
  raise AssertionError("expected ValueError")


def exercise_phil_model_add_scope_instance_grows_tree():
  """add_scope_instance grows the model tree and the extract list together."""
  from qttbx.phil import PhilModel
  m = PhilModel()
  master = libtbx.phil.parse('''
ncs_group {
  selection = "all"
    .type = str
}
'''.strip())
  master.objects[0].multiple = True
  m.initialize_model(master)
  master_item = m._root.child(0)
  assert master_item.childCount() == 1
  qpi = m.add_scope_instance("ncs_group")
  assert master_item.childCount() == 2, master_item.childCount()
  assert master_item.child(1).instance_index == 1
  assert qpi.isValid()
  assert len(m._extract.ncs_group) == 2
  print("exercise_phil_model_add_scope_instance_grows_tree OK")


def exercise_phil_model_remove_scope_instance_shrinks_tree():
  """remove_scope_instance shrinks both tree and extract list and re-numbers indices."""
  from qttbx.phil import PhilModel
  m = PhilModel()
  master = libtbx.phil.parse('''
ncs_group {
  selection = "all"
    .type = str
}
'''.strip())
  master.objects[0].multiple = True
  m.initialize_model(master)
  m.add_scope_instance("ncs_group")
  m.add_scope_instance("ncs_group")
  master_item = m._root.child(0)
  assert master_item.childCount() == 3
  # Set distinct values so re-numbering can be checked.
  m.set_phil_extract_value([("ncs_group", 0), "selection"], "A")
  m.set_phil_extract_value([("ncs_group", 1), "selection"], "B")
  m.set_phil_extract_value([("ncs_group", 2), "selection"], "C")
  m.remove_scope_instance("ncs_group", 1)
  assert master_item.childCount() == 2
  # The B instance is gone; A and C remain, re-numbered as 0 and 1.
  assert master_item.child(0).instance_index == 0
  assert master_item.child(1).instance_index == 1
  assert m.get_phil_extract_value([("ncs_group", 0), "selection"]) == "A"
  assert m.get_phil_extract_value([("ncs_group", 1), "selection"]) == "C"
  print("exercise_phil_model_remove_scope_instance_shrinks_tree OK")


def exercise_phil_model_add_scope_instance_treeview_sync():
  """An attached QTreeView reflects the new row after add_scope_instance."""
  from PySide2.QtWidgets import QTreeView
  from qttbx.phil import PhilModel
  _get_app()
  m = PhilModel()
  master = libtbx.phil.parse('''
ncs_group {
  selection = "all"
    .type = str
}
'''.strip())
  master.objects[0].multiple = True
  m.initialize_model(master)
  view = QTreeView()
  view.setModel(m)
  master_qmi = m.index(0, 0)
  assert m.rowCount(master_qmi) == 1
  m.add_scope_instance("ncs_group")
  assert m.rowCount(master_qmi) == 2
  print("exercise_phil_model_add_scope_instance_treeview_sync OK")


def exercise_phil_model_restore_defaults():
  """restore_defaults reverts edits and resets multi-scope instance counts."""
  from qttbx.phil import PhilModel
  m = PhilModel()
  master = libtbx.phil.parse('''
val = 5
  .type = int
ncs_group {
  selection = "all"
    .type = str
}
'''.strip())
  master.objects[1].multiple = True
  m.initialize_model(master)
  # Mutate.
  m.set_phil_extract_value("val", 99)
  m.add_scope_instance("ncs_group")
  assert m.get_phil_extract().val == 99
  assert len(m.get_phil_extract().ncs_group) == 2
  # Restore.
  m.restore_defaults()
  assert m.get_phil_extract().val == 5
  assert len(m.get_phil_extract().ncs_group) == 1
  print("exercise_phil_model_restore_defaults OK")


def exercise_phil_model_intrinsic_validation_int_bounds():
  """A row whose stored int is out of bounds reports a BackgroundRole pink + ToolTipRole error."""
  from PySide2.QtCore import Qt
  from qttbx.phil import PhilModel
  _get_app()
  m = PhilModel()
  m.initialize_model(libtbx.phil.parse('''
val = 5
  .type = int(value_min=1, value_max=10)
'''.strip()))
  idx = m.index_for_path("val")
  assert m.data(idx, Qt.BackgroundRole) is None
  assert m.data(idx, Qt.ToolTipRole) is None
  m.setData(idx, 99, Qt.EditRole)
  brush = m.data(idx, Qt.BackgroundRole)
  tip = m.data(idx, Qt.ToolTipRole)
  assert brush is not None
  assert "exceeds maximum" in tip, tip
  m.setData(idx, 5, Qt.EditRole)
  assert m.data(idx, Qt.BackgroundRole) is None
  assert m.data(idx, Qt.ToolTipRole) is None
  print("exercise_phil_model_intrinsic_validation_int_bounds OK")


def exercise_phil_model_intrinsic_validation_choice_membership():
  """A choice value not in the master's set is flagged."""
  from PySide2.QtCore import Qt
  from qttbx.phil import PhilModel
  _get_app()
  m = PhilModel()
  m.initialize_model(libtbx.phil.parse('''
mode = *fast slow exhaustive
  .type = choice
'''.strip()))
  idx = m.index_for_path("mode")
  m.setData(idx, "rocket", Qt.EditRole)
  tip = m.data(idx, Qt.ToolTipRole)
  assert tip is not None
  assert "unknown choice" in tip, tip
  print("exercise_phil_model_intrinsic_validation_choice_membership OK")


def exercise_phil_model_intrinsic_validation_at_initialize():
  """initialize_model runs the intrinsic check via _populate_tree.

  libtbx.phil rejects out-of-bounds defaults at parse/extract time, so this
  test directly mutates the item's stored value and re-runs the validator
  (the same path _populate_tree exercises) to confirm that initialize_model's
  per-leaf refresh would flag a bad value identically.
  """
  from PySide2.QtCore import Qt
  from qttbx.phil import PhilModel
  _get_app()
  m = PhilModel()
  m.initialize_model(libtbx.phil.parse('''
val = 5
  .type = int(value_min=1, value_max=10)
'''.strip()))
  idx = m.index_for_path("val")
  assert m.data(idx, Qt.BackgroundRole) is None
  # Drive the same code path _populate_tree uses on each leaf.
  item = idx.internalPointer()
  item.value = 99
  m._refresh_intrinsic_validity(item)
  assert m.data(idx, Qt.BackgroundRole) is not None
  assert "exceeds maximum" in m.data(idx, Qt.ToolTipRole)
  print("exercise_phil_model_intrinsic_validation_at_initialize OK")

def exercise_label_for_definition_helper():
  """label_for_definition returns short_caption when set, else a prettified name."""
  from qttbx.phil import label_for_definition
  scope = libtbx.phil.parse('''
atom_selection = "all"
  .type = str
labelled = 1
  .type = int
  .short_caption = "Custom Label"
'''.strip())
  d_no_caption = scope.objects[0]
  d_with_caption = scope.objects[1]
  assert label_for_definition(d_no_caption) == "Atom selection"
  assert label_for_definition(d_with_caption) == "Custom Label"
  print("exercise_label_for_definition_helper OK")


from qttbx.widgets.phil import PhilWidget

def _make_int_definition():
  scope = libtbx.phil.parse("""
n = 5
  .type = int(value_min=0, value_max=10)
""")
  return scope.objects[0]

def exercise_phil_widget_base_construct():
  d = _make_int_definition()
  w = PhilWidget(d)
  assert w.definition is d
  # Calling validate() on the abstract base hits the NotImplementedError in
  # value(), which populates _error and flips isValid() to False.
  w.validate()
  assert not w.isValid()
  print("exercise_phil_widget_base_construct OK")

def exercise_phil_widget_base_signals_exist():
  d = _make_int_definition()
  w = PhilWidget(d)
  # valueChanged and validityChanged are class attributes (PySide2 Signal)
  assert hasattr(type(w), "valueChanged")
  assert hasattr(type(w), "validityChanged")
  print("exercise_phil_widget_base_signals_exist OK")

from qttbx.widgets.phil.text_base import ValidatedLineEdit

def exercise_validated_line_edit_valid_to_invalid():
  def parse(s):
    n = int(s)              # raises ValueError on bad input
    if n < 0:
      raise ValueError("must be non-negative")
    return n
  le = ValidatedLineEdit(parse=parse)
  le.show()
  le.setText("3")
  le.validate()
  assert le.isValid()
  assert le.value() == 3
  le.setText("abc")
  le.validate()
  assert not le.isValid()
  assert "invalid literal" in le.errorString() or "abc" in le.errorString()
  le.setText("-1")
  le.validate()
  assert not le.isValid()
  assert "non-negative" in le.errorString()
  print("exercise_validated_line_edit_valid_to_invalid OK")

def exercise_validated_line_edit_emits_validity_changed():
  events = []
  def parse(s):
    return int(s)
  le = ValidatedLineEdit(parse=parse)
  le.validityChanged.connect(events.append)
  le.setText("5")
  le.validate()                # transitions to valid (was invalid before any text)
  le.setText("xx")
  le.validate()                # transitions to invalid
  assert events == [True, False], "events: " + repr(events)
  print("exercise_validated_line_edit_emits_validity_changed OK")

def exercise_validated_text_edit_valid_to_invalid():
  """ValidatedTextEdit transitions style and validity on bad input."""
  from qttbx.widgets.phil.text_base import ValidatedTextEdit
  _get_app()
  te = ValidatedTextEdit(parse=int)
  te.setPlainText("123")
  te.validate()
  assert te.isValid(), te.errorString()
  te.setPlainText("not-an-int")
  # textChanged auto-fires _on_text_changed -> validate
  assert not te.isValid()
  assert "invalid literal" in te.errorString() or "int" in te.errorString()
  print("exercise_validated_text_edit_valid_to_invalid OK")


def exercise_validated_text_edit_focus_out_commits():
  """focusOutEvent on a valid ValidatedTextEdit emits valueChanged."""
  from PySide2.QtGui import QFocusEvent
  from PySide2.QtCore import QEvent
  from qttbx.widgets.phil.text_base import ValidatedTextEdit
  _get_app()
  te = ValidatedTextEdit(parse=int)
  te.setPlainText("42")
  te.validate()
  emitted = []
  te.valueChanged.connect(lambda v: emitted.append(v))
  te.focusOutEvent(QFocusEvent(QEvent.FocusOut))
  assert emitted == [42], emitted
  print("exercise_validated_text_edit_focus_out_commits OK")

from qttbx.widgets.phil import (
  register_widget, widget_for_definition, PhilWidget)

class _DummyWidget(PhilWidget):
  pass

def exercise_registry_register_and_dispatch():
  register_widget("int", _DummyWidget)
  d = _make_int_definition()
  w = widget_for_definition(d)
  assert type(w) is _DummyWidget
  print("exercise_registry_register_and_dispatch OK")

def exercise_registry_unknown_type_raises():
  d = _make_int_definition()
  # Pretend it's an unknown type for this test only:
  saved = d.type.phil_type
  d.type.phil_type = "no_such_type_xyz"
  try:
    raised = False
    try:
      widget_for_definition(d)
    except ValueError:
      raised = True
    assert raised
  finally:
    d.type.phil_type = saved
  print("exercise_registry_unknown_type_raises OK")

def exercise_widget_for_definition_unregistered_raises_value_error():
  """Unregistered phil_type raises ValueError naming the type and the registered set."""
  from qttbx.widgets.phil import widget_for_definition
  _get_app()
  # Build a fake definition exposing an unregistered phil_type.
  class _FakeType:
    phil_type = "no_such_type_xyz"
    multi = False
  class _FakeDef:
    type = _FakeType()
    multiple = False
  try:
    widget_for_definition(_FakeDef())
  except ValueError as e:
    msg = str(e)
    assert "no_such_type_xyz" in msg
    assert "registered" in msg
    print("exercise_widget_for_definition_unregistered_raises_value_error OK")
    return
  raise AssertionError("expected ValueError")

def exercise_widget_for_definition_fallback_kwarg():
  """fallback= constructs the fallback class for unregistered phil_types."""
  from qttbx.widgets.phil import widget_for_definition
  from qttbx.widgets.phil.str_widget import StrWidget
  _get_app()
  class _FakeType:
    phil_type = "no_such_type_xyz"
    multi = False
  class _FakeDef:
    type = _FakeType()
    multiple = False
    def extract(self):
      return None
  w = widget_for_definition(_FakeDef(), fallback=StrWidget)
  assert type(w).__name__ == "StrWidget"
  print("exercise_widget_for_definition_fallback_kwarg OK")

from qttbx.widgets.phil.int_widget import IntWidget

def exercise_int_widget_round_trip():
  d = _make_int_definition()
  w = IntWidget(d)
  w.setValue(7)
  assert w.value() == 7
  assert w.isValid()
  w.setValue(None)
  assert w.value() is None
  print("exercise_int_widget_round_trip OK")

def exercise_int_widget_limits():
  d = _make_int_definition()             # value_min=0, value_max=10
  w = IntWidget(d)
  w.setValue(0)
  assert w.isValid()
  w.setValue(10)
  assert w.isValid()
  w._line_edit.setText("11")
  w._line_edit.validate()
  assert not w.isValid()
  assert "10" in w.errorString()
  w._line_edit.setText("-1")
  w._line_edit.validate()
  assert not w.isValid()
  assert "0" in w.errorString()
  print("exercise_int_widget_limits OK")

def exercise_int_widget_setvalue_idempotent():
  # setValue must NOT trigger valueChanged when the new value matches current.
  d = _make_int_definition()
  w = IntWidget(d)
  w.setValue(5)
  events = []
  w.valueChanged.connect(events.append)
  w.setValue(5)               # idempotent
  assert events == [], "spurious valueChanged: " + repr(events)
  print("exercise_int_widget_setvalue_idempotent OK")

from qttbx.widgets.phil.bool_widget import BoolWidget

def _make_bool_definition():
  return libtbx.phil.parse("flag = True\n  .type = bool").objects[0]

def exercise_bool_widget_round_trip():
  d = _make_bool_definition()
  w = BoolWidget(d)
  w.setValue(True)
  assert w.value() is True
  assert w.isValid()
  w.setValue(False)
  assert w.value() is False
  print("exercise_bool_widget_round_trip OK")

def exercise_bool_widget_tristate_for_none():
  d = _make_bool_definition()
  w = BoolWidget(d)
  w.setAllowNone(True)
  w.setValue(None)
  assert w.value() is None
  print("exercise_bool_widget_tristate_for_none OK")

from qttbx.widgets.phil.float_widget import FloatWidget

def _make_float_definition():
  return libtbx.phil.parse(
    "x = 1.5\n  .type = float(value_min=0.0, value_max=2.0)").objects[0]

def exercise_float_widget_round_trip():
  d = _make_float_definition()
  w = FloatWidget(d)
  w.setValue(1.25)
  assert approx_equal(w.value(), 1.25)
  assert w.isValid()
  print("exercise_float_widget_round_trip OK")

def exercise_float_widget_limits():
  d = _make_float_definition()           # value_min=0.0, value_max=2.0
  w = FloatWidget(d)
  w._line_edit.setText("2.5")
  w._line_edit.validate()
  assert not w.isValid()
  w._line_edit.setText("-0.1")
  w._line_edit.validate()
  assert not w.isValid()
  w._line_edit.setText("0.0")
  w._line_edit.validate()
  assert w.isValid()
  print("exercise_float_widget_limits OK")

from qttbx.widgets.phil.key_widget import KeyWidget

def _make_key_definition():
  return libtbx.phil.parse('output_label = "fobs"\n  .type = key').objects[0]

def exercise_key_widget_accepts_identifier():
  d = _make_key_definition()
  w = KeyWidget(d)
  w._line_edit.setText("abc_123")
  w._line_edit.validate()
  assert w.isValid()
  assert w.value() == "abc_123"
  print("exercise_key_widget_accepts_identifier OK")

def exercise_key_widget_rejects_non_identifier():
  d = _make_key_definition()
  w = KeyWidget(d)
  w._line_edit.setText("9starts_with_digit")
  w._line_edit.validate()
  assert not w.isValid()
  w._line_edit.setText("has space")
  w._line_edit.validate()
  assert not w.isValid()
  print("exercise_key_widget_rejects_non_identifier OK")

from qttbx.widgets.phil.choice_widget import ChoiceWidget

def _make_choice_definition():
  return libtbx.phil.parse(
    "mode = fast *thorough exhaustive\n  .type = choice").objects[0]

def exercise_choice_widget_round_trip():
  d = _make_choice_definition()
  w = ChoiceWidget(d)
  assert w.value() == "thorough"          # the * marks the default
  w.setValue("fast")
  assert w.value() == "fast"
  print("exercise_choice_widget_round_trip OK")

def exercise_choice_widget_rejects_unknown():
  d = _make_choice_definition()
  w = ChoiceWidget(d)
  raised = False
  try:
    w.setValue("not_a_choice")
  except ValueError:
    raised = True
  assert raised
  print("exercise_choice_widget_rejects_unknown OK")

def exercise_choice_widget_set_choices_dynamic():
  d = _make_choice_definition()
  w = ChoiceWidget(d)
  w.setValue("fast")
  w.setChoices(["fast", "extra", "more"])  # "fast" preserved
  assert w.value() == "fast"
  w.setChoices(["x", "y"])                  # current selection gone -> first item
  assert w.value() == "x"
  print("exercise_choice_widget_set_choices_dynamic OK")

def exercise_choice_widget_literal_none_label_does_not_collide():
  """A real choice whose text is "(none)" does NOT collide with the None sentinel."""
  from qttbx.widgets.phil.choice_widget import ChoiceWidget
  _get_app()
  scope = libtbx.phil.parse('''
mode = *fast "(none)" exhaustive
  .type = choice
'''.strip())
  d = scope.objects[0]
  w = ChoiceWidget(d)
  # The literal "(none)" choice is selectable AND distinct from None.
  w.setValue("(none)")
  assert w.value() == "(none)"          # the literal string, NOT None
  # Selecting None still works via the synthetic sentinel.
  w.setValue(None)
  assert w.value() is None
  # And a non-None setValue back recovers the literal again.
  w.setValue("(none)")
  assert w.value() == "(none)"
  print("exercise_choice_widget_literal_none_label_does_not_collide OK")

def _make_str_definition():
  scope = libtbx.phil.parse('''
title = "Untitled"
  .type = str
'''.strip())
  return scope.objects[0]


def exercise_str_widget_round_trip():
  """StrWidget round-trips strings and respects allow_none."""
  from qttbx.widgets.phil.str_widget import StrWidget
  _get_app()
  w = StrWidget(_make_str_definition())
  w.setValue("Hello world")
  assert w.value() == "Hello world"
  w.setValue(None)
  assert w.value() is None
  print("exercise_str_widget_round_trip OK")


def exercise_str_widget_min_max_length():
  """StrWidget rejects strings outside [min_length, max_length]."""
  from qttbx.widgets.phil.str_widget import StrWidget
  _get_app()
  d = _make_str_definition()
  w = StrWidget(d, min_length=3, max_length=8)
  w.setValue("abc")
  assert w.isValid(), w.errorString()
  w.setValue("ab")
  assert not w.isValid()
  w.setValue("abcdefghi")
  assert not w.isValid()
  print("exercise_str_widget_min_max_length OK")


def exercise_str_text_widget_round_trip_preserves_newlines():
  """StrTextWidget preserves embedded newlines."""
  from qttbx.widgets.phil.str_widget import StrTextWidget
  _get_app()
  w = StrTextWidget(_make_str_definition())
  multiline = "first line\nsecond line\nthird"
  w.setValue(multiline)
  assert w.value() == multiline
  print("exercise_str_text_widget_round_trip_preserves_newlines OK")


def _make_qstr_definition():
  scope = libtbx.phil.parse('''
label = "alpha"
  .type = qstr
'''.strip())
  return scope.objects[0]


def exercise_qstr_widget_round_trip():
  from qttbx.widgets.phil.qstr_widget import QstrWidget
  _get_app()
  w = QstrWidget(_make_qstr_definition())
  w.setValue("hello world")
  assert w.value() == "hello world"
  print("exercise_qstr_widget_round_trip OK")


def exercise_qstr_widget_rejects_dollar():
  """QstrWidget rejects '$' per wxtbx convention."""
  from qttbx.widgets.phil.qstr_widget import QstrWidget
  _get_app()
  w = QstrWidget(_make_qstr_definition())
  w.setValue("ok")
  assert w.isValid()
  # Bypass setValue's signal-block: set the text directly and validate.
  w._line_edit.setText("bad$value")
  w._line_edit.validate()
  assert not w.isValid()
  assert "$" in w.errorString()
  print("exercise_qstr_widget_rejects_dollar OK")


def exercise_qstr_text_widget_round_trip():
  from qttbx.widgets.phil.qstr_widget import QstrTextWidget
  _get_app()
  w = QstrTextWidget(_make_qstr_definition())
  w.setValue("first line\nsecond")
  assert w.value() == "first line\nsecond"
  print("exercise_qstr_text_widget_round_trip OK")


def _make_path_definition():
  scope = libtbx.phil.parse('''
input_file = None
  .type = path
'''.strip())
  return scope.objects[0]


def exercise_path_widget_round_trip():
  from qttbx.widgets.phil.path_widget import PathWidget
  _get_app()
  w = PathWidget(_make_path_definition())
  w.setValue("/tmp/example.pdb")
  assert w.value() == "/tmp/example.pdb"
  w.setValue(None)
  assert w.value() is None
  print("exercise_path_widget_round_trip OK")


def exercise_path_widget_must_exist():
  """PathWidget(must_exist=True) rejects non-existent paths."""
  import tempfile
  from qttbx.widgets.phil.path_widget import PathWidget
  _get_app()
  w = PathWidget(_make_path_definition(), must_exist=True)
  with tempfile.NamedTemporaryFile() as f:
    w.setValue(f.name)
    assert w.isValid(), w.errorString()
  # File now removed; validate again.
  w._line_edit.validate()
  assert not w.isValid()
  print("exercise_path_widget_must_exist OK")


def exercise_path_widget_browse_callback():
  """PathWidget._apply_browse_result emits valueChanged."""
  from qttbx.widgets.phil.path_widget import PathWidget
  _get_app()
  w = PathWidget(_make_path_definition())
  emitted = []
  w.valueChanged.connect(lambda v: emitted.append(v))
  w._apply_browse_result("/usr")
  assert emitted == ["/usr"]
  print("exercise_path_widget_browse_callback OK")


def exercise_paths_text_widget_round_trip():
  from qttbx.widgets.phil.path_widget import PathsTextWidget
  _get_app()
  w = PathsTextWidget(_make_path_definition())
  w.setValue(["/a", "/b", "/c"])
  assert w.value() == ["/a", "/b", "/c"]
  print("exercise_paths_text_widget_round_trip OK")


def _make_strings_definition():
  scope = libtbx.phil.parse('''
labels = None
  .type = strings
'''.strip())
  return scope.objects[0]


def exercise_strings_widget_round_trip():
  from qttbx.widgets.phil.strings_widget import StringsWidget
  _get_app()
  w = StringsWidget(_make_strings_definition())
  w.setValue(["a", "b", "c"])
  assert w.value() == ["a", "b", "c"]
  print("exercise_strings_widget_round_trip OK")


def exercise_strings_widget_size_limits():
  from qttbx.widgets.phil.strings_widget import StringsWidget
  _get_app()
  w = StringsWidget(_make_strings_definition(), size_min=2, size_max=4)
  w.setValue(["a", "b"])
  assert w.isValid(), w.errorString()
  w.setValue(["a"])
  assert not w.isValid()
  w.setValue(["a", "b", "c", "d", "e"])
  assert not w.isValid()
  print("exercise_strings_widget_size_limits OK")


def exercise_strings_text_widget_round_trip_strips_blank():
  from qttbx.widgets.phil.strings_widget import StringsTextWidget
  _get_app()
  w = StringsTextWidget(_make_strings_definition())
  w._text_edit.setPlainText("alpha\n\nbeta\n  \ngamma\n")
  w._text_edit.validate()
  assert w.value() == ["alpha", "beta", "gamma"]
  print("exercise_strings_text_widget_round_trip_strips_blank OK")


def _make_words_definition():
  scope = libtbx.phil.parse('''
tags = None
  .type = words
'''.strip())
  return scope.objects[0]


def exercise_words_widget_round_trip_quoted():
  """WordsWidget preserves quoted multi-word tokens."""
  from qttbx.widgets.phil.words_widget import WordsWidget
  _get_app()
  w = WordsWidget(_make_words_definition())
  w._line_edit.setText('alpha "two words" beta')
  w._line_edit.validate()
  values = [tw.value for tw in w.value()]
  assert values == ["alpha", "two words", "beta"], values
  print("exercise_words_widget_round_trip_quoted OK")


def exercise_words_text_widget_round_trip():
  from qttbx.widgets.phil.words_widget import WordsTextWidget
  _get_app()
  w = WordsTextWidget(_make_words_definition())
  w._text_edit.setPlainText("a b\nc d")
  w._text_edit.validate()
  values = [tw.value for tw in w.value()]
  assert values == ["a", "b", "c", "d"], values
  print("exercise_words_text_widget_round_trip OK")


def exercise_words_widget_round_trip_backslash():
  """WordsWidget round-trips tokens with backslash, quotes, and trailing backslash.

  The trailing-backslash case is the canonical demonstration: a value like
  ``ab\\`` formatted naively as ``"ab\\"`` makes the tokenizer treat the
  closing quote as an escaped literal and report a missing-closing-quote
  error. Using ``libtbx.phil.tokenizer.quote_python_str`` produces ``"ab\\\\"``,
  which round-trips correctly.
  """
  from libtbx.phil import tokenizer
  from qttbx.widgets.phil.words_widget import WordsWidget
  _get_app()
  w = WordsWidget(_make_words_definition())
  tokens_in = [
    tokenizer.word(value="plain"),
    tokenizer.word(value=r"a\b", quote_token='"'),
    tokenizer.word(value='has"quote', quote_token='"'),
    tokenizer.word(value="ab\\", quote_token='"'),     # trailing backslash
  ]
  w.setValue(tokens_in)
  values_out = [tw.value for tw in w.value()]
  assert w.isValid(), w.errorString()
  assert values_out == ["plain", r"a\b", 'has"quote', "ab\\"], values_out
  print("exercise_words_widget_round_trip_backslash OK")


def _make_ints_definition(value_min=None, value_max=None):
  if value_min is not None or value_max is not None:
    bounds = []
    if value_min is not None:
      bounds.append("value_min={vm}".format(vm=value_min))
    if value_max is not None:
      bounds.append("value_max={vm}".format(vm=value_max))
    type_str = "ints({b})".format(b=", ".join(bounds))
  else:
    type_str = "ints"
  scope = libtbx.phil.parse('''
nums = None
  .type = {t}
'''.strip().format(t=type_str))
  return scope.objects[0]


def exercise_ints_widget_round_trip():
  from qttbx.widgets.phil.ints_widget import IntsWidget
  _get_app()
  w = IntsWidget(_make_ints_definition())
  w.setValue([1, 2, 3])
  assert w.value() == [1, 2, 3]
  print("exercise_ints_widget_round_trip OK")


def exercise_ints_widget_per_element_bounds():
  from qttbx.widgets.phil.ints_widget import IntsWidget
  _get_app()
  w = IntsWidget(_make_ints_definition(value_min=0, value_max=10))
  w.setValue([1, 5, 9])
  assert w.isValid(), w.errorString()
  w.setValue([1, 11, 5])
  assert not w.isValid()
  print("exercise_ints_widget_per_element_bounds OK")


def exercise_ints_widget_size_bounds():
  from qttbx.widgets.phil.ints_widget import IntsWidget
  _get_app()
  w = IntsWidget(_make_ints_definition(), size_min=2, size_max=4)
  w.setValue([1])
  assert not w.isValid()
  w.setValue([1, 2, 3])
  assert w.isValid()
  print("exercise_ints_widget_size_bounds OK")


def exercise_ints_text_widget_round_trip():
  from qttbx.widgets.phil.ints_widget import IntsTextWidget
  _get_app()
  w = IntsTextWidget(_make_ints_definition())
  w.setValue([10, 20, 30])
  assert w.value() == [10, 20, 30]
  print("exercise_ints_text_widget_round_trip OK")


def _make_floats_definition(value_min=None, value_max=None):
  if value_min is not None or value_max is not None:
    bounds = []
    if value_min is not None:
      bounds.append("value_min={vm}".format(vm=value_min))
    if value_max is not None:
      bounds.append("value_max={vm}".format(vm=value_max))
    type_str = "floats({b})".format(b=", ".join(bounds))
  else:
    type_str = "floats"
  scope = libtbx.phil.parse('''
weights = None
  .type = {t}
'''.strip().format(t=type_str))
  return scope.objects[0]


def exercise_floats_widget_round_trip():
  from qttbx.widgets.phil.floats_widget import FloatsWidget
  _get_app()
  w = FloatsWidget(_make_floats_definition())
  w.setValue([1.5, 2.5, 3.5])
  out = w.value()
  for got, want in zip(out, [1.5, 2.5, 3.5]):
    assert approx_equal(got, want)
  print("exercise_floats_widget_round_trip OK")


def exercise_floats_widget_per_element_bounds():
  from qttbx.widgets.phil.floats_widget import FloatsWidget
  _get_app()
  w = FloatsWidget(_make_floats_definition(value_min=0.0, value_max=1.0))
  w.setValue([0.1, 0.5, 0.9])
  assert w.isValid(), w.errorString()
  w.setValue([0.1, 1.5, 0.5])
  assert not w.isValid()
  print("exercise_floats_widget_per_element_bounds OK")


def exercise_floats_text_widget_round_trip():
  from qttbx.widgets.phil.floats_widget import FloatsTextWidget
  _get_app()
  w = FloatsTextWidget(_make_floats_definition())
  w.setValue([10.0, 20.0])
  out = w.value()
  for got, want in zip(out, [10.0, 20.0]):
    assert approx_equal(got, want)
  print("exercise_floats_text_widget_round_trip OK")


def _make_choice_multi_definition(optional=True):
  opt_attr = "" if optional else "\n  .optional = False"
  scope = libtbx.phil.parse('''
methods = *fast slow exhaustive
  .type = choice(multi=True){opt}
'''.strip().format(opt=opt_attr))
  return scope.objects[0]


def exercise_choice_multi_widget_round_trip():
  from qttbx.widgets.phil.choice_widget import ChoiceMultiWidget
  _get_app()
  w = ChoiceMultiWidget(_make_choice_multi_definition())
  w.setValue(["fast", "exhaustive"])
  assert sorted(w.value()) == ["exhaustive", "fast"]
  print("exercise_choice_multi_widget_round_trip OK")


def exercise_choice_multi_widget_dispatch_via_registry():
  from qttbx.widgets.phil import widget_for_definition
  _get_app()
  d = _make_choice_multi_definition()
  w = widget_for_definition(d)
  # Duck-type check: it has the multi value() shape.
  w.setValue(["slow"])
  assert w.value() == ["slow"]
  # The class is the multi-variant.
  assert type(w).__name__ == "ChoiceMultiWidget"
  print("exercise_choice_multi_widget_dispatch_via_registry OK")


def exercise_choice_multi_widget_optional_false_requires_selection():
  from qttbx.widgets.phil.choice_widget import ChoiceMultiWidget
  _get_app()
  w = ChoiceMultiWidget(_make_choice_multi_definition(optional=False))
  w.setValue([])
  assert not w.isValid()
  w.setValue(["fast"])
  assert w.isValid()
  print("exercise_choice_multi_widget_optional_false_requires_selection OK")

from PySide2.QtCore import Qt, QCoreApplication
from qttbx.widgets.phil import PhilField

def exercise_phil_field_widget_to_model():
  m = PhilModel()
  m.initialize_model(_example_master())
  f = PhilField(m, "refinement.macro_cycles")
  inner = f.widget()
  inner.setValue(8)
  inner._line_edit._commit()                       # simulate focus loss
  assert m.get_phil_extract_value("refinement.macro_cycles") == 8
  print("exercise_phil_field_widget_to_model OK")

def exercise_phil_field_model_to_widget():
  m = PhilModel()
  m.initialize_model(_example_master())
  f = PhilField(m, "refinement.macro_cycles")
  idx = m.index_for_path("refinement.macro_cycles")
  m.setData(idx, 4, Qt.EditRole)
  assert f.widget().value() == 4
  print("exercise_phil_field_model_to_widget OK")

def exercise_phil_field_no_recursion():
  m = PhilModel()
  m.initialize_model(_example_master())
  f = PhilField(m, "refinement.macro_cycles")
  count = [0]
  orig_set = m.setData
  def counting_set(*a, **kw):
    count[0] += 1
    return orig_set(*a, **kw)
  m.setData = counting_set
  idx = m.index_for_path("refinement.macro_cycles")
  m.setData(idx, 9, Qt.EditRole)
  assert count[0] == 1, "recursion: setData called %d times" % count[0]
  print("exercise_phil_field_no_recursion OK")

def exercise_phil_field_lifecycle_no_signal_after_delete():
  m = PhilModel()
  m.initialize_model(_example_master())
  f = PhilField(m, "refinement.macro_cycles")
  f.deleteLater()
  QCoreApplication.processEvents()
  idx = m.index_for_path("refinement.macro_cycles")
  m.setData(idx, 6, Qt.EditRole)        # would crash if slot still bound
  print("exercise_phil_field_lifecycle_no_signal_after_delete OK")

def exercise_phil_field_error_icon_visibility():
  m = PhilModel()
  m.initialize_model(_example_master())
  f = PhilField(m, "refinement.macro_cycles")
  f.show()                                   # so isVisible() reflects show/hide
  inner = f.widget()
  inner._line_edit.setText("99")             # > value_max=20
  inner._line_edit.validate()
  assert not inner.isValid()
  assert f._error_icon.isVisible()
  assert "20" in f._error_icon.toolTip() or "20" in inner.toolTip()
  inner._line_edit.setText("5")              # back to valid
  inner._line_edit.validate()
  assert not f._error_icon.isVisible()
  print("exercise_phil_field_error_icon_visibility OK")


def exercise_phil_field_persistent_index_basic():
  """PhilField built with persistent_index= behaves the same as full_path=."""
  from qttbx.phil import PhilModel
  from qttbx.widgets.phil import PhilField
  _get_app()
  m = PhilModel()
  m.initialize_model(libtbx.phil.parse('''
refine {
  cycles = 3
    .type = int
}
'''.strip()))
  qpi = m.persistent_index_for_path("refine.cycles")
  field = PhilField(m, persistent_index=qpi)
  assert field.widget().value() == 3
  field.widget().setValue(7)
  field.commit()
  assert m.get_phil_extract().refine.cycles == 7
  print("exercise_phil_field_persistent_index_basic OK")


def exercise_phil_field_rejects_neither_or_both_kwargs():
  """PhilField raises if neither or both of full_path / persistent_index is given."""
  from qttbx.phil import PhilModel
  from qttbx.widgets.phil import PhilField
  _get_app()
  m = PhilModel()
  m.initialize_model(libtbx.phil.parse('''
refine {
  cycles = 3
    .type = int
}
'''.strip()))
  qpi = m.persistent_index_for_path("refine.cycles")
  # Neither kwarg.
  try:
    PhilField(m)
  except ValueError as e:
    assert "exactly one" in str(e)
  else:
    raise AssertionError("expected ValueError on neither kwarg")
  # Both kwargs.
  try:
    PhilField(m, full_path="refine.cycles", persistent_index=qpi)
  except ValueError as e:
    assert "exactly one" in str(e)
  else:
    raise AssertionError("expected ValueError on both kwargs")
  print("exercise_phil_field_rejects_neither_or_both_kwargs OK")


def exercise_phil_field_commit_returns_bool():
  """commit() returns True on valid, False on invalid."""
  from qttbx.phil import PhilModel
  from qttbx.widgets.phil import PhilField
  _get_app()
  m = PhilModel()
  m.initialize_model(libtbx.phil.parse('''
val = 5
  .type = int(value_min=1, value_max=10)
'''.strip()))
  f = PhilField(m, "val")
  # Valid input -> True.
  f.widget().setValue(7)
  assert f.commit() is True
  assert m.get_phil_extract().val == 7
  # Force the inner widget invalid by setting out-of-range text directly.
  f.widget()._line_edit.setText("99")
  f.widget()._line_edit.validate()
  assert not f.widget().isValid()
  assert f.commit() is False
  # Model still holds the previous valid value.
  assert m.get_phil_extract().val == 7
  print("exercise_phil_field_commit_returns_bool OK")


def exercise_phil_field_label_uses_short_caption_or_prettified():
  """PhilField shows short_caption when set; falls back to prettified name."""
  from qttbx.phil import PhilModel
  from qttbx.widgets.phil import PhilField
  _get_app()
  m = PhilModel()
  m.initialize_model(libtbx.phil.parse('''
atom_selection = "all"
  .type = str
labelled = 1
  .type = int
  .short_caption = "Custom Label"
'''.strip()))
  f_default = PhilField(m, "atom_selection")
  f_caption = PhilField(m, "labelled")
  assert f_default._label.text() == "Atom selection", f_default._label.text()
  assert f_caption._label.text() == "Custom Label", f_caption._label.text()
  print("exercise_phil_field_label_uses_short_caption_or_prettified OK")
from PySide2.QtWidgets import QTreeView, QStyleOptionViewItem
from qttbx.widgets.phil.delegate import PhilItemDelegate

def exercise_phil_item_delegate_round_trip():
  m = PhilModel()
  m.initialize_model(_example_master())
  view = QTreeView()
  view.setModel(m)
  d = PhilItemDelegate()
  view.setItemDelegate(d)
  idx = m.index_for_path("refinement.macro_cycles")
  editor = d.createEditor(view, QStyleOptionViewItem(), idx)
  d.setEditorData(editor, idx)
  assert editor.value() == 3
  # _example_master uses value_max=20, so 25 is definitely above the limit.
  editor.setValue(25)
  assert not editor.isValid()
  d.setModelData(editor, m, idx)      # invalid -> model unchanged
  assert m.get_phil_extract_value("refinement.macro_cycles") == 3
  editor.setValue(7)
  assert editor.isValid()
  d.setModelData(editor, m, idx)
  assert m.get_phil_extract_value("refinement.macro_cycles") == 7
  print("exercise_phil_item_delegate_round_trip OK")

from PySide2.QtWidgets import (
  QFormLayout, QDialog, QDialogButtonBox, QWidget as _QWidget)

def _v1_master():
  return libtbx.phil.parse("""
refinement {
  macro_cycles = 3
    .type = int(value_min=1, value_max=20)
  weight = 1.0
    .type = float(value_min=0.0, value_max=10.0)
  use_ncs = True
    .type = bool
  mode = fast *thorough exhaustive
    .type = choice
  output_label = result
    .type = key
}
""")

def exercise_v1_tree_and_form_sync():
  m = PhilModel()
  m.initialize_model(_v1_master())

  # Tree view side
  tree = QTreeView()
  tree.setModel(m)
  tree.setItemDelegate(PhilItemDelegate())

  # Form view side: one PhilField per parameter.
  form_host = _QWidget()
  form = QFormLayout(form_host)
  fields = {}
  for name in ("macro_cycles", "weight", "use_ncs", "mode", "output_label"):
    f = PhilField(m, "refinement.{}".format(name))
    form.addRow(f)
    fields[name] = f

  # Edit via the form -- verify model and tree see the change.
  fields["macro_cycles"].widget().setValue(7)
  fields["macro_cycles"].widget()._line_edit._commit()
  assert m.get_phil_extract_value("refinement.macro_cycles") == 7

  # Edit via setData (simulating tree-view delegate commit) -- verify form sees it.
  idx = m.index_for_path("refinement.weight")
  m.setData(idx, 5.5, Qt.EditRole)
  assert approx_equal(fields["weight"].widget().value(), 5.5)

  # Bool toggle through form.
  fields["use_ncs"].widget().setValue(False)
  fields["use_ncs"].commit()
  assert m.get_phil_extract_value("refinement.use_ncs") is False

  # Choice change through form.
  fields["mode"].widget().setValue("exhaustive")
  fields["mode"].commit()
  assert m.get_phil_extract_value("refinement.mode") == "exhaustive"

  # Validity gating on the macro_cycles field with invalid input.
  fields["macro_cycles"].widget()._line_edit.setText("99")  # > value_max=20
  fields["macro_cycles"].widget()._line_edit.validate()
  assert not fields["macro_cycles"].isValid()
  # commit() should be a no-op while invalid; model retains the prior value.
  fields["macro_cycles"].commit()
  assert m.get_phil_extract_value("refinement.macro_cycles") == 7
  print("exercise_v1_tree_and_form_sync OK")

def _v2_master():
  return libtbx.phil.parse("""
files {
  input_file = None
    .type = path
  labels = None
    .type = strings
  weights = None
    .type = floats(value_min=0.0)
}
notes = ""
  .type = str
methods = *fast slow exhaustive
  .type = choice(multi=True)
""")


def exercise_v2_tree_and_form_sync():
  """End-to-end sync between a tree view and a form for v2 types."""
  from PySide2.QtCore import Qt
  from PySide2.QtWidgets import QTreeView
  from qttbx.phil import PhilModel
  from qttbx.widgets.phil import PhilField
  from qttbx.widgets.phil.delegate import PhilItemDelegate
  from qttbx.widgets.phil.str_widget import StrTextWidget
  _get_app()

  m = PhilModel()
  m.initialize_model(_v2_master())

  view = QTreeView()
  view.setModel(m)
  view.setItemDelegate(PhilItemDelegate())

  fields = {
    "files.input_file": PhilField(m, "files.input_file"),
    "files.labels": PhilField(m, "files.labels"),
    "files.weights": PhilField(m, "files.weights"),
    "notes": PhilField(m, "notes", widget=StrTextWidget),
    "methods": PhilField(m, "methods"),
  }

  # Form -> model: edit each widget via setValue and confirm the model picks it up.
  fields["files.input_file"].widget().setValue("/tmp/x.pdb")
  fields["files.input_file"].commit()
  assert m.get_phil_extract().files.input_file == "/tmp/x.pdb"

  fields["files.labels"].widget().setValue(["F", "SIGF"])
  fields["files.labels"].commit()
  assert m.get_phil_extract().files.labels == ["F", "SIGF"]

  fields["files.weights"].widget().setValue([0.5, 1.0])
  fields["files.weights"].commit()
  for got, want in zip(m.get_phil_extract().files.weights, [0.5, 1.0]):
    assert approx_equal(got, want)

  fields["notes"].widget().setValue("First-pass refinement")
  fields["notes"].commit()
  assert m.get_phil_extract().notes == "First-pass refinement"

  fields["methods"].widget().setValue(["slow", "exhaustive"])
  fields["methods"].commit()
  assert sorted(m.get_phil_extract().methods) == ["exhaustive", "slow"]

  # Model -> form: drive the model directly and confirm the form reflects it.
  idx = m.index_for_path("notes")
  m.setData(idx, "Updated notes", Qt.EditRole)
  assert fields["notes"].widget().value() == "Updated notes"

  print("exercise_v2_tree_and_form_sync OK")


def run_all():
  _get_app()
  exercise_phil_model_get_definition()
  exercise_phil_model_thread_assertion_passes_in_gui_thread()
  exercise_phil_model_index_for_path()
  exercise_phil_model_populate_after_factory_removal()
  exercise_phil_item_extract_path_simple()
  exercise_phil_item_extract_path_with_instance_index()
  exercise_phil_model_multi_instance_tree()
  exercise_phil_model_indexed_get_set()
  exercise_phil_model_persistent_index_simple()
  exercise_phil_model_persistent_index_multi_scope()
  exercise_phil_model_persistent_index_requires_scope_indices()
  exercise_phil_model_add_scope_instance_grows_tree()
  exercise_phil_model_remove_scope_instance_shrinks_tree()
  exercise_phil_model_add_scope_instance_treeview_sync()
  exercise_phil_model_restore_defaults()
  exercise_phil_model_intrinsic_validation_int_bounds()
  exercise_phil_model_intrinsic_validation_choice_membership()
  exercise_phil_model_intrinsic_validation_at_initialize()
  exercise_label_for_definition_helper()
  exercise_phil_widget_base_construct()
  exercise_phil_widget_base_signals_exist()
  exercise_validated_line_edit_valid_to_invalid()
  exercise_validated_line_edit_emits_validity_changed()
  exercise_validated_text_edit_valid_to_invalid()
  exercise_validated_text_edit_focus_out_commits()
  exercise_registry_register_and_dispatch()
  exercise_registry_unknown_type_raises()
  exercise_widget_for_definition_unregistered_raises_value_error()
  exercise_widget_for_definition_fallback_kwarg()
  # Restore builtin widgets after exercise_registry_register_and_dispatch
  # replaced the "int" registration with _DummyWidget.
  from qttbx.widgets.phil import register_builtin_widgets
  register_builtin_widgets()
  exercise_int_widget_round_trip()
  exercise_int_widget_limits()
  exercise_int_widget_setvalue_idempotent()
  exercise_bool_widget_round_trip()
  exercise_bool_widget_tristate_for_none()
  exercise_float_widget_round_trip()
  exercise_float_widget_limits()
  exercise_key_widget_accepts_identifier()
  exercise_key_widget_rejects_non_identifier()
  exercise_choice_widget_round_trip()
  exercise_choice_widget_rejects_unknown()
  exercise_choice_widget_set_choices_dynamic()
  exercise_choice_widget_literal_none_label_does_not_collide()
  exercise_str_widget_round_trip()
  exercise_str_widget_min_max_length()
  exercise_str_text_widget_round_trip_preserves_newlines()
  exercise_qstr_widget_round_trip()
  exercise_qstr_widget_rejects_dollar()
  exercise_qstr_text_widget_round_trip()
  exercise_path_widget_round_trip()
  exercise_path_widget_must_exist()
  exercise_path_widget_browse_callback()
  exercise_paths_text_widget_round_trip()
  exercise_strings_widget_round_trip()
  exercise_strings_widget_size_limits()
  exercise_strings_text_widget_round_trip_strips_blank()
  exercise_words_widget_round_trip_quoted()
  exercise_words_text_widget_round_trip()
  exercise_words_widget_round_trip_backslash()
  exercise_ints_widget_round_trip()
  exercise_ints_widget_per_element_bounds()
  exercise_ints_widget_size_bounds()
  exercise_ints_text_widget_round_trip()
  exercise_floats_widget_round_trip()
  exercise_floats_widget_per_element_bounds()
  exercise_floats_text_widget_round_trip()
  exercise_choice_multi_widget_round_trip()
  exercise_choice_multi_widget_dispatch_via_registry()
  exercise_choice_multi_widget_optional_false_requires_selection()
  exercise_phil_field_widget_to_model()
  exercise_phil_field_model_to_widget()
  exercise_phil_field_no_recursion()
  exercise_phil_field_lifecycle_no_signal_after_delete()
  exercise_phil_field_error_icon_visibility()
  exercise_phil_field_persistent_index_basic()
  exercise_phil_field_rejects_neither_or_both_kwargs()
  exercise_phil_field_commit_returns_bool()
  exercise_phil_field_label_uses_short_caption_or_prettified()
  exercise_phil_item_delegate_round_trip()
  exercise_v1_tree_and_form_sync()
  exercise_v2_tree_and_form_sync()

if __name__ == "__main__":
  run_all()
  print(format_cpu_times())
  print("OK")
