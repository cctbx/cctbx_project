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

if __name__ == "__main__":
  run_all()
  print(format_cpu_times())
  print("OK")
