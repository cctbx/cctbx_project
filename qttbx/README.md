# qttbx

PySide2 (Qt5) GUI module for the cctbx project. The flagship subsystem is
`qttbx.widgets.phil` — a complete set of Qt widgets for editing
[PHIL parameters](https://phenix-online.org/phenix_html/static/papers/d-72-00018.pdf)
with built-in type and limit validation.

This README documents the PHIL widget package: its design, the catalog of
widget classes, and the patterns for using them.

---

## What this package provides

For every PHIL parameter type, `qttbx.widgets.phil` provides a Qt widget that:

- Validates user input against the PHIL type and any limits
  (`value_min`/`value_max`, `size_min`/`size_max`, choice membership, etc.)
- Renders inline (red border + tooltip) feedback when input is invalid —
  no modal dialogs.
- Round-trips Python values via `value()` / `setValue(value)`.
- Plugs into either a tree-view of the whole PHIL scope (one widget per leaf
  cell, via `PhilItemDelegate`) or a form-style dialog (one widget per field,
  via `PhilField`).

The package handles every PHIL type defined in `libtbx.phil` plus the
domain-specific types defined in `iotbx.phil` (`space_group`, `unit_cell`,
`atom_selection`). Custom site types can be registered with
`register_widget(phil_type, widget_cls)`.

---

## Architecture

Three layers, MVC-shaped:

```
                       ┌────────────────────┐
                       │   libtbx.phil      │   master scope (template)
                       │   .scope object    │   .extract() value tree
                       └──────────┬─────────┘
                                  │
                  ┌───────────────▼────────────────┐
   MODEL          │       qttbx.phil.PhilModel     │   QAbstractItemModel
                  │  • initialize_model(scope)     │   that wraps a PHIL
                  │  • get_phil_extract()          │   scope and its extract
                  │  • setData / get_phil_extract  │   together; emits
                  │  • add/remove_scope_instance   │   dataChanged on edits
                  │  • persistent_index_for_path   │
                  └────────┬────────────────┬──────┘
                           │                │
        ┌──────────────────▼────┐   ┌───────▼────────────────────┐
   CTRL │   PhilItemDelegate    │   │      PhilField             │
        │   (tree-view path)    │   │      (form-view path)      │
        │   QStyledItemDelegate │   │      QWidget container     │
        │   builds PhilWidget   │   │      builds PhilWidget +   │
        │   on createEditor     │   │      label + error icon    │
        └──────────────┬────────┘   └────────────┬───────────────┘
                       │                         │
                  ┌────▼─────────────────────────▼────┐
   VIEW           │           PhilWidget (base)       │   pure view; reads
                  │   IntWidget, BoolWidget, …        │   value from inner
                  │   StrWidget, StrTextWidget, …     │   Qt primitive on
                  │   MultipleWidget,                 │   demand (no cached
                  │   RepeatableScopeWidget,          │   _value field)
                  │   SpaceGroupWidget, …             │
                  └───────────────────────────────────┘
```

Two consumer paths — tree view (every leaf in one widget tree) or form view
(a hand-built dialog of `PhilField`s) — share the same `PhilModel`. Edits
made in one path propagate to the other automatically via Qt signals.

A `PhilWidget` never holds a separate `self._value` field — `value()` reads
from the underlying Qt primitive (`QLineEdit.text()`, `QSpinBox.value()`,
`QListWidget.itemWidget(...)`, etc.) so the widget state and the reported
value cannot diverge.

---

## Quick start: a minimal form

```python
from PySide2.QtWidgets import QApplication, QDialog, QFormLayout, QDialogButtonBox
import libtbx.phil
from qttbx.phil import PhilModel
from qttbx.widgets.phil import PhilField

master = libtbx.phil.parse("""
refinement {
  macro_cycles = 3
    .type = int(value_min=1, value_max=20)
  weight = 0.5
    .type = float(value_min=0.0)
  use_ncs = False
    .type = bool
}
""")

app = QApplication([])
model = PhilModel()
model.initialize_model(master)

dlg = QDialog()
form = QFormLayout(dlg)
fields = [
  PhilField(model, "refinement.macro_cycles"),
  PhilField(model, "refinement.weight"),
  PhilField(model, "refinement.use_ncs"),
]
for f in fields:
  form.addRow(f)

buttons = QDialogButtonBox(
  QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
form.addRow(buttons)
buttons.accepted.connect(dlg.accept)
buttons.rejected.connect(dlg.reject)

if dlg.exec_() == QDialog.Accepted:
  for f in fields:
    f.commit()                        # commit any pending edit
  params = model.get_phil_extract()
  print(params.refinement.macro_cycles, params.refinement.weight)
```

Each `PhilField` automatically picks the right widget class based on the
PHIL type at the given path. The `PhilModel` holds the single source of
truth; `PhilField`s synchronize bidirectionally with it.

---

## Widget catalog

All widgets are subclasses of `qttbx.widgets.phil.PhilWidget`. The
**registered default** is what `PhilField` / `PhilItemDelegate` produce for
each PHIL type unless overridden via `widget=` (form-view only).

### Scalar widgets (registered as defaults)

| PHIL type | Widget class | Module |
|---|---|---|
| `int` | `IntWidget` | `qttbx.widgets.phil.int_widget` |
| `bool` | `BoolWidget` (tri-state when `allow_none` or `use_auto`) | `qttbx.widgets.phil.bool_widget` |
| `float` | `FloatWidget` | `qttbx.widgets.phil.float_widget` |
| `key` | `KeyWidget` | `qttbx.widgets.phil.key_widget` |
| `choice` | `ChoiceWidget` (combobox; supports `setChoices()` for runtime lists) | `qttbx.widgets.phil.choice_widget` |
| `choice(multi=True)` | `ChoiceMultiWidget` (checkable list) — dispatched by `widget_for_definition` | `qttbx.widgets.phil.choice_widget` |
| `str` | `StrWidget` | `qttbx.widgets.phil.str_widget` |
| `qstr` | `QstrWidget` | `qttbx.widgets.phil.qstr_widget` |
| `path` | `PathWidget` (line edit + `[…]` browse button) | `qttbx.widgets.phil.path_widget` |

### List-of-scalars widgets (registered as defaults)

| PHIL type | Widget class (short) | Module |
|---|---|---|
| `strings` | `StringsWidget` | `qttbx.widgets.phil.strings_widget` |
| `words` | `WordsWidget` (uses `libtbx.phil.tokenizer`) | `qttbx.widgets.phil.words_widget` |
| `ints` | `IntsWidget` (per-element + size bounds, `allow_none_elements`) | `qttbx.widgets.phil.ints_widget` |
| `floats` | `FloatsWidget` (same; optional `decimals` kwarg) | `qttbx.widgets.phil.floats_widget` |

### Long-text variants (NOT registered — opt in via `widget=`)

These are used in form views when a single line is too cramped. The tree-view
delegate always uses the short variant.

| Long widget | Pair with | Module |
|---|---|---|
| `StrTextWidget` | `str` | `qttbx.widgets.phil.str_widget` |
| `QstrTextWidget` | `qstr` | `qttbx.widgets.phil.qstr_widget` |
| `PathsTextWidget` | multi-path entry | `qttbx.widgets.phil.path_widget` |
| `StringsTextWidget` | `strings` (one per line) | `qttbx.widgets.phil.strings_widget` |
| `WordsTextWidget` | `words` (one per line) | `qttbx.widgets.phil.words_widget` |
| `IntsTextWidget` | `ints` (one per line) | `qttbx.widgets.phil.ints_widget` |
| `FloatsTextWidget` | `floats` (one per line) | `qttbx.widgets.phil.floats_widget` |

Usage:

```python
from qttbx.widgets.phil.str_widget import StrTextWidget
PhilField(model, "refinement.notes", widget=StrTextWidget)
```

### Domain-specific widgets (self-register on import)

These cover the PHIL types defined in `iotbx.phil` (not plain `libtbx.phil`).
Tests and consumers MUST parse the master scope with
`iotbx.phil.parse(...)` (or pass `iotbx.phil.default_converter_registry` to
`libtbx.phil.parse`). Importing the module triggers `register_widget(...)`.

| PHIL type | Widget class | Validator |
|---|---|---|
| `space_group` | `SpaceGroupWidget` | `cctbx.sgtbx.space_group_info` (accepts SG number, full/short Hermann-Mauguin, Hall via fallback) |
| `unit_cell` | `UnitCellWidget` | `cctbx.uctbx.unit_cell` (six whitespace-separated floats) |
| `atom_selection` | `AtomSelectionWidget` (single line) / `AtomSelectionTextWidget` (multi-line, opt-in) | string-shape only in v4; semantic validation deferred |

Each module registers itself when imported:

```python
import qttbx.widgets.phil.space_group     # registers "space_group"
import qttbx.widgets.phil.unit_cell       # registers "unit_cell"
import qttbx.widgets.phil.atom_selection  # registers "atom_selection"
```

If your master scope uses any of these types, import the corresponding
module(s) **before** building any `PhilField` or `PhilItemDelegate`.

### Multiples

| Definition / scope | Widget class | Module |
|---|---|---|
| Definition with `.multiple = True` (and type not in `{ints, floats, strings, words}`) | `MultipleWidget` — wraps N inner widgets in `QListWidget`, plus `[+] [−] [↑] [↓]` toolbar | `qttbx.widgets.phil.multiple` |
| Scope with `.multiple = True` | `RepeatableScopeWidget` — `QTabWidget` with `[+]` / `[×]` corner buttons; one tab per scope-extract instance | `qttbx.widgets.phil.multiple` |

`MultipleWidget` is dispatched automatically by `widget_for_definition`.
`RepeatableScopeWidget` is instantiated explicitly by consumer code (see
[Scope multiples](#scope-multiples)).

---

## Tree-view consumer

A `QTreeView` over a `PhilModel` with `PhilItemDelegate` produces an
editor for every leaf parameter. Definition leaves use their registered
short widget; scope rows are read-only "headers" that expand to show
their children. Multi-scope rows expand to N children (one per instance).

```python
from PySide2.QtCore import QModelIndex
from PySide2.QtWidgets import QApplication, QPushButton, QTreeView, QVBoxLayout, QWidget
from qttbx.phil import PhilModel
from qttbx.widgets.phil.delegate import PhilItemDelegate
import libtbx.phil

master = libtbx.phil.parse("""
refinement {
  macro_cycles = 3
    .type = int(value_min=1, value_max=20)
  weight = 0.5
    .type = float
}
""")

app = QApplication([])
model = PhilModel()
model.initialize_model(master)

view = QTreeView()
view.setModel(model)
view.setItemDelegate(PhilItemDelegate())

# A "Run" button outside the tree view; flush any active cell editor before
# reading the model so the user's pending edit is captured.
def on_run_clicked():
  view.setCurrentIndex(QModelIndex())     # forces editor commit / discard
  params = model.get_phil_extract()
  print(params.refinement.macro_cycles, params.refinement.weight)

run_btn = QPushButton("Run")
run_btn.clicked.connect(on_run_clicked)

container = QWidget()
layout = QVBoxLayout(container)
layout.addWidget(view)
layout.addWidget(run_btn)
container.show()
app.exec_()
```

**Validation feedback in non-editing cells.** The model surfaces
definition-intrinsic validation (`value_min` / `value_max`, choice
membership, `size_min` / `size_max`) via `Qt.BackgroundRole` (faint pink)
and `Qt.ToolTipRole` (the error message) on column 1. This catches
out-of-range values written by external code (worker threads, scripts)
that never opened an editor. Widget-author-level constraints
(`must_exist=True` on PathWidget, `min_length` on StrWidget) are NOT
enforced at this layer; they remain editor-time-only.

---

## Form-view consumer

A `QFormLayout` of `PhilField`s gives a hand-laid-out dialog where each
field can be sized, ordered, and visually grouped independently of the
underlying PHIL tree.

```python
from PySide2.QtWidgets import QApplication, QDialog, QFormLayout, QDialogButtonBox
from qttbx.phil import PhilModel
from qttbx.widgets.phil import PhilField
from qttbx.widgets.phil.str_widget import StrTextWidget
import libtbx.phil


class RefinementDialog(QDialog):
  def __init__(self, model, parent=None):
    super().__init__(parent)
    layout = QFormLayout(self)

    self._fields = [
      PhilField(model, "refinement.macro_cycles"),
      PhilField(model, "refinement.weight"),
      # Override default short widget with the long variant:
      PhilField(model, "refinement.notes", widget=StrTextWidget),
    ]
    for f in self._fields:
      layout.addRow(f)

    self._buttons = QDialogButtonBox(
      QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    layout.addRow(self._buttons)
    self._buttons.accepted.connect(self.accept)
    self._buttons.rejected.connect(self.reject)

    # Wire OK enabled state to per-field validity.
    for f in self._fields:
      f.validityChanged.connect(self._update_ok)
    self._update_ok(True)

  def _update_ok(self, _):
    ok_btn = self._buttons.button(QDialogButtonBox.Ok)
    ok_btn.setEnabled(all(f.isValid() for f in self._fields))

  def accept(self):
    # commit() returns True on valid, False on invalid; refuse to close
    # if any field rejected its commit.
    if not all(f.commit() for f in self._fields):
      return
    super().accept()


master = libtbx.phil.parse("""
refinement {
  macro_cycles = 3
    .type = int(value_min=1, value_max=20)
  weight = 0.5
    .type = float(value_min=0.0)
  notes = ""
    .type = str
}
""")

app = QApplication([])
model = PhilModel()
model.initialize_model(master)

dlg = RefinementDialog(model)
if dlg.exec_() == QDialog.Accepted:
  params = model.get_phil_extract()
  print(params.refinement.macro_cycles)
  print(params.refinement.weight)
  print(repr(params.refinement.notes))
```

### Two views over the same model

```python
from PySide2.QtWidgets import (
  QApplication, QHBoxLayout, QTreeView, QWidget, QFormLayout,
)
from qttbx.phil import PhilModel
from qttbx.widgets.phil import PhilField
from qttbx.widgets.phil.delegate import PhilItemDelegate
import libtbx.phil

master = libtbx.phil.parse("""
refinement {
  macro_cycles = 3
    .type = int(value_min=1, value_max=20)
  weight = 0.5
    .type = float
}
""")

app = QApplication([])
model = PhilModel()
model.initialize_model(master)

# Tree view (left): every leaf editable inline.
tree = QTreeView()
tree.setModel(model)
tree.setItemDelegate(PhilItemDelegate())

# Form (right): a couple of hand-laid-out fields over the same model.
form_host = QWidget()
form = QFormLayout(form_host)
form.addRow(PhilField(model, "refinement.macro_cycles"))
form.addRow(PhilField(model, "refinement.weight"))

container = QWidget()
hbox = QHBoxLayout(container)
hbox.addWidget(tree)
hbox.addWidget(form_host)
container.show()
app.exec_()
# Edit a value in either pane — the other updates automatically via dataChanged.
```

---

## Validation contract

All editing widgets implement a uniform validity contract:

```python
widget.value()           # parse current Qt-primitive state; may raise
widget.setValue(v)       # idempotent, signal-blocked
widget.isValid()         # current parse is valid
widget.errorString()     # error message, or '' when valid
widget.setReadOnly(ro)   # propagates to the inner Qt primitive
widget.setAllowNone(b)   # empty input → None when True (default)
widget.setUseAuto(b)     # empty input → libtbx.Auto when True
```

**Signals:**

```python
widget.valueChanged       # Signal(object) — fires on commit (focus loss /
                          # Enter / choice change). Payload is the value.
widget.validityChanged    # Signal(bool)   — fires on transitions only
```

**Visual feedback:**
- Invalid input draws a red border on the inner line edit / text edit.
- The error message is in `errorString()`; `PhilField` shows it in a
  warning-icon tooltip.

**Theming.** Every validation color (invalid background, invalid border,
warning-icon foreground, low-emphasis side labels) derives from the
active `QApplication.palette()` via helpers in
`qttbx.widgets.phil._colors`. The widgets adapt to light and dark themes
without per-widget configuration: light palettes get a pale-pink invalid
background; dark palettes get a desaturated dark-red one. Tree-view
`Qt.BackgroundRole` re-evaluates per paint, so palette changes propagate
on the next repaint.

**Best-effort `valueChanged`:** the signal never fires for invalid
parses; consumers gate on `validityChanged` for "the widget is now in a
state where reads will succeed".

---

## Widget labels

All three label sources — `PhilItem.label_text` (tree-view column 0),
`PhilField`'s editor label, and `RepeatableScopeWidget` tab text — route
through a single helper:

```python
from qttbx.phil import label_for_definition

label_for_definition(definition)
```

The rule:

1. If the PHIL definition (or scope) declares `.short_caption = "..."`,
   that string is used verbatim.
2. Otherwise, the parameter name is prettified: underscores become spaces
   and the first letter is capitalized — e.g., `atom_selection` →
   `"Atom selection"`, `ncs_group` → `"Ncs group"`.

`RepeatableScopeWidget` appends ` <i>` (instance index, 1-based) so tabs
read `"NCS group 1"`, `"NCS group 2"` (or `"Ncs group 1"`, etc., when no
short_caption is set).

To customize a label, set `.short_caption` on the master PHIL:

```
ncs_group
  .short_caption = "NCS group"
{
  selection = "all"
    .type = str
}
```

---

## Multiple handling

PHIL `.multiple = True` applies to two distinct things and uses two
different UI strategies.

### Definition multiples (`MultipleWidget`)

A *definition* with `.multiple = True` produces a list value. The
registry automatically wraps such definitions in `MultipleWidget` (unless
the type is already list-native: `ints`, `floats`, `strings`, `words`).

```python
from PySide2.QtWidgets import QApplication, QDialog, QFormLayout, QDialogButtonBox
from qttbx.phil import PhilModel
from qttbx.widgets.phil import PhilField
import libtbx.phil

master = libtbx.phil.parse("""
selections = "all"
  .type = str
""")
master.objects[0].multiple = True

app = QApplication([])
model = PhilModel()
model.initialize_model(master)

dlg = QDialog()
form = QFormLayout(dlg)
field = PhilField(model, "selections")
# field.widget() is a MultipleWidget wrapping N StrWidgets, with
# [+] [−] [↑] [↓] toolbar buttons. Click [+] to add rows.
form.addRow(field)

# Optionally pre-seed three rows.
field.widget().setValue(["chain A", "chain B", "chain C"])

buttons = QDialogButtonBox(
  QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
form.addRow(buttons)
buttons.accepted.connect(dlg.accept)
buttons.rejected.connect(dlg.reject)

if dlg.exec_() == QDialog.Accepted:
  field.commit()
  print(model.get_phil_extract().selections)
```

### Scope multiples (`RepeatableScopeWidget`)

A *scope* with `.multiple = True` produces a list of sub-extract
instances. In a tree view, instances appear automatically as expandable
sibling subtrees. In a form view, use `RepeatableScopeWidget`:

```python
from PySide2.QtWidgets import QApplication, QDialog, QDialogButtonBox, QVBoxLayout
from qttbx.phil import PhilModel
from qttbx.widgets.phil.multiple import RepeatableScopeWidget
import libtbx.phil

master = libtbx.phil.parse("""
ncs_group {
  selection = "all"
    .type = str
  rotation = 0.0
    .type = float
}
""")
master.objects[0].multiple = True

app = QApplication([])
model = PhilModel()
model.initialize_model(master)

dlg = QDialog()
layout = QVBoxLayout(dlg)
ncs = RepeatableScopeWidget(model, "ncs_group")
# QTabWidget with one tab per ncs_group instance; corner buttons:
#   [+] adds an instance via model.add_scope_instance("ncs_group")
#   [×] removes the active tab via model.remove_scope_instance(...)
# Each tab's body is a QFormLayout of PhilFields constructed with a
# QPersistentModelIndex so they stay bound to the right instance even
# after sibling tabs are removed.
layout.addWidget(ncs)

buttons = QDialogButtonBox(
  QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
layout.addWidget(buttons)
buttons.accepted.connect(dlg.accept)
buttons.rejected.connect(dlg.reject)

if dlg.exec_() == QDialog.Accepted:
  groups = model.get_phil_extract().ncs_group
  for i, g in enumerate(groups):
    print(i, g.selection, g.rotation)
```

### Programmatic model mutation

```python
from PySide2.QtWidgets import QApplication
from qttbx.phil import PhilModel
import libtbx.phil

master = libtbx.phil.parse("""
ncs_group {
  selection = "all"
    .type = str
}
""")
master.objects[0].multiple = True

app = QApplication.instance() or QApplication([])
model = PhilModel()
model.initialize_model(master)

# Add a second instance; qpi survives sibling insertions/removals.
qpi = model.add_scope_instance("ncs_group")
print(qpi.isValid())                              # True
print(len(model.get_phil_extract().ncs_group))    # 2

# Add a third, then remove the middle one. Remaining siblings are
# re-numbered so instance_index stays contiguous (0..N-1).
model.add_scope_instance("ncs_group")
model.remove_scope_instance("ncs_group", instance_index=1)
print(len(model.get_phil_extract().ncs_group))    # 2
```

For **nested** multi-scopes, pass `scope_indices` to disambiguate
intermediate multi segments:

```python
from PySide2.QtWidgets import QApplication
from qttbx.phil import PhilModel
import libtbx.phil

master = libtbx.phil.parse("""
outer {
  ncs_group {
    selection = "all"
      .type = str
  }
}
""")
master.objects[0].multiple = True             # outer is multi
master.objects[0].objects[0].multiple = True  # ncs_group is also multi

app = QApplication.instance() or QApplication([])
model = PhilModel()
model.initialize_model(master)
model.add_scope_instance("outer")              # second outer instance

# Append a new ncs_group under outer instance 1.
model.add_scope_instance("outer.ncs_group", scope_indices=[1])
print(len(model.get_phil_extract().outer[1].ncs_group))  # 2
```

---

## Domain-specific widgets

### Required setup

The PHIL types `space_group`, `unit_cell`, `atom_selection` are registered
in `iotbx.phil.default_converter_registry`, NOT in plain `libtbx.phil`.
Two consequences:

1. Parse the master scope with `iotbx.phil.parse(...)`:

   ```python
   import iotbx.phil
   master = iotbx.phil.parse("""
   crystal {
     space_group = "P 1"
       .type = space_group
     unit_cell = None
       .type = unit_cell
   }
   selection = "all"
     .type = atom_selection
   """)
   ```

2. Import the corresponding widget module(s) before constructing any
   `PhilField` / `PhilItemDelegate` over them — importing the module
   triggers self-registration:

   ```python
   import qttbx.widgets.phil.space_group     # noqa: F401
   import qttbx.widgets.phil.unit_cell       # noqa: F401
   import qttbx.widgets.phil.atom_selection  # noqa: F401
   ```

### `SpaceGroupWidget`

Accepts three input forms:
- Space group **number** as a string: `"19"`, `"4"`.
- **Hermann-Mauguin**, full or short: `"P 21 21 21"`, `"P21"`,
  `"P 1 21 1"`, `"P_21_21_21"`.
- **Hall** symbol (via fallback): `"P 2yb"`, `"-P 2ac 2ab"`.

The widget tries `cctbx.sgtbx.space_group_info(symbol=text)` first; on
`RuntimeError` it falls back to
`cctbx.sgtbx.space_group_info(group=cctbx.sgtbx.space_group(hall_symbol=text))`.
The widget's `value()` is a `cctbx.sgtbx.space_group_info` object (or
`None` / `libtbx.Auto`).

A side label shows a canonical form when input is valid. By default it
shows the full Hermann-Mauguin symbol — the constructor kwarg
`side_label=SIDE_LABEL_HALL` switches it to the canonical Hall symbol.

**Default usage** — the script below opens a dialog with a single
SpaceGroupWidget. Try typing any of `"19"` / `"P21"` / `"P 2yb"` and
watch the side label show the canonical full Hermann-Mauguin form:

```python
from PySide2.QtWidgets import QApplication, QDialog, QFormLayout, QDialogButtonBox
import iotbx.phil
from qttbx.phil import PhilModel
from qttbx.widgets.phil import PhilField
import qttbx.widgets.phil.space_group  # noqa: F401 — triggers registration

master = iotbx.phil.parse("""
crystal {
  space_group = "P 1"
    .type = space_group
}
""")

app = QApplication([])
model = PhilModel()
model.initialize_model(master)

dlg = QDialog()
form = QFormLayout(dlg)
field = PhilField(model, "crystal.space_group")
form.addRow(field)

buttons = QDialogButtonBox(
  QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
form.addRow(buttons)
buttons.accepted.connect(dlg.accept)
buttons.rejected.connect(dlg.reject)

if dlg.exec_() == QDialog.Accepted:
  field.commit()
  sg = model.get_phil_extract().crystal.space_group
  print(sg.type().number())              # e.g. 19
  print(sg.type().lookup_symbol())       # e.g. 'P 21 21 21'
  print(sg.type().hall_symbol())         # e.g. ' P 2ac 2ab'
```

**Hall-labeled variant** — to switch the side label to the Hall symbol,
build `SpaceGroupWidget` directly (with `side_label=SIDE_LABEL_HALL`)
and wire it to the model manually. `PhilField` always uses the registry
default. Same dialog, different side-label content:

```python
from PySide2.QtCore import Qt
from PySide2.QtWidgets import (
  QApplication, QDialog, QFormLayout, QDialogButtonBox,
)
import iotbx.phil
from qttbx.phil import PhilModel
from qttbx.widgets.phil.space_group import (
  SpaceGroupWidget, SIDE_LABEL_HALL,
)

master = iotbx.phil.parse("""
crystal {
  space_group = "P 1"
    .type = space_group
}
""")

app = QApplication([])
model = PhilModel()
model.initialize_model(master)

dlg = QDialog()
form = QFormLayout(dlg)
sg_widget = SpaceGroupWidget(model.get_definition("crystal.space_group"),
                             side_label=SIDE_LABEL_HALL)
form.addRow("Space group:", sg_widget)
# Try typing "P 21" — side label shows " P 2yb" (canonical Hall) instead
# of the default "P 1 21 1" (full Hermann-Mauguin).

# Wire the widget to the model manually (PhilField bypassed).
sg_widget.valueChanged.connect(
  lambda v: model.setData(
    model.index_for_path("crystal.space_group"), v, Qt.EditRole))

buttons = QDialogButtonBox(
  QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
form.addRow(buttons)
buttons.accepted.connect(dlg.accept)
buttons.rejected.connect(dlg.reject)

if dlg.exec_() == QDialog.Accepted:
  sg = model.get_phil_extract().crystal.space_group
  print(sg.type().number(), sg.type().hall_symbol())
```

### `UnitCellWidget`

Accepts six whitespace-separated floats (`"a b c α β γ"`). Validated by
`cctbx.uctbx.unit_cell((...))`; rejects degenerate geometries (zero
axes, non-positive values, impossible angle combinations). Tooltip
shows the interpreted parameters labeled (one per line) when valid.
The widget's `value()` is a `cctbx.uctbx.unit_cell` object (or
`None` / `libtbx.Auto`).

**Default usage** — opens a dialog with a single UnitCellWidget.
Hover over a valid input to see the per-parameter tooltip; type
`"0 60 70 90 90 90"` (zero axis) or `"50 60 70"` (too few values) to
see invalid input drawn with a red border and a tooltip-surfaced
error message.

```python
from PySide2.QtWidgets import QApplication, QDialog, QFormLayout, QDialogButtonBox
import iotbx.phil
from qttbx.phil import PhilModel
from qttbx.widgets.phil import PhilField
import qttbx.widgets.phil.unit_cell  # noqa: F401 — triggers registration

master = iotbx.phil.parse("""
crystal {
  unit_cell = None
    .type = unit_cell
}
""")

app = QApplication([])
model = PhilModel()
model.initialize_model(master)

dlg = QDialog()
form = QFormLayout(dlg)
field = PhilField(model, "crystal.unit_cell")
form.addRow(field)

buttons = QDialogButtonBox(
  QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
form.addRow(buttons)
buttons.accepted.connect(dlg.accept)
buttons.rejected.connect(dlg.reject)

if dlg.exec_() == QDialog.Accepted:
  field.commit()
  uc = model.get_phil_extract().crystal.unit_cell
  if uc is not None:
    print(uc.parameters())               # e.g. (50.0, 60.0, 70.0, 90.0, 95.0, 90.0)
    print("%.4f" % uc.volume())          # e.g. 209200.8866
    print(uc.is_degenerate())            # False
```

### `AtomSelectionWidget` / `AtomSelectionTextWidget`

Currently does only string-shape validation (non-empty, no `$`). Real
syntactic validation against a model is deferred until a DataManager
widget exists. `value()` is a string.

```python
from PySide2.QtWidgets import QApplication, QDialog, QFormLayout, QDialogButtonBox
import iotbx.phil
from qttbx.phil import PhilModel
from qttbx.widgets.phil import PhilField
from qttbx.widgets.phil.atom_selection import AtomSelectionTextWidget
import qttbx.widgets.phil.atom_selection  # noqa: F401 — triggers registration

master = iotbx.phil.parse("""
selection = "all"
  .type = atom_selection
notes = "all"
  .type = atom_selection
""")

app = QApplication([])
model = PhilModel()
model.initialize_model(master)

dlg = QDialog()
form = QFormLayout(dlg)
form.addRow("Selection:", PhilField(model, "selection"))
# Long variant for free-form / multi-line selections (opt-in via widget=):
form.addRow("Notes:",
            PhilField(model, "notes", widget=AtomSelectionTextWidget))

buttons = QDialogButtonBox(
  QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
form.addRow(buttons)
buttons.accepted.connect(dlg.accept)
buttons.rejected.connect(dlg.reject)

if dlg.exec_() == QDialog.Accepted:
  ext = model.get_phil_extract()
  print(repr(ext.selection))
  print(repr(ext.notes))
```

---

## Custom widget registration

Any `PhilWidget` subclass can be registered for a PHIL type. The
example below registers a trivial uppercase-only string widget for
the standard `str` type, then opens a dialog over a master scope that
uses it via the registry. Run the script and notice the new widget
upper-cases everything you type:

```python
from PySide2.QtCore import QSignalBlocker
from PySide2.QtWidgets import (
  QApplication, QDialog, QDialogButtonBox, QFormLayout, QLineEdit, QHBoxLayout,
)
from qttbx.phil import PhilModel
from qttbx.widgets.phil import PhilField, PhilWidget, register_widget
import libtbx.phil


class UppercaseStrWidget(PhilWidget):
  """Trivial demo widget: a QLineEdit that upper-cases its content."""

  def __init__(self, definition, parent=None):
    super().__init__(definition, parent)
    self._edit = QLineEdit(self)
    layout = QHBoxLayout(self)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.addWidget(self._edit)
    self._edit.editingFinished.connect(
      lambda: self.valueChanged.emit(self.value()))

  def value(self):
    return self._edit.text().upper()

  def setValue(self, value):
    new = "" if value is None else str(value).upper()
    if new == self._edit.text():
      return
    blocker = QSignalBlocker(self._edit)
    self._edit.setText(new)
    del blocker

  def isValid(self):
    return True

  def errorString(self):
    return ""


# Override the registered "str" widget. Last-write-wins.
register_widget("str", UppercaseStrWidget)

master = libtbx.phil.parse("""
title = "untitled"
  .type = str
""")

app = QApplication([])
model = PhilModel()
model.initialize_model(master)

dlg = QDialog()
form = QFormLayout(dlg)
field = PhilField(model, "title")
form.addRow("Title:", field)

buttons = QDialogButtonBox(
  QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
form.addRow(buttons)
buttons.accepted.connect(dlg.accept)
buttons.rejected.connect(dlg.reject)

if dlg.exec_() == QDialog.Accepted:
  field.commit()
  print(model.get_phil_extract().title)   # always upper-case
```

After registration (typically at module bottom, like the v4 domain
widgets), every `PhilField` / `PhilItemDelegate` over a `my_exotic_type`
definition will use `MyExoticTypeWidget` automatically.

`register_widget` is last-write-wins, so site-specific code can override
a built-in by importing AFTER `qttbx.widgets.phil` has loaded.

For consumers building forms over uncertain PHIL trees, `widget_for_definition`
accepts an opt-in `fallback=` kwarg that constructs the given class for any
unregistered `phil_type` (rather than raising). The default behavior — raise
`ValueError` listing the registered types — is preferred when the master scope
is known, because silent fallback can mask typos and missing module imports.

---

## `PhilModel` API reference

Most consumers work with `PhilField` / `PhilItemDelegate` and never call
`PhilModel` directly. These are the methods that user code occasionally
needs:

| Method | Purpose |
|---|---|
| `initialize_model(scope)` | Required once after construction. Builds the tree and seeds the extract. |
| `get_phil_extract()` | Current value tree (a `libtbx.phil.scope_extract`). Read-only access for reading the params. |
| `get_master_phil()` | The original master scope (template). |
| `restore_defaults()` | Reset the working extract to the master defaults. Equivalent to `initialize_model(get_master_phil())`. Useful for "Reset to defaults" dialog buttons. |
| `get_working_phil(diff_only=True)` | A `libtbx.phil.scope` with the current values. `diff_only=True` returns only the non-default settings. |
| `get_phil_extract_value(full_path)` | Read a value by dotted path. For multi-scope paths, pass an indexed list (see below). |
| `set_phil_extract_value(full_path, value)` | Write a value by dotted path. |
| `index_for_path(full_path) -> QModelIndex` | Resolve a path to a Qt model index. Raises if the path traverses a multi-scope. |
| `persistent_index_for_path(full_path, scope_indices=None) -> QPersistentModelIndex` | Like above, but returns a persistent index and accepts scope indices to disambiguate multi-scope paths. |
| `add_scope_instance(full_path, scope_indices=None) -> QPersistentModelIndex` | Append a new instance to a multi-scope. Returns the new instance's persistent index. |
| `remove_scope_instance(full_path, instance_index, scope_indices=None)` | Remove the indexed instance. Re-numbers remaining siblings. |

### Indexed-path form for `get/set_phil_extract_value`

For paths that traverse a multi-scope, pass a list of segments where each
multi-scope segment is a `(name, instance_index)` tuple:

```python
from PySide2.QtWidgets import QApplication
from qttbx.phil import PhilModel
import libtbx.phil

master = libtbx.phil.parse("""
refinement {
  macro_cycles = 3
    .type = int
}
ncs_group {
  selection = "all"
    .type = str
}
""")
master.objects[1].multiple = True   # ncs_group is multi

app = QApplication.instance() or QApplication([])
model = PhilModel()
model.initialize_model(master)
model.add_scope_instance("ncs_group")    # second instance

# Single-instance dotted path:
print(model.get_phil_extract_value("refinement.macro_cycles"))    # 3

# Path under a multi-scope (the second ncs_group instance):
model.set_phil_extract_value([("ncs_group", 1), "selection"], "chain A")
print(model.get_phil_extract_value([("ncs_group", 1), "selection"]))  # chain A
```

(`PhilItem.extract_path()` produces exactly this form for a given tree
node.)

---

## `PhilField` API reference

```python
PhilField(model, full_path=None, widget=None, parent=None,
  persistent_index=None)
```

- `full_path` and `persistent_index` are **mutually exclusive**; exactly
  one must be given. For `RepeatableScopeWidget`-built sub-forms, pass
  `persistent_index=` so the field stays bound to the right scope instance
  across sibling insertions/removals.
- `widget=` — override the default class returned by the registry (e.g.
  to use a `*TextWidget` long variant in a form view).
- `parent=` — Qt parent widget.

| Method / signal | Purpose |
|---|---|
| `widget()` | The inner `PhilWidget`. |
| `isValid()` / `errorString()` | Forwarded from the inner widget. |
| `commit()` | Force a pending edit to land in the model **now** (no-op if invalid). Call before reading the model from a "Run" or "OK" handler. |
| `setEnabled(b)` / `setReadOnly(ro)` | Forwarded to the inner widget. |
| `validityChanged: Signal(bool)` | Forwarded from the inner widget. Wire to a `QDialogButtonBox` to gate the OK button. |

---

## Threading

`PhilModel`, `PhilField`, `PhilWidget`, `PhilItemDelegate`, `MultipleWidget`,
and `RepeatableScopeWidget` are **GUI-thread-only**. Worker threads that
need to mutate parameters must marshal the call onto the GUI thread:

```python
import threading
from PySide2.QtCore import QObject, Signal
from PySide2.QtWidgets import QApplication
from qttbx.phil import PhilModel
import libtbx.phil

master = libtbx.phil.parse("""
refinement {
  weight = 0.5
    .type = float
}
""")

app = QApplication([])
model = PhilModel()
model.initialize_model(master)


# A QObject whose signal carries the cross-thread payload. The signal is
# emitted from the worker thread; Qt's AutoConnection queues delivery to
# the slot on the bridge's home thread (the GUI thread).
class GuiBridge(QObject):
  apply_weight = Signal(float)


bridge = GuiBridge()


def _on_apply_weight(value):
  model.set_phil_extract_value("refinement.weight", value)
  print("weight ->", model.get_phil_extract().refinement.weight)
  app.quit()


bridge.apply_weight.connect(_on_apply_weight)


def worker():
  # Pretend we computed a new weight off-thread.
  new_value = 0.7
  # Emit a signal from any thread; Qt queues delivery to the GUI thread
  # because bridge lives there. (Direct PhilModel.setData calls from a
  # worker thread would trip the GUI-thread assertion.)
  bridge.apply_weight.emit(new_value)


threading.Thread(target=worker, daemon=True).start()
app.exec_()
```

`PhilModel.setData`, `add_scope_instance`, and `remove_scope_instance` carry
a development-time assertion that fires immediately if invoked off the GUI
thread.

---

## Patterns and gotchas

**Edits don't reach the model until `commit()`** when the user clicks
OK without leaving the focused widget. Always loop `f.commit()` over your
fields in `QDialog.accept()`, or call `view.setCurrentIndex(QModelIndex())`
on a tree view, before reading the model.

**`Auto` and `None` are runtime form decisions, not PHIL attributes.**
Use `widget.setAllowNone(True)` (default) and `widget.setUseAuto(True)`
to opt into either sentinel. The widgets emit `None` / `libtbx.Auto`
from `value()` accordingly when input is empty.

**Float comparisons:** use `libtbx.test_utils.approx_equal`, not `==`.

**Space group equality:** never compare via `str(...)`. Two HM/Hall
symbols can describe the same group. Compare with
`sg1.group() == sg2.group()`.

**Long-text widgets are NOT registered.** Pass them via `widget=` to
`PhilField` for the form-view path. The tree-view delegate always uses
the registered short variant.

**Domain widgets self-register on import.** A `PhilField` over an
`atom_selection` definition will fail with `ValueError` if the consumer
forgot to import `qttbx.widgets.phil.atom_selection` first. The error
message lists the registered phil_types and suggests the missing import.

**`MultipleWidget` initial value:** for definitions with a non-None
scalar default (e.g. `val = 1\n  .type = int\n  .multiple = True`),
the widget seeds one row showing the default. None default → zero rows.

**Tree views show one row per multi-scope instance.** No special widget
is needed; expanding the row reveals each instance's children.
`RepeatableScopeWidget` is only for form-view consumers.

---

## Testing your GUI

`qttbx/regression/tst_phil_widgets.py` is the canonical test file. The
end-to-end exercises (`exercise_v1_tree_and_form_sync`,
`exercise_v2_tree_and_form_sync`, `exercise_v3_definition_multi_and_scope_multi`,
`exercise_v4_domain_widgets_in_form`) are useful templates for your own
GUI-level integration tests. They run headless under `QT_QPA_PLATFORM=offscreen`
and are registered in `qttbx/run_tests.py`.

Run with:

```sh
libtbx.python qttbx/regression/tst_phil_widgets.py
# or:
libtbx.run_tests_parallel module=qttbx nproc=2 run_in_unique_dirs=True
```

---

## Module layout summary

```
qttbx/
├── phil.py                       # PhilModel + PhilItem (multi-instance aware)
├── widgets/
│   └── phil/
│       ├── __init__.py           # PhilWidget base, PhilField, registry, register_widget
│       ├── delegate.py           # PhilItemDelegate (tree-view controller)
│       ├── text_base.py          # ValidatedLineEdit, ValidatedTextEdit
│       ├── int_widget.py         # IntWidget
│       ├── bool_widget.py        # BoolWidget
│       ├── float_widget.py       # FloatWidget
│       ├── key_widget.py         # KeyWidget
│       ├── choice_widget.py      # ChoiceWidget, ChoiceMultiWidget
│       ├── str_widget.py         # StrWidget, StrTextWidget
│       ├── qstr_widget.py        # QstrWidget, QstrTextWidget
│       ├── path_widget.py        # PathWidget (browse), PathsTextWidget
│       ├── strings_widget.py     # StringsWidget, StringsTextWidget
│       ├── words_widget.py       # WordsWidget, WordsTextWidget
│       ├── ints_widget.py        # IntsWidget, IntsTextWidget
│       ├── floats_widget.py      # FloatsWidget, FloatsTextWidget
│       ├── multiple.py           # MultipleWidget, RepeatableScopeWidget
│       ├── space_group.py        # SpaceGroupWidget (self-registers)
│       ├── unit_cell.py          # UnitCellWidget (self-registers)
│       └── atom_selection.py     # AtomSelectionWidget, AtomSelectionTextWidget (self-registers)
└── regression/
    └── tst_phil_widgets.py       # ~115 exercise functions; integration tests
```
