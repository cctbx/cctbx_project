"""DataManagerItemDelegate -- custom paint for the Used-for and Delete
columns. Stale rows are rendered with a pale-pink background matching
existing qttbx convention (text_base.py uses #fff5f5)."""

from PySide2.QtCore import Qt, QSize, QRect, QEvent, Signal
from PySide2.QtGui import QColor, QBrush, QFontMetrics, QPalette
from PySide2.QtWidgets import QStyledItemDelegate

from qttbx.widgets.data_manager._table_model import (
  DataManagerTableModel, _trash_icon)


_CHIP_PAD_H = 8
_CHIP_PAD_V = 3
_CHIP_GAP_V = 4
_CHIP_RADIUS = 8
_ADD_LABEL = "+"


def _is_dark(palette):
  """True if the palette is a dark theme.

  Decision based on the Base role's lightness, matching the convention
  used by other qttbx widgets that need dark-mode-aware colors.
  """
  return palette.color(QPalette.Base).lightness() < 128


def _chip_background(palette):
  """Background brush for a normal Used-for chip.

  Uses the palette's Button role so the chip naturally tracks light/dark
  themes (light grey on a light theme, dark grey on a dark theme)."""
  return palette.color(QPalette.Button)


def _chip_text_color(palette):
  """Foreground color for chip text and the trailing ✕."""
  return palette.color(QPalette.ButtonText)


def _stale_background(palette):
  """Background color for stale rows -- pale pink in light mode,
  desaturated dark red in dark mode."""
  if _is_dark(palette):
    return QColor(64, 28, 28)
  return QColor("#fff5f5")


class DataManagerItemDelegate(QStyledItemDelegate):
  """Paint + interact for Used-for and Delete columns.

  Used-for column: each binding is a rounded chip on its own line with
  a trailing unbind glyph. A final "+ add" pseudo-chip opens the
  binding popup (signalled via ``addChipClicked``).

  Delete column: a centered delete glyph. Click emits ``deleteClicked``.

  Signals
  -------
  unbindClicked : (int, str)
    Row index and phil_path of the chip that was clicked to unbind.
  addChipClicked : (int, QRect)
    Row index and the anchor rect (in viewport coords) for the popup.
  deleteClicked : (int)
    Row index whose delete button was clicked.
  """

  unbindClicked = Signal(int, str)
  addChipClicked = Signal(int, QRect)
  deleteClicked = Signal(int)

  def __init__(self, parent=None):
    QStyledItemDelegate.__init__(self, parent)

  def sizeHint(self, option, index):
    """Return the size hint for ``index``.

    For the Used-for column, the height is computed as
    ``n_lines * line_h + (n_lines - 1) * gap`` where ``n_lines`` is the
    number of chips, with a minimum of one line so the ``+`` button
    is always visible at the top-right of the cell. When the row's
    data type has no compatible PHIL parameters (spec section 5: "if
    the file's data type has no compatible PHIL parameters, the cell
    is empty -- no chips, no ``+``"), the cell uses base sizing.
    """
    col = index.column()
    if col == DataManagerTableModel.COL_USED_FOR:
      model = index.model()
      chips = model.used_for_with_labels(index.row())
      has_add = model.has_compatible_params(index.row())
      n_chips = len(chips)
      if not has_add and n_chips == 0:
        return QStyledItemDelegate.sizeHint(self, option, index)
      # The ``+`` button shares the top line with the first chip (or
      # occupies the top line alone when there are no chips), so the
      # cell needs at least one line and never more than ``n_chips``
      # lines (or 1 if ``has_add`` and no chips).
      n_lines = max(n_chips, 1) if has_add else n_chips
      fm = QFontMetrics(option.font)
      line_h = fm.height() + 2 * _CHIP_PAD_V
      total = n_lines * line_h + max(0, n_lines - 1) * _CHIP_GAP_V + 4
      return QSize(option.rect.width() or 200, total)
    return QStyledItemDelegate.sizeHint(self, option, index)

  def paint(self, painter, option, index):
    """Paint ``index`` into ``painter``.

    Stale rows get a pale-pink background fill before the column-
    specific paint runs. The Used-for and Delete columns are painted
    fully here; other columns fall through to the default
    :class:`QStyledItemDelegate` paint.
    """
    model = index.model()
    is_stale = (model is not None
                and hasattr(model, "is_stale")
                and model.is_stale(index.row()))
    if is_stale:
      self._fill_stale_background(painter, option)
    col = index.column()
    if col == DataManagerTableModel.COL_USED_FOR:
      self._paint_used_for(painter, option, index)
      return
    if col == DataManagerTableModel.COL_DELETE:
      self._paint_delete(painter, option)
      return
    QStyledItemDelegate.paint(self, painter, option, index)

  def _fill_stale_background(self, painter, option):
    painter.fillRect(option.rect, QBrush(_stale_background(option.palette)))

  def _paint_used_for(self, painter, option, index):
    model = index.model()
    chips = model.used_for_with_labels(index.row())
    rect = option.rect.adjusted(2, 2, -2, -2)
    fm = QFontMetrics(option.font)
    line_h = fm.height() + 2 * _CHIP_PAD_V
    chip_bg = _chip_background(option.palette)
    chip_fg = _chip_text_color(option.palette)
    # Spec section 5: suppress the ``+`` button when the row's data
    # type has no compatible PHIL parameters.
    has_add = model.has_compatible_params(index.row())
    add_rect = None
    if has_add:
      add_w = fm.horizontalAdvance(_ADD_LABEL) + 2 * _CHIP_PAD_H
      # The ``+`` button always lives on the top row, anchored at the
      # right edge of the cell so it stays visible regardless of how
      # many chips are stacked below.
      add_rect = QRect(rect.right() - add_w + 1, rect.top(),
                       add_w, line_h)
    # Each chip drawn on its own line starting from the top. When the
    # ``+`` button is on the right, clip chip widths so they don't
    # overlap it.
    y = rect.top()
    available_w = rect.width()
    if add_rect is not None:
      available_w = max(0, add_rect.left() - rect.left() - _CHIP_GAP_V)
    for _phil_path, label in chips:
      label_text = u"%s ✕" % str(label)
      desired_w = fm.horizontalAdvance(label_text) + 2 * _CHIP_PAD_H
      chip_w = min(desired_w, available_w)
      chip_rect = QRect(rect.left(), y, chip_w, line_h)
      painter.save()
      painter.setBrush(QBrush(chip_bg))
      painter.setPen(Qt.NoPen)
      painter.drawRoundedRect(chip_rect, _CHIP_RADIUS, _CHIP_RADIUS)
      painter.setPen(chip_fg)
      painter.drawText(
        chip_rect.adjusted(_CHIP_PAD_H, 0, -_CHIP_PAD_H, 0),
        Qt.AlignVCenter | Qt.AlignLeft, label_text)
      painter.restore()
      y += line_h + _CHIP_GAP_V
    if add_rect is not None:
      painter.save()
      painter.setBrush(QBrush(chip_bg))
      painter.setPen(Qt.NoPen)
      painter.drawRoundedRect(add_rect, _CHIP_RADIUS, _CHIP_RADIUS)
      painter.setPen(chip_fg)
      painter.drawText(add_rect, Qt.AlignVCenter | Qt.AlignCenter,
                       _ADD_LABEL)
      painter.restore()

  def _paint_delete(self, painter, option):
    """Paint the system's standard trash icon centered in the cell.

    Uses :attr:`QStyle.SP_TrashIcon`, the native platform glyph (Finder
    trash on macOS, recycle bin on Windows, the active icon theme on
    Linux). Available since Qt 5.12.

    Falls back to a ``✕`` text glyph when the running Qt build doesn't
    expose ``SP_TrashIcon`` or returns a null icon (older Qt or a
    style without the standard pixmap).
    """
    icon = _trash_icon()
    if icon is None:
      painter.save()
      painter.setPen(option.palette.windowText().color())
      painter.drawText(option.rect, Qt.AlignCenter, u"✕")
      painter.restore()
      return
    # Size the icon to roughly the font's ascent so it tracks DPI and
    # font size automatically; clamp to a sensible minimum for tiny
    # fonts.
    fm = QFontMetrics(option.font)
    size = max(fm.ascent() + 2, 12)
    cx = option.rect.center().x()
    cy = option.rect.center().y()
    icon_rect = QRect(cx - size // 2, cy - size // 2, size, size)
    icon.paint(painter, icon_rect, Qt.AlignCenter)

  def editorEvent(self, event, model, option, index):
    """Translate mouse-release events into the delegate's signals.

    Clicks on a chip emit ``unbindClicked``; a click on the "+ add"
    pseudo-chip emits ``addChipClicked``; a click on the Delete column
    emits ``deleteClicked``. All other events fall through to the
    default :class:`QStyledItemDelegate.editorEvent` implementation.
    """
    if event.type() == QEvent.MouseButtonRelease:
      col = index.column()
      if col == DataManagerTableModel.COL_USED_FOR:
        return self._handle_used_for_click(event, option, index)
      if col == DataManagerTableModel.COL_DELETE:
        self.deleteClicked.emit(index.row())
        return True
    return QStyledItemDelegate.editorEvent(self, event, model, option, index)

  def _handle_used_for_click(self, event, option, index):
    model = index.model()
    chips = model.used_for_with_labels(index.row())
    rect = option.rect.adjusted(2, 2, -2, -2)
    fm = QFontMetrics(option.font)
    line_h = fm.height() + 2 * _CHIP_PAD_V
    # ``+`` button hit-test first: it sits at the top-right of the cell
    # and overlaps the first chip's row, so the click could fall on
    # either. Spec section 5: no ``+`` when no compatible PHIL params.
    has_add = model.has_compatible_params(index.row())
    add_rect = None
    if has_add:
      add_w = fm.horizontalAdvance(_ADD_LABEL) + 2 * _CHIP_PAD_H
      add_rect = QRect(rect.right() - add_w + 1, rect.top(),
                       add_w, line_h)
      if add_rect.contains(event.pos()):
        self.addChipClicked.emit(index.row(), add_rect)
        return True
    # Then chips on subsequent lines starting from the top, clipped on
    # the right so they don't overlap the ``+`` button.
    y = rect.top()
    available_w = rect.width()
    if add_rect is not None:
      available_w = max(0, add_rect.left() - rect.left() - _CHIP_GAP_V)
    for phil_path, label in chips:
      label_text = u"%s ✕" % str(label)
      desired_w = fm.horizontalAdvance(label_text) + 2 * _CHIP_PAD_H
      chip_w = min(desired_w, available_w)
      chip_rect = QRect(rect.left(), y, chip_w, line_h)
      if chip_rect.contains(event.pos()):
        self.unbindClicked.emit(index.row(), phil_path)
        return True
      y += line_h + _CHIP_GAP_V
    return False
