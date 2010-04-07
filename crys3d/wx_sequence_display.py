
from __future__ import division
import wx
import wx.lib.wordwrap
from scitbx.array_family import flex, shared
import math, sys, os

class SequencePanel (wx.PyPanel) :
  tooltip = "Double-click a residue to select it; hold down Shift to select \
multiple residues."
  __bg_color = (255,255,255)
  def __init__ (self, *args, **kwds) :
    wx.PyPanel.__init__(self, *args, **kwds)
    if self.__bg_color is not None :
      self.SetBackgroundColour(self.__bg_color)
    self.sequence = ""
    self.line_width = 50
    self.line_sep = 28
    self.start_offset = 0
    self.char_boxes = shared.stl_set_unsigned()
    self.flag_show_line_numbers = True
    self.flag_enable_selections = True
    self.flag_overwrite_selections = True
    self.flag_show_selections = True
    self.flag_show_tooltip = False
    self.highlights = []
    self.highlight_colors = []
    self.selected_residues = flex.bool()
    self.selection_color = (255, 255, 0)
    self._last_x = None
    self._last_y = None
    self.txt_font = wx.Font(14, wx.MODERN, wx.NORMAL, wx.NORMAL)
    self.Bind(wx.EVT_PAINT, self.OnPaint)
    self.Bind(wx.EVT_LEFT_DCLICK, self.OnDoubleClick)
    self.Bind(wx.EVT_LEFT_DOWN, self.OnLeftDown)
    self.Bind(wx.EVT_LEFT_UP, self.OnLeftUp)
    self.Bind(wx.EVT_RIGHT_DOWN, self.OnRightDown)
    self.Bind(wx.EVT_RIGHT_UP, self.OnRightUp)
    self.tip = wx.ToolTip(self.tooltip)
    self.SetToolTip(self.tip)
    self.tip.Enable(self.flag_show_tooltip)

  def enable_tooltip (self, enable=True) :
    self.tip.Enable(enable)

  def OnPaint (self, event) :
    dc = wx.AutoBufferedPaintDCFactory(self)
    gc = wx.GraphicsContext.Create(dc)
    self.paint(gc)

  def paint (self, gc) :
    self.paint_sequence(gc)

  def paint_sequence (self, gc) :
    gc.SetBrush(wx.TRANSPARENT_BRUSH)
    gc_txt_font = gc.CreateFont(self.txt_font, (0,0,0))
    gc.SetFont(gc_txt_font)
    i = 0
    ypos = self.line_sep
    xpos = 16
    black_pen = wx.Pen((0,0,0), 1)
    (char_w, char_h) = self.get_char_size(gc)
    if self.flag_show_selections :
      for j, (i_start, i_end) in enumerate(self.get_selected_ranges()) :
        gc.SetPen(wx.Pen((0,0,0), 1))
        gc.SetBrush(wx.Brush(self.selection_color))
        self.draw_box_around_residues(i_start, i_end, gc)
    gc.SetBrush(wx.TRANSPARENT_BRUSH)
    while i < len(self.sequence) :
      j = i
      line_start = self.start_offset + i + 1
      line_end = line_start + self.line_width - 1
      line = self.sequence[i:i+self.line_width]
      line_w = self.line_width * char_w
      label_w = gc.GetTextExtent("X"*8)[0]
      if self.flag_show_line_numbers :
        gc.DrawText("%6d  " % line_start, xpos, ypos)
        gc.DrawText("  %-6d" % line_end, xpos + label_w + line_w, ypos)
        xpos += label_w
      for j, char in enumerate(line) :
        i_seq = i+j
        if i_seq in self.highlights :
          color = self.highlight_colors[self.highlights.index(i_seq)]
          gc.SetFont(gc.CreateFont(self.txt_font, color))
        gc.DrawText(char, xpos, ypos)
        gc.SetFont(gc_txt_font)
        xpos += char_w
      i += self.line_width
      ypos += self.line_sep
      xpos = 16
    return gc

  def draw_box_around_residues (self, i_start, i_end, gc) :
    ranges = self.get_contiguous_ranges(i_start, i_end)
    char_w, char_h = self.get_char_size(gc)
    for (seg_start, seg_end) in ranges :
      (x1, y1, nx, ny) = self.get_char_position(seg_start)
      (nx, ny, x2, y2) = self.get_char_position(seg_end)
      line = gc.CreatePath()
      line.MoveToPoint(x1 - 1, y1 - 1)
      line.AddLineToPoint(x2 - 1, y1 - 1)
      line.AddLineToPoint(x2 - 1, y2 - 1)
      line.AddLineToPoint(x1 - 1, y2 - 1)
      line.CloseSubpath()
      gc.PushState()
      gc.FillPath(line)
      gc.StrokePath(line)
      gc.PopState()
      #gc.DrawRectangle(x1, y1, x2 - x1, y2 - y1)

  def DoGetBestSize (self) :
    dc = wx.ClientDC(self)
    dc.SetFont(self.txt_font)
    i = 0
    (panel_w, panel_h) = (32, 32)
    char_w, char_h = self.get_char_size(dc)
    line_w = char_w * self.line_width
    if self.flag_show_line_numbers :
      line_w += dc.GetTextExtent("X" * 16)[0]
    panel_w += line_w
    n_lines = int(math.ceil(len(self.sequence) / self.line_width))
    panel_h += n_lines * self.line_sep
    return (max(480, panel_w), max(240, panel_h))

  def set_sequence (self, seq) :
    self.sequence = "".join(seq.splitlines())
    self.build_boxes()
    self.clear_highlights()
    self.clear_selection()

  def get_char_size (self, dc=None) :
    if dc is None :
      dc = wx.ClientDC(self)
      dc.SetFont(self.txt_font)
    char_w, char_h = dc.GetTextExtent("X")
    char_w = max(10, char_w)
    char_h = max(16, char_h)
    return (char_w, char_h)

  def build_boxes (self) :
    dc = wx.ClientDC(self)
    dc.SetFont(self.txt_font)
    char_w, char_h = self.get_char_size(dc)
    x_start = 16
    y_start = self.line_sep
    if self.flag_show_line_numbers :
      x_start += dc.GetTextExtent("X" * 8)[0]
    self.char_boxes = []
    for i_seq in range(len(self.sequence)) :
      lines = int(math.floor((i_seq-self.start_offset) / self.line_width))
      y = y_start + (lines * self.line_sep)
      n_prev_chars = (i_seq - self.start_offset) % self.line_width
      x = x_start + (char_w * n_prev_chars)
      self.char_boxes.append((x, y, x + char_w - 1, y + char_h - 1))

  def get_char_position (self, i_seq) :
    return tuple(self.char_boxes[i_seq])

  def get_contiguous_ranges (self, i_start, i_end) :
    assert i_start < i_end or i_start == i_end
    (char_w, char_h) = self.get_char_size()
    line_start = int(math.floor((i_start-self.start_offset) / self.line_width))
    line_end = int(math.floor((i_end - self.start_offset) / self.line_width))
    ranges = []
    seg_start = i_start
    while line_start < line_end :
      seg_end = ((line_start + 1) * self.line_width) - 1
      ranges.append((seg_start, seg_end))
      line_start += 1
      seg_start = seg_end + 1
    ranges.append((seg_start, i_end))
    return ranges

  def get_selected_ranges (self) :
    ranges = []
    last_start = None
    last_end = None
    for i_seq in range(self.selected_residues.size()) :
      if self.selected_residues[i_seq] :
        if last_end is not None :
          if last_end < (i_seq - 1) :
            ranges.append((last_start, last_end))
            last_start = i_seq
        else :
          last_start = i_seq
        last_end = i_seq
      elif last_end is not None :
        ranges.append((last_start, last_end))
        last_start = None
        last_end = None
    return ranges

  def is_on_char (self, x, y) :
    return None # TODO

  def clear_highlights (self) :
    self.highlights = []
    self.highlight_colors = []

  def highlight_char (self, i_seq, color=(255,0,0)) :
    if i_seq in self.highlights :
      i = self.highlights.index(i_seq)
      self.highlights.pop(i)
      self.highlight_colors.pop(i)
    self.highlights.append(i_seq)
    self.highlight_colors.append(color)

  def clear_selection (self) :
    self.selected_residues = flex.bool(len(self.sequence), False)
    frame.statusbar.SetStatusText("")

  def select_chars (self, i_start, i_end, box=True) :
    assert i_end < len(self.sequence) and i_start != i_end
    for i_seq in range(i_start, i_end + 1) :
      self.selected_residues[i_seq] = box
    self.update_frame()

  def update_frame (self) :
    frame = self.GetParent()
    ranges = self.get_selected_ranges()
    if len(ranges) == 0 :
      txt = ""
    elif len(ranges) == 1 and ranges[0][0] == ranges[0][1] :
      txt = "SELECTED: residue %d" %  ranges[0][0]
    else :
      txt_ranges = []
      for (x, y) in ranges :
        if x == y : txt_ranges.append(str(x))
        else : txt_ranges.append("%d-%d" % (x, y))
      txt = "SELECTED: residues %s" % ", ".join(txt_ranges)
    frame.statusbar.SetStatusText(txt)

  def deselect_chars (self, i_start, i_end) :
    self.select_chars(i_start, i_end, False)

  def select_residue (self, x, y) :
    for i_seq, box in enumerate(self.char_boxes) :
      (x1, y1, x2, y2) = box
      if (x > x1 and x < x2) and (y > y1 and y < y2) :
        if self.selected_residues[i_seq] :
          self.selected_residues[i_seq] = False
        else :
          self.selected_residues[i_seq] = True
        self.update_frame()
        return True
    return False

  def OnClear (self, event) :
    self.clear_selection()
    self.Refresh()

  def OnSetMode (self, event) :
    self.flag_overwrite_selections = event.GetEventObject().GetValue()

  #---------------------------------------------------------------------
  # MOUSE EVENTS
  def clear_mouse (self) :
    self._last_x = None
    self._last_y = None

  def record_mouse (self, x, y) :
    self._last_x = None
    self._last_y = None

  def OnDoubleClick (self, event) :
    (x, y) = (event.GetX(), event.GetY())

  def OnLeftUp (self, event) :
    (x, y) = (event.GetX(), event.GetY())

  def OnLeftDown (self, event) :
    (x, y) = (event.GetX(), event.GetY())

  def OnRightDown (self, event) :
    (x, y) = (event.GetX(), event.GetY())
    if self.flag_show_tooltip :
      pass

  def OnRightUp (self, event) :
    (x, y) = (event.GetX(), event.GetY())

#-----------------------------------------------------------------------
class SequenceWithStructurePanel (SequencePanel) :
  tooltip = """\
Double-click on any residue or secondary-structure element to select the \
residue(s).  Holding down shift enables multiple selections."""
  def __init__ (self, *args, **kwds) :
    SequencePanel.__init__(self, *args, **kwds)
    self.line_sep = 64
    self.structure = ""
    self.selected_helices = []
    self.selected_strands = []
    self.selected_linkers = []
    self.flag_allow_select_structure = True

  def paint (self, gc) :
    self.paint_sequence(gc)
    self.paint_structure(gc)

  def paint_structure (self, gc) :
    helices = self.get_helices()
    helix_pen = wx.Pen('red', 2)
    h_helix_pen = wx.Pen((255,50,0), 2)
    for k, (i_start, i_end) in enumerate(helices) :
      color_start = (255, 50, 50)
      color_end = (255, 150, 150)
      bounds = self.get_region_bounds(i_start, i_end)
      if k in self.selected_helices :
        color_start = (200, 100, 100)
        color_end = (255, 200, 200)
        gc.SetPen(h_helix_pen)
      else :
        gc.SetPen(helix_pen)
      for box in bounds :
        (x1, y1, x2, y2) = box
        #gc.DrawRoundedRectangle(x1, y1, x2-x1, 16, 4)
        helix_brush = gc.CreateLinearGradientBrush(x1, y1, x1, y2, color_start,
          color_end)
        gc.SetBrush(helix_brush)
        gc.DrawRectangle(x1, y1, x2-x1, 16)
      gc.SetPen(helix_pen)
    strands = self.get_strands()
    strand_pen = wx.Pen((0, 100, 255), 2)
    h_strand_pen = wx.Pen('blue', 2)
    for k, (i_start, i_end) in enumerate(strands) :
      color_start = (0, 100, 255)
      color_end = (50, 200, 255)
      bounds = self.get_region_bounds(i_start, i_end)
      if k in self.selected_strands :
        color_start = (0, 200, 200)
        color_end = (50, 255, 255)
        gc.SetPen(h_strand_pen)
      else :
        gc.SetPen(strand_pen)
      for j, box in enumerate(bounds) :
        (x1, y1, x2, y2) = box
        strand_brush = gc.CreateLinearGradientBrush(x1, y1 + 4, x1, y2-4,
          color_start, color_end)
        gc.SetBrush(strand_brush)
        if j == (len(bounds) - 1) :
          segments = strand_as_arrow(x1, y1, x2, y2)
        else :
          segments = strand_as_box(x1, y1, x2, y2)
        line = gc.CreatePath()
        line.MoveToPoint(segments[0][0], segments[0][1])
        for segment in segments[1:] :
          line.AddLineToPoint(segment[0], segment[1])
        line.CloseSubpath()
        gc.PushState()
        gc.StrokePath(line)
        gc.FillPath(line)
        gc.PopState()
      gc.SetPen(strand_pen)
    missing = self.get_missing()
    missing_pen = wx.Pen((150, 150, 150), 4, style=wx.SHORT_DASH)
    for k, (i_start, i_end) in enumerate(missing) :
      bounds = self.get_region_bounds(i_start, i_end)
      for j, box in enumerate(bounds) :
        (x1, y1, x2, y2) = box
        line = gc.CreatePath()
        line.MoveToPoint(x1, y1 + 8)
        line.AddLineToPoint(x2, y1 + 8)
        line.CloseSubpath()
        gc.PushState()
        gc.SetPen(missing_pen)
        gc.StrokePath(line)
        gc.PopState()
    other = self.get_linkers()
    other_pen = wx.Pen((0,0,0), 4)
    h_other_pen = wx.Pen((255, 255, 0), 12)
    for k, (i_start, i_end) in enumerate(other) :
      bounds = self.get_region_bounds(i_start, i_end)
      for j, box in enumerate(bounds) :
        (x1, y1, x2, y2) = box
        line = gc.CreatePath()
        line.MoveToPoint(x1, y1 + 8)
        line.AddLineToPoint(x2, y1 + 8)
        line.CloseSubpath()
        if k in self.selected_linkers :
          gc.PushState()
          gc.SetPen(h_other_pen)
          gc.StrokePath(line)
          gc.PopState()
        gc.PushState()
        gc.SetPen(other_pen)
        gc.StrokePath(line)
        gc.PopState()

  def set_structure (self, ss) :
    self.structure = "".join(ss.splitlines())
    assert self.sequence == "" or len(self.sequence) == len(self.structure)
    self.clear_highlights()

  def apply_missing_residue_highlights (self) :
    ranges = self.get_missing()
    for i_start, i_end in ranges :
      for i_seq in range(i_start, i_end + 1) :
        self.highlights.append(i_seq)
        self.highlight_colors.append((150,150,150))

  def get_region_bounds (self, i_start, i_end) :
    (char_w, char_h) = self.get_char_size()
    ranges = self.get_contiguous_ranges(i_start, i_end)
    boxes = []
    for (seg_start, seg_end) in ranges :
      (x1, y1, nx, ny) = self.get_char_position(seg_start)
      y1 -= 24
      (nx, y2, x2, ny) = self.get_char_position(seg_end)
      y2 -= 8
      boxes.append((x1, y1, x2, y2))
    return boxes

  def get_sec_str_motif (self, code) :
    sections = []
    current_start = None
    for i, char in enumerate(self.structure) :
      if char == code :
        if current_start is None :
          current_start = i
      else :
        if current_start is not None :
          sections.append((current_start, i-1))
          current_start = None
    if current_start is not None :
      sections.append((current_start, len(self.structure) - 1))
    return sections

  def get_helices (self) :
    return self.get_sec_str_motif('H')

  def get_strands (self) :
    return self.get_sec_str_motif('S')

  def get_linkers (self) :
    return self.get_sec_str_motif('L')

  def get_missing (self) :
    return self.get_sec_str_motif('X')

  def get_region_by_id (self, code, idx) :
    sections = self.get_sec_str_motif(code)
    if idx < len(sections) :
      return sections[idx]

  def is_on_helix (self, x, y) :
    for k, (i_start, i_end) in enumerate(self.get_helices()) :
      bounds = self.get_region_bounds(i_start, i_end)
      for (x1, y1, x2, y2) in bounds :
        if (x > x1 and x < x2) and (y > y1 and y < y2) :
          return k
    return None

  def is_on_strand (self, x, y) :
    for k, (i_start, i_end) in enumerate(self.get_strands()) :
      bounds = self.get_region_bounds(i_start, i_end)
      for (x1, y1, x2, y2) in bounds :
        if (x > x1 and x < x2) and (y > (y1 + 2) and y < (y2 - 2)) :
          return k
    return None

  def is_on_linker (self, x, y) :
    for  k, (i_start, i_end) in enumerate(self.get_linkers()) :
      bounds = self.get_region_bounds(i_start, i_end)
      for (x1, y1, x2, y2) in bounds :
        if (x > x1 and x < x2) and (y > (y1 + 4) and y < (y2 - 4)) :
          return k
    return None

  def clear_structure_selections (self) :
    self.selected_helices = []
    self.selected_strands = []
    self.selected_linkers = []

  def select_helix (self, x, y) :
    idx = self.is_on_helix(x, y)
    if idx is not None :
      (helix_start, helix_end) = self.get_region_by_id('H', idx)
      if not idx in self.selected_helices :
        print "helix %d" % idx
        self.selected_helices.append(idx)
        self.select_chars(helix_start, helix_end)
      else :
        self.selected_helices.remove(idx)
        self.deselect_chars(helix_start, helix_end)
      return True
    return False

  def select_strand (self, x, y) :
    idx = self.is_on_strand(x, y)
    if idx is not None :
      (strand_start, strand_end) = self.get_region_by_id('S', idx)
      if not idx in self.selected_strands :
        print "strand %d" % idx
        self.selected_strands.append(idx)
        self.select_chars(strand_start, strand_end)
      else :
        self.selected_strands.remove(idx)
        self.deselect_chars(strand_start, strand_end)

  def select_linker (self, x, y) :
    idx = self.is_on_linker(x, y)
    if idx is not None :
      (linker_start, linker_end) = self.get_region_by_id('L', idx)
      if not idx in self.selected_linkers :
        print "linker %d" % idx
        self.selected_linkers.append(idx)
        self.select_chars(linker_start, linker_end)
      else :
        self.selected_linkers.remove(idx)
        self.deselect_chars(linker_start, linker_end)

  def OnDoubleClick (self, event) :
    (x, y) = (event.GetX(), event.GetY())
    if self.flag_overwrite_selections and not event.ShiftDown() :
      self.clear_structure_selections()
      self.clear_selection()
    if self.select_residue(x, y) :
      pass
    elif self.flag_allow_select_structure :
      if self.select_helix(x, y) :
        pass
      elif self.select_strand(x, y) :
        pass
      elif self.select_linker(x, y) :
        pass
      else :
        pass
    self.Refresh()

  def OnClear (self, event) :
    self.clear_structure_selections()
    SequencePanel.OnClear(self, event)

def strand_as_arrow (x1, y1, x2, y2) :
  segments = []
  segments.append((x1, y1 + 4))
  segments.append((x2 - 7, y1 + 4))
  segments.append((x2 - 7, y1))
  segments.append((x2, y1 + 8))
  segments.append((x2 - 7, y2))
  segments.append((x2 - 7, y2 - 4))
  segments.append((x1, y2 - 4))
  return segments

def strand_as_box (x1, y1, x2, y2) :
  segments = []
  segments.append((x1, y1 + 4))
  segments.append((x2, y1 + 4))
  segments.append((x2, y2 - 4))
  segments.append((x1, y2 - 4))
  return segments

########################################################################
class ControlPanel (wx.Panel) :
  def __init__ (self, *args, **kwds) :
    wx.Panel.__init__(self, *args, **kwds)
    szr = wx.BoxSizer(wx.HORIZONTAL)
    self.SetSizer(szr)
    self.multi_select_box = wx.CheckBox(self, -1,
      "Clicking overwrites selections")
    self.multi_select_box.SetValue(True)
    szr.Add(self.multi_select_box, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.clear_btn = wx.Button(self, -1, "Clear selection")
    szr.Add(self.clear_btn, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    szr.Fit(self)

  def bind_events (self, view) :
    self.Bind(wx.EVT_BUTTON, view.OnClear, self.clear_btn)
    self.Bind(wx.EVT_CHECKBOX, view.OnSetMode, self.multi_select_box)

class SequenceFrame (wx.Frame) :
  def __init__ (self, *args, **kwds) :
    wx.Frame.__init__(self, *args, **kwds)
    panel = SequenceWithStructurePanel(self, -1)
    szr = wx.BoxSizer(wx.VERTICAL)
    szr.Add(panel, 1, wx.EXPAND)
    self.SetSizer(szr)
    panel.enable_tooltip()
    panel2 = ControlPanel(self, -1, style=wx.SIMPLE_BORDER)
    panel2.bind_events(panel)
    szr.Add(panel2, 0, wx.EXPAND)
    self.statusbar = self.CreateStatusBar()
    self.Fit()
    self.panel = panel

  def set_sequence (self, seq) :
    self.panel.set_sequence(seq)

  def set_structure (self, sec_str) :
    self.panel.set_structure(sec_str)
    self.panel.apply_missing_residue_highlights()

if __name__ == "__main__" :
  app = wx.App(0)
  frame = SequenceFrame(None, -1, "Test sequence display")
  seq = """\
MFQAFPGDYDSGSRCSSSPSAESQYLSSVDSFGSPPTAAASQECAGLGEMPGSFVPTVTA
ITTSQDLQWLVQPTLISSMAQSQGQPLASQPPAVDPYDMPGTSYSTPGLSAYSTGGASGS
GGPSTSTTTSGPVSARPARARPRRPREETLTPEEEEKRRVRRERNKLAAAKCRNRRRELT
DRLQAETDQLEEEKAELESEIAELQKEKERLEFVLVAHKPGCKIPYEEGPGPGPLAEVRD
LPGSTSAKEDGFGWLLPPPPPPPLPFQSSRDAPPNLTASLFTHSEVQVLGDPFPVVSPSY
TSSFVLTCPEVSAFAGAQRTSGSEQPSDPLNSPSLLAL"""
  sec_str = """\
XXXXXXXXXXXLLLLLLLLLLLLLLLLSSSSSSSSSSSSSSSLLLLLLLHHHHHHHHLLS
SSSSSSSSLLLLLLLLLLLLSSSSSSSSSSLLLLLLLLSSSSSSSSSSLLLLLLLLHHHH
HHHHHHHHHHHHHHHHLLLLLHHHHHHHHHHHHHHHHHHHHLLLLLLLLLLLLLLLLLLL
LHHHHHHHHLLLLLLLLHHHHHHHHHHHHHHHHHLLLLLLLLLLLLHHHHHHHHHHHHHH
HHHHLLLLLLSSSSSSSSSSSSSSSSLLLLLLLLLLLLLSSSSSSSSSSSSSSSSSLLLL
LLLLLLLHHHHHHHHHHHHHHHHHHHHHHHHHLLLXXX"""
  frame.set_sequence(seq)
  frame.set_structure(sec_str)
  frame.Fit()
  frame.Show()
  app.MainLoop()
