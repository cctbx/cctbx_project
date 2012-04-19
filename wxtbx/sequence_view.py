
# Copyright 2010 University of California
#
# Sequence view and selection window, with optional secondary structure
# annotation.  The outer window and PDB processing methods depend on CCTBX,
# but the panels can be adapted to any other framework.

# FIXME use numbering/resids in PDB file

# TODO: show existing CCTBX atom selection?

from __future__ import division
import wx
import wx.lib.scrolledpanel
import cStringIO
import math, sys, os

WXTBX_SEQ_SELECT_NONE = 1
WXTBX_SEQ_SELECT_SINGLE = 2
WXTBX_SEQ_SELECT_MULTIPLE = 4
WXTBX_SEQ_SELECT_RANGE = 8
WXTBX_SEQ_SELECT_ANY = 16
WXTBX_SEQ_SHOW_LINE_NUMBERS = 32
WXTBX_SEQ_SHOW_LABELS = 64
WXTBX_SEQ_SHOW_TOOLTIP = 128
WXTBX_SEQ_SHOW_SELECTIONS = 256
WXTBX_SEQ_ENABLE_SELECT_STRUCTURE = 512
WXTBX_SEQ_ENABLE_SELECT_MISSING = 1024
WXTBX_SEQ_FANCY_HELICES = 2048
WXTBX_SEQ_SINGLE_CLICK_SELECTION = 4096
WXTBX_SEQ_MSA_SELECT_ALL = 8192

WXTBX_SEQ_DEFAULT_STYLE = WXTBX_SEQ_SHOW_LINE_NUMBERS | \
                          WXTBX_SEQ_SHOW_SELECTIONS | \
                          WXTBX_SEQ_ENABLE_SELECT_STRUCTURE | \
                          WXTBX_SEQ_ENABLE_SELECT_MISSING | \
                          WXTBX_SEQ_FANCY_HELICES

class sequence_panel (wx.PyPanel) :
  tooltip = "Double-click a residue to select it; hold down Shift to select \
multiple residues."
  __bg_color = (255,255,255)
  line_sep = 28
  line_width = 50
  line_height = 16
  line_indent = 16
  start_offset = 0

  def __init__ (self, *args, **kwds) :
    wx.PyPanel.__init__(self, *args, **kwds)
    if self.__bg_color is not None :
      self.SetBackgroundColour(self.__bg_color)
    #from scitbx.array_family import flex, shared
    self.sequences = []
    self.sequence_labels = []
    self.char_boxes = [] #shared.stl_set_unsigned()
    self.highlights = []
    self.highlight_colors = []
    self.selected_residues = [] #flex.bool()
    self.selection_color = (255, 255, 0)
    self._last_x = None
    self._last_y = None
    self._style = WXTBX_SEQ_DEFAULT_STYLE | WXTBX_SEQ_SELECT_ANY
    self.txt_font = wx.Font(14, wx.MODERN, wx.NORMAL, wx.NORMAL)
    self.Bind(wx.EVT_PAINT, self.OnPaint)
    self.Bind(wx.EVT_LEFT_DCLICK, self.OnDoubleClick)
    self.Bind(wx.EVT_LEFT_DOWN, self.OnLeftDown)
    self.Bind(wx.EVT_LEFT_UP, self.OnLeftUp)
    self.Bind(wx.EVT_RIGHT_DOWN, self.OnRightDown)
    self.Bind(wx.EVT_RIGHT_UP, self.OnRightUp)
    self.tip = wx.ToolTip(self.tooltip)
    self.SetToolTip(self.tip)
    self.tip.Enable(False)

  def SetStyle (self, style) :
    self._style = style
    self.enable_tooltip(style & WXTBX_SEQ_SHOW_TOOLTIP)

  def enable_tooltip (self, enable=True) :
    self.tip.Enable(enable)

  def OnPaint (self, event) :
    dc = wx.AutoBufferedPaintDCFactory(self)
    # XXX is there any reason not to do this on all systems?  test on Linux
    if (wx.Platform == "__WXMSW__") :
      dc.SetBackground(wx.WHITE_BRUSH)
      dc.SetBackgroundMode(wx.SOLID)
      dc.Clear()
    gc = wx.GraphicsContext.Create(dc)
    self.paint(gc)

  def paint (self, gc) :
    self.paint_sequence(gc)

  def paint_sequence (self, gc) :
    gc.SetBrush(wx.TRANSPARENT_BRUSH)
    gc_txt_font = gc.CreateFont(self.txt_font, (0,0,0))
    gc.SetFont(gc_txt_font)
    black_pen = wx.Pen((0,0,0), 1)
    (char_w, char_h) = self.get_char_size(gc)
    if (self._style & WXTBX_SEQ_SHOW_SELECTIONS) :
      if (self._style & WXTBX_SEQ_MSA_SELECT_ALL) :
        pass
      else :
        for k_seq in range(self.n_seqs()) :
          for j, (i_start, i_end) in enumerate(self.get_selected_ranges()) :
            gc.SetPen(wx.Pen((0,0,0), 1))
            gc.SetBrush(wx.Brush(self.selection_color))
            self.draw_box_around_residues(k_seq, i_start, i_end, gc)
    gc.SetBrush(wx.TRANSPARENT_BRUSH)
    for k_seq in range(self.n_seqs()) :
      xpos = self.line_indent
      ypos = self.get_line_sep()
      v_offset = k_seq * self.line_height
      i = 0
      while i < len(self.sequences[k_seq]) :
        j = i
        line_start = self.start_offset + i + 1
        line_end = line_start + self.line_width - 1
        line = self.sequences[k_seq][i:i+self.line_width]
        line_w = self.line_width * char_w
        label_start, label_end = self.get_line_labels(k_seq, i)
        label_w = self.get_label_width(gc)
        if (label_start is not None) :
          gc.DrawText(label_start, xpos, ypos + v_offset)
          xpos += label_w
        if (label_end is not None) :
          gc.DrawText(label_end, xpos + label_w + line_w, ypos + v_offset)
        for j, char in enumerate(line) :
          i_seq = i+j
          if i_seq in self.highlights :
            color = self.highlight_colors[self.highlights.index(i_seq)]
            gc.SetFont(gc.CreateFont(self.txt_font, color))
          gc.DrawText(char, xpos, ypos + v_offset)
          gc.SetFont(gc_txt_font)
          xpos += char_w
        i += self.line_width
        ypos += self.get_line_spacing()
        xpos = self.line_indent
    return gc

  def draw_box_around_residues (self, k_seq, i_start, i_end, gc) :
    ranges = self.get_contiguous_ranges(i_start, i_end)
    char_w, char_h = self.get_char_size(gc)
    for (seg_start, seg_end) in ranges :
      (x1, y1, nx, ny) = self.get_char_position(k_seq, seg_start)
      (nx, ny, x2, y2) = self.get_char_position(k_seq, seg_end)
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
    dc = wx.GraphicsContext.CreateMeasuringContext() #ClientDC(self)
    dc.SetFont(self.txt_font)
    i = 0
    (panel_w, panel_h) = (32, 32)
    char_w, char_h = self.get_char_size(dc)
    line_w = char_w * self.line_width
    if (self._style & WXTBX_SEQ_SHOW_LABELS) :
      line_w += self.get_label_width(dc)
    elif (self._style & WXTBX_SEQ_SHOW_LINE_NUMBERS) :
      line_w += dc.GetTextExtent("X" * 16)[0]
    panel_w += line_w
    n_lines = int(math.ceil(self.sequence_length() / self.line_width))
    panel_h += n_lines * self.get_line_spacing()
    return (max(480, panel_w), max(240, panel_h))

  def set_sequence (self, seq) :
    self.set_sequences([seq])

  def set_sequences (self, seqs, labels=()) :
    self.sequences = [ "".join(seq.splitlines()) for seq in seqs ]
    assert (len(set([ len(s) for s in self.sequences ])) == 1)
    assert (len(labels) == 0) or (len(labels) == len(seqs))
    if (len(self.sequence_labels) != 0) :
      assert (len(self.sequence_labels) == len(seqs))
    else :
      self.sequence_labels = labels
    self.build_boxes()
    self.clear_highlights()
    self.clear_selection()

  def set_sequence_labels (self, labels) :
    assert (len(labels) == len(self.sequences))
    self.sequence_labels = labels
    self.build_boxes()

  def n_seqs (self) :
    return len(self.sequences)

  def sequence_length (self) :
    if (self.n_seqs() > 0) :
      return len(self.sequences[0])
    return 0

  def get_line_sep (self) :
    return self.line_sep

  def get_line_spacing (self) :
    return ((self.n_seqs() - 1) * self.line_height) + self.get_line_sep()

  def get_line_labels (self, k_seq, i_seq) :
    line_start = self.start_offset + i_seq + 1
    line_end = line_start + self.line_width - 1
    if (self._style & WXTBX_SEQ_SHOW_LABELS) :
      label_start = self.sequence_labels[k_seq]
      label_end = None
    elif (self._style & WXTBX_SEQ_SHOW_LINE_NUMBERS) :
      label_start = "%6d  " % line_start
      label_end = "  %-6d" % line_end
    return (label_start, label_end)

  def get_label_width (self, dc) :
    if (self._style & WXTBX_SEQ_SHOW_LABELS) :
      if (len(self.sequence_labels) > 0) :
        n_chars = max([ len(l) for l in self.sequence_labels ])
      else :
        n_chars = 0
    elif (self._style & WXTBX_SEQ_SHOW_LINE_NUMBERS) :
      n_chars = 6
    else :
      n_chars = 0
    return dc.GetTextExtent("X" * (n_chars + 2))[0]

  def get_char_size (self, dc=None) :
    if dc is None :
      dc = wx.GraphicsContext.CreateMeasuringContext() #ClientDC(self)
      dc.SetFont(self.txt_font)
    line_w, char_h = dc.GetTextExtent("X" * 50)
    if wx.Platform == '__WXGTK__' :
      char_w = max(12, line_w / 50)
    elif wx.Platform == '__WXMAC__' :
      char_w = max(10, line_w / 50)
    elif (wx.Platform == '__WXMSW__') :
      char_w = max(10, line_w / 50)
    else :
      raise RuntimeError("Platform not supported!")
    char_h = max(16, char_h)
    return (char_w, char_h)

  def build_boxes (self) :
    dc = wx.GraphicsContext.CreateMeasuringContext() #ClientDC(self)
    dc.SetFont(self.txt_font)
    char_w, char_h = self.get_char_size(dc)
    x_start = 16
    x_start += self.get_label_width(dc)
    self.char_boxes = [ [] for k_seq in range(self.n_seqs()) ]
    for i_seq in range(self.sequence_length()) :
      n_lines = int(math.floor((i_seq-self.start_offset) / self.line_width))
      n_prev_chars = (i_seq - self.start_offset) % self.line_width
      x = x_start + (char_w * n_prev_chars)
      y_start = self.get_line_sep() + (n_lines * self.get_line_spacing())
      for k_seq in range(self.n_seqs()) :
        v_offset = k_seq * self.line_height
        y = y_start + v_offset
        self.char_boxes[k_seq].append((x, y, x + char_w - 1, y + char_h - 1))

  def get_char_position (self, k_seq, i_seq) :
    if i_seq >= len(self.char_boxes[k_seq]) :
      print i_seq, len(self.char_boxes[k_seq])
    return tuple(self.char_boxes[k_seq][i_seq])

  def get_contiguous_ranges (self, i_start, i_end) :
    assert (i_start < i_end) or (i_start == i_end)
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
    for i_seq, selected in enumerate(self.selected_residues) :
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
    if last_end is not None :
      ranges.append((last_start, last_end))
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

  def set_highlights (self, indices, **kwds) :
    for i_seq in indices :
      self.highlight_char(i_seq, **kwds)

  def clear_selection (self) :
    self.selected_residues = [False] * self.sequence_length()
    self.update_frame()

  def select_chars (self, i_start, i_end, box=True) :
    assert (i_end < self.sequence_length()) and (i_start <= i_end)
    for i_seq in range(i_start, i_end + 1) :
      self.selected_residues[i_seq] = box
    self.update_frame()

  def update_frame (self) :
    frame = self.GetParent().GetParent()
    frame.callback_on_select()

  def get_selection_info (self) :
    ranges = self.get_selected_ranges()
    if len(ranges) == 0 :
      txt = ""
    elif len(ranges) == 1 and ranges[0][0] == ranges[0][1] :
      txt = "SELECTED: residue %d" % (ranges[0][0] + 1)
    else :
      txt_ranges = []
      for (x, y) in ranges :
        if x == y : txt_ranges.append(str(x+1))
        else : txt_ranges.append("%d-%d" % (x+1, y+1))
      txt = "SELECTED: residues %s" % ", ".join(txt_ranges)
    return txt

  def get_atom_selection (self) :
    ranges = self.get_selected_ranges()
    if len(ranges) == 0 :
      return None
    elif len(ranges) == 1 and ranges[0][0] == ranges[0][1] :
      return "(resseq %d)" % (ranges[0][0] + 1)
    else :
      resseqs = []
      for (x, y) in ranges :
        if x == y : resseqs.append("resseq %d" % (x + 1))
        else : resseqs.append("resseq %d:%d" % (x+1, y+1))
      return "(" + " or ".join(resseqs) + ")"

  def deselect_chars (self, i_start, i_end) :
    self.select_chars(i_start, i_end, False)

  def select_residue (self, x, y) :
    for k_seq in range(self.n_seqs()) :
      for i_seq, box in enumerate(self.char_boxes[k_seq]) :
        if (self.sequences[k_seq][i_seq] == '-') :
          continue
        (x1, y1, x2, y2) = box
        if (x > x1 and x < x2) and (y > y1 and y < y2) :
          if self._style & WXTBX_SEQ_SELECT_RANGE :
            self.process_range_selection(i_seq)
          else :
            if self.selected_residues[i_seq] :
              self.selected_residues[i_seq] = False
            else :
              self.selected_residues[i_seq] = True
          self.update_frame()
          return True
    return False

  def process_range_selection (self, i_seq) :
    n_selected = self.selected_residues.count(True)
    if n_selected == 0 :
      self.selected_residues[i_seq] = True
    elif n_selected == 1 :
      j_seq = self.selected_residues.index(True)
      if j_seq != i_seq :
        for k_seq in range(min(i_seq, j_seq), max(i_seq+1,j_seq)) :
          self.selected_residues[k_seq] = True
    else :
      self.clear_selection()
      self.selected_residues[i_seq] = True

  def OnClear (self, event) :
    self.clear_selection()
    self.Refresh()

  def OnHelp (self, event) :
    self.enable_tooltip(True)

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

  def OnRightUp (self, event) :
    (x, y) = (event.GetX(), event.GetY())

#-----------------------------------------------------------------------
class sequence_with_structure_panel (sequence_panel) :
  tooltip = """\
Double-click on any residue or secondary-structure element to select the \
residue(s).  Holding down shift enables multiple selections."""
  line_sep = 64

  def __init__ (self, *args, **kwds) :
    sequence_panel.__init__(self, *args, **kwds)
    self.structure = ""
    self.selected_helices = []
    self.selected_strands = []
    self.selected_linkers = []
    self.selected_missing = []

  def paint (self, gc) :
    self.paint_sequence(gc)
    self.paint_structure(gc)

  def paint_structure (self, gc) :
    helices = self.get_helices()
    if self._style & WXTBX_SEQ_FANCY_HELICES :
      helix_pen = wx.Pen('red', 1)
      h_helix_pen = wx.Pen((255,50,50), 2)
      for k, (i_start, i_end) in enumerate(helices) :
        if k in self.selected_helices :
          color_inner_start = (255, 220, 220)
          color_inner_end = (255, 240, 240)
          color_outer_start = (255, 170, 170)
          color_outer_end = (255, 220, 220)
          gc.SetPen(h_helix_pen)
        else :
          color_inner_start = (255, 150, 150)
          color_inner_end = (255, 200, 200)
          color_outer_start = (255, 100, 100)
          color_outer_end = (255, 150, 150)
          gc.SetPen(helix_pen)
        bounds = self.get_region_bounds(i_start, i_end)
        for box in bounds :
          (x1, y1, x2, y2) = box
          strips = make_helix(x1 - 1, y1, x2, y2)
          def draw_strip (strip, brush) :
            gc.PushState()
            gc.SetBrush(brush) #wx.Brush(strip_color))
            path = gc.CreatePath()
            path.MoveToPoint(strip[0][0], strip[0][1])
            for x, y in strip[1:] :
              path.AddLineToPoint(x, y)
            path.CloseSubpath()
            gc.FillPath(path)
            gc.StrokePath(path)
            gc.PopState()
          inner_brush = gc.CreateLinearGradientBrush(x1-1, y1-2, x1-1, y2+2,
            color_inner_start, color_inner_end)
          for i, strip in enumerate(strips) :
            if (i % 2) == 0 :
              draw_strip(strip, inner_brush)
          outer_brush = gc.CreateLinearGradientBrush(x1-1, y1-2, x1-1, y2+2,
            color_outer_start, color_outer_end)
          for i, strip in enumerate(strips) :
            if (i % 2) == 1 :
              draw_strip(strip, outer_brush)
        gc.SetPen(helix_pen)
    else :
      helix_pen = wx.Pen('red', 1)
      h_helix_pen = wx.Pen((255,50,50), 2)
      for k, (i_start, i_end) in enumerate(helices) :
        bounds = self.get_region_bounds(i_start, i_end)
        if k in self.selected_helices :
          color_start = (200, 100, 100)
          color_end = (255, 200, 200)
          gc.SetPen(h_helix_pen)
        else :
          color_start = (255, 50, 50)
          color_end = (255, 150, 150)
          gc.SetPen(helix_pen)
        for box in bounds :
          (x1, y1, x2, y2) = box
          helix_brush = gc.CreateLinearGradientBrush(x1 - 1, y1, x1 - 1, y2,
            color_start, color_end)
          gc.SetBrush(helix_brush)
          gc.DrawRoundedRectangle(x1 - 1, y1, x2 - x1 + 1, 16, 4)
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
    if (wx.Platform == '__WXGTK__') : # dashed pen freezes on Linux
      missing_pen = wx.Pen((150, 150, 150), 4)
    elif (wx.Platform in ['__WXMAC__', '__WXMSW__']) :
      missing_pen = wx.Pen((150, 150, 150), 4, style=wx.SHORT_DASH)
    h_missing_pen = wx.Pen((255, 255, 0), 12)
    for k, (i_start, i_end) in enumerate(missing) :
      bounds = self.get_region_bounds(i_start, i_end)
      for j, box in enumerate(bounds) :
        (x1, y1, x2, y2) = box
        line = gc.CreatePath()
        line.MoveToPoint(x1, y1 + 8)
        line.AddLineToPoint(x2, y1 + 8)
        line.CloseSubpath()
        if k in self.selected_missing :
          gc.PushState()
          gc.SetPen(h_missing_pen)
          gc.StrokePath(line)
          gc.PopState()
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

  def get_line_sep (self) :
    #if (self.structure == "") :
    #  return sequence_panel.line_sep
    return self.line_sep

  def set_structure (self, ss) :
    self.structure = "".join(ss.splitlines())
    assert ((self.sequence_length() == 0) or
            (len(self.structure) == self.sequence_length()))
    self.clear_highlights()
    self.clear_structure_selections()

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
      (x1, y1, nx, ny) = self.get_char_position(0, seg_start)
      y1 -= 24
      (nx, y2, x2, ny) = self.get_char_position(0, seg_end)
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
    return self.is_on_region(x, y, self.get_helices(), tolerance=2)

  def is_on_strand (self, x, y) :
    return self.is_on_region(x, y, self.get_strands())

  def is_on_linker (self, x, y) :
    return self.is_on_region(x, y, self.get_linkers(), tolerance=-2)

  def is_on_missing (self, x, y) :
    return self.is_on_region(x, y, self.get_missing(), tolerance=-2)

  def is_on_region (self, x, y, regions, tolerance=0) :
    for k, (i_start, i_end) in enumerate(regions) :
      bounds = self.get_region_bounds(i_start, i_end)
      for (x1, y1, x2, y2) in bounds :
        if ((x > x1 and x < x2) and
            (y > (y1-tolerance) and y < (y2+tolerance))) :
          return k
    return None

  def clear_structure_selections (self) :
    self.selected_helices = []
    self.selected_strands = []
    self.selected_linkers = []
    self.selected_missing = []

  def select_helix (self, x, y) :
    idx = self.is_on_helix(x, y)
    if idx is not None :
      (helix_start, helix_end) = self.get_region_by_id('H', idx)
      if not idx in self.selected_helices :
        #print "helix %d" % idx
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
        #print "strand %d" % idx
        self.selected_strands.append(idx)
        self.select_chars(strand_start, strand_end)
      else :
        self.selected_strands.remove(idx)
        self.deselect_chars(strand_start, strand_end)
      return True
    return False

  def select_linker (self, x, y) :
    idx = self.is_on_linker(x, y)
    if idx is not None :
      (linker_start, linker_end) = self.get_region_by_id('L', idx)
      if not idx in self.selected_linkers :
        #print "linker %d" % idx
        self.selected_linkers.append(idx)
        self.select_chars(linker_start, linker_end)
      else :
        self.selected_linkers.remove(idx)
        self.deselect_chars(linker_start, linker_end)
      return True
    return False

  def select_missing (self, x, y) :
    if not (self._style & WXTBX_SEQ_ENABLE_SELECT_MISSING) :
      return False
    idx = self.is_on_missing(x, y)
    if idx is not None :
      (region_start, region_end) = self.get_region_by_id('X', idx)
      if not idx in self.selected_missing :
        #print "missing segment %d" % idx
        self.selected_missing.append(idx)
        self.select_chars(region_start, region_end)
      else :
        self.selected_missing.remove(idx)
        self.deselect_chars(region_start, region_end)
      return True
    return False

  def process_click (self, x, y, shift_down=False) :
    if not (self._style & WXTBX_SEQ_SELECT_MULTIPLE) :
      if (self._style & WXTBX_SEQ_SELECT_ANY) and shift_down :
        pass
      else :
        self.clear_structure_selections()
        if not (self._style & WXTBX_SEQ_SELECT_RANGE) :
          self.clear_selection()
    if self.select_residue(x, y) :
      pass
    elif (self._style & WXTBX_SEQ_ENABLE_SELECT_STRUCTURE) :
      if (self._style & WXTBX_SEQ_SELECT_RANGE) :
        self.clear_selection()
      if self.select_helix(x, y) :
        pass
      elif self.select_strand(x, y) :
        pass
      elif self.select_linker(x, y) :
        pass
      elif self.select_missing(x, y) :
        pass
      else :
        self.clear_selection()
    self.Refresh()

  def OnDoubleClick (self, event) :
    (x, y) = (event.GetX(), event.GetY())
    self.process_click(x, y, event.ShiftDown())

  def OnClear (self, event) :
    self.clear_structure_selections()
    sequence_panel.OnClear(self, event)

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

def make_helix (x1, y1, x2, y2, strip_width=10, strip_height=20) :
  y_center = y1 + ((y2 - y1) / 2)
  x_start = x1
  strips = []
  strips.append(((x_start, y_center),
                 (x_start + strip_width, y_center),
                 (x_start + (strip_width * 1.5), y_center + (strip_height/2)),
                 (x_start + (strip_width * 0.5), y_center + (strip_height/2))))
  x_start += (strip_width * 0.5)
  y_start = y_center + (strip_height / 2)
  y_end = y_center - (strip_height / 2)
  while x_start < (x2 - (strip_width * 1.5)) :
    strips.append(((x_start, y_start),
                   (x_start + strip_width, y_start),
                   (x_start + (strip_width * 2), y_end),
                   (x_start + strip_width, y_end)))
    x_start += strip_width
    y_tmp = y_start
    y_start = y_end
    y_end = y_tmp
  if x_start < (x2 - strip_width) :
    strips.append(((x_start, y_start),
                   (x_start + strip_width, y_start),
                   (x_start + (strip_width * 1.5), y_center),
                   (x_start + (strip_width * 0.5), y_center)))
  return strips

########################################################################
class control_panel (wx.Panel) :
  def __init__ (self, *args, **kwds) :
    wx.Panel.__init__(self, *args, **kwds)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    box = wx.BoxSizer(wx.HORIZONTAL)
    box.Add(wx.StaticText(self, -1, "Select chain:"), 0, wx.ALL, 5)
    self.chain_select = wx.Choice(self, -1)
    box.Add(self.chain_select, 0, wx.ALL, 5)
    self.clear_btn = wx.Button(self, -1, "Clear selection")
    box.Add(self.clear_btn, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.help_btn = wx.Button(self, wx.ID_HELP)
    box.Add(self.help_btn, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    szr.Add(box)
    szr.Fit(self)

  def bind_events (self, frame, view) :
    self.Bind(wx.EVT_BUTTON, view.OnClear, self.clear_btn)
    self.Bind(wx.EVT_BUTTON, view.OnHelp, self.help_btn)
    self.Bind(wx.EVT_CHOICE, frame.OnSelectChain, self.chain_select)
    self.Bind(wx.EVT_BUTTON, frame.OnHelp, self.help_btn)

# XXX keep this general so it can be embedded in other windows
class sequence_window (object) :
  cp_style = wx.NO_BORDER
  seq_panel_style = WXTBX_SEQ_DEFAULT_STYLE
  def create_main_panel (self, sizer=None) :
    if sizer is None :
      sizer = self.sizer
    outer_panel = wx.lib.scrolledpanel.ScrolledPanel(self, -1)
    szr2 = wx.BoxSizer(wx.VERTICAL)
    outer_panel.SetSizer(szr2)
    panel = sequence_with_structure_panel(outer_panel, -1)
    panel.SetStyle(self.seq_panel_style)
    szr2.Add(panel, 1, wx.EXPAND)
    szr2.Layout()
    szr2.Fit(outer_panel)
    sizer.Add(outer_panel, 1, wx.EXPAND)
    outer_panel.SetupScrolling()
    self.seq_panel = panel
    self.outer_panel = outer_panel
    if getattr(self, "control_panel", None) is not None :
      self.control_panel.bind_events(self, self.seq_panel)

  def create_control_panel (self, sizer=None) :
    if sizer is None :
      sizer = self.sizer
    cp = control_panel(self, -1, style=self.cp_style)
    sizer.Add(cp, 0, wx.EXPAND)
    self.control_panel = cp
    if getattr(self, "seq_panel", None) is not None :
      self.control_panel.bind_events(self, self.seq_panel)

  def load_pdb_file (self, file_name) :
    from iotbx import file_reader
    pdb_in = file_reader.any_file(file_name, force_type="pdb").file_object
    pdb_hierarchy = pdb_in.construct_hierarchy()
    pdb_hierarchy.atoms().reset_i_seq()
    xray_structure = pdb_in.xray_structure_simple()
    self.set_pdb_data(pdb_hierarchy, xray_structure)

  def set_pdb_data (self, pdb_hierarchy, xray_structure) :
    from mmtbx import secondary_structure
    sec_str = secondary_structure.manager(pdb_hierarchy, xray_structure)
    out = cStringIO.StringIO()
    sec_str.find_automatically(log=out)
    sec_str.show_summary()
    self.set_data(pdb_hierarchy, sec_str, auto_select=True)

  def set_data (self, pdb_hierarchy, sec_str, auto_select=True) :
    self.pdb_hierarchy = pdb_hierarchy
    self.sec_str = sec_str
    self._chain_cache = {}
    self._seq_cache = {}
    self._ss_cache = {}
    for chain in pdb_hierarchy.models()[0].chains() :
      conf = chain.conformers()[0]
      if not conf.is_protein() :
        print "Skipping non-protein chain %s" % chain.id
        continue
      if chain.id in self._chain_cache :
        n = 2
        while True :
          new_id = "%s (%d)" % (chain.id, n)
          if not (new_id in self._chain_cache) :
            self._chain_cache[new_id] = chain
            break
          n += 1
      else :
        self._chain_cache[chain.id] = chain
    self.control_panel.chain_select.SetItems(sorted(self._chain_cache.keys()))
    self.control_panel.Layout()
    if auto_select :
      self.control_panel.chain_select.SetSelection(0)
      chain_id = str(self.control_panel.chain_select.GetStringSelection())
      self.set_current_chain(chain_id)

  def set_current_chain (self, chain_id) :
    chain = self._chain_cache.get(chain_id, None)
    if chain is not None :
      chain_conf = chain.conformers()[0]
      if chain_id in self._seq_cache :
        seq = self._seq_cache[chain_id]
        ss = self._ss_cache[chain_id]
      else :
        helix_sele = self.sec_str.alpha_selection()
        sheet_sele = self.sec_str.beta_selection()
        seq = chain_conf.as_padded_sequence()
        ss = chain_conf.as_sec_str_sequence(helix_sele, sheet_sele)
        self._seq_cache[chain_id] = seq
        self._ss_cache[chain_id] = ss
      self.set_sequence(seq)
      self.set_structure(ss)
      self.reset_layout()
      self.seq_panel.Refresh()

  def set_chain_sequence_data (self, data) :
    assert isinstance(data, list)
    assert ((len(data) == 0) or (type(data[0]).__name__ == "chain"))
    self._seq_data = data
    ids = [ chain.chain_id for chain in data ]
    self.control_panel.chain_select.SetItems(ids)
    if (len(ids) > 0) :
      self.set_current_alignment(0)

  def set_current_alignment (self, index) :
    assert (index < len(self._seq_data))
    chain = self._seq_data[index]
    seqs = chain.get_alignment(include_sec_str=True)
    self.set_sequences([seqs[0], seqs[1]], labels=["PDB", "sequence"])
    self.set_structure(seqs[2])
    self.reset_layout()
    self.seq_panel.Refresh()

  def set_sequence (self, *args, **kwds) :
    self.seq_panel.set_sequence(*args, **kwds)

  def set_sequences (self, *args, **kwds) :
    self.seq_panel.set_sequences(*args, **kwds)

  def set_sequence_labels (self, *args, **kwds) :
    self.seq_panel.set_sequence_labels(*args, **kwds)

  def reset_layout (self) :
    pass

  def set_structure (self, sec_str) :
    self.seq_panel.set_structure(sec_str)
    self.seq_panel.apply_missing_residue_highlights()

  def get_selection_base (self) :
    chain_id = str(self.control_panel.chain_select.GetStringSelection())
    return "chain '%s'" % chain_id

  def callback_on_select (self) :
    pass

  def OnSelectChain (self, evt) :
    if (getattr(self, "_seq_data", None) is not None) :
      self.set_current_alignemnt(evt.GetEventObject().GetSelection())
    else :
      chain_id = str(evt.GetEventObject().GetStringSelection())
      self.set_current_chain(chain_id)

  def OnClose (self, evt) :
    self.Destroy()

  def OnDestroy (self, evt) :
    pass

  def OnHelp (self, evt) :
    pass

class sequence_frame_mixin (wx.Frame) :
  cp_style = wx.SIMPLE_BORDER
  def __init__ (self, *args, **kwds) :
    wx.Frame.__init__(self, *args, **kwds)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    self.sizer = szr
    self.create_main_panel()
    self.create_control_panel()
    self.statusbar = self.CreateStatusBar()
    self.Bind(wx.EVT_CLOSE, self.OnClose)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy)
    self._chain_cache = []
    self._selection_callback = None

  def reset_layout (self) :
    self.seq_panel.Layout()
    self.outer_panel.Layout()
    self.sizer.Layout()
    (w, h) = self.seq_panel.DoGetBestSize()
    (w2, h2) = self.control_panel.GetSize()
    if w <= 750 and h <= 550 :
      self.outer_panel.SetMinSize((w,h))
      self.SetSize((w+50, h+h2+50))
    else :
      self.SetSize((800,600))

class sequence_frame (sequence_frame_mixin, sequence_window) :
  def callback_on_select (self) :
    txt = self.seq_panel.get_selection_info()
    self.statusbar.SetStatusText(txt)

class msa_frame (sequence_frame_mixin, sequence_window) :
  seq_panel_style = WXTBX_SEQ_SHOW_SELECTIONS | \
                    WXTBX_SEQ_SHOW_LABELS | \
                    WXTBX_SEQ_ENABLE_SELECT_STRUCTURE | \
                    WXTBX_SEQ_ENABLE_SELECT_MISSING | \
                    WXTBX_SEQ_FANCY_HELICES | \
                    WXTBX_SEQ_SELECT_SINGLE

#-----------------------------------------------------------------------
def run (args) :
  pdb_file = args[-1]
  app = wx.App(0)
  frame = sequence_frame(None, -1, "Sequence display for %s" %
    os.path.basename(pdb_file))
  frame.load_pdb_file(pdb_file)
  if wx.Platform == '__WXMAC__' :
    frame.Fit()
  frame.Show()
  if "--range" in args :
    frame.seq_panel.SetStyle(WXTBX_SEQ_DEFAULT_STYLE|WXTBX_SEQ_SELECT_RANGE)
  app.MainLoop()

def run_test () :
  app = wx.App(0)
  exercise_simple_frame()
  app.MainLoop()

def exercise_simple_frame () :
  frame = sequence_frame(None, -1, "Demo of sequence display")
  frame.set_sequence('VLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQG')
  frame.set_structure('LLLHHHHHHHHHHHHHHHHLHHHHHHHHHHHHHHHLHHHHHLLHHHHLLLLHHHHHHLHHHHHHHHHHHHHHHHHHHHLLLLHHHHHHHHHHHHHLLLLLHHHHHHHHHHHHHHHHHHLLLLLLHHHHHHHHHHHHHHHHHHHHHHHHHLLLL')
  frame.reset_layout()
  frame.Show()

__test_seqs = [
  "---VLSEGEW---QLVLHVWAKVEADVAGHGQDILIRLFKSHPESTLEKFDRFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKK------GHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQGANQ",
  "------EGEWMRTQLVLHVWAKVEADVAGHGQDILIKLFKSHPE-TLEKFDRFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKMQVEHNGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQG---",
  "MAVVLSEGEWM--QLVLHVWAKVEADVAGHGQDILIRLFKSHPESTLEKFDRFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKK--VEH-GHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKEL-------"]
__test_labels = [
    "Protein A (H. sapiens)",
    "Protein B (M. musculus)",
    "Protein C (B. taurus)",]

def exercise_msa_frame () :
  frame = msa_frame(None, -1, "Demo of multiple sequence alignment")
  frame.set_sequences(__test_seqs, labels=__test_labels)
  frame.set_structure('------LLLHHHHHHHHHHHHHHHHLHHHHHHHHHHHHHHHLHH-HHHLLHHHHLLLLHHHHHHLHHHHHHHHHHHHHHHHHHHHLLLLLLLLLLHHHHHHHHHHHHHLLLLLHHHHHHHHHHHHHHHHHHLLLLLLHHHHHHHHHHHHHHHHHHHHHHHHHLLLL---')
  frame.reset_layout()
  frame.Show()

def exercise_msa_frame_large () :
  frame = msa_frame(None, -1, "Really big sequence alignment")
  seqs = [ s*10 for s in __test_seqs ]
  labels = __test_labels
  frame.set_sequences(seqs, labels=labels)
  frame.reset_layout()
  frame.Show()

def exercise () :
  app = wx.App(0)
  exercise_simple_frame()
  exercise_msa_frame()
  exercise_msa_frame_large()
  app.MainLoop()

if __name__ == "__main__" :
#  run_test()
  if "--test" in sys.argv :
    exercise()
  else :
    run(sys.argv[1:])
