from __future__ import absolute_import, division, print_function

import wx

__default_border__ = 5
class SizerContainer(object):
  __default_text_width__ = 480
  if (wx.Platform == '__WXMAC__'):
    __default_font_size__ = 12
    __default_label_size__ = 16
  elif (wx.Platform == '__WXGTK__'):
    __default_font_size__ = 11
    __default_label_size__ = 14
  else :
    __default_font_size__ = 11
    __default_label_size__ = 14

  def __init__(self):
    self._panel_list = []
    self._boxes = {}
    self._box_labels = []
    self._sizer_stack = []
    self._font_size = self.__default_font_size__
    self._label_size = self.__default_label_size__
    self._text_width = self.__default_text_width__
    self._current_sizer = None
    self._current_panel = None
    self.main_sizer = None
    self._label_font = wx.Font(self._label_size, wx.FONTFAMILY_DEFAULT,
      wx.NORMAL, wx.BOLD)

  def get_panel(self):
    return self._current_panel

  #--- properties
  def set_font_size(self, font_size):
    assert isinstance(font_size, int)
    self._font_size = font_size

  def reset_font_size(self):
    self._font_size = self.__default_font_size__

  def set_text_width(self, width):
    assert isinstance(width, int)
    self._text_width = width

  def set_panel(self, panel):
    self._current_panel = panel

  def set_main_sizer(self, sizer):
    assert (self.main_sizer is None) and (len(self._sizer_stack) == 0)
    self.main_sizer = sizer
    self._sizer_stack.append(sizer)
    self._current_sizer = sizer

  #--- sizer management
  def get_current_sizer(self):
    return self._current_sizer

  def _start_box(self, sizer_type, label="", border=0, indent=None,
      expand=False, proportion=0, center=False):
    assert (sizer_type in [wx.HORIZONTAL, wx.VERTICAL])
    self._sizer_stack.append(self._current_sizer)
    new_sizer = wx.BoxSizer(sizer_type)
    flags = wx.ALL
    if expand :
      flags |= wx.EXPAND
    if center :
      flags |= wx.ALIGN_CENTER
    self._current_sizer.Add(new_sizer, proportion, flags, border)
    self._current_sizer = new_sizer
    if (indent is not None):
      assert isinstance(indent, int)
      self._current_sizer.Add((indent, -1))

  def start_hbox(self, *args, **kwds):
    self._start_box(wx.HORIZONTAL, *args, **kwds)

  def start_vbox(self, *args, **kwds):
    self._start_box(wx.VERTICAL, *args, **kwds)

  def start_grid(self, rows=0, cols=2, proportion=0, flags=0):
    self._sizer_stack.append(self._current_sizer)
    grid_sizer = wx.FlexGridSizer(rows, cols, 0, 0)
    self._current_sizer.Add(grid_sizer, proportion, flags)
    self._current_sizer = grid_sizer

  def start_section(self, label="", proportion=0, label_font=None,
      label_size=None):
    box = wx.StaticBox(self.get_panel(), -1, label, style=wx.NO_BORDER)
    if (label_font is None):
      label_font = self._label_font
    if (label_size is None):
      label_size = self._label_size
    label_font.SetPointSize(label_size)
    box.SetFont(label_font)
    box_sizer = wx.StaticBoxSizer(box, wx.VERTICAL)
    self._current_sizer.Add(box_sizer, proportion, wx.ALL|wx.EXPAND, 5)
    self._sizer_stack.append(self._current_sizer)
    self._current_sizer = box_sizer
    self._boxes[label] = box
    self._box_labels.append(label)

  def end_current(self):
    self._current_sizer = self._sizer_stack.pop()

  def end_box(self):
    self.end_current()

  def end_grid(self):
    self.end_current()

  def end_section(self, pad=0):
    if (pad > 0):
      self._current_sizer.Add((10, pad))
    self.end_current()

  def start_toggle_section(self, label="Details. . .", label2=None,
      collapsed=True, proportion=0, expand=False):
    self.start_vbox(proportion=proportion, expand=expand)
    # TODO???

  #--- control layout
  def Add(self, control, proportion=0, flags=0, border=__default_border__):
    self.get_current_sizer().Add(control, proportion, flags, border)

  def Detach(self, control):
    control.GetSizer().Detach(control)

  def add_expanding_widget(self, control, proportion=0, border=5):
    self.Add(control, proportion, wx.ALL|wx.EXPAND, border)

  def add_centered_widget(self, control):
    self.Add(control, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

  def add_centered_widgets(self, *args):
    for control in args :
      self.add_centered_widget(control)

  def add_row_widget(self, control, proportion=0, border=10):
    self.Add(control, proportion, wx.ALL|wx.EXPAND, border)

  def add_spacer(self, width=10, height=-1):
    assert (width > 0) or (height > 0)
    self._current_sizer.Add((width, height))

  def place_widget_row(self, *args):
    self.start_box()
    for control in args :
      self.add_row_widget(control)
    self.add_spacer(10, 1)
    self.end_box()

  #--- text creation
  def _add_text(self, text, bold=False, mono=False, font_size=None,
      text_color=(0,0,0), width=None, set_explicit_width=False):
    if (font_size is None):
      font_size = self._font_size
    if (width is None):
      width = self._text_width
    size = wx.DefaultSize
    if (width is not None) and (set_explicit_width):
      size = (width, -1)
    text_widget = wx.StaticText(
      parent=self.get_panel(),
      id=-1,
      label=text,
      size=size)
    font = text_widget.GetFont()
    font.SetPointSize(font_size)
    if bold :
      font.SetWeight(wx.FONTWEIGHT_BOLD)
    if mono :
      font.SetFamily(wx.FONTFAMILY_MODERN)
    text_widget.SetFont(font)
    if (text_color != (0,0,0)):
      text_widget.SetForegroundColour((text_color))
    if (width is not None):
      text_widget.Wrap(width)
    self.add_centered_widget(text_widget)
    return text_widget

  def add_text(self, text, **kwds):
    return self._add_text(text, **kwds)

  def add_bold_text(self, text, **kwds):
    kwds['bold'] = True
    return self._add_text(text, **kwds)

  def add_mono_text(self, text, **kwds):
    kwds['mono'] = True
    return self._add_text(text, **kwds)

  def add_italic_text(self, text, **kwds) : # FIXME
    return self._add_text(text, **kwds)

  #--- other controls
  def add_bitmap(self, bmp):
    if isinstance(bmp, wx.Image):
      bmp = bmp.ConvertToBitmap()
    img = wx.StaticBitmap(self.get_panel(), -1, bmp)
    self.add_centered_widget(img)

  def add_line(self, width=None, expand=False):
    if (width is None):
      width = self._text_width
    line = wx.StaticLine(self.get_panel(), -1, size=(width,2),
      style=wx.LI_HORIZONTAL)
    flags = wx.ALIGN_CENTER
    if expand :
      flags |= wx.EXPAND
    self._current_sizer.Add(line, 0, flags, 5)

  #--- debugging
  def show_sizers(self):
    indent = 0
    for sizer in self._sizer_stack :
      print("%s%s" % (" "*indent, type(sizer).__name__))
      indent += 2

########################################################################
# TESTING
if (__name__ == "__main__"):
  class TestPanel(wx.Panel, SizerContainer):
    def __init__(self, *args, **kwds):
      wx.Panel.__init__(self, *args, **kwds)
      SizerContainer.__init__(self)
      self.set_panel(self)
      main_sizer = wx.BoxSizer(wx.VERTICAL)
      self.SetSizer(main_sizer)
      self.set_main_sizer(main_sizer)

  app = wx.App(0)
  frame = wx.Frame(None, -1, "Test frame")
  p = TestPanel(frame, -1)
  p.start_hbox()
  p.add_bold_text("Bold text")
  p.add_mono_text("Mono text")
  p.add_text("Plain text")
  p.end_box()
  p.add_bold_text("Bold text, font size = 14, color = red",
    text_color=(255,0,0),
    font_size=14)
  p.start_section("New section")
  p.start_hbox()
  try :
    import wxtbx.bitmaps
    p.add_bitmap(wxtbx.bitmaps.fetch_icon_bitmap("filesystems", "folder_home"))
    p.add_text("Caption with bitmap")
  except ImportError :
    p.add_text("wxtbx.bitmaps not found, disabled")
  p.end_box()
  p.end_section()
  p.start_section("Another section (proportion=1)", label_size=12, proportion=1)
  p.set_text_width(720)
  p.set_font_size(9)
  p.add_text("  ".join(["The quick brown fox jumped over the lazy dogs."]*10))
  p.reset_font_size()
  p.add_line(expand=True)
  p.add_italic_text("Normal size")
  p.end_section()
  p.main_sizer.Fit(p)
  frame.Fit()
  frame.Show()
  app.MainLoop()
