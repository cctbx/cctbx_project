from __future__ import division, print_function, absolute_import

'''
Author      : Lyubimov, A.Y.
Created     : 11/15/2018
Last Changed: 11/15/2018
Description : IOTA GUI base classes (with backwards compatibility for
              wxPython 3)
'''

import wx
from wxtbx import bitmaps

class IOTAFrameError(Exception):
  def __init__(self, msg):
    Exception.__init__(self, msg)

class IOTABaseFrame(wx.Frame):
  """ New frame that will show processing info """

  def __init__(self, parent, id, title, *args, **kwargs):
    wx.Frame.__init__(self, parent, id, title, *args, **kwargs)
    self.parent = parent

    self.main_sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.main_sizer)

  def initialize_toolbar(self):
    self.toolbar = self.CreateToolBar(style=wx.TB_3DBUTTONS | wx.TB_TEXT)

  def add_tool(self, id=wx.ID_ANY, label=None, bitmap=None,
               kind=wx.ITEM_NORMAL, shortHelp=None, longHelp=None):
    if type(bitmap) in (tuple, list):
      if len(bitmap) == 3:
        bsz = bitmap[2]
        bsc = (bitmap[2], bitmap[2])
      else:
        bsz = bsc = None
      if bitmap[0] == 'custom':
        bmp = bitmaps.fetch_custom_icon_bitmap(bitmap[1], size=bsz)
      else:
        bmp = bitmaps.fetch_icon_bitmap(bitmap[0], bitmap[1], scale=bsc)
    else:
      bmp = bitmap

    short_string = shortHelp if shortHelp else ''
    long_string = longHelp if longHelp else ''

    if 'classic' in wx.version():
      return self.toolbar.AddLabelTool(id=id,
                                       label=label,
                                       kind=kind,
                                       bitmap=bmp,
                                       shortHelp=short_string,
                                       longHelp=long_string)
    elif 'phoenix' in wx.version():
      return self.toolbar.AddTool(toolId=id,
                                  label=label,
                                  kind=kind,
                                  bitmap=bmp,
                                  shortHelp=short_string)
    else:
      raise IOTAFrameError(msg='Unrecognized wxPython version {}'.format(
        wx.version()))

  def add_toolbar_separator(self, stretch=False):
    if stretch and 'phoenix' in wx.version():
      self.toolbar.AddStretchableSpace()
    else:
      self.toolbar.AddSeparator()

  def realize_toolbar(self):
    self.toolbar.Realize()

  def set_tool_states(self, tools):
    for tool in tools:
      # toggle = tool[2] if len(tool) == 3 else None
      self.set_tool_state(tool=tool[0], enable=tool[1],
                          toggle=tool[2] if len(tool) == 3 else None)

  def set_tool_state(self, tool, enable=None, toggle=None):
    id = tool.GetId()
    if enable is not None:
      self.toolbar.EnableTool(id, enable)
    if toggle is not None:
      self.toolbar.ToggleTool(id, toggle)


  def place_and_size(self, set_size=False, set_by=None, center=False):
    """ Place and size the frame"""

    # Determine effective minimum size
    if set_size:
      self.SetMinSize(self.GetEffectiveMinSize())

    # Find mouse position
    if set_by == 'mouse':
      self.SetPosition(wx.GetMousePosition())
    elif set_by == 'parent':
      self.SetPosition(self.set_relative_position())
    else:
      self.SetPosition((0, 0))

    # Center on display
    if center:
      self.Center()

  def set_relative_position(self):
    """ Determines screen position w/ respect to parent window; will also
    detect if it goes beyond the display edge, and adjust """

    # Position proc window w/ respect to IOTA window
    mx, my = self.parent.GetPosition()
    px = mx + 50
    py = my + 50

    # Calculate if proc window is going out of bounds, and adjust
    disp_idx = wx.Display.GetFromWindow(self.parent)
    disp_geom = wx.Display(disp_idx).GetClientArea()
    dxmin = disp_geom[0]
    dxmax = disp_geom[0] + disp_geom[2]
    dymin = disp_geom[1]
    dymax = disp_geom[1] + disp_geom[3]

    pw, pl = self.GetSize()
    if not (px + pw * 1.1 in range(dxmin, dxmax)):
      px = dxmax - pw * 1.1
    if not (py + pl * 1.1 in range(dymin, dymax)):
      py = dymax - pl * 1.1

    return (px, py)


class IOTABasePanel(wx.Panel):
  def __init__(self, parent, *args, **kwargs):
    wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY, size=(800, 500),
                      *args, **kwargs)

    self.main_sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.main_sizer)
