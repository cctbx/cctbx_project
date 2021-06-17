from __future__ import absolute_import, division, print_function

'''
Author      : Lyubimov, A.Y.
Created     : 04/14/2014
Last Changed: 11/05/2018
Description : wxPython 3-4 compatibility tools

The context managers, classes, and other tools below can be used to make the
GUI code compatible with wxPython 3 and 4. Mostly, the tools convert the
functions, enumerations, and classes which have been renamed in wxPython 4;
the name mismatches result in exceptions.

Use case 1: subclassing wx.PyControl or wx.Control:

from wxtbx import wx4_compatibility as wx4c
WxCtrl = wx4c.get_wx_mod(wx, wx.Control)
class MyCustomControl(WxCtrl): ...


Use case 2: brush style (NOTE: you can do that with fonts as well, but it
doesn't seem to be necessary):

from wxtbx import wx4_compatibility as wx4c
bkgrd = self.GetBackgroundColour()
with wx4c.set_brush_style(wx.BRUSHSTYLE_SOLID) as bstyle:
    brush = wx.Brush(bkgrd, bstyle)


Use case 3: Toolbars

from wxtbx import wx4_compatibility as wx4c, bitmaps
class MyFrame(wx.Frame):
  def __init__(self, parent, id, title, *args, **kwargs):
    wx.Frame.__init__(self, parent, id, title, *args, **kwargs)

    self.toolbar = wx4c.ToolBar(self, style=wx.TB_TEXT)
    self.quit_button = self.toolbar.AddTool(toolId=wx.ID_ANY,
                                             label='Quit',
                                             kind=wx.ITEM_NORMAL,
                                             bitmap=bitmaps.fetch_icon_bitmap('actions', 'exit')
                                             shortHelp='Exit program')
    ...
    self.SetToolBar(self.toolbar)
    self.toolbar.Realize()


'''

import wx
from contextlib import contextmanager
import importlib

wx4 = wx.__version__[0] == '4'

modnames = [
  ('PyControl', 'Control'),
  ('PyDataObjectSimple', 'DataObjectSimple'),
  ('PyDropTarget', 'DropTarget'),
  ('PyEvtHandler', 'EvtHandler'),
  ('PyImageHandler', 'ImageHandler'),
  ('PyLocale', 'Locale'),
  ('PyLog', 'Log'),
  ('PyPanel', 'Panel'),
  ('PyPickerBase', 'PickerBase'),
  ('PyPreviewControlBar', 'PreviewControlBar'),
  ('PyPreviewFrame', 'PreviewFrame'),
  ('PyPrintPreview', 'PrintPreview'),
  ('PyScrolledWindow', 'ScrolledWindow'),
  ('PySimpleApp', 'App'),
  ('PyTextDataObject', 'TextDataObject'),
  ('PyTimer', 'Timer'),
  ('PyTipProvider', 'adv.TipProvider'),
  ('PyValidator', 'Validator'),
  ('PyWindow'', Window')
]

font_families = [
  (wx.DEFAULT, wx.FONTFAMILY_DEFAULT),
  (wx.DECORATIVE, wx.FONTFAMILY_DECORATIVE),
  (wx.ROMAN, wx.FONTFAMILY_ROMAN),
  (wx.SCRIPT, wx.FONTFAMILY_SCRIPT),
  (wx.SWISS, wx.FONTFAMILY_SWISS),
  (wx.MODERN, wx.FONTFAMILY_MODERN),
  (wx.TELETYPE, wx.FONTFAMILY_TELETYPE)
]

font_weights = [
  (wx.NORMAL, wx.FONTWEIGHT_NORMAL),
  (wx.LIGHT, wx.FONTWEIGHT_LIGHT),
  (wx.BOLD, wx.FONTWEIGHT_BOLD)
]

font_styles = [
  (wx.NORMAL, wx.FONTSTYLE_NORMAL),
  (wx.ITALIC, wx.FONTSTYLE_ITALIC),
  (wx.SLANT, wx.FONTSTYLE_SLANT)
]

pen_styles = [
  (wx.SOLID, wx.PENSTYLE_SOLID),
  (wx.DOT, wx.PENSTYLE_DOT),
  (wx.LONG_DASH, wx.PENSTYLE_LONG_DASH),
  (wx.SHORT_DASH, wx.PENSTYLE_SHORT_DASH),
  (wx.DOT_DASH, wx.PENSTYLE_DOT_DASH),
  (wx.USER_DASH, wx.PENSTYLE_USER_DASH),
  (wx.TRANSPARENT, wx.PENSTYLE_TRANSPARENT)
]

brush_styles = [
  (wx.SOLID, wx.BRUSHSTYLE_SOLID),
  (wx.TRANSPARENT, wx.BRUSHSTYLE_TRANSPARENT),
  (wx.STIPPLE_MASK_OPAQUE, wx.BRUSHSTYLE_STIPPLE_MASK_OPAQUE),
  (wx.STIPPLE_MASK, wx.BRUSHSTYLE_STIPPLE_MASK),
  (wx.STIPPLE, wx.BRUSHSTYLE_STIPPLE),
  (wx.BDIAGONAL_HATCH, wx.BRUSHSTYLE_BDIAGONAL_HATCH),
  (wx.CROSSDIAG_HATCH, wx.BRUSHSTYLE_CROSSDIAG_HATCH),
  (wx.FDIAGONAL_HATCH, wx.BRUSHSTYLE_FDIAGONAL_HATCH),
  (wx.CROSS_HATCH, wx.BRUSHSTYLE_CROSS_HATCH),
  (wx.HORIZONTAL_HATCH, wx.BRUSHSTYLE_HORIZONTAL_HATCH),
  (wx.VERTICAL_HATCH, wx.BRUSHSTYLE_VERTICAL_HATCH),
]

def find_module(module):
  for m in modnames:
    if module.__name__ in m:
      return m

def find_enum(enums, item):
  for en in enums:
    if item in en:
      value = en[1] if wx4 else en[0]
      return value

def get_wx_mod(base, module):
  mname = find_module(module)[1] if wx4 else find_module(module)[0]
  bname = base.__name__
  if '.' in mname:
    spl = [i for i in mname.split('.') if i != bname]
    modname = '.'.join(spl[:-1])
    mod = importlib.import_module('{}.{}'.format(bname, modname))
    return getattr(mod, spl[-1])
  else:
    return getattr(base, mname)

@contextmanager
def wx_mod(base, module):
  ''' Identify and import the appropriate wxPython module '''
  yield get_wx_mod(base, module)

@contextmanager
def set_font_style(style):
  yield find_enum(font_styles, style)

@contextmanager
def set_font_weight(weight):
  yield find_enum(font_weights, weight)

@contextmanager
def set_font_family(family):
  yield find_enum(font_families, family)

@contextmanager
def set_pen_style(style):
  yield find_enum(pen_styles, style)

@contextmanager
def set_brush_style(style):
  yield find_enum(brush_styles, style)

@contextmanager
def create_measuring_context():
  dc = wx.GraphicsContext.Create() if wx4 else \
    wx.GraphicsContext.CreateMeasuringContext()
  yield dc

class Wx3ToolBar(wx.ToolBar):
  ''' Special toolbar class that accepts wxPython 4-style AddTool command and
  converts it to a wxPython 3-style AddLabelTool command '''
  def __init__(self, parent, id=wx.ID_ANY, pos=wx.DefaultPosition,
               size=wx.DefaultSize, style=wx.TB_HORIZONTAL, name='toolbar'):
    wx.ToolBar.__init__(self, parent, id, pos, size, style, name)

  def AddTool(self, toolId, label, bitmap, bmpDisabled=wx.NullBitmap,
              kind=wx.ITEM_NORMAL, shortHelp='', longHelp='',
              clientData=None):
    ''' Override to make this a very thin wrapper for AddLabelTool, which in
    wxPython 3 is the same as AddTool in wxPython 4 '''
    return self.AddLabelTool(id=toolId, label=label, bitmap=bitmap,
                             bmpDisabled=bmpDisabled, kind=kind,
                             shortHelp=shortHelp, longHelp=longHelp,
                             clientData=clientData)

class Wx4ToolBar(wx.ToolBar):
  ''' Special toolbar class that accepts wxPython 3-style AddLabelTool command
  and converts it to a wxPython 4-style AddTool command '''
  def __init__(self, parent, id=wx.ID_ANY, pos=wx.DefaultPosition,
               size=wx.DefaultSize, style=wx.TB_HORIZONTAL, name='toolbar'):
    wx.ToolBar.__init__(self, parent, id, pos, size, style, name)

  def AddLabelTool(self, id, label, bitmap, bmpDisabled=wx.NullBitmap,
              kind=wx.ITEM_NORMAL, shortHelp='', longHelp='',
              clientData=None):
    ''' Override to make this a very thin wrapper for AddTool, which in
    wxPython 4 is the same as AddLabelTool in wxPython 3 '''
    return self.AddTool(toolId=id, label=label, bitmap=bitmap,
                        bmpDisabled=bmpDisabled, kind=kind,
                        shortHelp=shortHelp, longHelp=longHelp,
                        clientData=clientData)

# Use this ToolBar class to create toolbars in frames
ToolBar = Wx4ToolBar if wx4 else Wx3ToolBar
