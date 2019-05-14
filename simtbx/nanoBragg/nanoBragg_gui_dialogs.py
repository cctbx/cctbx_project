from __future__ import absolute_import, division, print_function

'''
Author      : Lyubimov, A.Y.
Created     : 12/12/2017
Last Changed: 12/12/2017
Description : IOTA GUI Dialogs
'''

import wx
from iota.components.iota_utils import UnicodeCharacters, WxFlags


# Platform-specific stuff
# TODO: Will need to test this on Windows at some point
if wx.Platform == '__WXGTK__':
  norm_font_size = 10
  button_font_size = 12
  LABEL_SIZE = 14
  CAPTION_SIZE = 12
elif wx.Platform == '__WXMAC__':
  norm_font_size = 12
  button_font_size = 14
  LABEL_SIZE = 14
  CAPTION_SIZE = 12
elif (wx.Platform == '__WXMSW__'):
  norm_font_size = 9
  button_font_size = 11
  LABEL_SIZE = 11
  CAPTION_SIZE = 9

# Initialize unicode font and wx flags
u = UnicodeCharacters()
f = WxFlags()


# ---------------------------------------------------------------------------- #

class BaseDialog(wx.Dialog):
  def __init__(self, parent, style=wx.DEFAULT_DIALOG_STYLE,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):
    wx.Dialog.__init__(self, parent, style=style, *args, **kwargs)

    self.envelope = wx.BoxSizer(wx.VERTICAL)
    self.main_sizer = wx.BoxSizer(wx.VERTICAL)
    self.envelope.Add(self.main_sizer, 1, flag=wx.EXPAND | wx.ALL, border=5)
    self.SetSizer(self.envelope)

    if label_style == 'normal':
      self.font = wx.Font(norm_font_size, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
    elif label_style == 'bold':
      self.font = wx.Font(norm_font_size, wx.DEFAULT, wx.NORMAL, wx.BOLD)
    elif label_style == 'italic':
      self.font = wx.Font(norm_font_size, wx.DEFAULT, wx.ITALIC, wx.NORMAL)
    elif label_style == 'italic_bold':
      self.font = wx.Font(norm_font_size, wx.DEFAULT, wx.ITALIC, wx.BOLD)

    if content_style == 'normal':
      self.cfont = wx.Font(norm_font_size, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
    elif content_style == 'bold':
      self.cfont = wx.Font(norm_font_size, wx.DEFAULT, wx.NORMAL, wx.BOLD)
    elif content_style == 'italic':
      self.cfont = wx.Font(norm_font_size, wx.DEFAULT, wx.ITALIC, wx.NORMAL)
    elif content_style == 'italic_bold':
      self.cfont = wx.Font(norm_font_size, wx.DEFAULT, wx.ITALIC, wx.BOLD)

class BeamlineOptions(BaseDialog):
  def __init__(self, parent, phil=None,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):
    dlg_style = wx.CAPTION | wx.CLOSE_BOX | wx.RESIZE_BORDER | wx.STAY_ON_TOP

    BaseDialog.__init__(self, parent, style=dlg_style,
                        label_style=label_style,
                        content_style=content_style,
                        size=(600, 500),
                        *args, **kwargs)
