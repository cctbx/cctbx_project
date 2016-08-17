from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 05/02/2014
Last Changed: 08/16/2016
Description : IOTA Custom Widgets
'''

import os
import wx
from wxtbx import metallicbutton as mb

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

icons = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'icons/')


class GradButton(mb.MetallicButton):
  def __init__(self, parent, label='', bmp=None, size=wx.DefaultSize,
               style=mb.MB_STYLE_BOLD_LABEL, handler_function=None,
               user_data=None, start_color=(218, 218, 218),
               gradient_percent=15.0, highlight_color=(230, 230, 230),
               label_size=LABEL_SIZE, caption_size=CAPTION_SIZE,
               button_margin=4, disable_after_click=0) :
    if isinstance(bmp, str) :
      bmp = self.StandardBitmap(bmp)
    bmp_size = bmp.GetSize()
    if bmp_size > size[1]:
      size = (size[0], 1.5 * bmp_size[1])
    mb.MetallicButton.__init__(self,
                               parent=parent,
                               label=label,
                               bmp=bmp,
                               size=size,
                               style=style,
                               name=str(user_data),
                               start_color=start_color,
                               gradient_percent=gradient_percent,
                               highlight_color=highlight_color,
                               label_size=label_size,
                               caption_size=caption_size,
                               button_margin=button_margin,
                               disable_after_click=disable_after_click
                               )
    if handler_function is not None:
     self.bind_event(wx.EVT_BUTTON, handler_function)

  def StandardBitmap(img_name, size=None):
    img_path = img_name
    img = wx.Image(img_path, type=wx.BITMAP_TYPE_ANY, index=-1)
    if size is not None:
     (w, h) = size
     img.Rescale(w, h)
    bmp = img.ConvertToBitmap()
    return bmp


# TODO: make these classes even more generic

class CtrlBase(wx.Panel):
  ''' Control panel base class '''
  def __init__(self,
               parent,
               label_style='normal',
               content_style='normal'):

    wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)
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

class InputCtrl(CtrlBase):
  ''' Generic panel that will place a text control, with a label and an
      optional Browse / magnifying-glass buttons into a window'''

  def __init__(self, parent,
               label='', label_size=(100, -1),
               label_style='normal',
               button=False, value=''):

    CtrlBase.__init__(self, parent=parent, label_style=label_style)

    output_box = wx.FlexGridSizer(1, 4, 0, 10)
    self.txt = wx.StaticText(self, label=label, size=label_size)
    self.txt.SetFont(self.font)
    output_box.Add(self.txt)

    self.ctr = wx.TextCtrl(self) #, size=ctr_size)
    self.ctr.SetValue(value)
    output_box.Add(self.ctr, flag=wx.EXPAND)

    self.btn_browse = wx.Button(self, label='Browse...')
    self.btn_mag = wx.BitmapButton(self,
                                   bitmap=wx.Bitmap('{}/16x16/viewmag.png'
                                                    ''.format(icons)))
    output_box.Add(self.btn_browse, flag=wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    output_box.Add(self.btn_mag, flag=wx.RESERVE_SPACE_EVEN_IF_HIDDEN)

    if not button:
      self.btn_browse.Hide()
      self.btn_mag.Hide()

    output_box.AddGrowableCol(1, 1)
    self.SetSizer(output_box)

class OptionCtrl(CtrlBase):
  ''' Generic panel will place a text control w/ label '''
  def __init__(self, parent, items,
               label='',
               label_size=(100, -1),
               label_style='normal',
               sub_labels=[],
               ctrl_size=(300, -1)):

    CtrlBase.__init__(self, parent=parent, label_style=label_style)

    if label != '':
      opt_box = wx.FlexGridSizer(1, len(items) * 2 + 1, 0, 10)
      self.txt = wx.StaticText(self, label=label, size=label_size)
      self.txt.SetFont(self.font)
      opt_box.Add(self.txt, flag=wx.ALIGN_CENTER_VERTICAL)
    else:
      opt_box = wx.FlexGridSizer(1, len(items) * 2, 0, 10)

    for key, value in items:
      if sub_labels != []:
        sub_label = sub_labels[items.index((key, value))].decode('utf-8')
      else:
        sub_label = key

      if len(items) > 1:
        opt_label = wx.StaticText(self, id=wx.ID_ANY, label=sub_label)
        opt_box.Add(opt_label, flag=wx.ALIGN_CENTER_VERTICAL)

      item = wx.TextCtrl(self, id=wx.ID_ANY, size=ctrl_size,
                         style=wx.TE_PROCESS_ENTER)
      item.SetValue(str(value))
      opt_box.Add(item, flag=wx.ALIGN_CENTER_VERTICAL)
      self.__setattr__(key, item)

    self.SetSizer(opt_box)


class SpinCtrl(CtrlBase):
  ''' Generic panel will place a spin control w/ label '''
  def __init__(self, parent,
               label='',
               label_size=(200, -1),
               label_style='normal',
               ctrl_size=(60, -1)):

    CtrlBase.__init__(self, parent=parent, label_style=label_style)

    ctr_box = wx.FlexGridSizer(1, 3, 0, 10)
    self.txt = wx.StaticText(self, label=label, size=label_size)
    self.txt.SetFont(self.font)
    self.ctr = wx.SpinCtrl(self, value='3', max=(20), min=(3), size=ctrl_size)
    ctr_box.Add(self.txt)
    ctr_box.Add(self.ctr)

    self.SetSizer(ctr_box)

class ChoiceCtrl(CtrlBase):
  ''' Generic panel will place a choice control w/ label '''

  def __init__(self, parent,
               choices,
               label='',
               label_size=(200, -1),
               label_style='normal',
               ctrl_size=(100, -1)):

    CtrlBase.__init__(self, parent=parent, label_style=label_style)
    ctr_box = wx.FlexGridSizer(1, 3, 0, 10)
    self.txt = wx.StaticText(self, label=label, size=label_size)
    self.txt.SetFont(self.font)
    self.ctr = wx.Choice(self, size=ctrl_size, choices=choices)
    ctr_box.Add(self.txt)
    ctr_box.Add(self.ctr)

    self.SetSizer(ctr_box)

class TableCtrl(CtrlBase):
  ''' Generic panel will place a table w/ x and y labels
      Data must be a list of lists for multi-column tables '''

  def __init__(self, parent,
               clabels=[],
               clabel_size=(200, -1),
               rlabels=[],
               rlabel_size=(200, -1),
               contents=[],
               label_style='normal',
               content_style='normal'):

    CtrlBase.__init__(self, parent=parent, label_style=label_style,
                      content_style=content_style)
    nrows = len(rlabels) + 1

    if len(clabels) == 0:
      ncols = 2
    else:
      ncols = len(clabels) + 1
    self.sizer = wx.FlexGridSizer(nrows, ncols, 10, 10)

    # add column labels (xlabels)
    if len(clabels) > 0:
      self.sizer.Add(wx.StaticText(self, label=''))
      for item in column_labels:
        clabel = wx.StaticText(self, label=i.decode('utf-8'), size=clabel_size)
        clabel.SetFont(self.font)
        self.sizer.Add(clabel)

    # add row labels and row contents
    for l in rlabels:
      row_label = wx.StaticText(self, label=l.decode('utf-8'), size=rlabel_size)
      row_label.SetFont(self.font)
      self.sizer.Add(row_label)

      # Add data to table
      c_index = rlabels.index(l)
      for item in contents[c_index]:
        cell = wx.StaticText(self, label=item.decode('utf-8'))
        cell.SetFont(self.cfont)
        self.sizer.Add(cell)

    self.SetSizer(self.sizer)

class PHILBox(CtrlBase):
  def __init__(self, parent,
               btn_import=True,
               btn_import_size=(120, -1),
               btn_import_label='Import PHIL',
               btn_export=False,
               btn_export_size=(120, -1),
               btn_export_label='Export PHIL',
               btn_default=True,
               btn_default_size=(120, -1),
               btn_default_label='Default PHIL',
               ctr_size=(-1, 125),
               ctr_value='',
               label_style='normal',
               content_style='normal'):

    CtrlBase.__init__(self, parent=parent, label_style=label_style,
                      content_style=content_style)

    self.sizer = wx.GridBagSizer(5, 5)
    self.SetSizer(self.sizer)
    self.sizer.AddGrowableCol(1)

    self.ctr = wx.richtext.RichTextCtrl(self,
                                        size=ctr_size,
                                        style=wx.VSCROLL,
                                        value=ctr_value)
    span_counter = 0
    if btn_import:
      self.btn_import = wx.Button(self,
                                  label=btn_import_label,
                                  size=btn_import_size)
      self.sizer.Add(self.btn_import, pos=(span_counter, 0))
      span_counter += 1
    if btn_export:
      self.btn_export = wx.Button(self,
                                  label=btn_export_label,
                                  size=btn_export_size)
      self.sizer.Add(self.btn_export, pos=(span_counter, 0))
      span_counter += 1
    if btn_default:
      self.btn_default = wx.Button(self,
                                   label=btn_default_label,
                                   size=btn_default_size)
      self.sizer.Add(self.btn_default, pos=(span_counter, 0))
      span_counter += 1

    if span_counter > 0:
      self.sizer.AddGrowableRow(span_counter)
      self.sizer.Add(self.ctr, pos=(0, 1), span=(span_counter + 1, 1),
                     flag=wx.EXPAND)
    elif span_counter == 0:
      self.sizer = wx.BoxSizer(wx.VERTICAL)
      self.sizer.Add(self.ctr, 1, flag=wx.EXPAND)
