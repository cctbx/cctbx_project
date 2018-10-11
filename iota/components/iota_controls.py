from __future__ import division, print_function, absolute_import
from past.builtins import range

'''
Author      : Lyubimov, A.Y.
Created     : 07/08/2016
Last Changed: 08/29/2018
Description : IOTA GUI controls
'''

import os
import wx
import wx.richtext
import wx.lib.agw.floatspin as fs
import wx.lib.agw.ultimatelistctrl as ulc
import wx.lib.agw.knobctrl as kc
from wx.lib.mixins.listctrl import ListCtrlAutoWidthMixin, ColumnSorterMixin
from wxtbx import metallicbutton as mb
from wxtbx import bitmaps
import wx.lib.buttons as btn

import numpy as np

try:  # for Py3 compatibility
    import itertools.izip as zip
except ImportError:
    pass

from iota.components.iota_misc import noneset

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

# --------------------------------- Widgets ---------------------------------- #

def StandardBitmap(img_name, size=None):
  img_path = img_name
  img = wx.Image(img_path, type=wx.BITMAP_TYPE_ANY, index=-1)
  if size is not None:
    (w, h) = size
    img.Rescale(w, h)
  bmp = img.ConvertToBitmap()
  return bmp

class GradButton(mb.MetallicButton):
  ''' Customized MetallicButton '''

  def __init__(self, parent, label='', bmp=None, size=wx.DefaultSize,
               style=mb.MB_STYLE_BOLD_LABEL,
               handler_function=None,
               user_data=None, start_color=(218, 218, 218),
               gradient_percent=15, highlight_color=(250, 250, 250),
               label_size=LABEL_SIZE, caption_size=CAPTION_SIZE,
               button_margin=4, disable_after_click=0):
    if isinstance(bmp, str):
      bmp = StandardBitmap(bmp)
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

class MiniButtonBox(wx.Panel):
  ''' A box with three mini buttons for IOTA panel '''

  def __init__(self, parent):
    wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)

    self.index = None
    self.btn_box = wx.BoxSizer(wx.HORIZONTAL)
    view_bmp = bitmaps.fetch_custom_icon_bitmap('image_viewer16')
    self.btn_viewer = btn.GenBitmapButton(self, bitmap=view_bmp)
    viewmag_bmp = bitmaps.fetch_icon_bitmap('actions', 'viewmag', size=16)
    self.btn_mag = btn.GenBitmapButton(self, bitmap=viewmag_bmp)
    star_on_bmp = bitmaps.fetch_icon_bitmap('actions', 'bookmark_gold',
                                            size=16)
    star_off_bmp = bitmaps.fetch_icon_bitmap('actions', 'bookmark_silver',
                                             size=16)
    self.btn_flag = btn.GenBitmapToggleButton(self, bitmap=star_off_bmp)
    self.btn_flag.SetBitmapSelected(star_on_bmp)
    self.btn_box.Add(self.btn_mag)
    self.btn_box.Add(self.btn_viewer, flag=wx.LEFT, border=5)
    self.btn_box.Add(self.btn_flag, flag=wx.LEFT, border=5)
    self.SetSizer(self.btn_box)
    self.SetBackgroundColour('white')
    self.Fit()

class MiniButtonBoxInput(wx.Panel):
  ''' A box with three mini buttons for Input file panel '''

  def __init__(self, parent):
    wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)

    self.index = None
    self.btn_box = wx.BoxSizer(wx.HORIZONTAL)

    viewmag_bmp = bitmaps.fetch_icon_bitmap('actions', 'viewmag', size=16)
    # self.btn_mag = btn.GenBitmapButton(self, bitmap=viewmag_bmp)
    self.btn_mag =GradButton(self, bmp=viewmag_bmp, size=(31, 30),
                             gradient_percent=0)

    del_bmp = bitmaps.fetch_icon_bitmap('actions', 'editdelete', size=16)
    # self.btn_delete = btn.GenBitmapButton(self, bitmap=del_bmp)
    self.btn_delete = GradButton(self, bmp=del_bmp, size=(31, 30),
                                 gradient_percent=0)

    info_bmp = bitmaps.fetch_icon_bitmap('actions', 'info', size=16)
    # self.btn_info = btn.GenBitmapButton(self, bitmap=info_bmp)
    self.btn_info = GradButton(self, bmp=info_bmp, size=(31, 30),
                               gradient_percent=0)

    self.btn_box.Add(self.btn_mag)
    self.btn_box.Add(self.btn_delete, flag=wx.LEFT, border=5)
    self.btn_box.Add(self.btn_info, flag=wx.LEFT, border=5)
    self.SetSizer(self.btn_box)
    self.SetBackgroundColour('white')
    self.Fit()


class DataTypeChoice(wx.Panel):
  def __init__(self, parent, choices):
    wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)
    self.index = None
    self.type = wx.Choice(self, -1, choices=choices)
    self.Fit()


# ---------------------------------- Inputs ---------------------------------- #

class InputListCtrl(ulc.UltimateListCtrl, ListCtrlAutoWidthMixin):
  ''' Customized UltimateListCtrl with auto-width mixin'''

  def __init__(self, parent, ID, n_cols=3, pos=wx.DefaultPosition,
               size=wx.DefaultSize, style=0):
    ulc.UltimateListCtrl.__init__(self, parent, ID, pos,
                                  size=size,
                                  agwStyle=style)
    ListCtrlAutoWidthMixin.__init__(self)


class VirtualInputListCtrl(ulc.UltimateListCtrl, ListCtrlAutoWidthMixin,
                           ColumnSorterMixin):
  ''' Customized Virtual UltimateListCtrl with auto-width mixin'''

  def __init__(self, parent, ID, n_cols=4, pos=wx.DefaultPosition,
               size=wx.DefaultSize, style=0):
    ulc.UltimateListCtrl.__init__(self, parent, ID, pos, size, agwStyle=style)
    ListCtrlAutoWidthMixin.__init__(self)
    ColumnSorterMixin.__init__(self, n_cols)

    self.data = {}

  def InitializeDataMap(self, data):
    self.data = data
    self.itemDataMap = self.data
    self.itemIndexMap = self.data.keys()
    self.SetItemCount(len(self.data))

  def GetListCtrl(self):
    return self

  def OnColClick(self, e):
    e.Skip()

  def OnGetItemToolTip(self, item, col):
    return None

  def OnGetItemTextColour(self, item, col):
    return None

  def OnGetItemText(self, item, col):
    index = self.itemIndexMap[item]
    s = str(self.itemDataMap[index][col])

    if os.path.isfile(s):
      s = os.path.basename(s)

    return s

  def OnGetItemAttr(self, item):
    return None

  def OnGetItemImage(self, item):
    return -1

  def getColumnText(self, index, col):
    item = self.GetItem(index, col)
    return item.GetText()

class InputListItem(object):
  ''' Class that will contain all the elements of an input list entry '''
  def __init__(self, path, type, buttons):
    self.id = None
    self.path = path
    self.type = type
    self.type_selection = 0
    self.buttons = buttons
    self.warning = False
    self.info = {'WARNING':'Unrecognized file format! \n '
                           'We made the best guess. Please select file format '
                           'from the drop-down menu before proceding or '
                           'choose a different file.',
                 'raw image folder':'This folder contains raw diffraction '
                                     'images.',
                 'image pickle folder':'This folder contains diffraction '
                                        'images converted into Python pickle ',
                 'image list':'This file contains a list of paths to '
                              'diffraction images (raw or pickles)',
                 'processed pickle folder':'This folder contains processed '
                                            'diffraction images in Python '
                                            'pickle format',
                 'processed pickle list':'This file contains paths to '
                                          'processed diffraction images in '
                                          'Python pickle format',
                 'image object folder':'This folder contains "image '
                                        'objects", which are Python pickle '
                                        'files containing information about '
                                        'individual raw images, which have '
                                        'already been processed. Use these if '
                                        'you only intend to run the "selection" '
                                        'part of the IOTA process.',
                 'raw image file':'This is a single raw diffraction image file',
                 'image pickle file':'This is a single diffraction image file '
                                'converted to Python pickle format',
                 'image object':'This is a Python pickle file containing '
                                'information about a single diffraction '
                                'image, which has already been processed.',
                 'processed image':'This is a Python pickle file containing '
                                   'integrated intensities',
                 'reference MTZ': 'Merged intensities in MTZ format used as '
                                  'an isomorphous reference for scaling, '
                                  'post-refinement and merging.',
                 'sequence': 'Sequence of the protein of interest, used for '
                             'pseudo-Wilson scaling',
                 'coordinates': 'Structure coordinates in PDB format',
                 'structure factors': 'Structure factors in MTZ format',
                 'background' : 'Radially averaged background (in sin(theta) '
                                '/ lambda'

    }

class FileListItem(object):
  ''' Class that will contain all the elements of a file list entry '''
  def __init__(self, path, items=None):
    ''' A generic control to generate an entry in the file list control
    :param path: absolute path to file
    :param items: a dictionary of items that will be converted to attributes
    '''
    self.id = None
    self.flag = None
    self.path = path
    self.warning = False

    if items is not None:
      for key, value in items.items():
        self.__setattr__(key, value)

# --------------------------------- Controls --------------------------------- #

class CtrlBase(wx.Panel):
  ''' Control panel base class '''
  def __init__(self,
               parent,
               label_style='normal',
               content_style='normal',
               size=wx.DefaultSize):

    wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY, size=size)
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

class DialogButtonsCtrl(CtrlBase):
  ''' Customizable "bottom of window" set of buttons for dialogs '''

  def __init__(self, parent,
               preset = None,
               buttons = None,
               button_size = wx.DefaultSize,
               choices = None,
               choice_label = None,
               choice_size = wx.DefaultSize):
    CtrlBase.__init__(self, parent=parent)

    main_sizer = wx.FlexGridSizer(1, 3, 0, 15)
    main_sizer.AddStretchSpacer()

    # Presets that make for convenient shorthand
    if preset == 'OK_CANCEL':
      buttons = [('OK', wx.ID_OK), ('Cancel', wx.ID_CANCEL)]
    elif preset == 'YES_NO':
      buttons = [('Yes', wx.ID_YES), ('No', wx.ID_NO)]
    elif preset == 'CLOSE':
      buttons = [('Close', wx.ID_CLOSE)]
    if preset == 'PROC_DIALOG':
      buttons = [('OK', wx.ID_OK), ('Cancel', wx.ID_CANCEL)]
      choices = ['Basic', 'Advanced', 'Developer']
      choice_label = 'Expert Level: '
      choice_size = (150, -1)

    if choices is not None:
      choice_sizer = wx.BoxSizer(wx.HORIZONTAL)
      self.choice_txt = wx.StaticText(self, label=choice_label,
                                      size=wx.DefaultSize)
      self.choice = wx.Choice(self, size=choice_size, choices=choices)
      self.choice.SetSelection(0)
      choice_sizer.Add(self.choice_txt, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL,
                       border=5)
      choice_sizer.Add(self.choice)
      main_sizer.Add(choice_sizer, flag=wx.ALL, border=10)

    if buttons is not None:
      btn_sizer = wx.FlexGridSizer(1, len(buttons), 0, 0)
      for key, id in buttons:
        btn_name = 'btn_{}'.format(str.lower(key))
        button = wx.Button(self, label=key, id=id, size=button_size)
        self.__setattr__(btn_name, button)
        btn_sizer.Add(button, flag=wx.RIGHT, border=5)
      main_sizer.Add(btn_sizer, flag=wx.ALL, border=10)

    main_sizer.AddGrowableCol(0)
    self.SetSizer(main_sizer)

class TextCtrlWithButtons(CtrlBase):
  ''' Text control with multiple buttons '''
  def __init__(self, parent,
               buttons = None,
               button_size = wx.DefaultSize,
               ctrl_value = '',
               ctrl_size = wx.DefaultSize,
               ctrl_label = ''):
    CtrlBase.__init__(self, parent=parent)

    main_sizer = wx.FlexGridSizer(1, 3, 0, 15)

    self.txt_label = wx.StaticText(self, label=ctrl_label, size=ctrl_size)
    self.txt_ctrl = wx.TextCtrl(self)
    self.txt_ctrl.SetValue(ctrl_value)
    main_sizer.Add(self.txt_label)
    main_sizer.Add(self.txt_ctrl, flag=wx.EXPAND)

    if buttons is not None:
      btn_sizer = wx.FlexGridSizer(1, len(buttons), 0, 0)
      for key, id in buttons:
        btn_name = 'btn_{}'.format(str.lower(key))
        button = wx.Button(self, label=key, id=id, size=button_size)
        self.__setattr__(btn_name, button)
        btn_sizer.Add(button, flag=wx.RIGHT, border=5)
      main_sizer.Add(btn_sizer)

    main_sizer.AddGrowableCol(1)
    self.SetSizer(main_sizer)

class InputCtrl(CtrlBase):
  ''' Generic panel that will place a text control with a label '''

  def __init__(self, parent,
               label='', label_size=(100, -1),
               label_style='normal',
               buttons=False, value=''):

    CtrlBase.__init__(self, parent=parent, label_style=label_style)

    self.value=value

    output_box = wx.FlexGridSizer(1, 4, 0, 10)
    self.txt = wx.StaticText(self, label=label, size=label_size)
    self.txt.SetFont(self.font)
    output_box.Add(self.txt)

    self.ctr = wx.TextCtrl(self)
    self.ctr.SetValue(self.value)
    output_box.Add(self.ctr, flag=wx.EXPAND)

    self.btn_browse = wx.Button(self, label='Browse...')
    viewmag_bmp = bitmaps.fetch_icon_bitmap('actions', 'viewmag', size=16)
    self.btn_mag = wx.BitmapButton(self, bitmap=viewmag_bmp)
    output_box.Add(self.btn_browse, flag=wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    output_box.Add(self.btn_mag, flag=wx.RESERVE_SPACE_EVEN_IF_HIDDEN)

    if not buttons:
      self.btn_browse.Hide()
      self.btn_mag.Hide()

    output_box.AddGrowableCol(1, 1)
    self.SetSizer(output_box)

  def reset_default(self):
    self.ctr.SetValue(self.value)

class ChoiceCtrl(CtrlBase):
  ''' Generic panel will place a choice control w/ label and a text control (
  optional) which will activate when a 'custom' setting is selected '''

  def __init__(self, parent,
               choices,
               custom_choices=[],
               label='',
               label_size=(200, -1),
               label_style='normal',
               ctrl_size=(100, -1),
               custom_ctrl_size=(200, -1)):

    CtrlBase.__init__(self, parent=parent, label_style=label_style)

    self.choices = choices
    self.custom_choices = custom_choices
    self.custom_value = None

    # Extend the sizer if a custom choice(s) specified
    if len(custom_choices) > 0:
      ctr_box = wx.FlexGridSizer(1, 4, 0, 10)
      ctr_box.AddGrowableCol(3)
    else:
      ctr_box = wx.FlexGridSizer(1, 3, 0, 10)
    self.SetSizer(ctr_box)

    self.txt = wx.StaticText(self, label=label, size=label_size)
    self.txt.SetFont(self.font)

    # Check if choices are tuples, extract data and assign to items if so
    if all(isinstance(i, tuple) for i in self.choices):
      items = [i[0] for i in choices]
      self.ctr = wx.Choice(self, size=ctrl_size, choices=items)
      for choice in self.choices:
        item_idx = self.ctr.FindString(choice[0])
        self.ctr.SetClientData(item_idx, choice[1])
    else:
      self.ctr = wx.Choice(self, size=ctrl_size, choices=self.choices)

    ctr_box.Add(self.txt, flag=wx.ALIGN_CENTER_VERTICAL)
    ctr_box.Add(self.ctr, flag=wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT)

    # Create binding for choice in case custom choice(s) specified
    if len(custom_choices) > 0:
      self.custom = wx.TextCtrl(self, size=custom_ctrl_size)
      self.custom.Disable()
      ctr_box.Add(self.custom, flag=wx.RIGHT | wx.LEFT | wx.EXPAND, border=10)
      self.Bind(wx.EVT_CHOICE, self.onCustomChoice, self.ctr)

  def onCustomChoice(self, e):
    self.set_choice(custom=self.custom_value)

  def set_choice(self, custom=None):
    if self.ctr.GetString(self.ctr.GetCurrentSelection()) in self.custom_choices:
      self.custom.Enable()
      if custom is not None:
        self.custom_value = custom
        self.custom.SetValue(self.custom_value)
      else:
        self.custom.SetValue('')
    else:
      self.custom.Disable()
      self.custom.SetValue('')

  def reset_default(self):
    self.ctr.SetSelection(0)
    if self.custom_choices != []:
      self.custom.Disable()
      self.custom.SetValue('')

class CheckListCtrl(CtrlBase):
  def __init__(self, parent,
               choices,
               label='',
               label_size=(200, -1),
               label_style='normal',
               ctrl_size=(150, -1),
               direction='horizontal'):

    CtrlBase.__init__(self, parent=parent, label_style=label_style)

    self.txt = wx.StaticText(self, label=label, size=label_size)
    self.txt.SetFont(self.font)
    self.ctr = wx.CheckListBox(self, size=ctrl_size, choices=choices)

    if label == '':
      ctr_box = wx.BoxSizer(wx.VERTICAL)
    else:
      if direction == 'horizontal':
        ctr_box = wx.FlexGridSizer(1, 2, 0, 10)
      elif direction == 'vertical':
        ctr_box = wx.FlexGridSizer(2, 1, 10, 0)
      ctr_box.Add(self.txt, flag=wx.ALIGN_CENTER_VERTICAL)

    ctr_box.Add(self.ctr, proportion=1,
                flag=wx.ALIGN_CENTER_VERTICAL | wx.EXPAND)

    self.SetSizer(ctr_box)

  def update(self, choices):
    self.ctr.Clear()
    for choice in choices:
      self.ctr.Append(choice)

class SpinCtrl(CtrlBase):
  ''' Generic panel will place a spin control w/ label '''
  def __init__(self, parent,
               label='',
               label_size=wx.DefaultSize,
               label_style='normal',
               checkbox=False,
               checkbox_state=False,
               checkbox_label='',
               ctrl_size=wx.DefaultSize,
               ctrl_value='3',
               ctrl_max=999999,
               ctrl_min=0,
               ctrl_step=1,
               ctrl_digits=0):

    CtrlBase.__init__(self, parent=parent, label_style=label_style)

    self.value = ctrl_value
    self.checkbox_state = checkbox_state
    self.toggle = None
    cols = 3

    if checkbox:
      assert checkbox_label != ''
      label = ''
      cols += 1
      self.toggle = wx.CheckBox(self, label=checkbox_label, size=label_size)
      self.toggle.SetValue(self.checkbox_state)

    ctr_box = wx.FlexGridSizer(1, cols, 0, 5)

    self.txt = wx.StaticText(self, label=label.decode('utf-8'),
                             size=label_size)
    self.txt.SetFont(self.font)
    self.ctr = fs.FloatSpin(self, value=ctrl_value, max_val=(ctrl_max),
                            min_val=(ctrl_min), increment=ctrl_step,
                            digits=ctrl_digits, size=ctrl_size)

    if checkbox:
      ctr_box.Add(self.toggle, flag=wx.ALIGN_CENTER_VERTICAL)
      self.toggle_boxes(flag_on=self.checkbox_state)
      self.Bind(wx.EVT_CHECKBOX, self.onToggle, self.toggle)

    ctr_box.Add(self.txt, flag=wx.ALIGN_CENTER_VERTICAL)
    ctr_box.Add(self.ctr, flag=wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT | wx.EXPAND)


    self.SetSizer(ctr_box)

  def onToggle(self, e):
    self.toggle_boxes(flag_on=self.toggle.GetValue())
    e.Skip()

  def toggle_boxes(self, flag_on=True):
    self.toggle.SetValue(flag_on)
    if flag_on:
      self.ctr.Enable()
      self.ctr.SetValue(int(self.value))
    else:
      self.value = self.ctr.GetValue()
      self.ctr.Disable()

  def reset_default(self):
    self.ctr.SetValue(int(self.value))

class OptionCtrl(CtrlBase):
  ''' Generic panel will place a text control w/ label
      - items: a list of tuples (can be one tuple) with 'key' and 'value'
   '''
  def __init__(self, parent, items,
               label='',
               label_size=(100, -1),
               label_style='normal',
               sub_labels=None,
               sub_label_justify=wx.ALIGN_LEFT,
               sub_label_vertical=wx.ALIGN_CENTER_VERTICAL,
               grid=None,
               checkbox=False,
               checkbox_label='',
               checkbox_state=False,
               ctrl_size=(300, -1),
               expand_rows=None,
               expand_cols=None):

    CtrlBase.__init__(self, parent=parent, label_style=label_style)

    self.items = items
    self.checkbox_state = checkbox_state
    self.toggle = None

    if grid is not None:
      rows = grid[0]
      cols = grid[1]
    else:
      rows = 1
      cols = len(self.items) * 2
    self.ctrl_box = wx.FlexGridSizer(rows, cols, 5, 10)

    for key, value in self.items:
      if sub_labels is None:
        sub_label = key
      else:
        sub_label = sub_labels[self.items.index((key, value))].decode('utf-8')


      if len(self.items) > 1:
        opt_label = wx.StaticText(self, id=wx.ID_ANY, label=sub_label)
        self.ctrl_box.Add(opt_label, flag=sub_label_vertical |
                                          sub_label_justify)

      if key == '':
        item = wx.StaticText(self, label='')
      else:
        item = wx.TextCtrl(self, id=wx.ID_ANY, size=ctrl_size,
                           style=wx.TE_PROCESS_ENTER)
        self.__setattr__(key, item)
        item.SetValue(str(value))
      self.ctrl_box.Add(item, proportion=1,
                        flag=wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT | wx.EXPAND)


    cols = 1
    if label != '':
      cols += 1
      self.txt = wx.StaticText(self, label=label, size=label_size)
      self.txt.SetFont(self.font)

    if checkbox:
      assert checkbox_label != ''
      cols += 1
      self.toggle = wx.CheckBox(self, label=checkbox_label, size=label_size)
      self.toggle.SetValue(self.checkbox_state)

    opt_box = wx.FlexGridSizer(1, cols, 5, 10)

    if label != '':
      opt_box.Add(self.txt, flag=sub_label_vertical)
    if checkbox:
      opt_box.Add(self.toggle, flag=sub_label_vertical)
      self.toggle_boxes(flag_on=False)

      self.Bind(wx.EVT_CHECKBOX, self.onToggle, self.toggle)

    if expand_rows is not None:
      assert isinstance(expand_rows, (list, tuple))
      for row in expand_rows:
        self.ctrl_box.AddGrowableRow(row)
      opt_box.AddGrowableRow(0)

    if expand_cols is not None:
      assert isinstance(expand_cols, (list, tuple))
      for col in expand_cols:
        self.ctrl_box.AddGrowableCol(col)
      opt_box.AddGrowableCol(cols-1)

    opt_box.Add(self.ctrl_box, flag=wx.EXPAND)
    self.SetSizer(opt_box)
    self.Layout()

  def onToggle(self, e):
    self.toggle_boxes(flag_on=self.toggle.GetValue())
    e.Skip()

  def toggle_boxes(self, flag_on=True):
    self.toggle.SetValue(flag_on)

    if self.items is not None:
      if flag_on:
        for item in self.items:
          widget = self.__getattribute__(item[0])
          widget.Enable()
          value = widget.GetValue()
          if noneset(value) is None:
            widget.SetValue(str(item[1]))
          else:
            widget.SetValue(str(value))
      else:
        for item in self.items:
          widget = self.__getattribute__(item[0])
          widget.Disable()
          # widget.Clear()

  def reset_default(self):
    if self.toggle is not None:
      self.toggle_boxes(flag_on=self.checkbox_state)
    else:
      for item in self.items:
        widget = self.__getattribute__(item[0])
        widget.Enable()
        widget.SetValue(item[1])


class TwoButtonCtrl(CtrlBase):
  ''' Generic panel that will place a text control, with a label and an
      optional large button, and an optional bitmap button'''

  def __init__(self, parent,
               label='', label_size=(100, -1),
               label_style='normal',
               text_style=wx.TE_LEFT,
               button1=False,
               button1_label='Browse...',
               button1_size=wx.DefaultSize,
               button2=False,
               button2_label='Default',
               button2_size=wx.DefaultSize,
               value=''):

    CtrlBase.__init__(self, parent=parent, label_style=label_style)

    output_box = wx.FlexGridSizer(1, 5, 0, 5)
    self.txt = wx.StaticText(self, label=label, size=label_size)
    self.txt.SetFont(self.font)
    output_box.Add(self.txt)

    self.ctr = wx.TextCtrl(self, style=text_style)
    self.ctr.SetValue(value)
    output_box.Add(self.ctr, flag=wx.EXPAND)

    if button1:
      self.button1 = wx.Button(self, label=button1_label, size=button1_size)
      output_box.Add(self.button1)

    if button2:
      self.button2 = wx.Button(self, label=button2_label, size=button2_size)
      output_box.Add(self.button2)

    output_box.AddGrowableCol(1, 1)
    self.SetSizer(output_box)

class KnobCtrl(CtrlBase):
  ''' From AGW KnobCtrl class, with attendant spin controls '''

  def __init__(self, parent,
               label='',
               label_size=(100, -1),
               label_style='normal',
               knob_size=(100, 100),
               spin_ctr_size=(60, -1),
               tags=True,
               tags_start=0,
               tags_end=360,
               tags_step=10,
               values_start=0,
               values_end=360,
               values_step=1,
               value=0):
    CtrlBase.__init__(self, parent=parent, label_style=label_style)

    self.knob_ctr = kc.KnobCtrl(self, -1, size=knob_size)
    if tags:
      self.knob_ctr.SetTags(range(tags_start, tags_end, tags_step))
    self.knob_ctr.SetAngularRange(values_start, values_end)
    self.knob_ctr.SetValue(value)
    self.knob_ctr.SetBoundingColour(wx.BLACK)
    self.knob_ctr.SetFirstGradientColour(self.GetBackgroundColour())
    self.knob_ctr.SetSecondGradientColour(self.GetBackgroundColour())

    self.value_ctr = SpinCtrl(self,
                              label=label,
                              label_size=label_size,
                              label_style=label_style,
                              ctrl_size=spin_ctr_size,
                              ctrl_value=str(value),
                              ctrl_max=values_end,
                              ctrl_min=values_start,
                              ctrl_step=values_step)

    ctrl_box = wx.FlexGridSizer(2, 1, 0, 10)

    ctrl_box.Add(self.knob_ctr)
    ctrl_box.Add(self.value_ctr)

    self.SetSizer(ctrl_box)

    self.Bind(kc.EVT_KC_ANGLE_CHANGING, self.onKC_Angle_Change, self.knob_ctr)
    self.Bind(fs.EVT_FLOATSPIN, self.onFS_Value_Change, self.value_ctr.ctr)

  def onKC_Angle_Change(self, e):
    print (self.knob_ctr.GetMaxValue())
    self.value_ctr.ctr.SetValue(self.knob_ctr.GetValue())
    e.Skip()

  def onFS_Value_Change(self, e):
    self.knob_ctr.SetValue(self.value_ctr.ctr.GetValue())


class CustomListCtrl(CtrlBase):
  def __init__(self, parent, size=wx.DefaultSize, content_style='normal'):
    CtrlBase.__init__(self, parent=parent, content_style=content_style,
                      size=size)

    self.sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.sizer)

    # Input List control
    self.control_sizer = wx.BoxSizer(wx.VERTICAL)
    self.ctr = InputListCtrl(self, ID=wx.ID_ANY, size=size,
                             style=ulc.ULC_REPORT |
                                   ulc.ULC_HRULES |
                                   ulc.ULC_VRULES |
                                   ulc.ULC_SINGLE_SEL |
                                   ulc.ULC_HAS_VARIABLE_ROW_HEIGHT |
                                   ulc.ULC_NO_HIGHLIGHT)
    self.control_sizer.Add(self.ctr, -1, flag=wx.EXPAND)
    self.sizer.Add(self.control_sizer, 1, flag=wx.EXPAND)


class CustomImageListCtrl(CtrlBase):
  def __init__(self, parent, size=wx.DefaultSize, content_style='normal'):
    CtrlBase.__init__(self, parent=parent, content_style=content_style)

    self.sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.sizer)

    # Input List control
    self.control_sizer = wx.BoxSizer(wx.VERTICAL)
    self.ctr = InputListCtrl(self, ID=wx.ID_ANY, size=size,
                             style=ulc.ULC_REPORT |
                                   ulc.ULC_HRULES |
                                   ulc.ULC_HAS_VARIABLE_ROW_HEIGHT |
                                   ulc.ULC_VRULES)
    self.control_sizer.Add(self.ctr, -1, flag=wx.EXPAND)
    self.sizer.Add(self.control_sizer, 1, flag=wx.EXPAND)

class VirtualImageListCtrl(CtrlBase):
  def __init__(self, parent, size=wx.DefaultSize, content_style='normal'):
    CtrlBase.__init__(self, parent=parent, content_style=content_style)

    self.sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.sizer)

    self.control_sizer = wx.BoxSizer(wx.VERTICAL)
    self.ctr = VirtualInputListCtrl(self, -1, size=size, n_cols=4,
                                    style=ulc.ULC_REPORT |
                                          ulc.ULC_VIRTUAL|
                                          ulc.ULC_HRULES |
                                          ulc.ULC_VRULES |
                                          ulc.ULC_HAS_VARIABLE_ROW_HEIGHT)

    self.control_sizer.Add(self.ctr, -1, flag=wx.EXPAND)
    self.sizer.Add(self.control_sizer, 1, flag=wx.EXPAND)


class PHILBox(CtrlBase):
  def __init__(self, parent,
               btn_clear_size=(120, -1),
               btn_clear_label='Clear PHIL',
               btn_import=True,
               btn_import_size=(120, -1),
               btn_import_label='Import PHIL',
               btn_export=False,
               btn_export_size=(120, -1),
               btn_export_label='Export PHIL',
               btn_default=True,
               btn_default_size=(120, -1),
               btn_default_label='Default PHIL',
               btn_pos = 'top',
               ctr_size=(-1, 125),
               ctr_value='',
               label_style='normal',
               content_style='normal'):

    CtrlBase.__init__(self, parent=parent, label_style=label_style,
                      content_style=content_style)

    self.main_sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.main_sizer)

    self.ctrl_sizer = wx.GridBagSizer(5, 5)

    self.ctr = wx.richtext.RichTextCtrl(self,
                                        size=ctr_size,
                                        style=wx.VSCROLL,
                                        value=ctr_value)
    span_counter = 1
    if btn_pos in ('left', 'right'):
      self.btn_sizer = wx.BoxSizer(wx.VERTICAL)
      b_flag = wx.BOTTOM
    elif btn_pos in ('top', 'bottom'):
      self.btn_sizer = wx.BoxSizer(wx.HORIZONTAL)
      b_flag = wx.RIGHT
    else:
      b_flag = wx.ALL

    if btn_import:
      self.btn_import = wx.Button(self,
                                  label=btn_import_label,
                                  size=btn_import_size)
      self.btn_sizer.Add(self.btn_import, flag=b_flag, border=5)
    if btn_export:
      self.btn_export = wx.Button(self,
                                  label=btn_export_label,
                                  size=btn_export_size)
      self.btn_sizer.Add(self.btn_export, flag=b_flag, border=5)
    if btn_default:
      self.btn_default = wx.Button(self,
                                   label=btn_default_label,
                                   size=btn_default_size)
      self.btn_sizer.Add(self.btn_default, flag=b_flag, border=5)

    self.btn_clear = wx.Button(self,
                               label=btn_clear_label,
                               size=btn_clear_size)
    self.btn_sizer.Add(self.btn_clear)
    self.Bind(wx.EVT_BUTTON, self.onClear, self.btn_clear)

    if btn_pos == 'left':
      self.ctrl_sizer.Add(self.btn_sizer, pos=(0, 0), flag=wx.ALL)
      self.ctrl_sizer.Add(self.ctr, pos=(0, 1), flag=wx.EXPAND)
      self.ctrl_sizer.AddGrowableCol(1)
      self.ctrl_sizer.AddGrowableRow(0)
    elif btn_pos == 'right':
      self.ctrl_sizer.Add(self.btn_sizer, pos=(0, 1), flag=wx.ALL)
      self.ctrl_sizer.Add(self.ctr, pos=(0, 0), flag=wx.EXPAND)
      self.ctrl_sizer.AddGrowableCol(0)
      self.ctrl_sizer.AddGrowableRow(0)
    elif btn_pos == 'top':
      self.ctrl_sizer.Add(self.btn_sizer, pos=(0, 0), flag=wx.ALL | wx.ALIGN_LEFT)
      self.ctrl_sizer.Add(self.ctr, pos=(1, 0), flag=wx.EXPAND)
      self.ctrl_sizer.AddGrowableCol(0)
      self.ctrl_sizer.AddGrowableRow(1)
    elif btn_pos == 'bottom':
      self.ctrl_sizer.Add(self.btn_sizer, pos=(1, 0), flag=wx.ALL | wx.ALIGN_RIGHT)
      self.ctrl_sizer.Add(self.ctr, pos=(0, 0), flag=wx.EXPAND)
      self.ctrl_sizer.AddGrowableCol(0)
      self.ctrl_sizer.AddGrowableRow(0)

    self.main_sizer.Add(self.ctrl_sizer, 1, flag=wx.EXPAND | wx.ALL, border=5)

  def onClear(self, e):
    self.reset_default()

  def reset_default(self):
    self.ctr.SetValue('')

class TableCtrl(CtrlBase):
  ''' Generic panel will place a table w/ x and y labels
      Data must be a list of lists for multi-column tables '''

  def __init__(self, parent,
               clabels=[],
               clabel_size=wx.DefaultSize,
               rlabels=[],
               rlabel_size=wx.DefaultSize,
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
      for item in clabels:
        clabel = wx.StaticText(self, label=item, size=clabel_size)
        clabel.SetFont(self.font)
        self.sizer.Add(clabel)

    # add row labels and row contents
    for l in rlabels:
      row_label = wx.StaticText(self, label=l, size=rlabel_size)
      row_label.SetFont(self.font)
      self.sizer.Add(row_label)

      # Add data to table
      c_index = rlabels.index(l)
      for item in contents[c_index]:
        cell = wx.StaticText(self, label=item)
        cell.SetFont(self.cfont)
        self.sizer.Add(cell)

    self.SetSizer(self.sizer)

# ------------------------------ Generic Plotter ----------------------------- #

class FastPlotter(CtrlBase):
  '''Panel with one or multiple charts in it with these characteristics:
        1. Chart is configurable (w/ axis sharing)
        2. Chart is fast-updated (i.e. new data are added to existing chart)
        3. Data are selectable (w/ up to three span selectors per chart)

     n_plots: number of plots, can be either int, signifying number of plots, or
              tuple, signifying how the plots are arranged
     axes: list of axis object names, e.g. ['ax1', 'ax2', 'ax3']; must be same
           length as number of plots
     legends: list of tuples; each tuple contains labels for x- and y-axes
     transparent: Boolean, set to True if want figure/plot to be transparent
  '''

  # NOTE: THIS DOES NOT YET WORK!

  # Constructor
  def __init__(self, parent,
               n_plots=1,
               axes=None,
               labels=None,
               chart_label='',
               sharex=False,
               sharey=False,
               transparent=True,
               slider=True):

    from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
    from matplotlib import pyplot as plt
    import numpy as np

    CtrlBase.__init__(self, parent)

    # Initialize the main sizer for the chart
    self.main_box = wx.StaticBox(self, label=chart_label)
    self.main_fig_sizer = wx.StaticBoxSizer(self.main_box, wx.VERTICAL)
    self.SetSizer(self.main_fig_sizer)

    # Create figure and axis array; flatten axis array to zip with other info
    self.track_figure, axarr = plt.subplots(n_plots, sharex=sharex)
    axarr_flat = [ax for sub_ax in axarr for ax in sub_ax]

    # Make sure the number of axes matches number of axis names / legends
    try:
      assert len(axarr_flat) == len(axes) == len(labels)
    except AssertionError:
      print ('FastPlotter Error: Match number of plots to number of ' \
            'axes and axis labels!')

    # Generate generic Nonetype axis object names & labels if none are provided
    if axes is None:
      axes = [None] * len(axes)
    if labels is None:
      labels = [None] * len(axes)
    axes_w_names = zip(axes, axarr_flat, labels)

    # Attach axis objects to class attributes (handling various cases)
    for axn in axes_w_names:
      ax_name = axn[0]
      if ax_name is None:
        ax_name = 'ax' + str(axes_w_names.index(axn))
      ax_obj = axn[1]
      ax_legend = axn[2]
      if transparent:
        ax_obj.patch.set_visible(False)
      if ax_legend is not None:
        if isinstance(ax_legend, (list, tuple)):
          ax_obj.set_xlabel(ax_legend[0])
          ax_obj.set_ylabel(ax_legend[1])
        else:
          ax_obj.set_ylabel(ax_legend)
      ax_obj.set_autoscaley_on(True)
      setattr(self, ax_name, ax_obj)

    # Data-handling attributes
    self.xdata = []
    self.ydata = []
    self.x_min = 0
    self.x_max = 1
    self.y_max = 1

    # Chart navigation attributes
    self.bracket_set = False
    self.button_hold = False
    self.plot_zoom = False
    self.chart_range = None
    self.selector = None
    self.max_lock = True
    self.patch_x = 0
    self.patch_x_last = 1
    self.patch_width = 1
    self.start_edge = 0
    self.end_edge = 1

    # Set general figure characteristics (default tight layout)
    if transparent:
      self.track_figure.patch.set_visible(False)
    self.track_figure.set_tight_layout(True)
    self.track_canvas = FigureCanvas(self, -1, self.track_figure)

    # Slider bar for the plot
    if slider:
      self.plot_sb = wx.Slider(self, minValue=0, maxValue=1)
      self.plot_sb.Hide()
      self.Bind(wx.EVT_SCROLL, self.onScroll, self.plot_sb)

    self.main_fig_sizer.Add(self.track_canvas, 1, wx.EXPAND)
    self.main_fig_sizer.Add(self.plot_sb, flag=wx.EXPAND)

    # Plot bindings
    self.track_figure.canvas.mpl_connect('button_press_event', self.onPress)

  def onScroll(self, e):
    sb_center = self.plot_sb.GetValue()
    half_span = (self.x_max - self.x_min) / 2
    if sb_center - half_span == 0:
      self.x_min = 0
      self.x_max = half_span * 2
    else:
      self.x_min = sb_center - half_span
      self.x_max = sb_center + half_span

    if self.plot_sb.GetValue() == self.plot_sb.GetMax():
      self.max_lock = True
    else:
      self.max_lock = False

    self.draw_plot()


  def draw_plot(self, new_x=None, new_y=None, new_data=False):

    if new_x is None:
      new_x = []
    if new_y is None:
      new_y = []

    nref_x = np.append(self.xdata, np.array(new_x).astype(np.double))
    nref_y = np.append(self.ydata, np.array(new_y).astype(np.double))
    self.xdata = nref_x
    self.ydata = nref_y
    nref_xy = zip(nref_x, nref_y)
    all_acc = [i[0] for i in nref_xy if i[1] >= min_bragg]
    all_rej = [i[0] for i in nref_xy if i[1] < min_bragg]

    if nref_x != [] and nref_y != []:
      if self.plot_zoom:
        if self.max_lock:
          self.x_max = np.max(nref_x)
          self.x_min = self.x_max - self.chart_range
      else:
        self.x_min = -1
        self.x_max = np.max(nref_x) + 1

      if min_bragg > np.max(nref_y):
        self.y_max = min_bragg + int(0.1 * min_bragg)
      else:
        self.y_max = np.max(nref_y) + int(0.1 * np.max(nref_y))

      self.track_axes.set_xlim(self.x_min, self.x_max)
      self.track_axes.set_ylim(0, self.y_max)

    else:
      self.x_min = -1
      self.x_max = 1

    acc = [i for i in all_acc if i > self.x_min and i < self.x_max]
    rej = [i for i in all_rej if i > self.x_min and i < self.x_max]

    self.acc_plot.set_xdata(nref_x)
    self.rej_plot.set_xdata(nref_x)
    self.acc_plot.set_ydata(nref_y)
    self.rej_plot.set_ydata(nref_y)
    self.acc_plot.set_markevery(acc)
    self.rej_plot.set_markevery(rej)

    self.Layout()

    count = '{}'.format(len([i for i in nref_xy if i[1] >= min_bragg]))
    self.main_window.tracker_panel.count_txt.SetLabel(count)
    self.main_window.tracker_panel.status_sizer.Layout()

    # Set up scroll bar
    if len(self.xdata) > 0:
      self.plot_sb.SetMax(np.max(nref_x))
      if self.max_lock:
        self.plot_sb.SetValue(self.plot_sb.GetMax())

    # Draw extended plots
    self.track_axes.draw_artist(self.acc_plot)
    self.track_axes.draw_artist(self.rej_plot)
