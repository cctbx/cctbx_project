from __future__ import division

from builtins import str
from builtins import object
'''
Author      : Lyubimov, A.Y.
Created     : 07/08/2016
Last Changed: 07/21/2017
Description : IOTA GUI controls
'''
from __future__ import print_function

import os
import wx
import wx.richtext
import wx.lib.agw.floatspin as fs
import wx.lib.agw.ultimatelistctrl as ulc
from wx.lib.mixins.listctrl import ListCtrlAutoWidthMixin, ColumnSorterMixin
from wxtbx import metallicbutton as mb
from wxtbx import bitmaps
import wx.lib.buttons as btn

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
    self.btn_mag = btn.GenBitmapButton(self, bitmap=viewmag_bmp)
    del_bmp = bitmaps.fetch_icon_bitmap('actions', 'editdelete', size=16)
    self.btn_delete = btn.GenBitmapButton(self, bitmap=del_bmp)
    info_bmp = bitmaps.fetch_icon_bitmap('actions', 'info', size=16)
    self.btn_info = btn.GenBitmapButton(self, bitmap=info_bmp)
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
    ulc.UltimateListCtrl.__init__(self, parent, ID, pos, size, agwStyle=style)
    ListCtrlAutoWidthMixin.__init__(self)

class VirtualInputListCtrl(ulc.UltimateListCtrl, ListCtrlAutoWidthMixin,
                           ColumnSorterMixin):
  ''' Customized Virtual UltimateListCtrl with auto-width mixin'''

  def __init__(self, parent, ID, n_cols=3, pos=wx.DefaultPosition,
               size=wx.DefaultSize, style=0):
    ulc.UltimateListCtrl.__init__(self, parent, ID, pos, size, agwStyle=style)
    ListCtrlAutoWidthMixin.__init__(self)
    ColumnSorterMixin.__init__(self, n_cols)

    self.data = {}

  def InitializeDataMap(self, data):
    self.data = data
    self.itemDataMap = self.data
    self.itemIndexMap = list(self.data.keys())
    self.SetItemCount(len(self.data))

  def GetListCtrl(self):
    return self

  def OnColClick(self, e):
    print("column clicked")
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
                 'raw image':'This is a single raw diffraction image file',
                 'image pickle':'This is a single diffraction image file '
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
                             'pseudo-Wilson scaling'
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
               label_size=(200, -1),
               label_style='normal',
               checkbox=False,
               checkbox_state=False,
               checkbox_label='',
               ctrl_size=(60, -1),
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

    ctr_box = wx.FlexGridSizer(1, cols, 0, 10)

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

  def toggle_boxes(self, flag_on):
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
               sub_labels=[],
               sub_label_justify=wx.ALIGN_LEFT,
               grid=None,
               checkbox=False,
               checkbox_label='',
               checkbox_state=False,
               ctrl_size=(300, -1)):

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
      if sub_labels != []:
        sub_label = sub_labels[self.items.index((key, value))].decode('utf-8')
      else:
        sub_label = key

      if len(self.items) > 1:
        opt_label = wx.StaticText(self, id=wx.ID_ANY, label=sub_label)
        self.ctrl_box.Add(opt_label, flag=wx.ALIGN_CENTER_VERTICAL | sub_label_justify)

      if key == '':
        item = wx.StaticText(self, label='')
      else:
        item = wx.TextCtrl(self, id=wx.ID_ANY, size=ctrl_size,
                           style=wx.TE_PROCESS_ENTER)
        self.__setattr__(key, item)
        item.SetValue(str(value))
      self.ctrl_box.Add(item, flag=wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT | wx.EXPAND)


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
      opt_box.Add(self.txt, flag=wx.ALIGN_CENTER_VERTICAL)
    if checkbox:
      opt_box.Add(self.toggle, flag=wx.ALIGN_CENTER_VERTICAL)
      self.toggle_boxes(flag_on=False)

      self.Bind(wx.EVT_CHECKBOX, self.onToggle, self.toggle)

    opt_box.Add(self.ctrl_box)
    self.SetSizer(opt_box)

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
          widget.SetValue(item[1])
      else:
        for item in self.items:
          widget = self.__getattribute__(item[0])
          widget.Disable()
          widget.Clear()

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

    output_box = wx.FlexGridSizer(1, 5, 0, 10)
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

class CustomListCtrl(CtrlBase):
  def __init__(self, parent, size=wx.DefaultSize, content_style='normal'):
    CtrlBase.__init__(self, parent=parent, content_style=content_style)

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
    self.ctr = VirtualInputListCtrl(self, -1, size=size, n_cols=3,
                                    style=ulc.ULC_REPORT |
                                          ulc.ULC_VIRTUAL|
                                          ulc.ULC_HRULES |
                                          ulc.ULC_VRULES |
                                          ulc.ULC_HAS_VARIABLE_ROW_HEIGHT)

    self.control_sizer.Add(self.ctr, -1, flag=wx.EXPAND)
    self.sizer.Add(self.control_sizer, 1, flag=wx.EXPAND)


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
      self.sizer.Add(self.ctr, pos=(0, 1), span=(span_counter + 1, 1),
                     flag=wx.EXPAND)
      self.sizer.AddGrowableCol(1)
    elif span_counter == 0:
      self.sizer.Add(self.ctr, pos=(0, 0), flag=wx.EXPAND)
      self.sizer.AddGrowableCol(0)
    self.sizer.AddGrowableRow(span_counter)

  def reset_default(self):
    self.ctr.SetValue('')

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
