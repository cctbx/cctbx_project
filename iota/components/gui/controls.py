from __future__ import absolute_import, division, print_function

'''
Author      : Lyubimov, A.Y.
Created     : 07/08/2016
Last Changed: 12/02/2019
Description : IOTA GUI controls
'''

import os
from glob import glob

import wx
import wx.richtext
import wx.lib.agw.floatspin as fs
import wx.lib.agw.ultimatelistctrl as ulc
import wx.lib.agw.knobctrl as kc
from wx.lib.mixins.listctrl import ListCtrlAutoWidthMixin, ColumnSorterMixin
import wx.lib.buttons as btn

from wxtbx import metallicbutton as mb
from wxtbx import bitmaps
from libtbx.utils import Sorry

from iota.components.iota_utils import noneset, InputFinder
from iota.components.iota_threads import ImageViewerThread

ginp = InputFinder()

# Platform-specific stuff
# TODO: Will need to test this on Windows at some point
if wx.Platform == '__WXGTK__':
  plot_font_size = 9
  norm_font_size = 9
  button_font_size = 10
  LABEL_SIZE = 12
  CAPTION_SIZE = 10
  python = 'python'
elif wx.Platform == '__WXMAC__':
  plot_font_size = 9
  norm_font_size = 12
  button_font_size = 14
  LABEL_SIZE = 14
  CAPTION_SIZE = 12
  python = "Python"
elif wx.Platform == '__WXMSW__':
  plot_font_size = 9
  norm_font_size = 9
  button_font_size = 11
  LABEL_SIZE = 11
  CAPTION_SIZE = 9
  python = "Python"  # TODO: make sure it's right!

# Metallicbutton globals (temporary!)
GRADIENT_NORMAL = 0
GRADIENT_PRESSED = 1
GRADIENT_HIGHLIGHT = 2

MB_STYLE_DEFAULT = 1
MB_STYLE_BOLD_LABEL = 2
MB_STYLE_DROPARROW = 4

# --------------------------- Master Widget Classes -------------------------- #

class IOTACtrl(object):
  ''' Mixin to add IOTA-specific functions to IOTA controls '''

  def bind_event(self, event_type, event_handler):
    ''' General function to bind an event to parent frame; absent that,
        to parent object; absent that, to current object
    :param event_type: e.g. wx.EVT_BUTTON
    :param event_handler: handler function for bound event
    '''

    try:
      self.parent.get_frame().Bind(event_type, event_handler, self)
    except Exception:
      try:
        self.parent.Bind(event_type, event_handler, self)
      except Exception:
        self.Bind(event_type, event_handler, self)

# --------------------------------- Widgets ---------------------------------- #

def StandardBitmap(img_name, size=None):
  img_path = img_name
  img = wx.Image(img_path, type=wx.BITMAP_TYPE_ANY, index=-1)
  if size is not None:
    (w, h) = size
    img.Rescale(w, h)
  bmp = img.ConvertToBitmap()
  return bmp

class GradButton(mb.MetallicButton, IOTACtrl):
  """ Customized MetallicButton """

  def __init__(self,
               parent,
               label='',
               bmp=None,
               size=wx.DefaultSize,
               style=mb.MB_STYLE_BOLD_LABEL,
               handler_function=None,
               user_data=None,
               start_color=(225, 225, 225),
               gradient_percent=-10,
               highlight_color=(250, 250, 250),
               label_size=LABEL_SIZE,
               caption_size=CAPTION_SIZE,
               button_margin=4,
               disable_after_click=0):
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

    self.user_data = user_data
    if handler_function is not None:
      self.bind_event(wx.EVT_BUTTON, handler_function)


class AddDeleteButtonBox(wx.Panel):
  """ A box with Add and Delete bitmap buttons; when Delete button is
      pressed, several delete options show up; an optional Undo button can be
      added """
  def __init__(self, parent, reset_button=False):
    wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)

    self.btn_sizer = wx.FlexGridSizer(1, 4, 0, 5)
    self.SetSizer(self.btn_sizer)

    bmp_add = bitmaps.fetch_icon_bitmap('actions', 'edit_add', size=16)
    bmp_del = bitmaps.fetch_icon_bitmap('actions', 'editdelete', size=16)
    bmp_rst = bitmaps.fetch_icon_bitmap('actions', 'recur', size=16)

    # Main buttons
    self.btn_add = GradButton(self, bmp=bmp_add, button_margin=2, size=(24, 24))
    self.btn_del = GradButton(self, bmp=bmp_del, button_margin=2, size=(24, 24))
    self.btn_rst = GradButton(self, bmp=bmp_rst, button_margin=2, size=(24, 24))
    self.btn_sizer.AddMany([self.btn_add,
                            self.btn_del,
                            self.btn_rst])
    if not reset_button:
      self.btn_rst.Hide()

    # Delete option buttons
    self.del_btn_sizer = wx.FlexGridSizer(1, 3, 0, 5)
    self.btn_del_lst = GradButton(self, label='Delete Last', label_size=10,
                                  style=mb.MB_STYLE_DEFAULT)
    self.btn_del_sel = GradButton(self, label='Delete Selected', label_size=10,
                                  style=mb.MB_STYLE_DEFAULT)
    self.btn_del_not = GradButton(self, label='Cancel', label_size=10,
                                  style=mb.MB_STYLE_DEFAULT)
    self.del_btn_sizer.AddMany([self.btn_del_sel,
                               self.btn_del_lst,
                               self.btn_del_not])
    self.btn_sizer.Add(self.del_btn_sizer, flag=wx.LEFT, border=10)
    self.Layout()

    self.btn_sizer.Hide(self.del_btn_sizer)

    # Bindings
    # Main Buttons
    self.Bind(wx.EVT_BUTTON, self.onDelete, self.btn_del)
    self.Bind(wx.EVT_BUTTON, self.onAdd, self.btn_add)
    self.Bind(wx.EVT_BUTTON, self.onReset, self.btn_rst)

    # Delete Buttons
    self.Bind(wx.EVT_BUTTON, self.onDeleteSelected, self.btn_del_sel)
    self.Bind(wx.EVT_BUTTON, self.onDeleteLast, self.btn_del_lst)
    self.Bind(wx.EVT_BUTTON, self.onCancelDelete, self.btn_del_not)

  def onDelete(self, e):
    self._show_delete_options()
    e.Skip()

  def onAdd(self, e):
    e.Skip()

  def onReset(self, e):
    e.Skip()

  def onDeleteSelected(self, e):
    self._hide_delete_options()
    e.Skip()

  def onDeleteLast(self, e):
    self._hide_delete_options()
    e.Skip()

  def onCancelDelete(self, e):
    self._hide_delete_options()
    e.Skip()

  def _show_delete_options(self):
    self.btn_del.Disable()
    self.btn_add.Disable()
    self.btn_rst.Disable()
    self.btn_sizer.Show(self.del_btn_sizer)
    self.Layout()

  def _hide_delete_options(self):
    self.btn_del.Enable()
    self.btn_add.Enable()
    self.btn_rst.Enable()
    self.btn_sizer.Hide(self.del_btn_sizer)
    self.Layout()



class MiniButtonBox(wx.Panel):
  """ A box with three mini buttons for IOTA panel """

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
  """ A box with three mini buttons for Input file panel """

  def __init__(self, parent):
    wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)

    self.index = None
    self.btn_box = wx.BoxSizer(wx.HORIZONTAL)

    viewmag_bmp = bitmaps.fetch_icon_bitmap('actions', 'viewmag', size=16)
    self.btn_mag = GradButton(self, bmp=viewmag_bmp, size=(24, 24),
                              button_margin=2, gradient_percent=-15)

    del_bmp = bitmaps.fetch_icon_bitmap('actions', 'editdelete', size=16)
    self.btn_delete = GradButton(self, bmp=del_bmp, size=(24, 24),
                                 button_margin=2, gradient_percent=-15)

    info_bmp = bitmaps.fetch_icon_bitmap('actions', 'info', size=16)
    self.btn_info = GradButton(self, bmp=info_bmp, size=(24, 24),
                               button_margin=2, gradient_percent=-15)

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


class IOTAButton(wx.Button, IOTACtrl):
  _id = -1
  ''' Master class for all basic IOTA buttons '''

  def __init__(self, parent, label=wx.EmptyString, handler_function=None):
    ''' Constructor
    :param label: button label
    :param handler_function: handler function for the button press event
    '''
    wx.Button.__init__(self, parent, label=label, id=self._id)

    self.window = getattr(parent, 'window', None)
    if not self.window:
      self.window = self.GetTopLevelParent()
    self.parent = parent

    if handler_function:
      self.bind_event(wx.EVT_BUTTON, handler_function)


class OKButton(IOTAButton):
  _id = wx.ID_OK


class CancelButton(IOTAButton):
  _id = wx.ID_CANCEL


# ---------------------------------- Inputs ---------------------------------- #

class InputListCtrl(ulc.UltimateListCtrl, ListCtrlAutoWidthMixin):
  """ Customized UltimateListCtrl with auto-width mixin"""

  def __init__(self, parent, ID, pos=wx.DefaultPosition, style=0,
               *args, **kwargs):
    ulc.UltimateListCtrl.__init__(self, parent, ID, pos, agwStyle=style,
                                  *args, **kwargs)
    ListCtrlAutoWidthMixin.__init__(self)


class VirtualInputListCtrl(ulc.UltimateListCtrl, ListCtrlAutoWidthMixin,
                           ColumnSorterMixin):
  """ Customized Virtual UltimateListCtrl with auto-width mixin"""

  def __init__(self, parent, ID, n_cols=4, pos=wx.DefaultPosition,
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
  """ Class that will contain all the elements of an input list entry """
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
  """ Class that will contain all the elements of a file list entry """
  def __init__(self, path, items=None):
    """ A generic control to generate an entry in the file list control
    :param path: absolute path to file
    :param items: a dictionary of items that will be converted to attributes
    """
    self.id = None
    self.flag = None
    self.path = path
    self.warning = False

    if items is not None:
      for key, value in items.items():
        self.__setattr__(key, value)

# --------------------------------- Controls --------------------------------- #

class CtrlBase(wx.Panel):
  """ Control panel base class
     @DynamicAttrs
  """
  def __init__(self,
               parent,
               label_style='normal',
               label_font_size=norm_font_size,
               content_style='normal',
               content_font_size=norm_font_size,
               size=wx.DefaultSize):

    wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY, size=size)

    self.window = self.GetTopLevelParent()

    # TODO: streamline this
    # Set control attributes
    self.expert_level = 0
    self.font = wx.Font(norm_font_size, wx.FONTFAMILY_DEFAULT,
                        wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL)
    self.cfont = wx.Font(norm_font_size, wx.FONTFAMILY_DEFAULT,
                         wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL)

    # Set font attributes for label
    if 'bold' in label_style:
      self.font.SetWeight(wx.FONTWEIGHT_BOLD)
    if 'italic' in label_style:
      self.font.SetStyle(wx.FONTSTYLE_ITALIC)
    if 'teletype' in label_style:
      self.font.SetFamily(wx.FONTFAMILY_TELETYPE)

    # Set font attributes for content
    if 'bold' in content_style:
      self.cfont.SetWeight(wx.FONTWEIGHT_BOLD)
    if 'italic' in content_style:
      self.cfont.SetStyle(wx.FONTSTYLE_ITALIC)
    if 'teletype' in content_style:
      self.cfont.SetFamily(wx.FONTFAMILY_TELETYPE)

    self.font.SetPointSize(label_font_size)
    self.cfont.SetPointSize(content_font_size)


class DialogButtonsCtrl(CtrlBase):
  """ Customizable "bottom of window" set of buttons for dialogs """

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
    elif preset == 'PROC_DIALOG':
      buttons = [('OK', wx.ID_OK), ('Cancel', wx.ID_CANCEL)]
      choices = ['Basic', 'Advanced', 'Developer']
      choice_label = 'Expert Level: '
      choice_size = (150, -1)
    elif preset == 'PHIL_DIALOG':
      buttons = [('OK', wx.ID_OK), ('Cancel', wx.ID_CANCEL)]
      choices = ['Basic', 'Intermediate', 'Advanced', 'Developer']
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
  """ Text control with multiple buttons """
  def __init__(self, parent,
               buttons = None,
               button_size = wx.DefaultSize,
               button_style = wx.BU_AUTODRAW,
               ctrl_value = '',
               ctrl_size = wx.DefaultSize,
               ctrl_label = '',
               ctrl_label_size=wx.DefaultSize):
    CtrlBase.__init__(self, parent=parent)

    main_sizer = wx.FlexGridSizer(1, 3, 0, 15)

    self.txt_label = wx.StaticText(self, label=ctrl_label, size=ctrl_label_size)
    self.txt_ctrl = wx.TextCtrl(self, size=ctrl_size,
                                style=wx.TE_PROCESS_ENTER | wx.TE_DONTWRAP)
    self.txt_ctrl.SetValue(ctrl_value)
    main_sizer.Add(self.txt_label)
    main_sizer.Add(self.txt_ctrl, flag=wx.EXPAND)

    if buttons is not None:
      btn_sizer = wx.FlexGridSizer(1, len(buttons), 0, 0)
      for key, id in buttons.items():
        btn_name = 'btn_{}'.format(str.lower(key))
        if type(id[1]) == str:
          button = wx.Button(self, label=id[1], id=id[0], size=button_size,
                             style=button_style)
        else:
          button = GradButton(self, bmp=id[1], start_color=(250, 250, 250),
                              size=(28, 28), gradient_percent=0)
        self.__setattr__(btn_name, button)
        btn_sizer.Add(button, flag=wx.RIGHT, border=5)
      main_sizer.Add(btn_sizer)

    main_sizer.AddGrowableCol(1)
    self.SetSizer(main_sizer)

class InputCtrl(CtrlBase):
  """ Generic panel that will place a text control with a label """

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
  """ Generic panel will place a choice control w/ label and a text control (
  optional) which will activate when a 'custom' setting is selected """

  def __init__(self, parent,
               choices,
               custom_choices=None,
               captions=None,
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
    if custom_choices:
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
      for c in self.choices:
        item_idx = self.ctr.FindString(c[0])
        self.ctr.SetClientData(item_idx, c[1])
    else:
      self.ctr = wx.Choice(self, size=ctrl_size, choices=self.choices)

    ctr_box.Add(self.txt, flag=wx.ALIGN_CENTER_VERTICAL)
    ctr_box.Add(self.ctr, flag=wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT)

    # Create binding for choice in case custom choice(s) specified
    if custom_choices:
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
    if self.custom_choices:
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
  """ Generic panel will place a spin control w/ label """
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
  """ Generic panel will place a text control w/ label
      - items: a list of tuples (can be one tuple) with 'key' and 'value'
   """
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
      assert isinstance(expand_cols, (list, tuple, int))
      if type(expand_cols) == int:
        self.ctrl_box.AddGrowableCol(expand_cols)
      else:
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
  """ Generic panel that will place a text control, with a label and an
      optional large button, and an optional bitmap button"""

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
  """ From AGW KnobCtrl class, with attendant spin controls """

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
      self.knob_ctr.SetTags(list(range(tags_start, tags_end, tags_step)))
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
  def __init__(self, parent, *args, **kwargs):
    custom_style = kwargs.pop('style', None)
    CtrlBase.__init__(self, parent=parent, *args, **kwargs)

    self.sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.sizer)

    style = ulc.ULC_REPORT | \
            ulc.ULC_HRULES | \
            ulc.ULC_VRULES | \
            ulc.ULC_SINGLE_SEL | \
            ulc.ULC_HAS_VARIABLE_ROW_HEIGHT \
            #| ulc.ULC_NO_HIGHLIGHT
    if custom_style:
      style |= custom_style

    # Input List control
    self.ctr = InputListCtrl(self, ID=wx.ID_ANY, style=style)
    self.sizer.Add(self.ctr, 1, flag=wx.EXPAND)
    self.ctr.SetFont(self.cfont)


class CustomImageListCtrl(CtrlBase):
  def __init__(self, parent, size=wx.DefaultSize, content_style='normal',
               *args, **kwargs):
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
    self.ctr = VirtualInputListCtrl(self, -1, size=size, style=ulc.ULC_REPORT |
                                                               ulc.ULC_VIRTUAL |
                                                               ulc.ULC_HRULES |
                                                               ulc.ULC_VRULES |
                                                               ulc.ULC_HAS_VARIABLE_ROW_HEIGHT)

    self.control_sizer.Add(self.ctr, -1, flag=wx.EXPAND)
    self.sizer.Add(self.control_sizer, 1, flag=wx.EXPAND)


class FileListCtrl(CustomListCtrl):
  """ File list window for the input tab """

  _file_types = [
    'text file',
    'binary file'
  ]

  _folder_types = [
    'text folder',
    'binary folder'
  ]
  _data_types = ['file', 'folder']

  _input_filter = None

  def __init__(self, parent, size=(-1, 400), file_types=None,
               folder_types=None, data_types=None, input_filter=None,
               *args, **kwargs):
    CustomListCtrl.__init__(self, parent=parent, size=size,
                            style=ulc.ULC_STICKY_HIGHLIGHT |
                                  ulc.ULC_STICKY_NOSELEVENT |
                                  ulc.ULC_BORDER_SELECT,
                            # style=ulc.ULC_HOT_TRACKING,
                            *args, **kwargs)

    self.parent = parent
    self.window = parent.window

    # Initialize dictionaries for imported data types
    self.all_data_images = {}
    self.all_img_objects = {}
    self.all_proc_pickles = {}
    self.image_count = 0

    # Add custom parameters
    if file_types:
      self._file_types = file_types
    if folder_types:
      self._folder_types = folder_types
    if data_types:
      self._data_types = data_types
    if input_filter:
      self._input_filter = input_filter

    # Generate columns
    self.ctr.InsertColumn(0, "")
    self.ctr.InsertColumn(1, "Input Path")
    self.ctr.InsertColumn(2, "Input Type")
    self.ctr.InsertColumn(3, "Action")
    self.ctr.setResizeColumn(2)

    # Add file / folder buttons
    self.button_sizer = wx.FlexGridSizer(1, 5, 0, 5)
    bmp_add = bitmaps.fetch_icon_bitmap('actions', 'edit_add', size=16)
    bmp_browse = bitmaps.fetch_icon_bitmap('actions', 'open', scale=(16, 16))
    self.inp_path = TextCtrlWithButtons(self, ctrl_label='Input path: ')
    self.btn_add_path = GradButton(self, bmp=bmp_add, label=' Add Path',
                                   label_size=norm_font_size)
    self.btn_add_path.Disable()
    self.btn_browse = GradButton(self, label=' Browse...',
                                 label_size=norm_font_size)
    self.txt_total_images = wx.StaticText(self, label='0 total images')
    self.button_sizer.Add(self.inp_path, flag=wx.EXPAND)
    self.button_sizer.Add(self.btn_add_path)
    self.button_sizer.Add(self.btn_browse)
    self.button_sizer.Add((20, 0))
    self.button_sizer.Add(self.txt_total_images, flag=wx.ALIGN_RIGHT)
    self.button_sizer.AddGrowableCol(0)

    self.sizer.Add(self.button_sizer, flag=wx.EXPAND | wx.TOP, border=10)

    # Event bindings
    self.Bind(wx.EVT_BUTTON, self.onAddPath, self.btn_add_path)
    self.Bind(wx.EVT_BUTTON, self.onBrowse, self.btn_browse)
    self.Bind(wx.EVT_TEXT_ENTER, self.onPath, self.inp_path.txt_ctrl)
    self.Bind(wx.EVT_TEXT, self.onPathTyped, self.inp_path.txt_ctrl)

  def onPath(self, e):
    self.add_path()
    self.inp_path.txt_ctrl.SetValue('')
    self.btn_add_path.Disable()

  def onPathTyped(self, e):
    txt = self.inp_path.txt_ctrl.GetValue()
    if txt and not txt.isspace():
      self.btn_add_path.Enable()
    else:
      self.btn_add_path.Disable()

  def onAddPath(self, e):
    self.add_path()
    self.btn_add_path.Disable()

  def add_path(self):
    path = self.inp_path.txt_ctrl.GetValue()
    self.inp_path.txt_ctrl.SetValue(os.path.abspath(path))

    # Unpack wildcards:
    paths = glob(pathname=path)

    # Add paths to input control
    for p in paths:
      if os.path.exists(p):
        self.add_item(os.path.abspath(p))

  def onBrowse(self, e):
    obj = e.GetEventObject()
    command_list = [('Browse folders...',
                     lambda evt: self.open_folder_dialog()),
                    ('Browse files...',
                     lambda evt: self.open_file_dialog())]
    browse_menu = Menu(self)
    browse_menu.add_commands(command_list)
    self.PopupMenu(browse_menu)
    browse_menu.Destroy()

  def open_file_dialog(self):
    wx.SystemOptions.SetOption("osx.openfiledialog.always-show-types", "1")
    file_dlg = wx.FileDialog(self,
                             message="Load File",
                             defaultDir=os.curdir,
                             defaultFile='*',
                             wildcard="All Files (*.*)|*.*|",
                             style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST |
                                   wx.FD_MULTIPLE)
    if file_dlg.ShowModal() == wx.ID_OK:
      files = file_dlg.GetPaths()
      for item in files:
        self.add_item(item)
    file_dlg.Destroy()

  def open_folder_dialog(self):
    dlg = wx.DirDialog(self, "Load Folder:",
                       style=wx.DD_DEFAULT_STYLE)
    if dlg.ShowModal() == wx.ID_OK:
      path = dlg.GetPath()
      self.add_item(path)
    dlg.Destroy()

  def set_type_choices(self, path):
    # Determine what type of input this is and present user with choices
    type_choices = ['[  SELECT INPUT TYPE  ]']
    preferred_selection = 0
    inputs, input_type, input_count = ginp.get_input(
      path,
      filter_results=(self._input_filter is not None),
      filter_type=self._input_filter)
    if os.path.isdir(path):
      type_choices.extend(self._folder_types)
      if input_type in type_choices:
        preferred_selection = type_choices.index(input_type)
    elif os.path.isfile(path):
      type_choices.extend(self._file_types)
      if input_type in type_choices:
        preferred_selection = type_choices.index(input_type)
    return inputs, input_count, type_choices, preferred_selection

  def add_item(self, path):
    # Generate item
    inputs, input_count, inp_choices, inp_sel = self.set_type_choices(path)
    type_choice = DataTypeChoice(self.ctr,
                                 choices=inp_choices)
    item = InputListItem(path=path,
                         type=type_choice,
                         buttons=MiniButtonBoxInput(self.ctr))

    self.Bind(wx.EVT_CHOICE, self.onTypeChoice, item.type.type)
    self.Bind(wx.EVT_BUTTON, self.onMagButton, item.buttons.btn_mag)
    self.Bind(wx.EVT_BUTTON, self.onDelButton, item.buttons.btn_delete)
    self.Bind(wx.EVT_BUTTON, self.onInfoButton, item.buttons.btn_info)

    # Insert list item
    idx = self.ctr.InsertStringItem(self.ctr.GetItemCount() + 1, '')
    self.ctr.SetStringItem(idx, 1, item.path)
    self.ctr.SetItemWindow(idx, 2, item.type, expand=True)
    self.ctr.SetItemWindow(idx, 3, item.buttons, expand=True)

    # Set drop-down selection, check it for data and open other tabs
    item.type.type.SetSelection(inp_sel)
    sel = item.type.type.GetString(inp_sel)
    have_data = True in [(t in sel) for t in self._data_types]
    if have_data:
      # Main Window may not have the run button yet
      if hasattr(self.window, 'btn_run'):
        self.window.btn_run.Enable()
      self.all_data_images[item.path] = inputs

      # Display input type
      input_type = item.type.type.GetString(inp_sel)

      # Display image count
      self.ctr.SetStringItem(idx, 0, str(input_count))

      # Place a "view" button
      if "image" in input_type:
        view_bmp = bitmaps.fetch_custom_icon_bitmap('hklview_2d',
                                                    scale=(16, 16))
        item.buttons.btn_mag.SetBitmapLabel(view_bmp)
    else:
      n_images = 0
      warn_bmp = bitmaps.fetch_icon_bitmap('actions', 'status_unknown',
                                           size=16)
      item.buttons.btn_info.SetBitmapLabel(warn_bmp)
      item.warning = True

    # Record index in all relevant places
    item.id = idx
    item.buttons.index = idx
    item.type.index = idx
    item.type_selection = inp_sel

    # Resize columns to fit content
    self.ctr.SetColumnWidth(0, width=-1)
    self.ctr.SetColumnWidth(2, width=-1)
    self.ctr.SetColumnWidth(3, width=-1)
    self.ctr.SetColumnWidth(1, width=-3)

    # Make sure all the choice lists are the same size
    if item.type.type.GetSize()[0] < self.ctr.GetColumnWidth(2) - 5:
       item.type.type.SetSize((self.ctr.GetColumnWidth(2) - 5, -1))

    # Attach data object to item
    self.ctr.SetItemData(item.id, item)

    if input_count > 0:
      self.image_count += input_count
      if self.image_count > 0:
        self.update_total_image_count()

    self.window.Layout()

  def onTypeChoice(self, e):
    type = e.GetEventObject().GetParent()
    item_data = self.ctr.GetItemData(type.index)
    item_data.type.type.SetSelection(type.type.GetSelection())
    item_data.type_selection = type.type.GetSelection()

    # Evaluate whether data folders / files are present
    data_items = 0
    for idx in range(self.ctr.GetItemCount()):
      if self.ctr.GetItemData(idx).type_selection != 0:
        data_items += 1
    self.window.set_tool_state(self.window.btn_run, (data_items > 0))
    e.Skip()

  def onMagButton(self, e):
    idx = e.GetEventObject().GetParent().index
    item_obj = self.ctr.GetItemData(idx)
    path = item_obj.path
    type = item_obj.type.type.GetString(item_obj.type_selection)

    if os.path.isfile(path):
      if 'image' in type:
        if 'file' in type:
          self.view_images([path], img_type=type)
        elif 'list' in type:
          with open(path, 'r') as f:
            file_list = [i.replace('\n', '') for i in f.readlines()]
            self.view_images(file_list, img_type=type)
      elif type == 'text' or 'list' in type:
        with open(path, 'r') as f:
          file_list = f.readlines()
          msg = ' '.join(file_list)
          textview = TextFileView(self, title=path, contents=msg)
          textview.ShowModal()
      else:
        wx.MessageBox('Unknown file format', 'Warning',
                      wx.OK | wx.ICON_EXCLAMATION)
    elif os.path.isdir(path):
      file_list, input_type, input_count = ginp.get_input(path)
      if 'image' in input_type:
        self.view_images(file_list, img_type=type)
      else:
        msg = ' '.join(file_list)
        textview = TextFileView(self, title=path, contents=msg)
        textview.ShowModal()

  def view_images(self, img_list, img_type=None):
    """ Launches image viewer (depending on backend) """
    # self.parent.input_phil.show()
    viewer = self.window.gparams.gui.image_viewer
    if viewer == 'cxi.view' and 'pickle' not in img_type:
        wx.MessageBox('cxi.view only accepts image pickles', 'Warning',
                      wx.OK | wx.ICON_EXCLAMATION)
    else:
      if len(img_list) > 10:
        view_warning = ViewerWarning(self, len(img_list))
        if view_warning.ShowModal() == wx.ID_OK:
          # parse 'other' entry
          img_no_string = str(view_warning.no_images).split(',')
          filenames = []
          for n in img_no_string:
            if '-' in n:
              img_limits = [int(i) for i in n.split('-')]
              start = min(img_limits)
              end = max(img_limits)
              if start <= len(img_list) and end <= len(img_list):
                filenames.extend(img_list[start:end])
            else:
              if int(n) <= len(img_list):
                filenames.append(img_list[int(n)])
          file_string = ' '.join(filenames)
        else:
          return
        view_warning.Close()
      elif viewer == 'distl.image_viewer' and len(img_list) > 1:
        wx.MessageBox('distl.image_viewer can show only one image', 'Warning',
                      wx.OK | wx.ICON_EXCLAMATION)
        file_string = img_list[0]
      else:
        file_string = ' '.join(img_list)

      viewer = ImageViewerThread(self,
                                 viewer=viewer,
                                 file_string=file_string,
                                 img_type=img_type)
      viewer.start()

  def onDelButton(self, e):
    item = e.GetEventObject().GetParent()
    self.delete_button(item.index)

  def delete_all(self):
    for idx in range(self.ctr.GetItemCount()):
      self.delete_button(index=0)
    self.window.btn_run.Disable()

  def delete_button(self, index):
    try:
      self.image_count -= int(self.ctr.GetItemText(index))
    except ValueError:
      self.image_count = 0
    self.all_data_images.pop(self.ctr.GetItemData(index).path, None)
    self.update_total_image_count()
    self.ctr.DeleteItem(index)

    # Refresh widget and list item indices
    for i in range(self.ctr.GetItemCount()):
      item_data = self.ctr.GetItemData(i)
      item_data.id = i
      item_data.buttons.index = i
      item_data.type.index = i
      type_choice = self.ctr.GetItemWindow(i, col=2)
      type_selection = item_data.type.type.GetSelection()
      type_choice.type.SetSelection(type_selection)
      self.ctr.SetItemData(i, item_data)

    if self.ctr.GetItemCount() == 0:
      self.window.btn_run.Disable()

  def onInfoButton(self, e):
    """ Info / alert / error button (will change depending on circumstance) """
    idx = e.GetEventObject().GetParent().index
    item_obj = self.ctr.GetItemData(idx)
    item_type = item_obj.type.type.GetString(item_obj.type_selection)

    if item_obj.warning:
      wx.MessageBox(item_obj.info['WARNING'], 'Warning', wx.OK |
                    wx.ICON_EXCLAMATION)
    else:
      wx.MessageBox(item_obj.info[item_type], 'Info', wx.OK |
                    wx.ICON_INFORMATION)

  def update_total_image_count(self):
    self.txt_total_images.SetLabel("{} total images".format(self.image_count))

class TextFileView(wx.Dialog):
  def __init__(self, parent,
               contents=None,
               *args, **kwargs):

    dlg_style = wx.CAPTION | wx.CLOSE_BOX | wx.RESIZE_BORDER | wx.STAY_ON_TOP

    wx.Dialog.__init__(self, parent, style=dlg_style,
                        size=(600, 500),
                        *args, **kwargs)

    self.main_sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.main_sizer)

    from wx.lib.scrolledpanel import ScrolledPanel

    self.txt_panel = ScrolledPanel(self)
    self.txt_sizer = wx.BoxSizer(wx.VERTICAL)
    self.txt_panel.SetSizer(self.txt_sizer)

    self.txt = wx.StaticText(self.txt_panel, label=contents)
    self.txt_sizer.Add(self.txt)

    self.txt_panel.SetupScrolling()
    self.main_sizer.Add(self.txt_panel, 1, flag=wx.EXPAND | wx.ALL, border=10)

    # Dialog control
    self.main_sizer.Add(self.CreateSeparatedButtonSizer(wx.OK),
                        flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL, border=10)


class ViewerWarning(wx.Dialog):
  def __init__(self, parent,
               img_list_length = None,
               *args, **kwargs):

    dlg_style = wx.CAPTION | wx.CLOSE_BOX | wx.RESIZE_BORDER | wx.STAY_ON_TOP
    wx.Dialog.__init__(self, parent, style=dlg_style,
                        size=(400, 400),
                        *args, **kwargs)

    self.main_sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.main_sizer)

    self.img_list_length = img_list_length
    self.no_images = 0

    self.opt_sizer = wx.FlexGridSizer(6, 3, 10, 10)

    self.rb1_img_view = wx.RadioButton(self, label='First 1 image',
                                       style=wx.RB_GROUP)
    self.rb2_img_view = wx.RadioButton(self, label='First 10 images')
    self.rb3_img_view = wx.RadioButton(self, label='First 50 images')
    self.rb4_img_view = wx.RadioButton(self, label='First 100 images')
    self.rb5_img_view = wx.RadioButton(self, label='All {} images'
                                      ''.format(self.img_list_length))
    self.rb_custom = wx.RadioButton(self, label='Other: ')
    self.opt_custom = wx.TextCtrl(self, size=(100, -1))
    self.opt_custom.Disable()
    self.txt_custom = wx.StaticText(self, label='images')
    self.txt_custom.Disable()

    self.opt_sizer.AddMany([self.rb1_img_view, (0, 0), (0, 0),
                            self.rb2_img_view, (0, 0), (0, 0),
                            self.rb3_img_view, (0, 0), (0, 0),
                            self.rb4_img_view, (0, 0), (0, 0),
                            self.rb5_img_view, (0, 0), (0, 0),
                            self.rb_custom, self.opt_custom, self.txt_custom])

    # Grey out irrelevant radio buttons
    if self.img_list_length < 100:
      self.rb4_img_view.Disable()
    if self.img_list_length < 50:
      self.rb3_img_view.Disable()
    if self.img_list_length < 10:
      self.rb3_img_view.Disable()

    self.main_sizer.Add(self.opt_sizer, flag=wx.ALL, border=10)

    # Dialog control
    self.main_sizer.Add(self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL),
                        flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL, border=10)

    self.rb1_img_view.Bind(wx.EVT_RADIOBUTTON, self.onCustom)
    self.rb2_img_view.Bind(wx.EVT_RADIOBUTTON, self.onCustom)
    self.rb3_img_view.Bind(wx.EVT_RADIOBUTTON, self.onCustom)
    self.rb4_img_view.Bind(wx.EVT_RADIOBUTTON, self.onCustom)
    self.rb5_img_view.Bind(wx.EVT_RADIOBUTTON, self.onCustom)
    self.rb_custom.Bind(wx.EVT_RADIOBUTTON, self.onCustom)

    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)

  def onCustom(self, e):
    if self.rb_custom.GetValue():
      self.txt_custom.Enable()
      self.opt_custom.Enable()
      if self.img_list_length < 25:
        value = str(self.img_list_length)
      else:
        value = '25'
      self.opt_custom.SetValue(value)
    else:
      self.txt_custom.Disable()
      self.opt_custom.Disable()
      self.opt_custom.SetValue('')

  def onOK(self, e):
    if self.rb1_img_view.GetValue():
      self.no_images = '1'
    elif self.rb2_img_view.GetValue():
      self.no_images = '1-10'
    elif self.rb3_img_view.GetValue():
      self.no_images = '1-50'
    elif self.rb4_img_view.GetValue():
      self.no_images = '1-100'
    elif self.rb5_img_view.GetValue():
      self.no_images = '1-{}'.format(self.img_list_length)
    elif self.rb_custom.GetValue():
      self.no_images = self.opt_custom.GetValue()
    self.EndModal(wx.ID_OK)

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
  """ Generic panel will place a table w/ x and y labels
      Data must be a list of lists for multi-column tables """

  def __init__(self, parent,
               clabels=None,
               clabel_size=wx.DefaultSize,
               rlabels=None,
               rlabel_size=wx.DefaultSize,
               contents=None,
               label_style='normal',
               content_style='normal'):

    CtrlBase.__init__(self, parent=parent, label_style=label_style,
                      content_style=content_style)
    if contents is None:
      contents = []
    if rlabels:
      nrows = len(rlabels) + 1
    else:
      nrows = 1

    if clabels:
      ncols = len(clabels) + 1
    else:
      ncols = 2
    self.sizer = wx.FlexGridSizer(nrows, ncols, 10, 10)

    # add column labels (xlabels)
    if clabels:
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
        if item is None:
          item = ''
        cell = wx.StaticText(self, label=item)
        cell.SetFont(self.cfont)
        self.sizer.Add(cell)


class RichTextTableCtrl(CtrlBase):
  """ Generic panel will place a table w/ x and y labels
      Data must be a list of lists for multi-column tables """

  def __init__(self, parent,
               clabels=None,
               rlabels=None,
               contents=None,
               label_style='normal',
               content_style='teletype bold'):

    CtrlBase.__init__(self, parent=parent, label_style=label_style,
                      content_style=content_style)

    # Generate column widths
    if clabels:
      col_w = [len(l)+3 for l in clabels]
    else:
      col_w = [len(i)+3 for i in contents[0]]
    for row in contents:
      for item in row:
        idx = row.index(item)
        item_width = len(item) + 3
        col_w[idx] = item_width if item_width > col_w[idx] else col_w[idx]

    # Generate row label width
    row_w = 0
    for rlabel in rlabels:
      rlabel_width = len(rlabel) + 3
      row_w = rlabel_width if rlabel_width > row_w else row_w

    # Generate table
    if clabels:
      table_txt = ' ' * row_w
      for l in clabels:
        idx = clabels.index(l)
        spacer = ' ' * (col_w[idx] - len(l))
        label = l + spacer
        table_txt += label
      table_txt += '\n'
    else:
      table_txt = ''

    lines = []
    for l in rlabels:
      # Set row label
      spacer = ' ' * (row_w - len(l))
      line = l + spacer

      # Add data to table
      c_index = rlabels.index(l)
      row_contents = contents[c_index]
      for item in row_contents:
        i_idx = row_contents.index(item)
        spacer = ' ' * (col_w[i_idx] - len(item))
        if item is None:
          line += spacer
        else:
          line += item + spacer

      lines.append(line)

    table_txt += '\n'.join(lines)

    # Generate RichTextCtrl
    self.ctr = wx.richtext.RichTextCtrl(self, style=wx.NO_BORDER |
                                                    wx.TE_DONTWRAP |
                                                    wx.richtext.RE_READONLY)
    self.ctr.BeginLineSpacing(lineSpacing=15)
    self.cfont = wx.Font(norm_font_size, wx.FONTFAMILY_TELETYPE,
                         wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL, False)

    self.ctr.BeginFont(self.cfont)
    self.ctr.SetValue(value=table_txt)
    self.ctr.SetBackgroundColour(self.GetBackgroundColour())

    # Create sizer and add control to it
    self.sizer = wx.BoxSizer()
    self.sizer.Add(self.ctr, 1, flag=wx.EXPAND)
    self.SetSizer(self.sizer)

    self.Layout()


class Menu(wx.Menu):
  """ Customizable context menu for IOTA GUI (based on the base Menu class
      for Phenix GUI) """

  def __init__(self, frame, *args, **kwargs):
    """ Constructor
    :param frame: parent window
    """
    wx.Menu.__init__(self, *args, **kwargs)
    self.frame = frame

  def add_commands(self, command_list=None):
    if not command_list:
      raise Sorry('Cannot initialize context menu: no commands found!')
    for (label, function) in command_list:
      self.add_command(label, function)

  def add_command(self, label, function):
    if label is None:
      self.AppendSeparator()
    else:
      menu_item = self.Append(wx.ID_ANY, label)
      if function is not None:
        self.frame.Bind(wx.EVT_MENU, function, menu_item)

  def get_event_item_label (self, event):
    item = self.FindItemById(event.GetId())
    return item.GetText()

# ---end
