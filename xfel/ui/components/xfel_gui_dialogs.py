from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 06/04/2016
Last Changed: 06/04/2016
Description : XFEL UI Custom Dialogs
'''

import wx
import os
from wx.lib.mixins.listctrl import TextEditMixin, getListCtrlSelection
from wx.lib.scrolledpanel import ScrolledPanel

import xfel.ui.components.xfel_gui_controls as gctr

icons = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'icons/')

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


class EdListCtrl(wx.ListCtrl, TextEditMixin):
  ''' TextEditMixin allows any column to be edited. '''

  # ----------------------------------------------------------------------
  def __init__(self, parent, ID=wx.ID_ANY, pos=wx.DefaultPosition,
               size=wx.DefaultSize, style=0):
    """Constructor"""
    wx.ListCtrl.__init__(self, parent, ID, pos, size, style)
    TextEditMixin.__init__(self)

class BaseDialog(wx.Dialog):
  def __init__(self, parent,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):
    wx.Dialog.__init__(self, parent, *args, **kwargs)

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


# --------------------------------- Dialogs ---------------------------------- #

class SettingsDialog(BaseDialog):
  ''' Initial settings for cctbx.xfel; accepts DB credentials (may be
  separate dialog), populates experiment name / tag, separate dialog for
  multiprocessing and other settings; starts sentinels on OK '''

  def __init__(self, parent,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):


    BaseDialog.__init__(self, parent,
                        label_style=label_style,
                        content_style=content_style,
                        size=(600, 200),
                        style=wx.NO_BORDER,
                        *args, **kwargs)

    self.main_sizer = wx.BoxSizer(wx.VERTICAL)

    # Experiment tag and DB Credentials button
    self.db_cred = gctr.TextButtonCtrl(self,
                                       label='Experiment Tag',
                                       label_style='bold',
                                       label_size=(150, -1),
                                       big_button=True,
                                       big_button_label='DB Credentials...',
                                       big_button_size=(130, -1),
                                       value='')
    self.main_sizer.Add(self.db_cred,
                        flag=wx.EXPAND | wx.ALL,
                        border=10)

    # Experiment name control
    self.experiment = gctr.TextButtonCtrl(self,
                                          label='Experiment',
                                          label_style='bold',
                                          label_size=(150, -1),
                                          big_button_size=(130, -1),
                                          value='')
    self.main_sizer.Add(self.experiment,
                        flag=wx.EXPAND | wx.ALL,
                        border=10)

    # Output folder text control w/ Browse / magnifying glass button
    current_folder = os.path.abspath(os.curdir)
    self.output = gctr.TextButtonCtrl(self,
                                      label='Output',
                                      label_style='bold',
                                      label_size=(150, -1),
                                      big_button=True,
                                      big_button_label='Browse...',
                                      big_button_size=(120, -1),
                                      value=current_folder)
    self.main_sizer.Add(self.output,
                        flag=wx.EXPAND | wx.ALL,
                        border=10)

    self.btn_mp = wx.Button(self, label='Multiprocessing...')
    self.btn_op = wx.Button(self, label='Advanced Settings...')
    self.btn_OK = wx.Button(self, label="OK", id=wx.ID_OK)
    self.btn_cancel = wx.Button(self, label="Cancel", id=wx.ID_CANCEL)

    button_sizer = wx.FlexGridSizer(1, 5, 0, 10)
    button_sizer.AddMany([(self.btn_mp),
                          (self.btn_op),
                          (0,0),
                          (self.btn_OK),
                          (self.btn_cancel)])

    button_sizer.AddGrowableCol(2)
    self.main_sizer.Add(button_sizer,
                        flag=wx.EXPAND | wx.ALL,
                        border=10)
    self.SetSizer(self.main_sizer)


class TagDialog(BaseDialog):
  def __init__(self, parent,
               label_style='bold',
               content_style='normal',
               tags=[],
               *args, **kwargs):
    BaseDialog.__init__(self, parent, label_style=label_style,
                        content_style=content_style, *args, **kwargs)

    self.db_tags = tags
    self.deleted_tags = []
    self.new_tags = []
    self.edited_tags =[]

    self.main_sizer = wx.BoxSizer(wx.VERTICAL)
    self.top_sizer = wx.BoxSizer(wx.HORIZONTAL)

    self.button_panel = wx.Panel(self)
    self.button_sizer = wx.BoxSizer(wx.VERTICAL)
    self.btn_add = wx.Button(self.button_panel, size=(120, -1),
                             label='Add Tag')
    self.btn_rmv = wx.Button(self.button_panel, size=(120, -1),
                             label='Remove Tags')
    self.btn_clr = wx.Button(self.button_panel, size=(120, -1),
                             label='Clear All')
    self.button_sizer.Add(self.btn_add)
    self.button_sizer.Add(self.btn_rmv)
    self.button_sizer.Add(self.btn_clr)
    self.button_panel.SetSizer(self.button_sizer)

    self.tag_panel = ScrolledPanel(self, size=(500, 400))
    self.tag_list = EdListCtrl(self.tag_panel,
                               style=wx.LC_REPORT | wx.SUNKEN_BORDER)
    self.tag_sizer = wx.BoxSizer(wx.VERTICAL)
    self.tag_panel.SetSizer(self.tag_sizer)

    self.tag_list.InsertColumn(0, 'Sample Tag', width=200)
    self.tag_list.InsertColumn(1, 'Comments', width=300)

    # Populate tags with current values from db
    if len(self.db_tags) > 0:
      self.index = 0
      for tag in self.db_tags:
        self.tag_list.InsertStringItem(self.index, tag.name)
        self.tag_list.SetStringItem(self.index, 1, tag.comment)
        self.tag_list.SetItemData(self.index, tag.tag_id)
        self.index += 1

    self.tag_sizer.Add(self.tag_list, 1, flag=wx.EXPAND)

    # Add panels to main sizer
    self.top_sizer.Add(self.button_panel,
                       flag=wx.LEFT, border=10)
    self.top_sizer.Add(self.tag_panel,
                       flag=wx.EXPAND | wx.RIGHT | wx.LEFT, border=10)
    self.main_sizer.Add(self.top_sizer,
                        flag=wx.EXPAND| wx.TOP | wx.BOTTOM, border=10)
    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
    self.main_sizer.Add(dialog_box,
                   flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                   border=10)

    self.SetSizer(self.main_sizer)
    self.Layout()

    # Button bindings
    self.Bind(wx.EVT_BUTTON, self.onAdd, self.btn_add)
    self.Bind(wx.EVT_BUTTON, self.onRemove, self.btn_rmv)
    self.Bind(wx.EVT_BUTTON, self.onClear, self.btn_clr)
    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)

  def onAdd(self, e):
    ''' Add a string item to list; focus on item & provide default tag name'''
    new_tag = ('default tag {}'.format(self.index), '')
    self.new_tags.append(new_tag)
    self.tag_list.InsertStringItem(self.index, new_tag[0])
    self.tag_list.SetStringItem(self.index, 1, new_tag[1])
    self.tag_list.SetItemData(self.index, None)
    self.tag_list.Select(self.index)
    self.tag_list.Focus(self.index)
    self.index += 1

  def onRemove(self, e):
    selected_indices = getListCtrlSelection(self.tag_list)
    tag_ids = [self.tag_list.GetItemData(i) for i in selected_indices]
    print tag_ids
    self.deleted_tags = [i for i in self.db_tags if i.tag_id in tag_ids]

    for i in range(len(selected_indices)):
      to_delete = getListCtrlSelection(self.tag_list)
      self.tag_list.DeleteItem(to_delete[0])

  def onClear(self, e):
    warning = wx.MessageDialog(self,
                               message='Are you sure you want to delete all '
                                       'tags? \n This cannot be reversed!',
                               caption='Warning',
                               style=wx.YES_NO | wx.ICON_EXCLAMATION)
    if (warning.ShowModal() == wx.ID_YES):
      self.tag_list.DeleteAllItems()
      self.index = 0
      self.db_tags = []

  def onOK(self, e):
    count = self.tag_list.GetItemCount()
    if count == 0:
      warning = wx.MessageDialog(self,
                                 message='Must have at least one tag!',
                                 caption='Warning',
                                 style=wx.OK | wx.ICON_EXCLAMATION)
      warning.ShowModal()
    else:
      try:
        # Update names for edited tags
        # TODO: Actually connect to actual tags in DB
        all_items = [(self.tag_list.GetItemData(i),
                    self.tag_list.GetItem(itemId=i, col=0),
                    self.tag_list.GetItem(itemId=i, col=1))
                   for i in range(count)]
        edited_items = [i for i in all_items if i[0] is not None]
        for tag in self.db_tags:
          for item in edited_items:
            if tag.tag_id == item[0]:
              if tag.name != item[1]:
                tag.name = item[1]
              if tag.comment != item[2]:
                tag.comment = item[2]

        # Delete tags from DB
        for tag in self.deleted_tags:
          self.parent.parent.parent.db.delete_tag(tag)

        # Add new tags to DB
        for tag in self.new_tags:
          self.parent.parent.parent.db.create_tags(name=tag[0], comment=tag[1])
      except Exception:
        pass

      e.Skip()
