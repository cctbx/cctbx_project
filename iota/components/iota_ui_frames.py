from __future__ import division, print_function, absolute_import

from past.builtins import range

import iota.components.iota_ui_dialogs

'''
Author      : Lyubimov, A.Y.
Created     : 01/17/2017
Last Changed: 11/05/2018
Description : IOTA GUI Windows / frames
'''

import os
import wx
from wxtbx import bitmaps
import wx.lib.buttons as btn
from wx import richtext as rt
from wx.aui import AuiNotebook

import numpy as np
import time
import warnings
import multiprocessing

import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.patches as mpatches
import matplotlib.path as path
from mpl_toolkits.mplot3d import Axes3D

assert Axes3D

from libtbx import easy_run
from libtbx import easy_pickle as ep
from libtbx.utils import to_unicode

from prime.postrefine.mod_mx import mx_handler
import prime.postrefine.mod_plotter as ppl

from iota.components.iota_analysis import Analyzer, Plotter
import iota.components.iota_ui_controls as ct
import iota.components.iota_threads as thr
import iota.components.iota_ui_dialogs as d
import iota.components.iota_utils as ut

f = ut.WxFlags()
ginp = ut.InputFinder()

# Platform-specific stuff
# TODO: Will need to test this on Windows at some point
if wx.Platform == '__WXGTK__':
  plot_font_size = 10
  norm_font_size = 10
  button_font_size = 12
  LABEL_SIZE = 14
  CAPTION_SIZE = 12
  python = 'python'
elif wx.Platform == '__WXMAC__':
  plot_font_size = 9
  norm_font_size = 12
  button_font_size = 14
  LABEL_SIZE = 14
  CAPTION_SIZE = 12
  python = "Python"
elif (wx.Platform == '__WXMSW__'):
  plot_font_size = 9
  norm_font_size = 9
  button_font_size = 11
  LABEL_SIZE = 11
  CAPTION_SIZE = 9
  python = "Python"  # TODO: make sure it's right!

# ------------------------------ Input Window -------------------------------- #

class InputWindow(iota.components.iota_ui_dialogs.BasePanel):
  """ Input window - data input, description of project """

  def __init__(self, parent, phil):
    iota.components.iota_ui_dialogs.BasePanel.__init__(self, parent=parent)

    # Generate default parameters
    self.input_phil = phil
    self.target_phil = None
    self.gparams = self.input_phil.extract()
    self.imageset = None

    if str(self.gparams.advanced.processing_backend).lower() == 'cctbx':
      target = self.gparams.cctbx_ha14.target
    elif str(self.gparams.advanced.processing_backend).lower() == 'dials':
      target = self.gparams.cctbx_xfel.target
    else:
      target = None

    if str(target).lower != 'none':
      try:
        with open(target, 'r') as pf:
          self.target_phil = pf.read()
      except Exception:
        self.target_phil = None

    self.project_folder = ct.InputCtrl(self, label='Project Folder: ',
                                       label_size=(150, -1),
                                       label_style='bold',
                                       value=os.path.abspath(os.curdir),
                                       buttons=True)

    self.project_title = ct.InputCtrl(self, label='Description',
                                      label_size=(150, -1))

    # List control to add / manage input items
    self.input = FileListCtrl(self)

    # Put everything into main sizer
    self.main_sizer.Add(self.project_title, flag=f.expand, border=10)
    self.main_sizer.Add(self.project_folder, flag=f.expand, border=10)
    self.main_sizer.Add(self.input, 1, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(wx.StaticLine(self), flag=wx.EXPAND)

    # Options
    opt_box = wx.FlexGridSizer(1, 5, 0, 15)
    total_procs = multiprocessing.cpu_count()
    self.opt_spc_nprocs = ct.SpinCtrl(self,
                                      label='No. Processors: ',
                                      label_size=(120, -1),
                                      ctrl_min = 1,
                                      ctrl_value = str(int(total_procs / 2)))
    self.opt_btn_import = wx.Button(self, label='Import options...')
    self.opt_btn_process = wx.Button(self, label='Processing options...')
    self.opt_btn_analysis = wx.Button(self, label='Analysis options...')
    opt_box.AddMany([(self.opt_spc_nprocs),
                     (0, 0),
                     (self.opt_btn_import),
                     (self.opt_btn_process),
                     (self.opt_btn_analysis)])
    opt_box.AddGrowableCol(1)
    self.main_sizer.Add(opt_box, flag=f.expand, border=10)

    # Button bindings
    self.project_folder.btn_browse.Bind(wx.EVT_BUTTON, self.onOutputBrowse)
    self.project_folder.btn_mag.Bind(wx.EVT_BUTTON, self.onMagButton)
    self.Bind(wx.EVT_BUTTON, self.onImportOptions, self.opt_btn_import)
    self.Bind(wx.EVT_BUTTON, self.onProcessOptions, self.opt_btn_process)
    self.Bind(wx.EVT_BUTTON, self.onAnalysisOptions, self.opt_btn_analysis)


  def onMagButton(self, e):
    dirview = d.DirView(self, title='Current Folder')
    if dirview.ShowModal() == wx.ID_OK:
      dirview.Destroy()


  def onImportOptions(self, e):
    e.Skip()

  def onProcessOptions(self, e):
    e.Skip()

  def onAnalysisOptions(self, e):
    e.Skip()


  def onInfo(self, e):
    """ On clicking the info button """
    info_txt = '''Input diffraction images here. IOTA accepts either raw images (mccd, cbf, img, etc.) or image pickles. Input can be either a folder with images, or a text file with a list of images.'''
    info = wx.MessageDialog(None, info_txt, 'Info', wx.OK)
    info.ShowModal()


  def onOutputBrowse(self, e):
    """ On clicking the Browse button: show the DirDialog and populate 'Output'
        box w/ selection """
    dlg = wx.DirDialog(self, "Choose the output directory:",
                       style=wx.DD_DEFAULT_STYLE)
    if dlg.ShowModal() == wx.ID_OK:
      self.project_folder.ctr.SetValue(dlg.GetPath())
    dlg.Destroy()
    e.Skip()

class FileListCtrl(ct.CustomListCtrl):
  """ File list window for the input tab """

  def __init__(self, parent, size=(-1, 300)):
    ct.CustomListCtrl.__init__(self, parent=parent, size=size)

    self.parent = parent
    self.main_window = parent.GetParent()

    # Initialize dictionaries for imported data types
    self.all_data_images = {}
    self.all_img_objects = {}
    self.all_proc_pickles = {}
    self.image_count = 0

    # Generate columns
    self.ctr.InsertColumn(0, "")
    self.ctr.InsertColumn(1, "Input Path")
    self.ctr.InsertColumn(2, "Input Type")
    self.ctr.InsertColumn(3, "Action")
    self.ctr.setResizeColumn(2)

    # Add file / folder buttons
    # self.button_sizer = wx.BoxSizer(wx.HORIZONTAL)
    self.button_sizer = wx.FlexGridSizer(1, 3, 0, 10)
    self.btn_add_file = wx.Button(self, label='Add File...')
    self.btn_add_dir = wx.Button(self, label='Add Folder...')
    self.txt_total_images = wx.StaticText(self, label='')
    self.button_sizer.Add(self.btn_add_file)
    self.button_sizer.Add(self.btn_add_dir)
    self.button_sizer.Add(self.txt_total_images)
    self.button_sizer.AddGrowableCol(1)

    self.sizer.Add(self.button_sizer, flag=wx.EXPAND | wx.TOP | wx.BOTTOM,
                   border=10)

    # Event bindings
    self.Bind(wx.EVT_BUTTON, self.onAddFile, self.btn_add_file)
    self.Bind(wx.EVT_BUTTON, self.onAddFolder, self.btn_add_dir)

  def onAddFile(self, e):
    file_dlg = wx.FileDialog(self,
                             message="Load File",
                             defaultDir=os.curdir,
                             defaultFile="*",
                             wildcard="*",
                             style=wx.OPEN | wx.FD_FILE_MUST_EXIST |
                                   wx.FD_MULTIPLE)
    if file_dlg.ShowModal() == wx.ID_OK:
      files = file_dlg.GetPaths()
      for item in files:
        self.add_item(item)
    file_dlg.Destroy()
    e.Skip()

  def onAddFolder(self, e):
    dlg = wx.DirDialog(self, "Load Folder:",
                       style=wx.DD_DEFAULT_STYLE)
    if dlg.ShowModal() == wx.ID_OK:
      self.add_item(dlg.GetPath())
    dlg.Destroy()
    e.Skip()

  def set_type_choices(self, path):
    # Determine what type of input this is and present user with choices
    # (this so far works for images ONLY)
    type_choices = ['[  SELECT INPUT TYPE  ]']
    preferred_selection = 0
    inputs, input_type = ginp.get_input(path)

    if os.path.isdir(path):
      type_choices.extend(['raw image folder', 'image pickle folder'])
      if input_type in type_choices:
        preferred_selection = type_choices.index(input_type)
    elif os.path.isfile(path):
      if input_type in ('image pickle file', 'raw image file'):
        type_choices.extend(['raw image file', 'image pickle file'])
        if input_type in type_choices:
          preferred_selection = type_choices.index(input_type)
      elif input_type in ('raw image list', 'image pickle list'):
        type_choices.extend(['raw image list', 'image pickle list'])
        if input_type in type_choices:
          preferred_selection = type_choices.index(input_type)
    return inputs, type_choices, preferred_selection

  def add_item(self, path):
    # Generate item
    inputs, inp_choices, inp_sel = self.set_type_choices(path)
    type_choice = ct.DataTypeChoice(self.ctr,
                                    choices=inp_choices)
    item = ct.InputListItem(path=path,
                            type=type_choice,
                            buttons=ct.MiniButtonBoxInput(self.ctr))

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
    if item.type.type.GetString(inp_sel) in ['raw image folder',
                                             'image pickle folder',
                                             'image pickle file',
                                             'raw image file',
                                             'raw image list',
                                             'image pickle list']:
      self.main_window.toolbar.EnableTool(self.main_window.tb_btn_run.GetId(), True)
      self.all_data_images[item.path] = inputs

      # Calculate # of images and display w/ item
      self.ctr.SetStringItem(idx, 0, str(len(inputs)))

      if "image" in item.type.type.GetString(inp_sel):
        view_bmp = bitmaps.fetch_custom_icon_bitmap('image_viewer16')
        item.buttons.btn_mag.SetBitmapLabel(view_bmp)
    else:
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

    if len(inputs) > 0:
      self.image_count += len(inputs)
      if self.image_count > 0:
        self.update_total_image_count()

    self.main_window.Layout()

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
    if data_items > 0:
      self.main_window.toolbar.EnableTool(self.main_window.tb_btn_run.GetId(),
                                          True)
    else:
      self.main_window.toolbar.EnableTool(self.main_window.tb_btn_run.GetId(),
                                          False)
    e.Skip()


  def onMagButton(self, e):
    idx = e.GetEventObject().GetParent().index
    item_obj = self.ctr.GetItemData(idx)
    path = item_obj.path
    type = item_obj.type.type.GetString(item_obj.type_selection)

    if os.path.isfile(path):
      if type in ('raw image file', 'image pickle file'):
        self.view_images([path], img_type=type)
      elif type in ('raw image list', 'image pickle list'):
        with open(path, 'r') as f:
          file_list = [i.replace('\n', '') for i in f.readlines()]
          self.view_images(file_list, img_type=type)
      elif type == 'text':
        with open(path, 'r') as f:
          file_list = f.readlines()
          msg = ' '.join(file_list)
          textview = d.TextFileView(self, title=path, contents=msg)
          textview.ShowModal()
      else:
        wx.MessageBox('Unknown file format', 'Warning',
                      wx.OK | wx.ICON_EXCLAMATION)
    elif os.path.isdir(path):
      file_list, _ = ginp.get_input(path)
      self.view_images(file_list, img_type=type)

  def view_images(self, img_list, img_type=None):
    """ Launches image viewer (depending on backend) """
    viewer = self.parent.gparams.gui.image_viewer
    if viewer == 'cxi.view' and 'pickle' not in img_type:
        wx.MessageBox('cxi.view only accepts image pickles', 'Warning',
                      wx.OK | wx.ICON_EXCLAMATION)
    else:
      if len(img_list) > 10:
        view_warning = d.ViewerWarning(self, len(img_list))
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

      viewer = thr.ImageViewerThread(self,
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

# ----------------------------  Processing Window ---------------------------  #

class LogTab(wx.Panel):
  def __init__(self, parent):
    wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)

    self.log_sizer = wx.BoxSizer(wx.VERTICAL)
    self.log_window = rt.RichTextCtrl(self,
                                      style=rt.RE_MULTILINE |
                                            rt.RE_READONLY |
                                            wx.TE_DONTWRAP)
    self.log_window.SetFont(wx.Font(9, wx.TELETYPE, wx.NORMAL, wx.NORMAL, False))
    self.log_sizer.Add(self.log_window, proportion=1, flag= wx.EXPAND | wx.ALL, border=10)

    self.find_string = ct.TextCtrlWithButtons(self,
                                              buttons=[('Forward', -1),
                                                       ('Reverse', -1)],
                                              ctrl_label='Find String:')
    self.log_sizer.Add(self.find_string, flag=wx.EXPAND | wx.ALL, border=10)
    self.SetSizer(self.log_sizer)

    self.Bind(wx.EVT_BUTTON, self.onSearchForward, self.find_string.btn_forward)
    self.Bind(wx.EVT_BUTTON, self.onSearchReverse, self.find_string.btn_reverse)

  def onSearchForward(self, e):
    if self.log_window.GetCaretPosition() == -1:
      self.log_window.SetCaretPosition(0)
    pos = self.log_window.GetCaretPosition()
    search_string = self.find_string.txt_ctrl.GetValue().lower()
    log_string = self.log_window.GetValue()[pos:-1].lower()
    if search_string.replace(' ', '') not in ('', None):
      found_pos = log_string.find(search_string)
      if found_pos == -1:
        if pos > 0:
          msg_text = 'String Not Found! Search From the Top?'
          msg = wx.MessageDialog(None, msg_text, 'Not Found!',
                                 wx.YES_NO | wx.ICON_QUESTION)
          if msg.ShowModal() == wx.ID_YES:
            log_string = self.log_window.GetValue()[0:pos].lower()
            found_pos = log_string.find(search_string)
            if found_pos == -1:
              wx.MessageBox('String Not Found!', 'Not Found!',
                            wx.OK | wx.ICON_EXCLAMATION)
              return
          else:
            return
        else:
          wx.MessageBox('String Not Found!', 'Not Found!',
                        wx.OK | wx.ICON_EXCLAMATION)
          return
      else:
        found_pos += pos
      sel_range = (found_pos, found_pos + len(search_string))
      self.log_window.SetSelectionRange(sel_range)
    else:
      found_pos = 0
    self.log_window.SetCaretPosition(found_pos + len(search_string))
    if not self.log_window.IsPositionVisible(found_pos):
      self.log_window.ShowPosition(found_pos)

  def onSearchReverse(self, e):
    if self.log_window.GetCaretPosition() == -1:
      self.log_window.SetCaretPosition(0)
    pos = self.log_window.GetCaretPosition()
    search_string = self.find_string.txt_ctrl.GetValue().lower()
    full_log = self.log_window.GetValue()
    log_string = full_log[0:pos].lower()
    log_end = len(full_log)
    if search_string.replace(' ', '') not in ('', None):
      found_pos = log_string.rfind(search_string)
      if found_pos == -1:
        if pos < log_end:
          msg_text = 'String Not Found! Search From the Bottom?'
          msg = wx.MessageDialog(None, msg_text, 'Not Found!',
                                 wx.YES_NO | wx.ICON_QUESTION)
          if msg.ShowModal() == wx.ID_YES:
            log_string = full_log[pos:-1].lower()
            found_pos = log_string.rfind(search_string)
            if found_pos == -1:
              wx.MessageBox('String Not Found!', 'Not Found!',
                            wx.OK | wx.ICON_EXCLAMATION)
              return
            else:
              found_pos += pos
          else:
            return
        else:
          wx.MessageBox('String Not Found!', 'Not Found!',
                        wx.OK | wx.ICON_EXCLAMATION)
          return
      sel_range = (found_pos, found_pos + len(search_string))
      self.log_window.SetSelectionRange(sel_range)
    else:
      found_pos = 0
    self.log_window.SetCaretPosition(found_pos - len(search_string))
    if not self.log_window.IsPositionVisible(found_pos):
      self.log_window.ShowPosition(found_pos)


# class ProcessingTab(d.ScrolledPanel):
class ProcessingTab(wx.Panel):
  def __init__(self, parent):
    wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)
    self.gparams = None
    self.init = None
    self.parent = parent

    self.info = None

    self.hkl_view_axis = 'l'
    self.pick = {'image':None, 'index':0, 'axis':None, 'picked':False}
    self.proc_fnames = None

    self.main_fig_sizer = wx.GridBagSizer(0, 0)

    # Set regular font
    plt.rc('font', family='sans-serif', size=plot_font_size)
    plt.rc('mathtext', default='regular')

    # Integration figure (resolution & reflections / frame)
    self.int_panel = wx.Panel(self)
    int_sizer = wx.BoxSizer(wx.VERTICAL)
    self.int_panel.SetSizer(int_sizer)

    # Image info sizer
    self.info_sizer = wx.GridBagSizer(0, 5)
    self.info_txt = wx.TextCtrl(self.int_panel, style=wx.TE_READONLY)
    view_bmp = bitmaps.fetch_custom_icon_bitmap('image_viewer16')
    r_bmp = bitmaps.fetch_icon_bitmap('actions', '1rightarrow', size=16)
    l_bmp = bitmaps.fetch_icon_bitmap('actions', '1leftarrow', size=16)
    self.btn_right = btn.GenBitmapButton(self, bitmap=r_bmp)
    self.btn_left = btn.GenBitmapButton(self, bitmap=l_bmp)
    self.btn_viewer = btn.GenBitmapButton(self, bitmap=view_bmp)
    self.info_sizer.Add(self.info_txt, pos=(0, 1), flag=wx.EXPAND)
    self.info_sizer.Add(self.btn_left, pos=(0, 2))
    self.info_sizer.Add(self.btn_right, pos=(0, 3))
    self.info_sizer.Add(self.btn_viewer, pos=(0, 4))
    self.info_sizer.AddGrowableCol(1)
    int_sizer.Add(self.info_sizer, flag=wx.EXPAND)
    self.info_txt.Hide()
    self.btn_right.Hide()
    self.btn_left.Hide()
    self.btn_viewer.Hide()

    self.btn_right.Disable()
    self.btn_left.Disable()
    self.btn_viewer.Disable()

    self.Bind(wx.EVT_BUTTON, self.onImageView, self.btn_viewer)
    self.Bind(wx.EVT_BUTTON, self.onArrow, self.btn_right)
    self.Bind(wx.EVT_BUTTON, self.onArrow, self.btn_left)

    # Charts
    self.int_figure = Figure(figsize=(1, 2.5))
    self.int_figure.patch.set_visible(False)    # create transparent background
    int_gsp = gridspec.GridSpec(2, 1, wspace=0, hspace=0)

    # Resolution / No. strong reflections chart
    self.res_axes = self.int_figure.add_subplot(int_gsp[0])
    self.res_axes.set_ylabel('Resolution')
    self.res_axes.yaxis.get_major_ticks()[0].label1.set_visible(False)
    self.res_axes.yaxis.get_major_ticks()[-1].label1.set_visible(False)
    plt.setp(self.res_axes.get_xticklabels(), visible=False)
    self.nsref_axes = self.int_figure.add_subplot(int_gsp[1])
    self.nsref_axes.set_xlabel('Frame')
    self.nsref_axes.set_ylabel('Strong Spots')
    self.nsref_axes.yaxis.get_major_ticks()[0].label1.set_visible(False)
    self.nsref_axes.yaxis.get_major_ticks()[-1].label1.set_visible(False)

    self.info_txt.Show()
    self.btn_right.Show()
    self.btn_left.Show()
    self.btn_viewer.Show()

    self.int_figure.set_tight_layout(True)
    self.int_canvas = FigureCanvas(self.int_panel, -1, self.int_figure)
    int_sizer.Add(self.int_canvas, 1, flag=wx.EXPAND)

    # Wilson (<I> vs. res) plot
    self.wp_panel = wx.Panel(self)
    wp_sizer = wx.BoxSizer(wx.VERTICAL)
    self.wp_panel.SetSizer(wp_sizer)

    self.wp_figure = Figure(figsize=(0.3, 0.6))
    self.wp_figure.patch.set_visible(False)
    self.wp_axes = self.wp_figure.add_subplot(111)
    self.wp_axes.set_ylabel("<I>")

    self.wp_figure.set_tight_layout(True)
    self.wp_canvas = FigureCanvas(self.wp_panel, -1, self.wp_figure)
    wp_sizer.Add(self.wp_canvas, 1, flag=wx.EXPAND)

    # HKL (or slice) plot
    self.hkl_panel = wx.Panel(self)
    hkl_sizer = wx.BoxSizer(wx.VERTICAL)
    self.hkl_panel.SetSizer(hkl_sizer)

    self.hkl_figure = Figure(figsize=(0.3, 0.3))
    self.hkl_figure.patch.set_visible(False)    # create transparent background

    self.hkl_axes = self.hkl_figure.add_subplot(111, frameon=False)
    self.hkl_axes.set_xticks([])
    self.hkl_axes.set_yticks([])

    self.hkl_canvas = FigureCanvas(self.hkl_panel, -1, self.hkl_figure)
    hkl_sizer.Add(self.hkl_canvas, 1,flag=wx.EXPAND)

    self.hkl_sg = ct.OptionCtrl(self.hkl_panel, items=[('sg', 'P1')],
                                checkbox=True, checkbox_label='Space Group: ',
                                label_size=wx.DefaultSize,
                                ctrl_size=wx.DefaultSize)
    hkl_sizer.Add(self.hkl_sg, flag=wx.ALIGN_CENTER)

    # Processing summary figure
    self.proc_panel = wx.Panel(self, size=(-1, 120))
    proc_sizer = wx.BoxSizer(wx.VERTICAL)
    self.proc_panel.SetSizer(proc_sizer)

    self.proc_figure = Figure(figsize=(0.5, 2.5))
    self.proc_figure.patch.set_visible(False)
    self.sum_axes = self.proc_figure.add_subplot(111)
    self.sum_axes.axis('off')
    self.proc_figure.set_tight_layout(True)
    self.proc_canvas = FigureCanvas(self.proc_panel, -1, self.proc_figure)
    proc_sizer.Add(self.proc_canvas, flag=wx.EXPAND | wx.BOTTOM, border=10)

    self.main_fig_sizer.Add(self.int_panel, pos=(0, 0), span=(2, 6),
                            flag=wx.EXPAND)
    self.main_fig_sizer.Add(self.wp_panel, pos=(2, 0), span=(2, 4),
                            flag=wx.EXPAND)
    self.main_fig_sizer.Add(self.hkl_panel, pos=(2, 4), span=(2, 2),
                            flag=wx.EXPAND)
    self.main_fig_sizer.Add(self.proc_panel, pos=(4, 0), span=(1, 6),
                            flag=wx.EXPAND)

    self.main_fig_sizer.AddGrowableCol(0)
    self.main_fig_sizer.AddGrowableCol(1)
    self.main_fig_sizer.AddGrowableCol(2)
    self.main_fig_sizer.AddGrowableCol(3)
    self.main_fig_sizer.AddGrowableCol(4)
    self.main_fig_sizer.AddGrowableCol(5)
    self.main_fig_sizer.AddGrowableRow(1)
    self.main_fig_sizer.AddGrowableRow(2)
    self.main_fig_sizer.AddGrowableRow(3)

    cid = self.int_canvas.mpl_connect('pick_event', self.on_pick)
    sid = self.proc_canvas.mpl_connect('pick_event', self.on_bar_pick)
    xid = self.int_canvas.mpl_connect('button_press_event', self.on_button_press)
    xid = self.hkl_canvas.mpl_connect('button_press_event', self.on_hkl_press)
    xid = self.proc_canvas.mpl_connect('button_press_event', self.on_button_press)
    xid = self.proc_canvas.mpl_connect('button_release_event',
                                       self.on_button_release)

    self.SetSizer(self.main_fig_sizer)

  def draw_summary(self):
    try:
      nonzero = self.info.status_summary['nonzero']
      names = self.info.status_summary['names']
      patches = self.info.status_summary['patches']

      self.sum_axes.clear()
      self.sum_axes.axis([0, 1.01*len(self.info.img_list), -0.5, 0.75])
      self.sum_axes.axis('off')

      for i in range(len(nonzero)):
        percent = np.round(nonzero[i][2] / len(self.info.img_list) * 100)

        previous = [j[2] for j in nonzero[:i]]
        lf = np.sum(previous)
        barh = self.sum_axes.barh(0, nonzero[i][2],
                                  left=lf,
                                  color=nonzero[i][1],
                                  align='center',
                                  label='{}%'.format(percent),
                                  picker=True)
        patch = barh[0]
        bl = patch.get_xy()
        patches[i] = ((bl[0], bl[0] + patch.get_width()),
                      (bl[1], bl[1] + patch.get_height()))
        x = 0.5 * patch.get_width() + bl[0]
        y = 0.5 * patch.get_height() + bl[1]
        self.sum_axes.text(x, y, "{}%".format(percent),
                           ha='center', va='center')
      self.sum_axes.legend(names, ncol=len(nonzero),
                           bbox_to_anchor=(0, -1, 1, 1),
                           loc='upper center', fontsize=9,
                           frameon=False, handlelength=1)

      self.bracket = mpatches.Rectangle((0, 0), len(self.info.img_list), 1,
                                        facecolor='white', edgecolor='black',
                                        zorder=2, lw=2, visible=False)
      self.sum_axes.add_patch(self.bracket)

      # If bar is clicked, find category and files for display
      if self.pick['picked'] and self.pick['axis'] == 'summary':
        idx = self.pick['index']
        p_idx = [patches.index(i) for i in patches if
                 i[0][0] <= idx <= i[0][1]][0]
        cat = names[p_idx]
        if self.gparams.advanced.processing_backend == 'ha14':
          if cat == 'failed indexing / integration':
            cat = 'failed grid search'
          elif cat == 'failed filter':
            cat = 'failed prefilter'

        if cat not in ('integrated', 'not processed'):
          self.proc_fnames = [i.img_path for i in self.info.finished_objects
                              if i.fail == cat]
        elif cat == 'integrated':
          self.proc_fnames = [i.img_path for i in self.info.finished_objects if
                              i.fail is None and i.final['final'] is not None]
        elif cat == 'not processed':
          self.proc_fnames = [i.img_path for i in self.info.finished_objects if
                              i.fail is None and i.final['final'] is None]
        else:
          self.proc_fnames = ''

        px = patches[p_idx][0][0]
        pw = patches[p_idx][0][1] - patches[p_idx][0][0]
        py = patches[p_idx][1][0]
        ph = patches[p_idx][1][1] - patches[p_idx][1][0]
        self.bracket.set_bounds(px, py, pw, ph)
        self.bracket.set_visible(True)
        self.proc_panel.Refresh()

      self.proc_canvas.draw()

    except ValueError as e:
      print ('SUMMARY PLOT ERROR: ', e)

  def draw_plots(self):
    if sum(self.info.nref_list) > 0 and sum(self.info.res_list) > 0:
      self.draw_integration_plots()
      if self.info.b_factors:
        self.draw_b_factors()
      if self.info.merged_indices:
        self.draw_measured_indices()
      self.int_canvas.draw()
    self.Layout()

  def draw_integration_plots(self):
    try:
      # Strong reflections per frame
      self.nsref_axes.clear()
      nsref_ylabel = 'Reflections (I/{0}(I) > {1})' \
                     ''.format(r'$\sigma$', self.gparams.cctbx_ha14.selection.min_sigma)
      self.nsref = self.nsref_axes.scatter(self.info.nsref_x,
                                           self.info.nsref_y,
                                           s=45, marker='o',
                                           edgecolors='black',
                                           color='#ca0020',
                                           picker=True)
      nsref_med = self.nsref_axes.axhline(self.info.nsref_median,
                                          c='#ca0020', ls='--')

      self.nsref_axes.set_xlim(0, np.nanmax(self.info.nsref_x) + 2)
      nsref_ymax = np.nanmax(self.info.nsref_y) * 1.25 + 10
      if nsref_ymax == 0:
        nsref_ymax = 100
      self.nsref_axes.set_ylim(ymin=0, ymax=nsref_ymax)
      self.nsref_axes.set_ylabel(nsref_ylabel, fontsize=10)
      self.nsref_axes.set_xlabel('Frame')
      self.nsref_axes.yaxis.get_major_ticks()[0].label1.set_visible(False)
      self.nsref_axes.yaxis.get_major_ticks()[-1].label1.set_visible(False)

    except ValueError as e:
      print ("NSREF ERROR: ", e)

    try:
      # Resolution per frame
      self.res_axes.clear()
      res_m = np.isfinite(self.info.res_y)

      self.res = self.res_axes.scatter(self.info.res_x[res_m],
                                       self.info.res_y[res_m],
                                       s=45, marker='o', edgecolors='black',
                                       color='#0571b0', picker=True)

      with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        res_med = self.res_axes.axhline(self.info.res_median,
                                        c='#0571b0', ls='--')

      self.res_axes.set_xlim(0, np.nanmax(self.info.res_x) + 2)
      res_ymax = np.nanmax(self.info.res_y) * 1.1
      res_ymin = np.nanmin(self.info.res_y) * 0.9
      if res_ymin == res_ymax:
        res_ymax = res_ymin + 1
      self.res_axes.set_ylim(ymin=res_ymin, ymax=res_ymax)
      res_ylabel = 'Resolution ({})'.format(r'$\AA$')
      self.res_axes.set_ylabel(res_ylabel, fontsize=10)
      self.res_axes.yaxis.get_major_ticks()[0].label1.set_visible(False)
      self.res_axes.yaxis.get_major_ticks()[-1].label1.set_visible(False)
      plt.setp(self.res_axes.get_xticklabels(), visible=False)

      self.nsref_pick, = self.nsref_axes.plot(self.info.nsref_x[0],
                                              self.info.nsref_y[0],
                                              marker='o',
                                              ms=12, alpha=0.5,
                                              zorder=2,
                                              color='yellow', visible=False)
      self.res_pick, = self.res_axes.plot(self.info.res_x[0],
                                          self.info.res_x[0],
                                          marker='o',
                                          zorder=2,
                                          ms=12, alpha=0.5,
                                          color='yellow', visible=False)
      if self.pick['picked'] and self.pick['axis'] != 'summary':
        img = self.pick['image']
        idx = self.pick['index']
        axis = self.pick['axis']
        if axis == 'nsref':
          if not np.isnan(self.info.nsref_y[idx]):
            self.nsref_pick.set_visible(True)
            self.nsref_pick.set_data(self.info.nsref_x[idx],
                                     self.info.nsref_y[idx])
        elif axis == 'res':
          if not np.isnan(self.info.res_y[idx]):
            self.res_pick.set_visible(True)
            self.res_pick.set_data(self.info.res_x[idx], self.info.res_y[idx])
        self.info_txt.SetValue(img)
        self.btn_left.Enable()
        self.btn_right.Enable()
        self.btn_viewer.Enable()

    except ValueError as e:
      print ("RES ERROR: ", e)


  def draw_b_factors(self):
    self.wp_axes.clear()
    self.wp_axes.set_xlabel('B-factor')
    self.wp_axes.set_ylabel('Count')
    self.wp_axes.set_title('Wilson B-factor Histogram')
    if self.info.b_factors:
      self.wp_axes.hist(self.info.b_factors, 50, normed=False,
                        facecolor='#4575b4',
                        histtype='stepfilled')

  def draw_measured_indices(self):
    # Draw a h0, k0, or l0 slice of merged data so far
    self.hkl_axes.clear()
    try:
      self.hkl_colorbar.remove()
    except Exception:
      pass

    slice = self.info.merged_indices.slice(axis=self.hkl_view_axis,
                                           slice_start=0, slice_end=0)

    hkl = [i[0] for i in slice]
    freq = [i[1] for i in slice]

    if self.hkl_view_axis == 'l':
      x = [i[0] for i in hkl]
      y = [i[1] for i in hkl]
      self.hkl_axes.set_xlabel('h', weight='bold')
      self.hkl_axes.set_ylabel('k', weight='bold', rotation='horizontal')
    elif self.hkl_view_axis == 'k':
      x = [i[0] for i in hkl]
      y = [i[2] for i in hkl]
      self.hkl_axes.set_xlabel('h', weight='bold')
      self.hkl_axes.set_ylabel('l', weight='bold', rotation='horizontal')
    elif self.hkl_view_axis == 'h':
      x = [i[1] for i in hkl]
      y = [i[2] for i in hkl]
      self.hkl_axes.set_xlabel('k', weight='bold')
      self.hkl_axes.set_ylabel('l', weight='bold', rotation='horizontal')
    else:
      # if for some reason this goes to plotting without any indices
      x = np.zeros(500)
      y = np.zeros(500)
      freq = np.zeros(500)

    # Format plot
    pt_size = int(max(self.hkl_panel.GetSize()) / 250)
    if pt_size == 0:
      pt_size = 1

    hkl_scatter = self.hkl_axes.scatter(x, y, c=freq, cmap="jet", s=pt_size,
                                        edgecolor='none')
    self.hkl_axes.axhline(0, lw=0.5, c='black', ls='-')
    self.hkl_axes.axvline(0, lw=0.5, c='black', ls='-')
    self.hkl_axes.set_xticks([])
    self.hkl_axes.set_yticks([])
    self.hkl_axes.xaxis.set_label_coords(x=1, y=0.5)
    self.hkl_axes.yaxis.set_label_coords(x=0.5, y=1)

    try:
      xmax = abs(max(x, key=abs))
      ymax = abs(max(y, key=abs))
      self.hkl_axes.set_xlim(xmin=-xmax, xmax=xmax)
      self.hkl_axes.set_ylim(ymin=-ymax, ymax=ymax)
    except ValueError:
      pass

    norm = colors.Normalize(vmin=0, vmax=np.max(freq))
    self.hkl_colorbar = self.hkl_figure.colorbar(hkl_scatter,
                                                 ax=self.hkl_axes, cmap='jet',
                                                 norm=norm,
                                                 orientation='vertical',
                                                 aspect=40)

  def onImageView(self, e):
    filepath = self.info_txt.GetValue()
    viewer = self.gparams.gui.image_viewer
    if os.path.isfile(filepath):
      viewer = thr.ImageViewerThread(self,
                                     viewer=viewer,
                                     file_string=filepath)
      viewer.start()

  def view_proc_images(self):
    if self.proc_fnames is not None:
      view_warning = d.ViewerWarning(self.parent, len(self.proc_fnames))
      if view_warning.ShowModal() == wx.ID_OK:
        # parse 'other' entry
        img_no_string = str(view_warning.no_images).split(',')
        filenames = []
        for n in img_no_string:
          if '-' in n:
            img_limits = n.split('-')
            start = int(min(img_limits))
            end = int(max(img_limits))
            if start <= len(self.proc_fnames) and end <= len(self.proc_fnames):
              filenames.extend(self.proc_fnames[start:end])
          else:
            if int(n) <= len(self.proc_fnames):
              filenames.append(self.proc_fnames[int(n)])
        file_string = ' '.join(filenames)
        viewer = self.gparams.gui.image_viewer
        viewer = thr.ImageViewerThread(self,
                                       viewer=viewer,
                                       file_string=file_string)
        viewer.start()

  def onArrow(self, e):
    idx = self.pick['index']
    self.pick['image'] = None
    if e.GetId() == self.btn_right.GetId():
      direction = 1
    elif e.GetId() == self.btn_left.GetId():
      direction = -1

    search = True
    while search:
      idx += direction
      try:
        if self.pick['axis'] == 'nsref':
          if not np.isnan(self.nsref_y[idx]):
            self.pick['index'] = idx
            self.pick['image'] = self.info.finished_objects[idx].img_path
            search = False
          else:
            search = True
        if self.pick['axis'] == 'res':
          if not np.isnan(self.res_y[idx]):
            self.pick['index'] = idx
            self.pick['image'] = self.info.finished_objects[idx].img_path
            search = False
          else:
            search = True
      except IndexError as e:
        search = False
        self.pick['index'] = idx
        self.pick['image'] = None

    self.draw_plots()

  def on_pick(self, event):
    self.pick['picked'] = True
    idx = int(round(event.mouseevent.xdata)) - 1
    obj_i = [i for i in self.info.finished_objects if i.img_index == idx + 1]
    img = obj_i[0].img_path

    print ('{}: {}'.format(idx+1, img))
    self.pick['image'] = img
    self.pick['index'] = idx

    if event.mouseevent.inaxes == self.nsref_axes:
      self.pick['axis'] = 'nsref'
    elif event.mouseevent.inaxes == self.res_axes:
      self.pick['axis'] = 'res'
    self.draw_plots()

  def on_bar_pick(self, event):
    self.show_image_group(e=event.mouseevent)

  def show_image_group(self, e):
    self.pick['picked'] = True
    if e.inaxes == self.sum_axes:
      self.pick['axis'] = 'summary'
      self.pick['index'] = e.xdata
    self.draw_summary()

  def on_hkl_press(self, event):
    if event.inaxes == self.hkl_axes:
      if self.hkl_view_axis == 'h':
        self.hkl_view_axis = 'k'
      elif self.hkl_view_axis == 'k':
        self.hkl_view_axis = 'l'
      elif self.hkl_view_axis == 'l':
        self.hkl_view_axis = 'h'
      self.draw_plots()

  def on_button_press(self, event):
    if event.button != 1:
      self.pick['picked'] = False
      if event.inaxes == self.sum_axes:
        self.bracket.set_visible(False)
        self.draw_summary()
      elif event.inaxes == self.nsref_axes:
        self.nsref_pick.set_visible(False)
      elif event.inaxes == self.res_axes:
        self.res_pick.set_visible(False)

      self.draw_plots()

    if event.dblclick:
      self.dblclick = True
    else:
      self.dblclick = False

  def on_button_release(self, event):
    if event.button == 1 and self.dblclick:
      self.show_image_group(e=event)
      self.view_proc_images()

class LiveAnalysisTab(d.ScrolledPanel):
  def __init__(self,
               parent,
               init=None,
               gparams=None,
               finished_objects=None):
    self.parent = parent
    self.init = init
    self.gparams = gparams
    self.finished_objects = finished_objects
    self.cluster_info = None
    self.prime_info = None
    self.tb1 = None

    d.ScrolledPanel.__init__(self, parent)
    self.main_fig_sizer = wx.GridBagSizer(0, 0)

    # Set regular font
    plt.rc('font', family='sans-serif', size=plot_font_size)
    plt.rc('mathtext', default='regular')

    # UC Histogram / cluster figure
    self.uc_panel = wx.Panel(self)
    uc_box = wx.StaticBox(self.uc_panel, label='Unit Cell Histograms')
    uc_sizer = wx.StaticBoxSizer(uc_box, wx.VERTICAL)
    self.uc_panel.SetSizer(uc_sizer)
    self.uc_figure = Figure(figsize=(1, 2.5))
    self.uc_figure.patch.set_visible(False)  # create transparent background

    uc_gsub = gridspec.GridSpec(2, 3, wspace=0, hspace=0)
    self.a_axes = self.uc_figure.add_subplot(uc_gsub[0])
    self.a_axes.xaxis.get_major_ticks()[0].label1.set_visible(False)
    self.a_axes.xaxis.get_major_ticks()[-1].label1.set_visible(False)
    self.b_axes = self.uc_figure.add_subplot(uc_gsub[1], sharey=self.a_axes)
    self.b_axes.xaxis.get_major_ticks()[0].label1.set_visible(False)
    self.b_axes.xaxis.get_major_ticks()[-1].label1.set_visible(False)
    plt.setp(self.b_axes.get_yticklabels(), visible=False)
    self.c_axes = self.uc_figure.add_subplot(uc_gsub[2], sharey=self.a_axes)
    self.c_axes.xaxis.get_major_ticks()[0].label1.set_visible(False)
    self.c_axes.xaxis.get_major_ticks()[-1].label1.set_visible(False)
    plt.setp(self.c_axes.get_yticklabels(), visible=False)
    self.alpha_axes = self.uc_figure.add_subplot(uc_gsub[3])
    self.alpha_axes.xaxis.get_major_ticks()[0].label1.set_visible(False)
    self.alpha_axes.xaxis.get_major_ticks()[-1].label1.set_visible(False)
    self.beta_axes = self.uc_figure.add_subplot(uc_gsub[4],
                                                sharey=self.alpha_axes)
    self.beta_axes.xaxis.get_major_ticks()[0].label1.set_visible(False)
    self.beta_axes.xaxis.get_major_ticks()[-1].label1.set_visible(False)
    plt.setp(self.beta_axes.get_yticklabels(), visible=False)
    self.gamma_axes = self.uc_figure.add_subplot(uc_gsub[5],
                                                 sharey=self.alpha_axes)
    self.gamma_axes.xaxis.get_major_ticks()[0].label1.set_visible(False)
    self.gamma_axes.xaxis.get_major_ticks()[-1].label1.set_visible(False)
    plt.setp(self.gamma_axes.get_yticklabels(), visible=False)

    self.uc_figure.set_tight_layout(True)
    self.uc_canvas = FigureCanvas(self.uc_panel, -1, self.uc_figure)
    uc_sizer.Add(self.uc_canvas, 1, flag=wx.EXPAND)

    # UC Clustering Result
    self.cluster_panel = wx.Panel(self)
    cluster_box = wx.StaticBox(self.cluster_panel, label='Unit Cell Clustering')
    cluster_box_sizer = wx.StaticBoxSizer(cluster_box, wx.VERTICAL)
    self.cluster_panel.SetSizer(cluster_box_sizer)
    self.cluster_list = ct.CustomListCtrl(self.cluster_panel, size=(-1, 100))
    self.cluster_list.ctr.InsertColumn(0, "#")
    self.cluster_list.ctr.InsertColumn(1, "Lattice")
    self.cluster_list.ctr.InsertColumn(2, "Unit Cell", width=200)
    self.cluster_list.ctr.setResizeColumn(3)
    cluster_box_sizer.Add(self.cluster_list, proportion=1, flag=wx.EXPAND)

    # PRIME result
    self.tb1_panel = wx.Panel(self)
    tb1_box = wx.StaticBox(self.tb1_panel, label='PRIME Merging Statistics ('
                                                 'no postref)')
    self.tb1_box_sizer = wx.StaticBoxSizer(tb1_box, wx.VERTICAL)
    self.tb1_panel.SetSizer(self.tb1_box_sizer)

    # Analysis Options
    self.opt_panel = wx.Panel(self)
    self.opt_sizer = wx.BoxSizer(wx.VERTICAL)
    self.opt_panel.SetSizer(self.opt_sizer)
    self.pg_uc = ct.OptionCtrl(self.opt_panel,
                               items=[('pg', ''), ('uc', '')],
                               sub_labels=('symmetry', 'cell'),
                               grid=(2, 2),
                               expand_cols=(1),
                               label_size=wx.DefaultSize,
                               ctrl_size=wx.DefaultSize)
    self.opt_sizer.Add(self.pg_uc, flag=wx.EXPAND | wx.ALL, border=5)

    self.btn_run_analysis = wx.Button(self.opt_panel, label='Run Analysis')
    self.opt_sizer.Add(self.btn_run_analysis, wx.TOP, border=10)

    self.main_fig_sizer.Add(self.uc_panel, pos=(0, 0), span=(2, 4),
                            flag=wx.EXPAND)
    self.main_fig_sizer.Add(self.cluster_panel, pos=(2, 0), span=(2, 3),
                            flag=wx.EXPAND)
    self.main_fig_sizer.Add(self.tb1_panel, pos=(2, 3), span=(1, 1),
                            flag=wx.EXPAND)
    self.main_fig_sizer.Add(self.opt_panel, pos=(3, 3), flag=wx.EXPAND)

    self.main_fig_sizer.AddGrowableCol(0)
    self.main_fig_sizer.AddGrowableCol(1)
    self.main_fig_sizer.AddGrowableCol(2)
    self.main_fig_sizer.AddGrowableRow(1)
    self.main_fig_sizer.AddGrowableRow(2)
    self.main_fig_sizer.AddGrowableRow(3)

    self.SetSizer(self.main_fig_sizer)

  def draw_plots(self):
    if self.cluster_info is not None:
      self.report_clustering_results()
    if self.prime_info is not None:
      self.report_prime_results()
    if self.finished_objects:
      self.draw_uc_histograms()

    self.Layout()
    self.SetupScrolling()

  def report_clustering_results(self):
    self.cluster_list.ctr.DeleteAllItems()
    clusters = sorted(self.cluster_info, key=lambda i:i['number'], reverse=True)
    for c in clusters:
      i = clusters.index(c)
      idx = self.cluster_list.ctr.InsertStringItem(i, str(c['number']))
      self.cluster_list.ctr.SetStringItem(idx, 1, str(c['pg']))
      self.cluster_list.ctr.SetStringItem(idx, 2, str(c['uc']))

  def report_prime_results(self):
    # Remove previous data if exists
    if self.tb1 is not None:
      self.tb1.Destroy()

    try:
      self.plot = ppl.Plotter(info=self.prime_info)
      self.tb1_labels, self.tb1_data = self.plot.table_one()
      self.tb1 = ct.TableCtrl(self.tb1_panel,
                              rlabels=self.tb1_labels,
                              contents=self.tb1_data,
                              label_style='bold')
      self.tb1_box_sizer.Add(self.tb1, 1, flag=wx.EXPAND | wx.ALL, border=10)
    except Exception as e:
      print ('PRIME PLOTTER ERROR: ', e)

  def calculate_uc_histogram(self, a, axes, xticks_loc='top', set_ylim=False):
    n, bins = np.histogram(a, 50)
    left = np.array(bins[:-1])
    right = np.array(bins[1:])
    bottom = np.zeros(len(left))
    top = bottom + n
    XY = np.array(
      [[left, left, right, right], [bottom, top, top, bottom]]).T
    barpath = path.Path.make_compound_path_from_polys(XY)
    patch = mpatches.PathPatch(barpath, fc='#4575b4', lw=0, alpha=0.75)
    axes.add_patch(patch)

    axes.set_xlim(left[0], right[-1])
    if set_ylim:
      axes.set_ylim(bottom.min(), 1.05 * top.max())

    # axes.hist(a, 50, normed=False, facecolor='#4575b4',
    #                  histtype='stepfilled')
    axes.xaxis.get_major_ticks()[0].label1.set_visible(False)
    axes.xaxis.get_major_ticks()[-1].label1.set_visible(False)

    if xticks_loc == 'top':
      axes.xaxis.tick_top()
    elif xticks_loc == 'bottom':
      axes.xaxis.tick_bottom()


  def draw_uc_histograms(self):
    try:
      # Unit cell histograms
      self.a_axes.clear()
      self.b_axes.clear()
      self.c_axes.clear()
      self.alpha_axes.clear()
      self.beta_axes.clear()
      self.gamma_axes.clear()

      a = [i.final['a'] for i in self.finished_objects]
      b = [i.final['b'] for i in self.finished_objects]
      c = [i.final['c'] for i in self.finished_objects]
      alpha = [i.final['alpha'] for i in self.finished_objects]
      beta = [i.final['beta'] for i in self.finished_objects]
      gamma = [i.final['gamma'] for i in self.finished_objects]

      self.calculate_uc_histogram(a, self.a_axes, set_ylim=True)
      edge_ylabel = 'a, b, c ({})'.format(r'$\AA$')
      self.a_axes.set_ylabel(edge_ylabel)

      self.calculate_uc_histogram(b, self.b_axes)
      plt.setp(self.b_axes.get_yticklabels(), visible=False)

      self.calculate_uc_histogram(c, self.c_axes)
      plt.setp(self.c_axes.get_yticklabels(), visible=False)

      self.calculate_uc_histogram(alpha, self.alpha_axes,
                                  xticks_loc='bottom', set_ylim=True)
      ang_ylabel = '{}, {}, {} ({})'.format(r'$\alpha$', r'$\beta$',
                                            r'$\gamma$', r'$^\circ$')
      self.alpha_axes.set_ylabel(ang_ylabel)

      self.calculate_uc_histogram(beta, self.beta_axes, xticks_loc='bottom')
      plt.setp(self.beta_axes.get_yticklabels(), visible=False)

      self.calculate_uc_histogram(gamma, self.gamma_axes, xticks_loc='bottom')
      plt.setp(self.gamma_axes.get_yticklabels(), visible=False)

    except ValueError as e:
      print ('UC HISTOGRAM ERROR: ', e)


class SummaryTab(d.ScrolledPanel):
  def __init__(self,
               parent,
               init=None,
               gparams=None,
               final_objects=None,
               out_dir=None,
               plot=None):
    d.ScrolledPanel.__init__(self, parent)

    self.parent = parent
    self.final_objects = final_objects
    self.gparams = gparams
    self.out_dir = out_dir
    self.plot = plot
    self.init = init

    summary_sizer = wx.BoxSizer(wx.VERTICAL)

    sfont = wx.Font(norm_font_size, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
    bfont = wx.Font(norm_font_size, wx.DEFAULT, wx.NORMAL, wx.BOLD)
    self.SetFont(bfont)

    # Run information
    run_box = wx.StaticBox(self, label='Run Information')
    run_box.SetFont(sfont)
    run_box_sizer = wx.StaticBoxSizer(run_box, wx.VERTICAL)
    run_box_grid = wx.FlexGridSizer(3, 2, 5, 20)
    self.title_txt = wx.StaticText(self, label='')
    self.title_txt.SetFont(sfont)
    self.folder_txt = wx.StaticText(self, label='')
    self.folder_txt.SetFont(sfont)

    run_box_grid.AddMany([(wx.StaticText(self, label='Title')),
                          (self.title_txt, 1, wx.EXPAND),
                          (wx.StaticText(self, label='Directory')),
                          (self.folder_txt, 1, wx.EXPAND)])

    run_box_grid.AddGrowableCol(1, 1)
    run_box_sizer.Add(run_box_grid, flag=wx.EXPAND | wx.ALL, border=10)

    summary_sizer.Add(run_box_sizer, flag=wx.EXPAND | wx.ALL, border=10)

    # Integration summary
    if self.gparams.advanced.processing_backend == 'ha14':
      int_box = wx.StaticBox(self, label='Analysis of Integration')
      int_box.SetFont(sfont)
      int_box_sizer = wx.StaticBoxSizer(int_box, wx.HORIZONTAL)
      int_box_grid = wx.FlexGridSizer(4, 5, 5, 20)
      int_btn_sizer = wx.BoxSizer(wx.VERTICAL)

      # Grid search summary
      self.sih_min = wx.StaticText(self, label='4.0')
      self.sih_min.SetFont(sfont)
      self.sih_max = wx.StaticText(self, label='8.0')
      self.sih_max.SetFont(sfont)
      self.sih_avg = wx.StaticText(self, label='6.0')
      self.sih_avg.SetFont(sfont)
      self.sih_std = wx.StaticText(self, label='0.05')
      self.sih_std.SetFont(sfont)
      self.sph_min = wx.StaticText(self, label='5.0')
      self.sph_min.SetFont(sfont)
      self.sph_max = wx.StaticText(self, label='12.0')
      self.sph_max.SetFont(sfont)
      self.sph_avg = wx.StaticText(self, label='8.5')
      self.sph_avg.SetFont(sfont)
      self.sph_std = wx.StaticText(self, label='0.15')
      self.sph_std.SetFont(sfont)
      self.spa_min = wx.StaticText(self, label='10.0')
      self.spa_min.SetFont(sfont)
      self.spa_max = wx.StaticText(self, label='20.0')
      self.spa_max.SetFont(sfont)
      self.spa_avg = wx.StaticText(self, label='15.0')
      self.spa_avg.SetFont(sfont)
      self.spa_std = wx.StaticText(self, label='0.01')
      self.spa_std.SetFont(sfont)

      int_box_grid.AddMany([(wx.StaticText(self, label='')),
                            (wx.StaticText(self, label='min')),
                            (wx.StaticText(self, label='max')),
                            (wx.StaticText(self, label='avg')),
                            (wx.StaticText(self, label='std')),
                            (
                            wx.StaticText(self, label='minimum signal height')),
                            (self.sih_min), (self.sih_max),
                            (self.sih_avg), (self.sih_std),
                            (wx.StaticText(self, label='minimum spot height')),
                            (self.sph_min), (self.sph_max),
                            (self.sph_avg), (self.sph_std),
                            (wx.StaticText(self, label='minimum spot area')),
                            (self.spa_min), (self.spa_max),
                            (self.spa_avg), (self.spa_std)
                            ])

      # Button & binding for heatmap display
      heatmap_bmp = bitmaps.fetch_custom_icon_bitmap('heatmap24')
      self.int_heatmap = ct.GradButton(self,
                                       bmp = heatmap_bmp,
                                       label='  Spotfinding Heatmap',
                                       size=(250, -1))
      int_btn_sizer.Add(self.int_heatmap)
      self.Bind(wx.EVT_BUTTON, self.onPlotHeatmap, self.int_heatmap)

      # Insert into sizers
      int_box_sizer.Add(int_box_grid, flag=wx.ALL, border=10)
      int_box_sizer.AddStretchSpacer()
      int_box_sizer.Add(int_btn_sizer, flag=wx.ALL, border=10)
      summary_sizer.Add(int_box_sizer, flag=wx.EXPAND | wx.ALL, border=10)


    # Dataset Info
    dat_box = wx.StaticBox(self, label='Dataset Information')
    dat_box.SetFont(sfont)
    dat_box_sizer = wx.StaticBoxSizer(dat_box, wx.HORIZONTAL)
    dat_box_grid = wx.FlexGridSizer(4, 2, 5, 20)
    dat_btn_sizer = wx.BoxSizer(wx.VERTICAL)

    self.pg_txt = wx.StaticText(self, label='P4')
    self.pg_txt.SetFont(sfont)
    self.uc_txt = wx.StaticText(self, label='79 79 38 90 90 90')
    self.uc_txt.SetFont(sfont)
    res = to_unicode(u'50.0 - 1.53 {}'.format(u'\u212B'))
    self.rs_txt = wx.StaticText(self, label=res)
    self.rs_txt.SetFont(sfont)
    self.xy_txt = wx.StaticText(self, label='X = 224.90 mm, Y = 225.08 mm')
    self.xy_txt.SetFont(sfont)

    dat_box_grid.AddMany([(wx.StaticText(self, label='Bravais lattice: ')),
                          (self.pg_txt),
                          (wx.StaticText(self, label='Unit cell: ')),
                          (self.uc_txt),
                          (wx.StaticText(self, label='Resolution: ')),
                          (self.rs_txt),
                          (wx.StaticText(self, label='Beam XY: ')),
                          (self.xy_txt)
                          ])

    # Buttons for res. histogram and beam xy plot
    hist_bmp = bitmaps.fetch_icon_bitmap('mimetypes', 'spreadsheet', size=32,
                                         scale=(24, 24))
    self.dat_reshist = ct.GradButton(self,
                                     bmp=hist_bmp,
                                     label='  Resolution Histogram', size=(250, -1))
    beamXY_bmp = bitmaps.fetch_custom_icon_bitmap('scatter_plot_24')
    self.dat_beamxy = ct.GradButton(self,
                                    bmp=beamXY_bmp,
                                    label='  Beam XY Plot', size=(250, -1))
    self.dat_beam3D = ct.GradButton(self,
                                    bmp=beamXY_bmp,
                                    label='  Beam XYZ Plot', size=(250, -1))
    dat_btn_sizer.Add(self.dat_reshist)
    dat_btn_sizer.Add(self.dat_beamxy, flag=wx.TOP, border=5)
    dat_btn_sizer.Add(self.dat_beam3D, flag=wx.TOP, border=5)
    self.Bind(wx.EVT_BUTTON, self.onPlotBeamXY, self.dat_beamxy)
    self.Bind(wx.EVT_BUTTON, self.onPlotBeam3D, self.dat_beam3D)
    self.Bind(wx.EVT_BUTTON, self.onPlotResHist, self.dat_reshist)

    # Insert into sizers
    dat_box_sizer.Add(dat_box_grid, flag=wx.ALL, border=10)
    dat_box_sizer.AddStretchSpacer()
    dat_box_sizer.Add(dat_btn_sizer, flag=wx.ALL, border=10)
    summary_sizer.Add(dat_box_sizer, flag=wx.EXPAND | wx.ALL, border=10)

    # Clustering info
    self.cluster_panel = wx.Panel(self)
    cluster_box = wx.StaticBox(self.cluster_panel, label='Unit Cell Clustering')
    cluster_box_sizer = wx.StaticBoxSizer(cluster_box, wx.VERTICAL)
    self.cluster_panel.SetSizer(cluster_box_sizer)
    self.cluster_info = ct.CustomListCtrl(self.cluster_panel, size=(-1, 100))
    self.cluster_info.ctr.InsertColumn(0, "#")
    self.cluster_info.ctr.InsertColumn(1, "Lattice")
    self.cluster_info.ctr.InsertColumn(2, "Unit Cell")
    self.cluster_info.ctr.InsertColumn(3, "Filename")
    self.cluster_info.ctr.setResizeColumn(3)
    cluster_box_sizer.Add(self.cluster_info, proportion=1, flag=wx.EXPAND)
    summary_sizer.Add(self.cluster_panel, flag=wx.EXPAND | wx.ALL, border=10)

    # Hide if not done
    if not self.gparams.analysis.run_clustering:
      self.cluster_panel.Hide()

    # Summary
    smr_box = wx.StaticBox(self, label='Run Summary')
    smr_box.SetFont(sfont)
    smr_box_sizer = wx.StaticBoxSizer(smr_box, wx.HORIZONTAL)
    smr_box_grid = wx.FlexGridSizer(8, 2, 5, 20)
    smr_btn_sizer = wx.BoxSizer(wx.VERTICAL)

    self.readin_txt = wx.StaticText(self, label='250')
    self.readin_txt.SetFont(sfont)
    self.nodiff_txt = wx.StaticText(self, label='100')
    self.nodiff_txt.SetFont(sfont)
    self.w_diff_txt = wx.StaticText(self, label='150')
    self.w_diff_txt.SetFont(sfont)
    self.noint_txt = wx.StaticText(self, label='30')
    self.noint_txt.SetFont(sfont)
    self.final_txt = wx.StaticText(self, label='100')
    self.final_txt.SetFont(sfont)

    smr_box_grid.AddMany([(wx.StaticText(self, label='Read in: ')),
                          (self.readin_txt),
                          (wx.StaticText(self, label='No diffraction:')),
                          (self.nodiff_txt),
                          (wx.StaticText(self, label='Have diffraction: ')),
                          (self.w_diff_txt)])

    prime_bmp = bitmaps.fetch_custom_icon_bitmap('prime32', scale=(24, 24))
    self.smr_runprime = ct.GradButton(self,
                                      bmp=prime_bmp,
                                      label='  Run PRIME', size=(250, -1))
    cluster_bmp = bitmaps.fetch_custom_icon_bitmap('distance_difference',
                                                   scale=(24, 24))
    self.smr_runcluster = ct.GradButton(self,
                                      bmp=cluster_bmp,
                                      label='  Run CLUSTER', size=(250, -1))

    smr_btn_sizer.Add(self.smr_runcluster)
    smr_btn_sizer.Add(self.smr_runprime, flag=wx.TOP, border=5)
    self.Bind(wx.EVT_BUTTON, self.onPRIME, self.smr_runprime)
    self.Bind(wx.EVT_BUTTON, self.onCLUSTER, self.smr_runcluster)

    if self.gparams.advanced.processing_backend == 'ha14':
      self.noprf_txt = wx.StaticText(self, label='20')
      self.noprf_txt.SetFont(sfont)
      smr_box_grid.AddMany([(wx.StaticText(self,
                                           label='Failed indexing / integration')),
                              (self.noint_txt),
                              (wx.StaticText(self, label='Failed filter')),
                              (self.noprf_txt)])
    elif self.gparams.advanced.processing_backend == 'cctbx.xfel':
      self.nospf_txt = wx.StaticText(self, label='10')
      self.nospf_txt.SetFont(sfont)
      self.noidx_txt = wx.StaticText(self, label='20')
      self.noidx_txt.SetFont(sfont)
      self.noflt_txt = wx.StaticText(self, label='20')
      self.noflt_txt.SetFont(sfont)
      smr_box_grid.AddMany([(wx.StaticText(self,
                                           label='Failed spotfinding')),
                            (self.nospf_txt),
                            (wx.StaticText(self, label='Failed indexing')),
                            (self.noidx_txt),
                            (wx.StaticText(self, label='Failed integration')),
                            (self.noint_txt),
                            (wx.StaticText(self, label='Failed filter')),
                            (self.noflt_txt)])
    smr_box_grid.AddMany([(wx.StaticText(self,
                                         label='Final integrated pickles')),
                          (self.final_txt)])


    smr_box_sizer.Add(smr_box_grid, flag=wx.ALL, border=10)
    smr_box_sizer.AddStretchSpacer()
    smr_box_sizer.Add(smr_btn_sizer, flag=wx.ALL,  border=10)
    summary_sizer.Add(smr_box_sizer, flag=wx.EXPAND | wx.ALL, border=10)

    self.SetFont(sfont)
    self.SetSizer(summary_sizer)
    self.SetupScrolling()

  def onPRIME(self, e):
    from prime.postrefine.mod_gui_init import PRIMEWindow
    self.prime_window = PRIMEWindow(self, -1, title='PRIME',
                                    prefix=self.gparams.advanced.prime_prefix)
    self.prime_window.load_script(out_dir=self.out_dir)
    self.prime_window.SetMinSize(self.prime_window.GetEffectiveMinSize())
    self.prime_window.Show(True)

  def onCLUSTER(self, e):
    cluster_dlg = d.ClusterDialog(self)
    cluster_dlg.write_files.SetValue(self.gparams.analysis.cluster_write_files)
    cluster_dlg.cluster_threshold.ctr.SetValue(self.gparams.analysis.cluster_threshold)
    cluster_dlg.cluster_limit.ctr.SetValue(self.gparams.analysis.cluster_limit)
    if self.gparams.analysis.cluster_n_images > 0:
      cluster_dlg.cluster_n_images.ctr.SetValue(self.gparams.analysis.cluster_n_images)

    if (cluster_dlg.ShowModal() == wx.ID_OK):
      self.cluster_panel.Show()
      self.Layout()
      self.gparams.analysis.run_clustering = True
      self.gparams.analysis.cluster_write_files = cluster_dlg.write_files.GetValue()
      self.gparams.analysis.cluster_threshold = cluster_dlg.cluster_threshold.ctr.GetValue()
      self.gparams.analysis.cluster_limit = cluster_dlg.cluster_limit.ctr.GetValue()
      if cluster_dlg.cluster_n_images.toggle.GetValue():
        self.gparams.analysis.cluster_n_images = int(cluster_dlg.cluster_n_images.ctr.GetValue())
      else:
        self.gparams.analysis.cluster_n_images = 0
      self.init.params = self.gparams

      analysis = Analyzer(init=self.init,
                          all_objects=self.final_objects,
                          gui_mode=True)
      clusters = analysis.unit_cell_analysis()
      if clusters:
        self.report_clustering_results(clusters=clusters)

  def report_clustering_results(self, clusters):
    self.cluster_info.ctr.DeleteAllItems()
    clusters = sorted(clusters, key=lambda i:i['number'], reverse=True)
    for c in clusters:
      i = clusters.index(c)
      idx = self.cluster_info.ctr.InsertStringItem(i, str(c['number']))
      self.cluster_info.ctr.SetStringItem(idx, 1, str(c['pg']))
      self.cluster_info.ctr.SetStringItem(idx, 2, str(c['uc']))
      if c['filename'] not in ('*', None):
        self.cluster_info.ctr.SetStringItem(idx, 3, c['filename'])

    self.Refresh()
    self.SetupScrolling()

  def onPlotHeatmap(self, e):
    if self.final_objects is not None:
      self.plot.plot_spotfinding_heatmap()

  def onPlotBeamXY(self, e):
    if self.final_objects is not None:
      self.plot.plot_beam_xy()

  def onPlotBeam3D(self, e):
    if self.final_objects is not None:
      self.plot.plot_beam_xy(threeD=True)

  def onPlotResHist(self, e):
    if self.final_objects is not None:
      self.plot.plot_res_histogram()

class ProcWindow(wx.Frame):
  """ New frame that will show processing info """

  def __init__(self, parent, id, title, phil, target_phil=None):
    wx.Frame.__init__(self, parent, id, title,
                      size=(800, 900),
                      style= wx.SYSTEM_MENU | wx.CAPTION |
                             wx.CLOSE_BOX | wx.RESIZE_BORDER)

    self.parent = parent
    self.logtext = ''
    self.obj_counter = 0
    self.bookmark = 0
    self.iota_phil = phil
    self.gparams = phil.extract()

    self.target_phil = target_phil
    self.state = 'process'
    self.recovery = False

    self.abort_initiated = False
    self.run_aborted = False
    self.monitor_mode = False
    self.monitor_mode_timeout = None
    self.timeout_start = None
    self.find_new_images = self.monitor_mode
    self.start_object_finder = True
    self.plotter_time = []
    self.obj_reader_time = []

    self.running_cluster = False
    self.running_prime = False
    self.running_manual_analysis = False
    self.job_id = None

    self.finished_objects = []
    self.read_object_files = []
    self.new_images = []
    self.indices = []
    self.b_factors = []
    self.cluster_info = None
    self.pparams = None
    self.prime_info = None

    # Initialize mx handler
    self.mxh = mx_handler()

    self.main_panel = wx.Panel(self)
    self.main_sizer = wx.BoxSizer(wx.VERTICAL)

    # Toolbar
    self.proc_toolbar = self.CreateToolBar(style=wx.TB_3DBUTTONS | wx.TB_TEXT)
    abort_bmp = bitmaps.fetch_icon_bitmap('actions', 'stop')
    self.tb_btn_abort = self.proc_toolbar.AddLabelTool(wx.ID_ANY, label='Abort',
                                                bitmap=abort_bmp,
                                                shortHelp='Abort')
    resume_bmp = bitmaps.fetch_icon_bitmap('actions', 'quick_restart')
    self.tb_btn_resume = self.proc_toolbar.AddLabelTool(wx.ID_ANY,
                                                        label='Resume',
                                                        bitmap=resume_bmp,
                                                        shortHelp='Resume aborted run')
    self.proc_toolbar.EnableTool(self.tb_btn_resume.GetId(), False)
    self.proc_toolbar.AddSeparator()
    watch_bmp = bitmaps.fetch_icon_bitmap('apps', 'search')
    self.tb_btn_monitor = self.proc_toolbar.AddCheckLabelTool(wx.ID_ANY,
                                                label='Monitor',
                                                bitmap=watch_bmp,
                                                shortHelp='Monitor Mode')
    hist_bmp = bitmaps.fetch_icon_bitmap('mimetypes', 'spreadsheet', size=32)
    self.tb_btn_analysis = self.proc_toolbar.AddCheckLabelTool(wx.ID_ANY,
                                              label='Analysis',
                                              bitmap=hist_bmp,
                                              shortHelp='Toggle Runtime Analysis Tab')
    self.proc_toolbar.Realize()

    # Status box
    self.status_panel = wx.Panel(self.main_panel)
    self.status_sizer = wx.BoxSizer(wx.VERTICAL)
    self.status_box = wx.StaticBox(self.status_panel, label='Status')
    self.status_box_sizer = wx.StaticBoxSizer(self.status_box, wx.HORIZONTAL)
    self.status_txt = wx.StaticText(self.status_panel, label='')
    self.status_box_sizer.Add(self.status_txt, flag=wx.ALL | wx.ALIGN_CENTER,
                              border=10)
    self.status_sizer.Add(self.status_box_sizer,
                          flag=wx.EXPAND | wx.ALL, border=3)
    self.status_panel.SetSizer(self.status_sizer)

    # Tabbed output window(s)
    self.proc_panel = wx.Panel(self.main_panel)
    # self.proc_nb = wx.Notebook(self.proc_panel, style=0)
    self.proc_nb = AuiNotebook(self.proc_panel, style=wx.aui.AUI_NB_TOP)
    self.proc_tab = ProcessingTab(self.proc_nb)
    self.log_tab = LogTab(self.proc_nb)
    self.chart_tab = LiveAnalysisTab(self.proc_nb)
    self.proc_nb.AddPage(self.log_tab, 'Log')
    self.proc_nb.AddPage(self.proc_tab, 'Processing')
    self.proc_nb.SetSelection(1)
    self.proc_sizer = wx.BoxSizer(wx.VERTICAL)
    self.proc_sizer.Add(self.proc_nb, 1, flag=wx.EXPAND | wx.ALL, border=3)
    self.proc_panel.SetSizer(self.proc_sizer)

    self.main_sizer.Add(self.status_panel, flag=wx.EXPAND | wx.ALL, border=3)
    self.main_sizer.Add(self.proc_panel, 1, flag=wx.EXPAND | wx.ALL, border=3)
    self.main_panel.SetSizer(self.main_sizer)

    #Processing status bar
    self.sb = self.CreateStatusBar()
    self.sb.SetFieldsCount(2)
    self.sb.SetStatusWidths([-1, -2])

    # Output gauge in status bar
    self.gauge_process = wx.Gauge(self.sb, -1, style=wx.GA_HORIZONTAL | wx.GA_SMOOTH)
    rect = self.sb.GetFieldRect(0)
    self.gauge_process.SetPosition((rect.x + 2, rect.y + 2))
    self.gauge_process.SetSize((rect.width - 4, rect.height - 4))
    self.gauge_process.Hide()

    # Output polling timer
    self.proc_timer = wx.Timer(self)
    self.chart_timer = wx.Timer(self)

    # PostEvent bindings
    self.Bind(thr.EVT_ALLDONE, self.onFinishedProcess)
    self.Bind(thr.EVT_IMGDONE, self.onFinishedImageFinder)
    self.Bind(thr.EVT_OBJDONE, self.onFinishedObjectReader)
    self.Bind(thr.EVT_CLUSTERDONE, self.onFinishedCluster)
    self.Bind(thr.EVT_PRIMEDONE, self.onFinishedPRIME)

    # Event bindings
    self.sb.Bind(wx.EVT_SIZE, self.onStatusBarResize)
    self.Bind(wx.EVT_TIMER, self.onProcTimer, id=self.proc_timer.GetId())
    self.Bind(wx.EVT_TIMER, self.onChartTimer, id=self.chart_timer.GetId())
    self.Bind(wx.aui.EVT_AUINOTEBOOK_PAGE_CHANGED, self.onTabChange,
              self.proc_nb)

    # Button bindings
    self.Bind(wx.EVT_TOOL, self.onAbort, self.tb_btn_abort)
    self.Bind(wx.EVT_TOOL, self.onResume, self.tb_btn_resume)
    self.Bind(wx.EVT_TOOL, self.onMonitor, self.tb_btn_monitor)
    self.Bind(wx.EVT_TOOL, self.onAnalysis, self.tb_btn_analysis)
    self.Bind(wx.EVT_BUTTON, self.onAnalysisManualRun,
              self.chart_tab.btn_run_analysis)

    # Proc Tab bindings
    self.Bind(wx.EVT_TEXT_ENTER, self.onSGTextEnter, self.proc_tab.hkl_sg.sg)
    self.Bind(wx.EVT_CHECKBOX, self.onSGCheckbox, self.proc_tab.hkl_sg.toggle)

    # Determine if monitor mode was previously selected
    if self.gparams.gui.monitor_mode:
      self.proc_toolbar.ToggleTool(self.tb_btn_monitor.GetId(), True)
      self.monitor_mode = True
      if self.gparams.gui.monitor_mode_timeout:
        if self.gparams.gui.monitor_mode_timeout_length is None:
          self.monitor_mode_timeout = 30
        else:
          self.monitor_mode_timeout = self.gparams.gui.monitor_mode_timeout_length

  def set_position(self):
    """ Determines screen position w/ respect to parent window; will also
    detect if it goes beyond the display edge, and adjust """

    # Position proc window w/ respect to IOTA window
    mx, my = self.parent.GetPosition()
    px = mx + 50
    py = my + 50

    # Calculate if proc window is going out of bounds, and adjust
    disp_idx = wx.Display.GetFromWindow(window=self.parent)
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
    self.SetPosition((px, py))

  def onSGTextEnter(self, e):
    self.info.user_sg = str(self.proc_tab.hkl_sg.sg.GetValue())

  def onSGCheckbox(self, e):
    # TODO: find a way to re-draw just this plot on the fly (maybe put all
    # TODO: plot generation in yet another thread?)
    if e.IsChecked():
      if self.gparams.cctbx_xfel.target_space_group:
        self.proc_tab.hkl_sg.sg.SetValue(str(self.gparams.cctbx_xfel.target_space_group))
      else:
        self.proc_tab.hkl_sg.sg.SetValue('P1')
      self.info.user_sg = str(self.proc_tab.hkl_sg.sg.GetValue())
    else:
      self.info.user_sg = 'P1'
      self.proc_tab.hkl_sg.sg.SetValue('')

  def onTabChange(self, e):
  #   """ Only update displays if user changes to the tab """
  #   tab = self.proc_nb.GetSelection()
  #   if tab == 2:
  #     self.chart_timer.Start(15000)
  #   else:
  #     self.chart_timer.Stop()
    pass

  def onAnalysisManualRun(self, e):
    # Run analysis calculations when Run Analysis button is pressed; this
    # will enable analysis run even when IOTA isn't running anymore
    if not (self.running_cluster or self.running_prime):
      self.running_manual_analysis = True
      self.chart_tab.btn_run_analysis.Disable()
      self.run_clustering_thread()

  def onAnalysis(self, e):
    if self.proc_toolbar.GetToolState(self.tb_btn_analysis.GetId()):
      if hasattr(self, 'proc_thread'):
        self.chart_timer.Start(15000)
      self.proc_nb.InsertPage(n=2, page=self.chart_tab, text='Analysis',
                              select=True)
      if self.proc_nb.GetSelection() != 2:
        self.proc_nb.SetSelection(2)
      self.plot_live_analysis(force_plot=True)
    else:
      if self.chart_timer.IsRunning():
        self.chart_timer.Stop()
      if self.proc_nb.GetSelection() == 2:
        self.proc_nb.SetSelection(1)
      self.proc_nb.RemovePage(2)

  def onMonitor(self, e):
    if self.proc_toolbar.GetToolState(self.tb_btn_monitor.GetId()):
      self.monitor_mode = True
      if self.gparams.gui.monitor_mode_timeout:
        if self.gparams.gui.monitor_mode_timeout_length is None:
          self.monitor_mode_timeout = 30
        else:
          self.monitor_mode_timeout = self.gparams.gui.monitor_mode_timeout_length
    elif not self.proc_toolbar.GetToolState(self.tb_btn_monitor.GetId()):
      self.monitor_mode = False
      self.monitor_mode_timeout = None
    self.find_new_images = self.monitor_mode

  def onStatusBarResize(self, e):
    rect = self.sb.GetFieldRect(0)
    self.gauge_process.SetPosition((rect.x + 2, rect.y + 2))
    self.gauge_process.SetSize((rect.width - 4, rect.height - 4))
    self.Refresh()
    self.Layout()

  def onAbort(self, e):
    if self.gparams.mp.method == 'lsf':
      kill_command = 'bkill -J {}'.format(self.job_id)
      easy_run.fully_buffered(kill_command)
    with open(self.tmp_abort_file, 'w') as af:
      af.write('')
    self.abort_initiated = True
    self.status_txt.SetForegroundColour('red')
    self.status_txt.SetLabel('Aborting...')
    self.proc_toolbar.EnableTool(self.tb_btn_abort.GetId(), False)

  def onResume(self, e):
    """ Restarts an aborted run if the processing window is still open.
    Basically goes through self.finished_objects, extracts the raw image
    names and regenerates the self.img_list to only have those image paths;
    then finds any 'new' images (which includes unprocessed images as well as
    any images that may have been added during the abort pause) and runs
    processing """

    # Remove abort signal file(s)
    if os.path.isfile(self.tmp_abort_file):
      os.remove(self.tmp_abort_file)
    if os.path.isfile(self.tmp_aborted_file):
        os.remove(self.tmp_aborted_file)
    self.abort_initiated = False
    self.run_aborted = False

    # Reset tmp folder
    if not os.path.isdir(self.init.tmp_base):
      os.makedirs(self.init.tmp_base)

    # Reset toolbar buttons
    self.proc_toolbar.EnableTool(self.tb_btn_abort.GetId(), True)
    self.proc_toolbar.EnableTool(self.tb_btn_resume.GetId(), False)
    self.proc_toolbar.EnableTool(self.tb_btn_monitor.GetId(), True)

    # Run processing, etc.
    self.proc_timer.Start(1000)
    self.chart_timer.Start(15000)
    self.state = 'resume'
    self.process_images()

  def recover(self, int_path, status, init, params):
    self.recovery = True
    self.init = init
    self.gparams = params
    self.tmp_abort_file = os.path.join(int_path, '.abort.tmp')
    self.tmp_aborted_file = os.path.join(int_path, '.aborted.tmp')

    self.info = ep.load(self.init.info_file)

    self.status_txt.SetLabel('Searching in {} ...'.format(int_path))
    self.state = status
    self.finish_process()

  def run(self, init):
    # Initialize IOTA parameters and log
    self.init = init
    self.init.params = self.gparams
    self.init.target_phil = self.target_phil
    good_init, msg = self.init.run()

    # Start process
    if good_init:
      self.tmp_abort_file = os.path.join(self.init.int_base, '.abort.tmp')
      self.tmp_aborted_file = os.path.join(self.init.int_base, '.aborted.tmp')
      self.status_txt.SetForegroundColour('black')
      self.status_txt.SetLabel('Running...')

      self.info = self.init.info
      self.process_images()
      self.good_to_go = True
      self.proc_timer.Start(3000)
    else:
      print (msg)
      self.good_to_go = False

  def extend_image_lists(self, img_list=None, img_only=False):
    ''' Extends image or input lists, as well as info lists, w/ new entries

    :param img_list: list of newly-found images
    :param img_only: do not extend info lists
    '''

    if img_list:
      self.init.input_list.extend(img_list)
      self.info.img_list.extend(img_list)
      if not img_only:
        self.info.nref_list.extend([0] * len(img_list))
        self.info.nref_xaxis.extend([i[0] for i in img_list])
        self.info.res_list.extend([0] * len(img_list))
    else:
      self.info.img_list = [[i, len(self.init.input_list) + 1, j] for
                            i, j in enumerate(self.init.input_list, 1)]
      self.info.nref_list = [0] * len(self.info.img_list)
      self.info.nref_xaxis = [i[0] for i in self.info.img_list]
      self.info.res_list = [0] * len(self.info.img_list)


  def set_proc_iterable(self):
    ''' Sets the following iterables for image processing

    :return: IOTA image iterable
    '''

    if self.state == 'new images':
      self.extend_image_lists(self.new_images)
      msg = 'Processing additional {} images ({} total)...'\
            ''.format(len(self.new_images), len(self.info.img_list))
      return self.new_images, msg

    elif self.state == 'resume':
      subset = ( self.gparams.advanced.image_range.flag_on or
                 self.gparams.advanced.random_sample.flag_on )
      if self.new_images and not subset:
        self.extend_image_lists(self.new_images)
      done_img_list = [i.raw_img for i in self.info.finished_objects]
      iterable = [i for i in self.info.img_list if i[2] not in done_img_list]

      if iterable:
        msg = 'Processing {} remaining images ({} total)...' \
              ''.format(len(iterable), len(self.info.img_list))
      else:
        msg = 'All {} images processed!'.format(len(self.info.img_list))
        self.state = 'finished'
      return iterable, msg

    else:  # New run
      self.extend_image_lists()
      msg = 'Processing {} images...'.format(len(self.info.img_list))
      return self.info.img_list, msg


  def process_images(self):
    """ One-fell-swoop importing / triaging / integration of images """

    # Generate iterable as per run status
    iterable, msg = self.set_proc_iterable()
    self.info.unprocessed_images = len(iterable)

    if self.state == 'finished':
        self.finish_process()

    # Set font properties for status window
    # font = self.sb.GetFont()
    # font.SetWeight(wx.NORMAL)
    # self.status_txt.SetFont(font)
    self.status_txt.SetForegroundColour('black')
    self.status_txt.SetLabel(msg)

    # TODO: Find a way to manage the iterable correctly
    iter_path = os.path.join(self.init.int_base, 'iter.cfg')
    ep.dump(iter_path, iterable)
    init_path = os.path.join(self.init.int_base, 'init.cfg')
    ep.dump(init_path, self.init)

    self.proc_thread = thr.JobSubmitThread(self, params=self.gparams,
                                           output_folder=self.init.int_base,
                                           abort_file=self.tmp_abort_file)
    self.proc_thread.start()


  def analyze_results(self, analysis=None):
    if len(self.final_objects) == 0:
      self.status_txt.SetForegroundColour('red')
      self.status_txt.SetLabel('No images successfully integrated')

    else:
      self.status_txt.SetForegroundColour('black')
      self.status_txt.SetLabel('Analyzing results...')

      # Do analysis
      if analysis is None:
        self.recovery = False
        analysis = Analyzer(self.init,
                            self.info.finished_objects,
                            gui_mode=True)
      plot = Plotter(self.gparams,
                     self.final_objects,
                     self.init.viz_base)

      # Initialize summary tab
      prime_file = os.path.join(self.init.int_base,
                                '{}.phil'.format(self.gparams.advanced.prime_prefix))
      self.summary_tab = SummaryTab(self.proc_nb,
                                    init=self.init,
                                    gparams=self.gparams,
                                    final_objects=self.final_objects,
                                    out_dir=os.path.dirname(prime_file),
                                    plot=plot)

      # Run information
      self.summary_tab.title_txt.SetLabel(ut.noneset(self.gparams.description))
      self.summary_tab.folder_txt.SetLabel(self.gparams.output)

      # Analysis of integration
      if self.gparams.advanced.processing_backend == 'ha14':
        self.summary_tab.sih_min.SetLabel("{:4.0f}".format(np.min(analysis.s)))
        self.summary_tab.sih_max.SetLabel("{:4.0f}".format(np.max(analysis.s)))
        self.summary_tab.sih_avg.SetLabel("{:4.2f}".format(np.mean(analysis.s)))
        self.summary_tab.sih_std.SetLabel("{:4.2f}".format(np.std(analysis.s)))
        self.summary_tab.sph_min.SetLabel("{:4.0f}".format(np.min(analysis.h)))
        self.summary_tab.sph_max.SetLabel("{:4.0f}".format(np.max(analysis.h)))
        self.summary_tab.sph_avg.SetLabel("{:4.2f}".format(np.mean(analysis.h)))
        self.summary_tab.sph_std.SetLabel("{:4.2f}".format(np.std(analysis.h)))
        self.summary_tab.spa_min.SetLabel("{:4.0f}".format(np.min(analysis.a)))
        self.summary_tab.spa_max.SetLabel("{:4.0f}".format(np.max(analysis.a)))
        self.summary_tab.spa_avg.SetLabel("{:4.2f}".format(np.mean(analysis.a)))
        self.summary_tab.spa_std.SetLabel("{:4.2f}".format(np.std(analysis.a)))

      # Dataset information
      if self.recovery:
        if hasattr(analysis, 'clusters'):
          clusters = analysis.clusters
          pg = clusters[0]['pg']
          uc = clusters[0]['uc']
        else:
          clusters = []
          if hasattr(analysis, 'cons_pg'):
            pg = analysis.cons_pg
          else:
            pg = None
          if hasattr(analysis, 'cons_uc'):
            uc = analysis.cons_uc
          else:
            uc = None
      else:
        analysis.print_results()
        clusters = analysis.unit_cell_analysis()
        pg = clusters[0]['pg']
        uc = clusters[0]['uc']

      if clusters:
        self.summary_tab.report_clustering_results(clusters=clusters)

      if pg is not None:
        self.summary_tab.pg_txt.SetLabel(str(pg))
      if uc is not None:
        if type(uc) is str:
          unit_cell = uc
        else:
          unit_cell = " ".join(['{:4.1f}'.format(i) for i in uc])
      else:
        unit_cell = ''

      self.summary_tab.uc_txt.SetLabel(unit_cell)
      res = to_unicode(u"{:4.2f} - {:4.2f} {}".format(np.mean(analysis.lres),
                                np.mean(analysis.hres),
                                u'\u212B'))
      self.summary_tab.rs_txt.SetLabel(res)

      if self.recovery and (
              hasattr(analysis, 'beamX_mm') and
              hasattr(analysis, 'beamY_mm') and
              hasattr(analysis, 'pixel_size')
              ):
        beamX_mm = analysis.beamX_mm
        beamY_mm = analysis.beamY_mm
        pixel_size = analysis.pixel_size
      else:
        with warnings.catch_warnings():
          # To catch any 'mean of empty slice' runtime warnings
          warnings.simplefilter("ignore", category=RuntimeWarning)
          beamxy_calc = plot.calculate_beam_xy()
          beamX, beamY = beamxy_calc[:2]
          pixel_size = beamxy_calc[-1]
          beamX_mm = np.median(beamX)
          beamY_mm = np.median(beamY)
      beamX_px = beamX_mm / pixel_size
      beamY_px = beamY_mm / pixel_size
      beamXY = "X = {:4.1f} mm / {:4.0f} px\n" \
               "Y = {:4.1f} mm / {:4.0f} px" \
               "".format(beamX_mm, beamX_px, beamY_mm, beamY_px)
      self.summary_tab.xy_txt.SetLabel(beamXY)

      # Summary
      if self.recovery:
        self.summary_tab.readin_txt.SetLabel(str(analysis.n_all_objects))
        self.summary_tab.nodiff_txt.SetLabel(str(analysis.n_no_diff_objects))
        self.summary_tab.w_diff_txt.SetLabel(str(analysis.n_diff_objects))
        if self.gparams.advanced.processing_backend == 'ha14':
          self.summary_tab.noint_txt.SetLabel(str(analysis.n_not_int_objects))
          self.summary_tab.noprf_txt.SetLabel(str(analysis.n_filter_fail_objects))
        elif self.gparams.advanced.processing_backend == 'cctbx.xfel':
          self.summary_tab.nospf_txt.SetLabel(str(analysis.n_not_spf_objects))
          self.summary_tab.noidx_txt.SetLabel(str(analysis.n_not_idx_objects))
          self.summary_tab.noint_txt.SetLabel(str(analysis.n_not_int_objects))
          self.summary_tab.noflt_txt.SetLabel(str(analysis.n_filter_fail_objects))
        self.summary_tab.final_txt.SetLabel(str(analysis.n_final_objects))
      else:
        self.summary_tab.readin_txt.SetLabel(str(len(analysis.all_objects)))
        self.summary_tab.nodiff_txt.SetLabel(str(len(analysis.no_diff_objects)))
        self.summary_tab.w_diff_txt.SetLabel(str(len(analysis.diff_objects)))
        if self.gparams.advanced.processing_backend == 'ha14':
          self.summary_tab.noint_txt.SetLabel(str(len(analysis.not_int_objects)))
          self.summary_tab.noprf_txt.SetLabel(str(len(analysis.filter_fail_objects)))
        elif self.gparams.advanced.processing_backend == 'cctbx.xfel':
          self.summary_tab.nospf_txt.SetLabel(str(len(analysis.not_spf_objects)))
          self.summary_tab.noidx_txt.SetLabel(str(len(analysis.not_idx_objects)))
          self.summary_tab.noint_txt.SetLabel(str(len(analysis.not_int_objects)))
          self.summary_tab.noflt_txt.SetLabel(str(len(analysis.filter_fail_objects)))
        self.summary_tab.final_txt.SetLabel(str(len(analysis.final_objects)))

        # Generate input file for PRIME
        analysis.print_summary()
        analysis.make_prime_input(filename=prime_file)

      # Display summary
      self.proc_nb.AddPage(self.summary_tab, 'Summary', select=True)
      self.proc_nb.SetSelection(self.proc_nb.GetPageCount())

      # Signal end of run
      font = self.sb.GetFont()
      font.SetWeight(wx.BOLD)
      self.status_txt.SetFont(font)
      self.status_txt.SetForegroundColour('blue')
      self.status_txt.SetLabel('DONE')

    # Finish up
    self.display_log()
    self.plot_integration(force_plot=True)
    self.plot_live_analysis(force_plot=True)

    # Stop timer
    self.proc_timer.Stop()
    self.chart_timer.Stop()

  def display_log(self):
    if os.path.isfile(self.init.logfile):
      with open(self.init.logfile, 'r') as out:
        out.seek(self.bookmark)
        output = out.read()
        self.bookmark = out.tell()

      ins_pt = self.log_tab.log_window.GetInsertionPoint()
      self.log_tab.log_window.AppendText(output)
      self.log_tab.log_window.SetInsertionPoint(ins_pt)

  def plot_integration(self, force_plot=False):
    """ This function will plot fast-drawing runtime processing charts on the
        "Processing" tab """
    sw = wx.StopWatch()
    if (self.info.nref_list and self.info.res_list) :
      self.proc_tab.init = self.init
      self.proc_tab.gparams = self.gparams
      self.proc_tab.info = self.info
      if self.proc_nb.GetSelection() == 1 or force_plot:
        s = time.time()
        self.proc_tab.draw_plots()
        self.proc_tab.draw_summary()
    self.plotter_time.append(sw.Time())

  def plot_live_analysis(self, force_plot=False):
    """ This function will plot in-depth analysis that will (importantly)
        involve expensive and slow-drawing charts on the Live Analysis tab """

    if self.proc_nb.GetSelection() == 2 or force_plot:
      self.chart_tab.init = self.init
      self.chart_tab.gparams = self.gparams
      self.chart_tab.finished_objects = self.info.finished_objects

      if self.cluster_info:
        self.chart_tab.cluster_info = self.cluster_info

      if self.prime_info:
        self.chart_tab.prime_info = self.prime_info

      self.chart_tab.draw_plots()

  def onChartTimer(self, e):
    self.plot_live_analysis()
    if not (self.running_cluster or self.running_prime):
      self.run_clustering_thread()

  def run_clustering_thread(self):
    self.chart_tab.btn_run_analysis.Disable()
    self.running_cluster = True
    self.cluster_thread = thr.ClusterThread(self,
                                            iterable=self.info.cluster_iterable)
    self.cluster_thread.start()

  def onFinishedCluster(self, e):
    self.cluster_info = e.GetValue()
    self.running_cluster = False

    if self.cluster_info is not None:
      cl_sorted = sorted(self.cluster_info, key=lambda i: i['number'],
                         reverse=True)
      self.info.best_pg = cl_sorted[0]['pg'].split('/')[0]
      self.info.best_uc = cl_sorted[0]['uc']

    # Output cluster results
    ep.dump(self.info.cluster_info_file, obj=self.cluster_info)

    if not (hasattr(self, 'prime_thread') and self.prime_thread.isAlive()):
      self.pparams = None
      self.run_prime_thread()

  def run_prime_thread(self):
    # Run PRIME (basic merge only)
    self.chart_tab.btn_run_analysis.Disable()
    self.running_prime = True
    pg = ut.makenone(self.chart_tab.pg_uc.pg.GetValue())
    uc = ut.makenone(self.chart_tab.pg_uc.uc.GetValue())

    self.prime_thread = thr.PRIMEThread(self,
                                        info=self.info,
                                        init=self.init,
                                        best_pg=pg, best_uc=uc)

    self.prime_thread.start()

  def onFinishedPRIME(self, e):
    self.prime_info = e.GetValue()
    self.running_prime = False

    if self.running_manual_analysis:
      self.plot_live_analysis(force_plot=True)
      self.running_manual_analysis = False
    self.chart_tab.btn_run_analysis.Enable()

  def adjust_timer(self, new_interval):
    self.proc_timer.Stop()
    self.proc_timer.Start(new_interval)
    print ("IOTA: TIMER SET TO {} SECONDS".format(new_interval/1000))


  def onProcTimer(self, e):
    sw_means = np.mean(self.plotter_time) + np.mean(self.obj_reader_time)
    timer_adjustment = np.ceil(sw_means*5/1000)*1000
    if not (
            np.isnan(timer_adjustment) or
            timer_adjustment == self.proc_timer.GetInterval()
    ):
      self.adjust_timer(new_interval=timer_adjustment)

    # Info object handling
    if self.info:
      # Save info object (could be trouble)
      info_filepath = os.path.join(self.init.int_base, 'proc.info')
      ep.dump(info_filepath, self.info)

      # Plot integration and display lob
      self.plot_integration()
      self.display_log()

    if self.abort_initiated:
      if self.gparams.mp.method == 'lsf':
        info_command = 'bjobs -J {}'.format(self.job_id)
        lsf_info = easy_run.fully_buffered(info_command).stdout_lines
        self.run_aborted = (lsf_info == [])
      else:
        self.run_aborted = (not self.proc_thread.is_alive())
        if not self.run_aborted:
          self.proc_thread.abort()

      if self.run_aborted:
        print ("IOTA: RUN ABORTED!")
        self.finish_process()
    else:
      self.run_aborted = False

    # If all current images processed, check for new ones and/or finish process
    if self.info.unprocessed_images <= 0:
      self.monitor_filesystem()
    else:
      # Instantiate and start processing thread
      if not (hasattr(self, 'object_reader') and self.object_reader.is_alive()):
        self.obj_sw = wx.StopWatch()
        self.object_reader = thr.ObjectReaderThread(self, info=self.info,
                                          source=self.info.obj_list_file)
        self.object_reader.start()

  def monitor_filesystem(self):
    # Check if all images have been looked at; if yes, finish process
    if self.monitor_mode:
      if self.find_new_images:
        self.get_images_from_filesystem()
      else:
        self.find_new_images = (not self.new_images)

      if self.new_images:
        self.status_txt.SetLabel(
          'Found {} new images'.format(len(self.new_images)))
        self.timeout_start = None
        self.state = 'new images'
        self.process_images()
      else:
        if self.monitor_mode_timeout:
          if self.timeout_start is None:
            self.timeout_start = time.time()
          else:
            interval = time.time() - self.timeout_start
            if interval >= self.monitor_mode_timeout:
              self.status_txt.SetLabel('Timed out. Finishing...')
              self.finish_process()
            else:
              timeout_msg = 'No images found! Timing out in {} seconds' \
                            ''.format(
                int(self.monitor_mode_timeout - interval))
              self.status_txt.SetLabel(timeout_msg)
        else:
          self.status_txt.SetLabel('No new images found! Waiting ...')
    else:
      self.status_txt.SetLabel('Wrapping up ...')
      self.finish_process()

  def get_images_from_filesystem(self):
    img_finder = thr.ImageFinderThread(self,
                                       image_paths=self.gparams.input,
                                       image_list=self.info.img_list)
    img_finder.start()

  def onFinishedProcess(self, e):
    if not self.monitor_mode:
      self.finish_process()
    else:
      self.find_new_images = True

  def onFinishedImageFinder(self, e):
    self.new_images = e.GetValue()

  def onFinishedObjectReader(self, info):
    # Read info object
    self.info = info.GetValue()

    # Update w/ latest Cluster and PRIME results
    self.info.prime_info = self.prime_info
    self.info.cluster_info = self.cluster_info

    if not self.recovery:
      self.obj_reader_time.append(self.obj_sw.Time())


  def finish_process(self):
    import shutil
    self.proc_timer.Stop()
    self.chart_timer.Stop()

    if self.finished_objects is None:
      font = self.sb.GetFont()
      font.SetWeight(wx.BOLD)
      self.status_txt.SetFont(font)
      self.status_txt.SetForegroundColour('blue')
      self.status_txt.SetLabel('OBJECT READ-IN ERROR! NONE IMPORTED')
      return

    if self.run_aborted:
      self.gauge_process.Hide()
      font = self.sb.GetFont()
      font.SetWeight(wx.BOLD)
      self.status_txt.SetFont(font)
      self.status_txt.SetForegroundColour('red')
      self.status_txt.SetLabel('ABORTED BY USER')
      self.proc_toolbar.EnableTool(self.tb_btn_resume.GetId(), True)
      try:
        shutil.rmtree(self.init.tmp_base)
      except Exception:
        pass
      print('JOB TERMINATED!')
      return
    elif self.recovery:
      self.gauge_process.Hide()
      font = self.sb.GetFont()
      font.SetWeight(wx.BOLD)
      self.status_txt.SetFont(font)
      run_no = int(self.init.int_base.split('/')[-1])
      self.status_txt.SetLabel('Run #{} Loaded!'.format(run_no))
      self.proc_toolbar.EnableTool(self.tb_btn_abort.GetId(), False)
      self.proc_toolbar.EnableTool(self.tb_btn_monitor.GetId(), False)
      self.proc_toolbar.ToggleTool(self.tb_btn_monitor.GetId(), False)

      self.display_log()

      analysis_file = os.path.join(self.init.int_base, 'analysis.pickle')
      if os.path.isfile(analysis_file):
        analysis = ep.load(analysis_file)
      else:
        analysis = None

      # TODO: check if you even need this for backwards compatibility
      from multiprocessing import cpu_count
      n_proc = int(cpu_count() * 0.75)
      if not self.info.finished_objects:  # Find finished objects if none loaded
        source = self.info.obj_list_file  # default to tmp/finished_objects.lst
        if not os.path.isfile(source):    # if not, find object folder
          source = self.init.obj_base
        self.object_reader = thr.ObjectReaderThread(self,
                                                    source=source,
                                                    info=self.info,
                                                    n_proc=n_proc)
        self.object_reader.start()

      # Check for and recover clustering data
      if not self.info.cluster_info:
        cluster_info_file = os.path.join(self.init.int_base,
                                         'cluster_info.pickle')
        if os.path.isfile(cluster_info_file):
          self.cluster_info = ep.load(cluster_info_file)
      else:
        self.cluster_info = self.info.cluster_info

        # Check for and recover live PRIME results
      if not self.info.prime_info:
        live_prime_phil = os.path.join(self.init.int_base, 'live_prime.phil')
        if os.path.isfile(live_prime_phil):
          import iotbx.phil as ip
          from prime.postrefine.mod_input import master_phil as prime_master_phil
          with open(live_prime_phil, 'r') as lpf:
            contents = lpf.read()
          prime_phil = ip.parse(contents)
          prime_phil = prime_master_phil.fetch(prime_phil)
          self.pparams = prime_phil.extract()
          live_prime_info_file = os.path.join(self.init.int_base,
                                              'life_prime_info.pickle')
          if os.path.isfile(live_prime_info_file):
            self.prime_info = ep.load(live_prime_info_file)
      else:
        self.prime_info = self.info.prime_info

      if (
              str(self.state).lower() == 'finished' or
              len( self.info.finished_objects) >= len(self.info.img_list)
      ):
        self.final_objects = [i for i in self.info.finished_objects if not i.fail]
        self.analyze_results(analysis=analysis)

      else:
        if len(self.info.finished_objects) > 0:
          self.plot_integration(force_plot=True)
          self.plot_live_analysis(force_plot=True)
        if os.path.isfile(os.path.join(self.init.int_base, 'init.cfg')):
          self.proc_toolbar.EnableTool(self.tb_btn_resume.GetId(), True)
      return
    else:
      self.final_objects = [i for i in self.info.finished_objects if not i.fail]
      self.gauge_process.Hide()
      self.proc_toolbar.EnableTool(self.tb_btn_abort.GetId(), False)
      self.proc_toolbar.EnableTool(self.tb_btn_monitor.GetId(), False)
      self.proc_toolbar.ToggleTool(self.tb_btn_monitor.GetId(), False)

      if len(self.final_objects) > 0:
        # Signal end of run
        self.plot_integration(force_plot=True)
        self.plot_live_analysis(force_plot=True)
        self.analyze_results()
      else:
        font = self.sb.GetFont()
        font.SetWeight(wx.BOLD)
        self.status_txt.SetFont(font)
        self.status_txt.SetForegroundColour('blue')
        self.status_txt.SetLabel('NO IMAGES PROCESSED')

      print ('JOB FINISHED!')

      try:
        shutil.rmtree(self.init.tmp_base)
      except Exception:
        pass
