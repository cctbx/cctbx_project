from __future__ import division, print_function

from six.moves import range
'''
Author      : Lyubimov, A.Y.
Created     : 12/12/2017
Last Changed: 01/29/2018
Description : SIMTBX (nanoBragg) GUI Windows / frames
'''

import os
import wx
import numpy as np

from wxtbx import bitmaps

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure

from iotbx import phil as ip
from simtbx.nanoBragg import nanoBragg_gui_dialogs as dlg
from simtbx.nanoBragg import nanoBragg_threads as thr
from iota.components import iota_ui_controls as ct
from iota.components.iota_utils import InputFinder, WxFlags, noneset

ginp = InputFinder()
import time

f = WxFlags()

# ------------------------------ Input Window -------------------------------- #

class BasePanel(wx.Panel):
  def __init__(self, parent):
    wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY, size=(800, 800))

    self.main_sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.main_sizer)


class TopPanel(BasePanel):
  def __init__(self, parent):
    BasePanel.__init__(self, parent=parent)

    self.parent = parent
    self.img_filename = None
    self.input_phil = None
    self.sparams = None
    self.coord_filename = None
    self.mtz_filename = None
    self.stol_filename = None
    self.ref_img_filename = None
    self.dblclick = False


    self.project_folder = ct.InputCtrl(self,
                                       label='Project Folder: ',
                                       label_size=(150, -1),
                                       label_style='bold',
                                       value=os.path.abspath(os.curdir),
                                       buttons=True)

    self.project_title = ct.InputCtrl(self,
                                      label='Description',
                                      label_size=(150, -1),
                                      label_style='normal')

    self.splitter = wx.SplitterWindow(self, style=wx.SP_LIVE_UPDATE |
                                                  wx.SP_3DSASH |
                                                  wx.SP_NOBORDER)
    self.file_panel = wx.Panel(self.splitter, size=(-1, 200))
    self.file_sizer = wx.BoxSizer(wx.VERTICAL)
    self.file_panel.SetSizer(self.file_sizer)
    self.preview_panel = wx.Panel(self.splitter, size=(-1, 500))
    self.preview_sizer = wx.BoxSizer(wx.VERTICAL)
    self.preview_panel.SetSizer(self.preview_sizer)
    self.splitter.SplitHorizontally(self.file_panel, self.preview_panel, 200)

    self.input = FileListCtrl(self.file_panel)
    self.file_sizer.Add(self.input, 1, flag=wx.ALL | wx.EXPAND, border=10)

    # Image preview box w/ options
    prev_box = wx.GridBagSizer(5, 15)

    # self.opt_knb_start = ct.KnobCtrl(self.preview_panel,
    #                                  label='start',
    #                                  label_size=(40, -1),
    #                                  spin_ctr_size=(50, -1),
    #                                  knob_size=(120, 120),
    #                                  values_start=0,
    #                                  values_end=360,
    #                                  values_step=1,
    #                                  value=0)
    #
    # self.opt_knb_end = ct.KnobCtrl(self.preview_panel,
    #                                label='end',
    #                                label_size=(40, -1),
    #                                spin_ctr_size=(50, -1),
    #                                knob_size=(120, 120),
    #                                values_start=0,
    #                                values_end=360,
    #                                values_step=1,
    #                                value=360)

    self.opt_spc_start = ct.SpinCtrl(self.preview_panel,
                                     label='start:',
                                     label_size=(50, -1),
                                     ctrl_size=(60, -1),
                                     ctrl_value='0',
                                     ctrl_max=360,
                                     ctrl_min=0,
                                   ctrl_step=0.1,
                                   ctrl_digits=1)

    self.opt_spc_end = ct.SpinCtrl(self.preview_panel,
                                   label='finish:',
                                   label_size=(50, -1),
                                   ctrl_size=(60, -1),
                                   ctrl_value='360',
                                   ctrl_max=360,
                                   ctrl_min=0,
                                   ctrl_step=0.1,
                                   ctrl_digits=1)

    self.opt_spc_osc = ct.SpinCtrl(self.preview_panel,
                                   label='step:',
                                   label_size=(50, -1),
                                   ctrl_size=(60, -1),
                                   ctrl_value='1.0',
                                   ctrl_max=360,
                                   ctrl_min=0,
                                   ctrl_step=0.1,
                                   ctrl_digits=2)

    self.opt_btn_prev = wx.Button(self.preview_panel,
                                  label='PREVIEW IMAGE')

    self.opt_chk_bkg = wx.CheckBox(self.preview_panel, label='Add Background')
    self.opt_chk_bkg.SetValue(True)
    self.opt_chk_noise = wx.CheckBox(self.preview_panel, label='Add Noise')
    self.opt_chk_noise.SetValue(True)
    self.opt_chk_rand = wx.CheckBox(self.preview_panel,
                                    label='Randomize Orientation')
    self.opt_chk_rand.SetValue(True)

    self.opt_spc_scale = ct.SpinCtrl(self.preview_panel,
                                     label='Crystal size (um): ',
                                     label_size=(120, -1),
                                     ctrl_size=(100, -1),
                                     ctrl_min = 1,
                                     ctrl_max = 1000,
                                     ctrl_step = 1,
                                     ctrl_value = 30)

    self.img_figure = Figure(figsize=(3, 3))
    self.img_axes = self.img_figure.add_subplot(111)

    self.img_axes.set_frame_on(False)
    self.img_axes.axis('off')
    self.img_axes.set_aspect('equal')

    self.img_figure.patch.set_visible(False)
    self.img_canvas = FigureCanvas(self.preview_panel, -1, self.img_figure)

    prev_box.Add(self.opt_spc_start, pos=(0, 0))
    prev_box.Add(self.opt_spc_end, pos=(0, 1))
    prev_box.Add(self.opt_spc_osc, pos=(1, 0))
    prev_box.Add(self.opt_chk_bkg, flag=wx.EXPAND, pos=(4, 0), span=(1, 2))
    prev_box.Add(self.opt_chk_noise, flag=wx.EXPAND, pos=(5, 0), span=(1, 2))
    prev_box.Add(self.opt_chk_rand, flag=wx.EXPAND, pos=(6, 0), span=(1, 2))
    prev_box.Add(self.opt_spc_scale, flag=wx.EXPAND, pos=(7, 0), span=(1, 2))
    prev_box.Add(self.opt_btn_prev, flag=wx.EXPAND, pos=(8,0), span=(1, 2))
    prev_box.Add(self.img_canvas, pos=(0, 2), span=(9, 1),
                 flag=wx.EXPAND)
    prev_box.AddGrowableCol(2)
    prev_box.AddGrowableRow(8)

    self.preview_sizer.Add(prev_box, 1, flag=wx.EXPAND | wx.ALL, border=10)

    self.main_sizer.Add(self.project_title, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.project_folder, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.splitter, 1, flag=wx.EXPAND)

    # Button bindings
    self.Bind(wx.EVT_BUTTON, self.onRunPreview, self.opt_btn_prev)

    # Thread bindings
    self.Bind(thr.EVT_NBDONE, self.onFinishedSimThread)

    # Image bindings
    xid = self.img_canvas.mpl_connect('button_press_event', self.on_button_press)
    xid = self.img_canvas.mpl_connect('button_release_event',
                                      self.on_button_release)

  def onRunPreview(self, e):
    e.Skip()

  def generate_phil(self):

    coord_path = None
    FCalc_path = None
    bkg_path = None
    img_path = None

    # Grab inputs (if any) from file list control
    idxs = self.input.ctr.GetItemCount()
    inputs = [self.input.ctr.GetItemData(i) for i in range(idxs)]

    for idx in range(idxs):
      item = self.input.ctr.GetItemData(idx)
      item_type_sel = item.type_selection
      item_type = item.type.type.GetString(item_type_sel)

      if item_type == 'coordinates':
        coord_path = item.path
      elif item_type == 'structure factors':
        FCalc_path = item.path
      elif item_type == 'background':
        bkg_path = item.path
      elif item_type == 'raw image file':
        img_path = item.path

    simtbx_phil_string = '\n'.join([
    'description = {}'.format(noneset(self.project_title.ctr.GetValue())),
    'output = {}'.format(noneset(self.project_folder.ctr.GetValue())),
    'reference_coordinates = {}'.format(coord_path),
    'reference_FCalc = {}'.format(FCalc_path),
    'reference_image = {}'.format(img_path),
    'radial_average_background = {}'.format(bkg_path),
    'dataset',
    '{',
    '  start_phi = {}'.format(str(self.opt_spc_start.ctr.GetValue())),
    '  finish_phi = {}'.format(str(self.opt_spc_end.ctr.GetValue())),
    '  oscillation = {}'.format(str(self.opt_spc_osc.ctr.GetValue())),
    '  }'
    ])

    self.input_phil = ip.parse(simtbx_phil_string)

  def run_simulator(self, init=None):
    pass

  def run_preview(self):
    img_filename = 'test_image.{}'.format(self.sparams.image_format)
    self.img_filename = os.path.join(self.sparams.output, img_filename)
    self.start_timer = time.time()
    self.sim = thr.nanoBraggThread(self,
                                   params=self.sparams,
                                   add_background=self.opt_chk_bkg.GetValue(),
                                   add_noise=self.opt_chk_noise.GetValue(),
                                   randomize=self.opt_chk_rand.GetValue(),
                                   pixel_scale=self.opt_spc_scale.ctr.GetValue(),
                                   preview=True)
    self.sim.start()

  def onFinishedSimThread(self, e):
    pixels = e.GetValue()
    self.display_image(pixels=pixels)
    print('TOTAL TIME = ', time.time() - self.start_timer)

  def display_image(self, pixels=None):
    if pixels is None:
      pixels = np.random.randint(low=0, high=1500000, size=(2576, 2576))
    else:
      pixels = pixels.as_numpy_array()

    clim = (pixels.min(), np.percentile(pixels, 99))

    self.img_axes.imshow(pixels,
                         #interpolation='nearest',
                         cmap='gray_r',
                         clim=clim)
    self.img_figure.subplots_adjust(left=0, bottom=0, right=1, top=1)

    self.preview_panel.Layout()
    print('DEBUG: AVERAGE PIXEL VALUE = ', np.mean(pixels))
    print('DONE!')

  def on_button_press(self, e):
    if e.button == 1 and e.dblclick:
      self.dblclick = True
    else:
      self.dblclick = False

  def on_button_release(self, e):
    if e.button == 1 and self.dblclick:
      self.view_image()

  def view_image(self):
    viewer = thr.ImageViewerThread(self,
                                   file_string=self.img_filename)
    viewer.start()

class FileListCtrl(ct.CustomListCtrl):
  ''' File list window for the input tab '''

  def __init__(self, parent, size=(-1, 250)):
    ct.CustomListCtrl.__init__(self, parent=parent, size=size)

    self.parent = parent
    self.main_window = parent.GetParent()

    # Generate columns
    self.ctr.InsertColumn(0, "Input Path")
    self.ctr.InsertColumn(1, "Input Type")
    self.ctr.InsertColumn(2, "Action")
    self.ctr.setResizeColumn(1)

    # Add file / folder buttons
    self.button_sizer = wx.BoxSizer(wx.HORIZONTAL)
    self.btn_add_file = wx.Button(self, label='Add File...')
    #self.btn_add_dir = wx.Button(self, label='Add Folder...')
    self.button_sizer.Add(self.btn_add_file)
    #self.button_sizer.Add(self.btn_add_dir, flag=wx.LEFT, border=10)

    self.sizer.Add(self.button_sizer, flag=wx.TOP | wx.BOTTOM, border=10)

    # Event bindings
    self.Bind(wx.EVT_BUTTON, self.onAddFile, self.btn_add_file)

    self.Layout()

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

  def set_type_choices(self, path):
    # Determine what type of input this is and present user with choices
    # (this so far works for images ONLY)
    type_choices = ['[  SELECT INPUT TYPE  ]',
                    'coordinates',
                    'structure factors',
                    'background',
                    'raw image file']
    preferred_selection = 0
    inputs, input_type = ginp.get_input(path, filter=False)

    if input_type == 'data (MTZ)':
      input_type = 'structure factors'
    elif input_type == 'text' and path.endswith('stol'):
      input_type = 'background'

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
    # item.buttons.btn_mag.Bind(wx.EVT_BUTTON, self.onMagButton)
    # item.buttons.btn_delete.Bind(wx.EVT_BUTTON, self.onDelButton)
    # item.buttons.btn_info.Bind(wx.EVT_BUTTON, self.onInfoButton)
    self.Bind(wx.EVT_BUTTON, self.onMagButton, item.buttons.btn_mag)
    self.Bind(wx.EVT_BUTTON, self.onDelButton, item.buttons.btn_delete)
    self.Bind(wx.EVT_BUTTON, self.onInfoButton, item.buttons.btn_info)

    # Insert list item
    idx = self.ctr.InsertStringItem(self.ctr.GetItemCount() + 1, item.path)
    self.ctr.SetItemWindow(idx, 1, item.type, expand=True)
    self.ctr.SetItemWindow(idx, 2, item.buttons, expand=True)

    # Set drop-down selection, check it for data and open other tabs
    item.type.type.SetSelection(inp_sel)
    if item.type.type.GetString(inp_sel) in ('coordinates',
                                             'structure factors',
                                             'background',
                                             'raw image file'):
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
    self.ctr.SetColumnWidth(1, width=-1)
    self.ctr.SetColumnWidth(2, width=-1)
    self.ctr.SetColumnWidth(0, width=-3)

    # Make sure all the choice lists are the same size
    if item.type.type.GetSize()[0] < self.ctr.GetColumnWidth(2) - 5:
       item.type.type.SetSize((self.ctr.GetColumnWidth(2) - 5, -1))

    # Attach data object to item
    self.ctr.SetItemData(item.id, item)

    self.Layout()

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
      elif type == 'text':
        with open(path, 'r') as f:
          file_list = f.readlines()
          msg = ' '.join(file_list)
          textview = dlg.TextFileView(self, title=path, contents=msg)
          textview.ShowModal()
      else:
        wx.MessageBox('Unknown file format', 'Warning',
                      wx.OK | wx.ICON_EXCLAMATION)

  def view_images(self, img_list, img_type=None):
    ''' Launches image viewer (depending on backend) '''
    viewer = self.parent.gparams.advanced.image_viewer
    if viewer == 'cxi.view' and 'pickle' not in img_type:
        wx.MessageBox('cxi.view only accepts image pickles', 'Warning',
                      wx.OK | wx.ICON_EXCLAMATION)
    else:
      if len(img_list) > 10:
        view_warning = dlg.ViewerWarning(self, len(img_list))
        if view_warning.ShowModal() == wx.ID_OK:
          # parse 'other' entry
          img_no_string = str(view_warning.no_images).split(',')
          filenames = []
          for n in img_no_string:
            if '-' in n:
              img_limits = n.split('-')
              start = int(min(img_limits))
              end = int(max(img_limits))
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
    self.ctr.DeleteItem(index)

    # Refresh widget and list item indices
    if self.ctr.GetItemCount() > 0:
      for i in range(self.ctr.GetItemCount()):
        item_data = self.ctr.GetItemData(i)
        item_data.id = i
        item_data.buttons.index = i
        item_data.type.index = i
        type_choice = self.ctr.GetItemWindow(i, col=1)
        type_selection = item_data.type.type.GetSelection()
        type_choice.type.SetSelection(type_selection)
        self.ctr.SetItemData(i, item_data)

  def onInfoButton(self, e):
    ''' Info / alert / error button (will change depending on circumstance) '''
    idx = e.GetEventObject().GetParent().index
    item_obj = self.ctr.GetItemData(idx)
    item_type = item_obj.type.type.GetString(item_obj.type_selection)

    if item_obj.warning:
      wx.MessageBox(item_obj.info['WARNING'], 'Warning', wx.OK |
                    wx.ICON_EXCLAMATION)
    else:
      wx.MessageBox(item_obj.info[item_type], 'Info', wx.OK |
                    wx.ICON_INFORMATION)
