from __future__ import division

# LIBTBX_SET_DISPATCHER_NAME iota.test2

'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 04/14/2016
Description : IOTA GComps
'''

import wx
import os
from cctbx.uctbx import unit_cell

icons = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ic/')

# ----------------------------------- Window Panels ------------------------------------ #

class InputWindow(wx.Panel):
  ''' Input window - data input, description of project '''

  def __init__(self, parent):
    super(InputWindow, self).__init__(parent)

    # Generate default parameters
    from iota.components.iota_input import master_phil
    self.gparams = master_phil.extract()

    main_box = wx.StaticBox(self, label='Input', size=(770, -1))
    vbox = wx.StaticBoxSizer(main_box, wx.VERTICAL)

    # Integration software choice
    int_choice_box = wx.BoxSizer(wx.HORIZONTAL)
    progs = ['cctbx.xfel', 'DIALS']
    self.int_txt_choice = wx.StaticText(self, label='Integrate with: ', size=(120, -1))
    self.int_cbx_choice = wx.Choice(self, choices=progs)
    int_choice_box.Add(self.int_txt_choice)
    int_choice_box.Add(self.int_cbx_choice, flag=wx.LEFT, border=5)
    vbox.Add(int_choice_box, flag=wx.LEFT | wx.TOP, border=10)

    # Input box and Browse button
    input_box = wx.BoxSizer(wx.HORIZONTAL)
    self.inp_txt = wx.StaticText(self, label='Input path: ', size=(120, -1))
    self.inp_btn_browse = wx.Button(self, label='Browse...')
    self.inp_btn_mag = wx.BitmapButton(self, bitmap=wx.Bitmap('{}/16x16/viewmag.png'.format(icons)))
    ctr_length = main_box.GetSize()[0] - self.inp_txt.GetSize()[0] - self.inp_btn_browse.GetSize()[0] - \
                 self.inp_btn_mag.GetSize()[0] - 55
    self.inp_ctr = wx.TextCtrl(self, size=(ctr_length, -1))
    input_box.Add(self.inp_txt)
    input_box.Add(self.inp_ctr, flag=wx.LEFT, border=5)
    input_box.Add(self.inp_btn_browse, flag=wx.LEFT, border=10)
    input_box.Add(self.inp_btn_mag, flag=wx.LEFT | wx.RIGHT, border=10)
    vbox.Add(input_box, flag=wx.LEFT | wx.TOP, border=10)

    # Output box and Browse button
    output_box = wx.BoxSizer(wx.HORIZONTAL)
    self.out_txt = wx.StaticText(self, label='Output path: ', size=(120, -1))
    self.out_btn_browse = wx.Button(self, label='Browse...')
    self.out_btn_mag = wx.BitmapButton(self, bitmap=wx.Bitmap('{}/16x16/viewmag.png'.format(icons)))
    self.out_ctr = wx.TextCtrl(self, size=(ctr_length, -1))
    self.out_ctr.SetValue(os.path.abspath(os.curdir))
    output_box.Add(self.out_txt)
    output_box.Add(self.out_ctr, flag=wx.LEFT, border=5)
    output_box.Add(self.out_btn_browse, flag=wx.LEFT, border=10)
    output_box.Add(self.out_btn_mag, flag=wx.LEFT | wx.RIGHT | wx.ALIGN_CENTER_HORIZONTAL, border=10)
    vbox.Add(output_box, flag=wx.LEFT | wx.TOP, border=10)

    # Script input box and browse button
    # script_box = wx.BoxSizer(wx.HORIZONTAL)
    # self.script_txt = wx.StaticText(self, label='IOTA Script: ', size=(120, -1))
    # self.script_btn_browse = wx.Button(self, label='Browse...')
    # self.script_btn_mag = wx.BitmapButton(self, bitmap=wx.Bitmap('{}/16x16/viewmag.png'.format(icons)))
    # self.script_ctr = wx.TextCtrl(self, size=(ctr_length, -1))
    # script_box.Add(self.script_txt)
    # script_box.Add(self.script_ctr, flag=wx.LEFT, border=5)
    # script_box.Add(self.script_btn_browse, flag=wx.LEFT, border=10)
    # script_box.Add(self.script_btn_mag, flag=wx.LEFT | wx.RIGHT | wx.ALIGN_CENTER_HORIZONTAL, border=10)
    # vbox.Add(script_box, flag=wx.LEFT | wx.TOP, border=10)

    # Title box
    title_box = wx.BoxSizer(wx.HORIZONTAL)
    self.inp_des_txt = wx.StaticText(self, label='Job title: ', size=(120, -1))
    self.inp_des_ctr = wx.TextCtrl(self, size=(ctr_length, -1))
    title_box.Add(self.inp_des_txt)
    title_box.Add(self.inp_des_ctr, flag=wx.LEFT, border=5)
    title_box.Add((10, -1))
    vbox.Add(title_box, flag=wx.LEFT | wx.TOP, border=10)

    # Input / Main Options
    opt_box_1 = wx.BoxSizer(wx.HORIZONTAL)
    self.opt_chk_random = wx.CheckBox(self, label='Random subset')
    self.opt_chk_random.SetValue(False)
    self.opt_txt_random = wx.StaticText(self, label='Images in subset:')
    self.opt_txt_random.Disable()
    self.opt_spc_random = wx.SpinCtrl(self, value='5', size=(80, -1))
    self.opt_spc_random.Disable()
    self.opt_txt_nprocs = wx.StaticText(self, label='Number of processors: ')
    self.opt_spc_nprocs = wx.SpinCtrl(self, value='8', size=(80, -1))
    nproc_spacer = main_box.GetSize()[0] - self.opt_chk_random.GetSize()[0] - self.opt_txt_random.GetSize()[0] - \
                   self.opt_spc_random.GetSize()[0] - self.opt_txt_nprocs.GetSize()[0] - \
                   self.opt_spc_nprocs.GetSize()[0] - self.inp_btn_browse.GetSize()[0] - 60
    opt_box_1.Add(self.opt_chk_random)
    opt_box_1.Add((10, -1))
    opt_box_1.Add(self.opt_txt_random, flag=wx.LEFT, border=10)
    opt_box_1.Add(self.opt_spc_random, flag=wx.LEFT, border=10)
    opt_box_1.Add(self.opt_txt_nprocs, flag=wx.LEFT, border=nproc_spacer)
    opt_box_1.Add(self.opt_spc_nprocs, flag=wx.LEFT, border=10)
    vbox.Add(opt_box_1, flag=wx.LEFT | wx.TOP, border=10)

    # Buttons for Other Options
    opt_box_2 = wx.BoxSizer(wx.HORIZONTAL)
    self.prime_txt_prefix = wx.StaticText(self, label='PRIME input prefix:')
    self.prime_ctr_prefix = wx.TextCtrl(self, size=(150, -1))
    self.prime_ctr_prefix.SetValue('prime')
    self.opt_btn_import = wx.Button(self, label='Import options...')
    self.opt_btn_process = wx.Button(self, label='Processing options...')
    self.opt_btn_analysis = wx.Button(self, label='Analysis options...')
    opt_box_2.Add(self.prime_txt_prefix)
    opt_box_2.Add(self.prime_ctr_prefix, flag=wx.LEFT, border=5)
    opt_box_2.Add(self.opt_btn_import, flag=wx.LEFT, border=10)
    opt_box_2.Add(self.opt_btn_process, flag=wx.LEFT, border=10)
    opt_box_2.Add(self.opt_btn_analysis, flag=wx.LEFT, border=10)
    vbox.Add((-1, 20))
    vbox.Add(opt_box_2, flag=wx.LEFT | wx.BOTTOM, border=10)
    vbox.Add((-1, 10))

    self.SetSizer(vbox)

    # Button bindings
    self.inp_btn_browse.Bind(wx.EVT_BUTTON, self.onInputBrowse)
    self.out_btn_browse.Bind(wx.EVT_BUTTON, self.onOutputBrowse)
    self.opt_chk_random.Bind(wx.EVT_CHECKBOX, self.onRandomCheck)
    self.opt_btn_import.Bind(wx.EVT_BUTTON, self.onImportOptions)
    self.opt_btn_process.Bind(wx.EVT_BUTTON, self.onProcessOptions)
    self.opt_btn_analysis.Bind(wx.EVT_BUTTON, self.onAnalysisOptions)


  def onImportOptions(self, e):
    ''' On clicking the Import Options button opens the Import Options window '''
    imp_dialog = ImportWindow(self,
                              title='Import Options',
                              style=wx.DEFAULT_DIALOG_STYLE | wx.STAY_ON_TOP)
    imp_dialog.Fit()

    # Set values to defaults or previously-selected params
    if str(self.gparams.image_conversion.rename_pickle_prefix).lower() == 'none':
      imp_dialog.conv_rb1_prefix.SetValue(True)
    elif str(self.gparams.image_conversion.rename_pickle_prefix).lower() == 'auto':
      imp_dialog.conv_rb2_prefix.SetValue(True)
    else:
      imp_dialog.conv_rb3_prefix.SetValue(True)
    imp_dialog.rbPrefix(wx.EVT_RADIOBUTTON)
    imp_dialog.conv_ctr_prefix.SetValue(str(self.gparams.image_conversion.rename_pickle_prefix))
    if str(self.gparams.image_conversion.square_mode).lower() == 'none':
      imp_dialog.mod_rb1_square.SetValue(True)
    elif str(self.gparams.image_conversion.square_mode).lower() == 'crop':
      imp_dialog.mod_rb2_square.SetValue(True)
    elif str(self.gparams.image_conversion.square_mode).lower() == 'pad':
      imp_dialog.mod_rb3_square.SetValue(True)
    imp_dialog.mod_ctr_beamstop.SetValue(str(self.gparams.image_conversion.beamstop))
    imp_dialog.mod_ctr_dist.SetValue(str(self.gparams.image_conversion.distance))
    imp_dialog.mod_ctr_beamx.SetValue(str(self.gparams.image_conversion.beam_center.x))
    imp_dialog.mod_ctr_beamy.SetValue(str(self.gparams.image_conversion.beam_center.y))
    if str(self.gparams.image_triage.type).lower() == 'none':
      imp_dialog.trg_rb1_mode.SetValue(True)
    elif str(self.gparams.image_triage.type).lower() == 'simple':
      imp_dialog.trg_rb2_mode.SetValue(True)
    elif str(self.gparams.image_triage.type).lower() == 'grid_search':
      imp_dialog.trg_rb3_mode.SetValue(True)
    imp_dialog.rbTriageMode(wx.EVT_RADIOBUTTON)
    imp_dialog.trg_ctr_bragg.SetValue(str(self.gparams.image_triage.min_Bragg_peaks))
    imp_dialog.trg_gs_hmin.SetValue(str(self.gparams.image_triage.grid_search.height_min))
    imp_dialog.trg_gs_hmax.SetValue(str(self.gparams.image_triage.grid_search.height_max))
    imp_dialog.trg_gs_amin.SetValue(str(self.gparams.image_triage.grid_search.area_min))
    imp_dialog.trg_gs_amax.SetValue(str(self.gparams.image_triage.grid_search.area_max))

    # Get values and change params:
    if (imp_dialog.ShowModal() == wx.ID_OK):
      if imp_dialog.conv_rb1_prefix.GetValue():
        self.gparams.image_conversion.rename_pickle_prefix = 'None'
      elif imp_dialog.conv_rb2_prefix.GetValue():
        self.gparams.image_conversion.rename_pickle_prefix = 'Auto'
      elif imp_dialog.conv_rb3_prefix.GetValue():
        self.gparams.image_conversion.rename_pickle_prefix = imp_dialog.conv_ctr_prefix.GetValue()
      if imp_dialog.mod_rb1_square.GetValue():
        self.gparams.image_conversion.square_mode = 'None'
      elif imp_dialog.mod_rb2_square.GetValue():
        self.gparams.image_conversion.square_mode = 'crop'
      elif imp_dialog.mod_rb3_square.GetValue():
        self.gparams.image_conversion.square_mode = 'pad'
      self.gparams.image_conversion.beamstop = float(imp_dialog.mod_ctr_beamstop.GetValue())
      self.gparams.image_conversion.distance = float(imp_dialog.mod_ctr_dist.GetValue())
      self.gparams.image_conversion.beam_center.x = float(imp_dialog.mod_ctr_beamx.GetValue())
      self.gparams.image_conversion.beam_center.y = float(imp_dialog.mod_ctr_beamy.GetValue())
      if imp_dialog.trg_rb1_mode.GetValue():
        self.gparams.image_triage.type = 'None'
      elif imp_dialog.trg_rb2_mode.GetValue():
        self.gparams.image_triage.type = 'simple'
      elif imp_dialog.trg_rb3_mode.GetValue():
        self.gparams.image_triage.type = 'grid_search'
      self.gparams.image_triage.min_Bragg_peaks= int(imp_dialog.trg_ctr_bragg.GetValue())
      self.gparams.image_triage.grid_search.area_min = int(imp_dialog.trg_gs_amin.GetValue())
      self.gparams.image_triage.grid_search.area_max = int(imp_dialog.trg_gs_amax.GetValue())
      self.gparams.image_triage.grid_search.height_min = int(imp_dialog.trg_gs_hmin.GetValue())
      self.gparams.image_triage.grid_search.height_max = int(imp_dialog.trg_gs_hmax.GetValue())

    imp_dialog.Destroy()

  def onProcessOptions(self, e):
    # For cctbx.xfel options
    if self.int_cbx_choice.GetCurrentSelection() == 0:
      int_dialog = CCTBXOptions(self, title='cctbx.xfel Options', style=wx.DEFAULT_DIALOG_STYLE | wx.STAY_ON_TOP)
      int_dialog.Fit()

      # Set values to defaults or previously-selected params
      if str(self.gparams.cctbx.target).lower() == 'none':
        int_dialog.cctbx_ctr_target.SetValue('')
      else:
        int_dialog.cctbx_ctr_target.SetValue(self.gparams.cctbx.target)
      if self.gparams.cctbx.grid_search.type == 'None':
        int_dialog.gs_rb1_type.SetValue(True)
      elif self.gparams.cctbx.grid_search.type == 'brute_force':
        int_dialog.gs_rb2_type.SetValue(True)
      elif self.gparams.cctbx.grid_search.type == 'smart':
        int_dialog.gs_rb3_type.SetValue(True)
      int_dialog.onGSType(wx.EVT_RADIOBUTTON)
      int_dialog.gs_ctr_height.SetValue(str(self.gparams.cctbx.grid_search.height_median))
      int_dialog.gs_ctr_hrange.SetValue(str(self.gparams.cctbx.grid_search.height_range))
      int_dialog.gs_ctr_area.SetValue(str(self.gparams.cctbx.grid_search.area_median))
      int_dialog.gs_ctr_arange.SetValue(str(self.gparams.cctbx.grid_search.area_range))
      int_dialog.gs_chk_sih.SetValue(self.gparams.cctbx.grid_search.sig_height_search)
      int_dialog.sel_chk_only.SetValue(self.gparams.cctbx.selection.select_only.flag_on)
      int_dialog.onSelOnly(wx.EVT_CHECKBOX)
      if str(self.gparams.cctbx.selection.select_only.grid_search_path).lower() == "none":
        int_dialog.sel_ctr_objpath.SetValue('')
      else:
        int_dialog.sel_ctr_objpath.SetValue(str(self.gparams.cctbx.selection.select_only.grid_search_path))
      int_dialog.sel_ctr_minsigma.SetValue(str(self.gparams.cctbx.selection.min_sigma))
      if self.gparams.cctbx.selection.select_by == 'epv':
        int_dialog.sel_rb1_selby.SetValue(True)
      elif self.gparams.cctbx.selection.select_by == 'mosaicity':
        int_dialog.sel_rb1_selby.SetValue(True)
      if self.gparams.cctbx.selection.prefilter.flag_on:
        if str(self.gparams.cctbx.selection.prefilter.target_pointgroup).lower() != 'none':
          int_dialog.sel_chk_lattice.SetValue(True)
          int_dialog.sel_ctr_lattice.Enable()
          int_dialog.sel_ctr_lattice.SetValue(str(self.gparams.cctbx.selection.prefilter.target_pointgroup))
        if str(self.gparams.cctbx.selection.prefilter.target_unit_cell).lower() != 'none':
            t_uc = self.gparams.cctbx.selection.prefilter.target_unit_cell.parameters()
            int_dialog.sel_chk_unitcell.SetValue(True)
            int_dialog.sel_txt_uc_a.Enable()
            int_dialog.sel_ctr_uc_a.Enable()
            int_dialog.sel_ctr_uc_a.SetValue(str(t_uc[0]))
            int_dialog.sel_txt_uc_b.Enable()
            int_dialog.sel_ctr_uc_b.Enable()
            int_dialog.sel_ctr_uc_b.SetValue(str(t_uc[1]))
            int_dialog.sel_txt_uc_c.Enable()
            int_dialog.sel_ctr_uc_c.Enable()
            int_dialog.sel_ctr_uc_c.SetValue(str(t_uc[2]))
            int_dialog.sel_txt_uc_alpha.Enable()
            int_dialog.sel_ctr_uc_alpha.Enable()
            int_dialog.sel_ctr_uc_alpha.SetValue(str(t_uc[3]))
            int_dialog.sel_txt_uc_beta.Enable()
            int_dialog.sel_ctr_uc_beta.Enable()
            int_dialog.sel_ctr_uc_beta.SetValue(str(t_uc[4]))
            int_dialog.sel_txt_uc_gamma.Enable()
            int_dialog.sel_ctr_uc_gamma.Enable()
            int_dialog.sel_ctr_uc_gamma.SetValue(str(t_uc[5]))
            int_dialog.sel_txt_uc_tol.Enable()
            int_dialog.sel_ctr_uc_tol.Enable()
            int_dialog.sel_ctr_uc_tol.SetValue(str(self.gparams.cctbx.selection.prefilter.target_uc_tolerance))
        if str(self.gparams.cctbx.selection.prefilter.min_reflections).lower() != 'none':
          int_dialog.sel_chk_minref.SetValue(True)
          int_dialog.sel_ctr_minref.Enable()
          int_dialog.sel_ctr_minref.SetValue(str(self.gparams.cctbx.selection.prefilter.min_reflections))
        if str(self.gparams.cctbx.selection.prefilter.min_resolution).lower() != 'none':
          int_dialog.sel_chk_res.SetValue(True)
          int_dialog.sel_ctr_res.Enable()
          int_dialog.sel_ctr_res.SetValue(str(self.gparams.cctbx.selection.prefilter.min_resolution))

      # Get values and set parameters
      if (int_dialog.ShowModal() == wx.ID_OK):
        if int_dialog.cctbx_ctr_target.GetValue() == '':
          self.gparams.cctbx.target = None
        else:
          self.gparams.cctbx.target = int_dialog.cctbx_ctr_target.GetValue()
        if int_dialog.gs_rb1_type.GetValue():
          self.gparams.cctbx.grid_search.type = 'None'
        elif int_dialog.gs_rb2_type.GetValue():
          self.gparams.cctbx.grid_search.type = 'brute_force'
        elif int_dialog.gs_rb3_type.GetValue():
          self.gparams.cctbx.grid_search.type = 'smart'
        if int_dialog.gs_ctr_area.GetValue() != '':
          self.gparams.cctbx.grid_search.area_median = int(int_dialog.gs_ctr_area.GetValue())
        if int_dialog.gs_ctr_arange.GetValue() != '':
          self.gparams.cctbx.grid_search.area_range = int(int_dialog.gs_ctr_arange.GetValue())
        if int_dialog.gs_ctr_height.GetValue() != '':
          self.gparams.cctbx.grid_search.height_median = int(int_dialog.gs_ctr_height.GetValue())
        if int_dialog.gs_ctr_hrange.GetValue() != '':
          self.gparams.cctbx.grid_search.height_range = int(int_dialog.gs_ctr_hrange.GetValue())
        self.gparams.cctbx.grid_search.sig_height_search = int_dialog.gs_chk_sih.GetValue()
        self.gparams.cctbx.selection.select_only.flag_on = int_dialog.sel_chk_only.GetValue()
        if int_dialog.sel_ctr_objpath.GetValue() == '':
          self.gparams.cctbx.selection.select_only.grid_search_path = None
        else:
          self.gparams.cctbx.selection.select_only.grid_search_path = int_dialog.sel_ctr_objpath.GetValue()
        if int_dialog.sel_ctr_minsigma.GetValue() != '':
          self.gparams.cctbx.selection.min_sigma = float(int_dialog.sel_ctr_minsigma.GetValue())
        if int_dialog.sel_rb1_selby.GetValue():
          self.gparams.cctbx.selection.select_by = 'epv'
        elif int_dialog.sel_rb2_selby.GetValue():
          self.gparams.cctbx.selection.select_by = 'mosaicity'
        if int_dialog.sel_chk_lattice.GetValue():
          self.gparams.cctbx.selection.prefilter.flag_on = True
          self.gparams.cctbx.selection.prefilter.target_pointgroup = int_dialog.sel_ctr_lattice.GetValue()
        else:
          self.gparams.cctbx.selection.prefilter.target_pointgroup = None
        if (  int_dialog.sel_chk_unitcell.GetValue() and
                int_dialog.sel_ctr_uc_a.GetValue() != '' and
                int_dialog.sel_ctr_uc_b.GetValue() != '' and
                int_dialog.sel_ctr_uc_c.GetValue() != '' and
                int_dialog.sel_ctr_uc_alpha.GetValue() != '' and
                int_dialog.sel_ctr_uc_beta.GetValue() != '' and
                int_dialog.sel_ctr_uc_gamma.GetValue() != ''
             ):
          self.gparams.cctbx.selection.prefilter.flag_on = True
          uc_tuple = (float(int_dialog.sel_ctr_uc_a.GetValue()), float(int_dialog.sel_ctr_uc_b.GetValue()),
                      float(int_dialog.sel_ctr_uc_c.GetValue()), float(int_dialog.sel_ctr_uc_alpha.GetValue()),
                      float(int_dialog.sel_ctr_uc_beta.GetValue()), float(int_dialog.sel_ctr_uc_gamma.GetValue()))
          self.gparams.cctbx.selection.prefilter.target_unit_cell = unit_cell(uc_tuple)
          if int_dialog.sel_ctr_uc_tol.GetValue() != '':
            self.gparams.cctbx.selection.prefilter.target_uc_tolerance = float(int_dialog.sel_ctr_uc_tol.GetValue())
          else:
            self.gparams.cctbx.selection.prefilter.target_uc_tolerance = 0.05
        else:
          self.gparams.cctbx.selection.prefilter.target_unit_cell = None
        if int_dialog.sel_chk_minref.GetValue():
          self.gparams.cctbx.selection.prefilter.flag_on = True
          self.gparams.cctbx.selection.prefilter.min_reflections = int(int_dialog.sel_ctr_minref.GetValue())
        else:
          self.gparams.cctbx.selection.prefilter.min_reflections = None
        if int_dialog.sel_chk_res.GetValue():
          self.gparams.cctbx.selection.prefilter.flag_on = True
          self.gparams.cctbx.selection.prefilter.min_resolution = float(int_dialog.sel_ctr_res.GetValue())
        else:
          self.gparams.cctbx.selection.prefilter.min_resolution = None
        if ( int_dialog.sel_chk_lattice.GetValue() == False and
             int_dialog.sel_chk_unitcell.GetValue() == False and
             int_dialog.sel_chk_minref.GetValue() == False and
             int_dialog.sel_chk_res.GetValue() == False
            ):
          self.gparams.cctbx.selection.prefilter.flag_on = False

    # For DIALS options
    elif self.int_cbx_choice.GetCurrentSelection() == 1:
      int_dialog = DIALSOptions(self, title='DIALS Options', style=wx.DEFAULT_DIALOG_STYLE | wx.STAY_ON_TOP)
      int_dialog.Fit()

      # Set values to defaults or previously-selected params
      int_dialog.d_ctr_smin.SetValue(str(self.gparams.dials.min_spot_size))
      int_dialog.d_ctr_gthr.SetValue(str(self.gparams.dials.global_threshold))

      # Get values and set parameters
      if (int_dialog.ShowModal() == wx.ID_OK):
        if int_dialog.d_ctr_smin.GetValue() != '':
          self.gparams.dials.min_spot_size = int(int_dialog.d_ctr_smin.GetValue())
        if int_dialog.d_ctr_gthr.GetValue() != '':
          self.gparams.dials.global_threshold = int(int_dialog.d_ctr_gthr.GetValue())

    int_dialog.Destroy()

  def onAnalysisOptions(self, e):
    an_dialog = AnalysisWindow(self, title='Dataset Analysis Options', style=wx.DEFAULT_DIALOG_STYLE | wx.STAY_ON_TOP)
    an_dialog.Fit()

    # Set values to defaults or previously-selected params
    an_dialog.an_chk_cluster.SetValue(self.gparams.analysis.run_clustering)
    if an_dialog.an_chk_cluster.GetValue():
      an_dialog.an_ctr_cluster.Enable()
      an_dialog.an_ctr_cluster.SetValue(str(self.gparams.analysis.cluster_threshold))
      an_dialog.an_txt_cluster.Enable()
    if str(self.gparams.analysis.viz).lower() == 'none':
      an_dialog.an_cbx_viz.SetSelection(0)
    elif str(self.gparams.analysis.viz).lower() == 'integration':
      an_dialog.an_cbx_viz.SetSelection(1)
    elif str(self.gparams.analysis.viz).lower() == 'cv_vectors':
      an_dialog.an_cbx_viz.SetSelection(2)
    an_dialog.an_chk_charts.SetValue(self.gparams.analysis.charts)

    # Get values and set parameters
    if (an_dialog.ShowModal() == wx.ID_OK):
      self.gparams.analysis.run_clustering = an_dialog.an_chk_cluster.GetValue()
      self.gparams.analysis.cluster_threshold = float(an_dialog.an_ctr_cluster.GetValue())
      if an_dialog.an_cbx_viz.GetSelection() == 0:
        self.gparams.analysis.viz = 'None'
      elif an_dialog.an_cbx_viz.GetSelection() == 1:
        self.gparams.analysis.viz = 'integration'
      elif an_dialog.an_cbx_viz.GetSelection() == 2:
        self.gparams.analysis.viz = 'cv_vectors'
      self.gparams.analysis.charts = an_dialog.an_chk_charts.GetValue()

    an_dialog.Destroy()

  def onRandomCheck(self, e):
    ''' On random subset check, enable spin control to select size of subset (default = 5) '''
    ischecked = e.GetEventObject().GetValue()
    if ischecked:
      self.opt_spc_random.Enable()
      self.opt_txt_random.Enable()
    else:
      self.opt_spc_random.Disable()
      self.opt_txt_random.Disable()

  def onInfo(self, e):
    ''' On clicking the info button '''
    info_txt = '''Input diffraction images here. IOTA accepts either raw images (mccd, cbf, img, etc.) or image pickles. Input can be either a folder with images, or a text file with a list of images.'''
    info = wx.MessageDialog(None, info_txt, 'Info', wx.OK)
    info.ShowModal()

  def onInputBrowse(self, e):
    """ On clincking the Browse button: show the DirDialog and populate 'Input' box w/ selection """
    dlg = wx.DirDialog(self, "Choose the input directory:", style=wx.DD_DEFAULT_STYLE)
    if dlg.ShowModal() == wx.ID_OK:
      self.inp_ctr.SetValue(dlg.GetPath())
    dlg.Destroy()
    e.Skip()

  def onOutputBrowse(self, e):
    """ On clicking the Browse button: show the DirDialog and populate 'Output' box w/ selection """
    dlg = wx.DirDialog(self, "Choose the output directory:", style=wx.DD_DEFAULT_STYLE)
    if dlg.ShowModal() == wx.ID_OK:
      self.out_ctr.SetValue(dlg.GetPath())
    dlg.Destroy()
    e.Skip()


class ImportWindow(wx.Dialog):
  # Import window - image import, modification and triage

  def __init__(self, *args, **kwargs):
    super(ImportWindow, self).__init__(*args, **kwargs)

    main_sizer = wx.BoxSizer(wx.VERTICAL)

    main_box = wx.StaticBox(self, label='Import')
    vbox = wx.StaticBoxSizer(main_box, wx.VERTICAL)

    # Conversion prefix options
    conv_box_1 = wx.BoxSizer(wx.HORIZONTAL)
    self.conv_txt_prefix = wx.StaticText(self, label='Rename converted images:')
    self.conv_rb1_prefix = wx.RadioButton(self, label='None', style=wx.RB_GROUP)
    self.conv_rb1_prefix.SetValue(False)
    self.conv_rb2_prefix = wx.RadioButton(self, label='Auto')
    self.conv_rb2_prefix.SetValue(True)
    self.conv_rb3_prefix = wx.RadioButton(self, label='Custom:')
    self.conv_rb3_prefix.SetValue(False)
    self.conv_ctr_prefix = wx.TextCtrl(self)
    self.conv_ctr_prefix.Disable()
    conv_box_1.Add(self.conv_txt_prefix, flag=wx.LEFT, border=5)
    conv_box_1.Add(self.conv_rb1_prefix, flag=wx.LEFT, border=5)
    conv_box_1.Add(self.conv_rb2_prefix, flag=wx.LEFT, border=5)
    conv_box_1.Add(self.conv_rb3_prefix, flag=wx.LEFT, border=5)
    conv_box_1.Add(self.conv_ctr_prefix, flag=wx.LEFT, border=5)
    vbox.Add(conv_box_1, flag=wx.RIGHT | wx.LEFT | wx.TOP, border=10)

    # Image modification separator
    mod_box_1 = wx.BoxSizer(wx.HORIZONTAL)
    # self.mod_txt_title = wx.StaticText(self, label='Image modification')
    # line_len = 525 - self.mod_txt_title.GetSize()[0]
    self.mod_stl = wx.StaticLine(self, size=(550, -1))
    # mod_box_1.Add(self.mod_txt_title, flag=wx.LEFT, border=5)
    # mod_box_1.Add((10, -1))
    mod_box_1.Add(self.mod_stl, flag=wx.ALIGN_BOTTOM)
    vbox.Add((-1, 10))
    vbox.Add(mod_box_1, flag=wx.RIGHT | wx.LEFT | wx.TOP | wx.BOTTOM, border=10)

    # Image modification options
    conv_box_2 = wx.BoxSizer(wx.HORIZONTAL)
    self.mod_txt_square = wx.StaticText(self, label='Square image:')
    self.mod_rb1_square = wx.RadioButton(self, label='None', style=wx.RB_GROUP)
    self.mod_rb1_square.SetValue(False)
    self.mod_rb2_square = wx.RadioButton(self, label='Crop')
    self.mod_rb2_square.SetValue(False)
    self.mod_rb3_square = wx.RadioButton(self, label='Pad')
    self.mod_rb3_square.SetValue(True)
    conv_box_2.Add(self.mod_txt_square, flag=wx.LEFT, border=5)
    conv_box_2.Add(self.mod_rb1_square, flag=wx.LEFT, border=5)
    conv_box_2.Add(self.mod_rb2_square, flag=wx.LEFT, border=5)
    conv_box_2.Add(self.mod_rb3_square, flag=wx.LEFT, border=5)
    self.mod_txt_beamstop = wx.StaticText(self, label='Beamstop shadow threshold:')
    self.mod_ctr_beamstop = wx.TextCtrl(self, size=(40, -1))
    conv_box_2.Add(self.mod_txt_beamstop, flag=wx.LEFT, border=25)
    conv_box_2.Add(self.mod_ctr_beamstop, flag=wx.LEFT, border=5)
    vbox.Add(conv_box_2, flag=wx.RIGHT | wx.LEFT | wx.TOP, border=10)

    conv_box_3 = wx.BoxSizer(wx.HORIZONTAL)
    self.mod_txt_mm1 = wx.StaticText(self, label='mm')
    self.mod_txt_dist = wx.StaticText(self, label='Detector distance:')
    self.mod_ctr_dist = wx.TextCtrl(self, size=(50, -1))
    conv_box_3.Add(self.mod_txt_dist, flag=wx.LEFT, border=5)
    conv_box_3.Add(self.mod_ctr_dist, flag=wx.LEFT, border=5)
    conv_box_3.Add(self.mod_txt_mm1, flag=wx.LEFT, border=1)
    self.mod_txt_beamxy = wx.StaticText(self, label='Direct beam')
    self.mod_txt_beamx = wx.StaticText(self, label='X =')
    self.mod_txt_mm2 = wx.StaticText(self, label='mm, ')
    self.mod_ctr_beamx = wx.TextCtrl(self, size=(50, -1))
    self.mod_txt_beamy = wx.StaticText(self, label='Y =')
    self.mod_txt_mm3 = wx.StaticText(self, label='mm')
    self.mod_ctr_beamy = wx.TextCtrl(self, size=(50, -1))
    conv_box_3.Add(self.mod_txt_beamxy, flag=wx.LEFT, border=35)
    conv_box_3.Add(self.mod_txt_beamx, flag=wx.LEFT, border=5)
    conv_box_3.Add(self.mod_ctr_beamx)
    conv_box_3.Add(self.mod_txt_mm2, flag=wx.LEFT, border=1)
    conv_box_3.Add(self.mod_txt_beamy, flag=wx.LEFT, border=5)
    conv_box_3.Add(self.mod_ctr_beamy)
    conv_box_3.Add(self.mod_txt_mm3, flag=wx.LEFT, border=1)
    vbox.Add(conv_box_3, flag=wx.RIGHT | wx.LEFT | wx.TOP, border=10)

    # Image triage separator
    trg_box = wx.BoxSizer(wx.HORIZONTAL)
    # self.trg_txt_title = wx.StaticText(self, label='Image triage')
    # line_len = 525 - self.trg_txt_title.GetSize()[0]
    self.trg_stl = wx.StaticLine(self, size=(550, -1))
    # trg_box.Add(self.trg_txt_title, flag=wx.LEFT, border=5)
    # trg_box.Add((10, -1))
    trg_box.Add(self.trg_stl, flag=wx.ALIGN_BOTTOM)
    vbox.Add((-1, 10))
    vbox.Add(trg_box, flag=wx.RIGHT | wx.LEFT | wx.TOP | wx.BOTTOM, border=10)

    # Image triage options
    conv_box_4 = wx.BoxSizer(wx.HORIZONTAL)
    self.trg_txt_mode = wx.StaticText(self, label='Triage mode:')
    self.trg_rb1_mode = wx.RadioButton(self, label='None', style=wx.RB_GROUP)
    self.trg_rb1_mode.SetValue(False)
    self.trg_rb2_mode = wx.RadioButton(self, label='Simple')
    self.trg_rb2_mode.SetValue(True)
    self.trg_rb3_mode = wx.RadioButton(self, label='Grid Search')
    self.trg_rb3_mode.SetValue(False)
    self.trg_txt_bragg = wx.StaticText(self, label='Minimum Bragg peaks:')
    self.trg_ctr_bragg = wx.TextCtrl(self, size=(40, -1))
    self.trg_ctr_bragg.SetValue('10')
    conv_box_4.Add(self.trg_txt_mode, flag=wx.LEFT, border=5)
    conv_box_4.Add(self.trg_rb1_mode, flag=wx.LEFT, border=5)
    conv_box_4.Add(self.trg_rb2_mode, flag=wx.LEFT, border=5)
    conv_box_4.Add(self.trg_rb3_mode, flag=wx.LEFT, border=5)
    conv_box_4.Add(self.trg_txt_bragg, flag=wx.LEFT, border=15)
    conv_box_4.Add(self.trg_ctr_bragg, flag=wx.LEFT | wx.ALIGN_TOP, border=5)
    vbox.Add(conv_box_4, flag=wx.RIGHT | wx.LEFT | wx.TOP, border=10)

    # Triage grid search box
    conv_box_5 = wx.BoxSizer(wx.HORIZONTAL)
    self.trg_txt_gs = wx.StaticText(self, label='Triage grid search: ')
    self.trg_txt_H = wx.StaticText(self, label='spot height =')
    self.trg_gs_hmin = wx.TextCtrl(self, size=(40, -1))
    self.trg_txt_Hdash = wx.StaticText(self, label='-')
    self.trg_gs_hmax = wx.TextCtrl(self, size=(40, -1))
    self.trg_txt_A = wx.StaticText(self, label=',   spot area =')
    self.trg_gs_amin = wx.TextCtrl(self, size=(40, -1))
    self.trg_txt_Adash = wx.StaticText(self, label='-')
    self.trg_gs_amax = wx.TextCtrl(self, size=(40, -1))
    self.trg_txt_gs.Disable()
    self.trg_txt_H.Disable()
    self.trg_gs_hmin.Disable()
    self.trg_txt_Hdash.Disable()
    self.trg_gs_hmax.Disable()
    self.trg_txt_A.Disable()
    self.trg_gs_amin.Disable()
    self.trg_txt_Adash.Disable()
    self.trg_gs_amax.Disable()
    conv_box_5.Add(self.trg_txt_gs, flag=wx.LEFT, border=5)
    conv_box_5.Add(self.trg_txt_H, flag=wx.LEFT, border=5)
    conv_box_5.Add(self.trg_gs_hmin, flag=wx.LEFT, border=1)
    conv_box_5.Add(self.trg_txt_Hdash, flag=wx.LEFT, border=1)
    conv_box_5.Add(self.trg_gs_hmax, flag=wx.LEFT, border=1)
    conv_box_5.Add(self.trg_txt_A, flag=wx.LEFT, border=1)
    conv_box_5.Add(self.trg_gs_amin, flag=wx.LEFT, border=1)
    conv_box_5.Add(self.trg_txt_Adash, flag=wx.LEFT, border=1)
    conv_box_5.Add(self.trg_gs_amax, flag=wx.LEFT, border=1)
    vbox.Add(conv_box_5, flag=wx.RIGHT | wx.LEFT, border=5)
    vbox.Add((-1, 10))

    # Button bindings
    # Prefix radio buttons
    self.conv_rb1_prefix.Bind(wx.EVT_RADIOBUTTON, self.rbPrefix)
    self.conv_rb2_prefix.Bind(wx.EVT_RADIOBUTTON, self.rbPrefix)
    self.conv_rb3_prefix.Bind(wx.EVT_RADIOBUTTON, self.rbPrefix)

    # Triage radio buttons
    self.trg_rb1_mode.Bind(wx.EVT_RADIOBUTTON, self.rbTriageMode)
    self.trg_rb2_mode.Bind(wx.EVT_RADIOBUTTON, self.rbTriageMode)
    self.trg_rb3_mode.Bind(wx.EVT_RADIOBUTTON, self.rbTriageMode)

    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
    main_sizer.Add(vbox, flag=wx.ALL, border=15)
    main_sizer.Add(dialog_box, flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.LEFT | wx.RIGHT | wx.TOP | wx.BOTTOM, border=10)

    self.SetSizer(main_sizer)

  def rbPrefix(self, e):
    ''' Handles converted pickle prefix radio button '''

    if self.conv_rb1_prefix.GetValue():
      self.conv_ctr_prefix.Disable()
    if self.conv_rb2_prefix.GetValue():
      self.conv_ctr_prefix.Disable()
    if self.conv_rb3_prefix.GetValue():
      self.conv_ctr_prefix.Enable()

  def rbTriageMode(self, e):
    ''' Handles triage mode radio button '''

    if self.trg_rb1_mode.GetValue():
      self.trg_txt_gs.Disable()
      self.trg_txt_H.Disable()
      self.trg_gs_hmin.Disable()
      self.trg_txt_Hdash.Disable()
      self.trg_gs_hmax.Disable()
      self.trg_txt_A.Disable()
      self.trg_gs_amin.Disable()
      self.trg_txt_Adash.Disable()
      self.trg_gs_amax.Disable()
      self.trg_txt_bragg.Disable()
      self.trg_ctr_bragg.Disable()

    if self.trg_rb2_mode.GetValue():
      self.trg_txt_gs.Disable()
      self.trg_txt_H.Disable()
      self.trg_gs_hmin.Disable()
      self.trg_txt_Hdash.Disable()
      self.trg_gs_hmax.Disable()
      self.trg_txt_A.Disable()
      self.trg_gs_amin.Disable()
      self.trg_txt_Adash.Disable()
      self.trg_gs_amax.Disable()
      self.trg_txt_bragg.Enable()
      self.trg_ctr_bragg.Enable()

    if self.trg_rb3_mode.GetValue():
      self.trg_txt_gs.Enable()
      self.trg_txt_H.Enable()
      self.trg_gs_hmin.Enable()
      self.trg_txt_Hdash.Enable()
      self.trg_gs_hmax.Enable()
      self.trg_txt_A.Enable()
      self.trg_gs_amin.Enable()
      self.trg_txt_Adash.Enable()
      self.trg_gs_amax.Enable()
      self.trg_txt_bragg.Enable()
      self.trg_ctr_bragg.Enable()

class CCTBXOptions(wx.Dialog):
  # CCTBX.XFEL options

  def __init__(self, *args, **kwargs):
    super(CCTBXOptions, self).__init__(*args, **kwargs)

    main_sizer = wx.BoxSizer(wx.VERTICAL)

    main_box = wx.StaticBox(self, label='cctbx.xfel Options')
    cctbx_box = wx.StaticBoxSizer(main_box, wx.VERTICAL)

    # Target file input
    target_box = wx.BoxSizer(wx.HORIZONTAL)
    self.cctbx_txt_target = wx.StaticText(self, label='CCTBX.XFEL target file:')
    self.cctbx_btn_target = wx.Button(self, label='Browse...')
    ctr_length = 540 - self.cctbx_btn_target.GetSize()[0] - self.cctbx_txt_target.GetSize()[0]
    self.cctbx_ctr_target = wx.TextCtrl(self, size=(ctr_length, -1))
    target_box.Add(self.cctbx_txt_target)
    target_box.Add(self.cctbx_ctr_target, flag=wx.LEFT, border=5)
    target_box.Add(self.cctbx_btn_target, flag=wx.LEFT, border=5)
    cctbx_box.Add((-1, 25))
    cctbx_box.Add(target_box, flag=wx.LEFT | wx.RIGHT, border=5)

    # Grid search separator
    gs_sep_box = wx.BoxSizer(wx.HORIZONTAL)
    self.gs_stl = wx.StaticLine(self, size=(550, -1))
    gs_sep_box.Add(self.gs_stl, flag=wx.ALIGN_BOTTOM)
    cctbx_box.Add((-1, 10))
    cctbx_box.Add(gs_sep_box, flag=wx.LEFT | wx.TOP | wx.BOTTOM, border=5)

    # Grid search options
    # Type selection
    gs_type_box = wx.BoxSizer(wx.HORIZONTAL)
    gs_txt_type = wx.StaticText(self, label='Grid Search:')
    self.gs_rb1_type = wx.RadioButton(self, label='None', style=wx.RB_GROUP)
    self.gs_rb1_type.SetValue(False)
    self.gs_rb2_type = wx.RadioButton(self, label='Brute Force')
    self.gs_rb2_type.SetValue(True)
    self.gs_rb3_type = wx.RadioButton(self, label='Smart')
    self.gs_rb3_type.SetValue(False)
    self.gs_chk_sih = wx.CheckBox(self, label='Signal height search')
    self.gs_chk_sih.SetValue(False)
    gs_type_box.Add(gs_txt_type)
    gs_type_box.Add(self.gs_rb1_type, flag=wx.LEFT, border=5)
    gs_type_box.Add(self.gs_rb2_type, flag=wx.LEFT, border=5)
    gs_type_box.Add(self.gs_rb3_type, flag=wx.LEFT, border=5)
    gs_type_box.Add(self.gs_chk_sih, flag=wx.LEFT, border=45)
    cctbx_box.Add((-1, 10))
    cctbx_box.Add(gs_type_box, flag=wx.LEFT | wx.BOTTOM, border=5)

    # Grid search parameter options
    gs_prm_box = wx.BoxSizer(wx.HORIZONTAL)
    self.gs_txt_height = wx.StaticText(self, label="Minimum spot height")
    self.gs_ctr_height = wx.TextCtrl(self, size=(50, -1))
    self.gs_txt_hrange = wx.StaticText(self, label="+/-")
    self.gs_ctr_hrange = wx.TextCtrl(self, size=(30, -1))
    self.gs_txt_area = wx.StaticText(self, label="Minimum spot area")
    self.gs_ctr_area = wx.TextCtrl(self, size=(50, -1))
    self.gs_txt_arange = wx.StaticText(self, label="+/-")
    self.gs_ctr_arange = wx.TextCtrl(self, size=(30, -1))
    gs_prm_box.Add(self.gs_txt_height)
    gs_prm_box.Add(self.gs_ctr_height, flag=wx.LEFT, border=5)
    gs_prm_box.Add(self.gs_txt_hrange, flag=wx.LEFT, border=5)
    gs_prm_box.Add(self.gs_ctr_hrange, flag=wx.LEFT, border=1)
    gs_prm_box.Add(self.gs_txt_area, flag=wx.LEFT, border=35)
    gs_prm_box.Add(self.gs_ctr_area, flag=wx.LEFT, border=5)
    gs_prm_box.Add(self.gs_txt_arange, flag=wx.LEFT, border=5)
    gs_prm_box.Add(self.gs_ctr_arange, flag=wx.LEFT, border=1)
    cctbx_box.Add((-1, 10))
    cctbx_box.Add(gs_prm_box, flag=wx.LEFT | wx.TOP | wx.BOTTOM, border=5)

    # Selection separator
    sel_sep_box = wx.BoxSizer(wx.HORIZONTAL)
    self.sel_stl = wx.StaticLine(self, size=(550, -1))
    sel_sep_box.Add(self.sel_stl, flag=wx.ALIGN_BOTTOM)
    cctbx_box.Add((-1, 10))
    cctbx_box.Add(sel_sep_box, flag=wx.LEFT | wx.TOP | wx.BOTTOM, border=5)

    # Selection options
    sel_box = wx.BoxSizer(wx.HORIZONTAL)
    self.sel_chk_only = wx.CheckBox(self, label="Select only")
    self.sel_chk_only.SetValue(False)
    self.sel_txt_objpath = wx.StaticText(self, label='Image objects:')
    self.sel_txt_objpath.Disable()
    self.sel_btn_objpath = wx.Button(self, label='Browse...')
    self.sel_btn_objpath.Disable()
    ctr_length = 525 - self.sel_btn_objpath.GetSize()[0] - self.sel_txt_objpath.GetSize()[0] - \
                 self.sel_chk_only.GetSize()[0]
    self.sel_ctr_objpath = wx.TextCtrl(self, size=(ctr_length, -1))
    self.sel_ctr_objpath.Disable()
    sel_box.Add(self.sel_chk_only)
    sel_box.Add(self.sel_txt_objpath, flag=wx.LEFT, border=15)
    sel_box.Add(self.sel_ctr_objpath, flag=wx.LEFT, border=5)
    sel_box.Add(self.sel_btn_objpath, flag=wx.LEFT, border=5)
    cctbx_box.Add((-1, 10))
    cctbx_box.Add(sel_box, flag=wx.LEFT | wx.RIGHT, border=5)

    sel_box_2 = wx.BoxSizer(wx.HORIZONTAL)
    self.sel_txt_selby = wx.StaticText(self, label='Select by')
    self.sel_rb1_selby = wx.RadioButton(self, label='Ewald proximal volume', style=wx.RB_GROUP)
    self.sel_rb1_selby.SetValue(True)
    self.sel_rb2_selby = wx.RadioButton(self, label='mosaicity')
    self.sel_rb2_selby.SetValue(False)
    self.sel_txt_minsigma = wx.StaticText(self, label='I/SigI cutoff:')
    self.sel_ctr_minsigma = wx.TextCtrl(self, size=(50, -1))
    self.sel_ctr_minsigma.SetValue('5')
    sel_box_2.Add(self.sel_txt_selby)
    sel_box_2.Add(self.sel_rb1_selby, flag=wx.LEFT, border=5)
    sel_box_2.Add(self.sel_rb2_selby, flag=wx.LEFT, border=5)
    sel_box_2.Add(self.sel_txt_minsigma, flag=wx.LEFT, border=35)
    sel_box_2.Add(self.sel_ctr_minsigma, flag=wx.LEFT, border=5)
    cctbx_box.Add((-1, 10))
    cctbx_box.Add(sel_box_2, flag=wx.LEFT | wx.RIGHT, border=5)

    sel_box_3 = wx.BoxSizer(wx.HORIZONTAL)
    self.sel_txt_filter = wx.StaticText(self, label='Filter by:')
    self.sel_chk_lattice = wx.CheckBox(self, label='Bravais lattice: ')
    self.sel_ctr_lattice = wx.TextCtrl(self, size=(100, -1))
    self.sel_ctr_lattice.Disable()
    sel_box_3.Add(self.sel_txt_filter)
    sel_box_3.Add(self.sel_chk_lattice, flag=wx.LEFT, border=5)
    sel_box_3.Add(self.sel_ctr_lattice, flag=wx.LEFT, border=5)
    cctbx_box.Add((-1, 10))
    cctbx_box.Add(sel_box_3, flag=wx.LEFT | wx.RIGHT, border=5)

    sel_box_4 = wx.BoxSizer(wx.HORIZONTAL)
    uc_spacer_1 = self.sel_txt_filter.GetSize()[0] + 5
    self.sel_chk_unitcell = wx.CheckBox(self, label='Unit cell: ')
    self.sel_txt_uc_a = wx.StaticText(self, label='a =')
    self.sel_ctr_uc_a = wx.TextCtrl(self, size=(60, -1))
    self.sel_txt_uc_b = wx.StaticText(self, label='b =')
    self.sel_ctr_uc_b = wx.TextCtrl(self, size=(60, -1))
    self.sel_txt_uc_c = wx.StaticText(self, label='c =')
    self.sel_ctr_uc_c = wx.TextCtrl(self, size=(60, -1))
    self.sel_txt_uc_a.Disable()
    self.sel_ctr_uc_a.Disable()
    self.sel_txt_uc_b.Disable()
    self.sel_ctr_uc_b.Disable()
    self.sel_txt_uc_c.Disable()
    self.sel_ctr_uc_c.Disable()
    sel_box_4.Add(self.sel_chk_unitcell, flag=wx.LEFT, border=uc_spacer_1)
    sel_box_4.Add(self.sel_txt_uc_a, flag=wx.LEFT, border=5)
    sel_box_4.Add(self.sel_ctr_uc_a, flag=wx.LEFT, border=5)
    sel_box_4.Add(self.sel_txt_uc_b, flag=wx.LEFT, border=5)
    sel_box_4.Add(self.sel_ctr_uc_b, flag=wx.LEFT, border=5)
    sel_box_4.Add(self.sel_txt_uc_c, flag=wx.LEFT, border=5)
    sel_box_4.Add(self.sel_ctr_uc_c, flag=wx.LEFT, border=5)

    cctbx_box.Add((-1, 10))
    cctbx_box.Add(sel_box_4, flag=wx.LEFT | wx.RIGHT, border=5)

    sel_box_4a = wx.BoxSizer(wx.HORIZONTAL)
    uc_spacer_2 = uc_spacer_1 + self.sel_chk_unitcell.GetSize()[0] + 5
    self.sel_txt_uc_alpha = wx.StaticText(self, label=u'\N{GREEK SMALL LETTER ALPHA} =')
    self.sel_ctr_uc_alpha = wx.TextCtrl(self, size=(60, -1))
    self.sel_txt_uc_beta = wx.StaticText(self, label=u'\N{GREEK SMALL LETTER BETA} =')
    self.sel_ctr_uc_beta = wx.TextCtrl(self, size=(60, -1))
    self.sel_txt_uc_gamma = wx.StaticText(self, label=u'\N{GREEK SMALL LETTER GAMMA} =')
    self.sel_ctr_uc_gamma = wx.TextCtrl(self, size=(60, -1))
    self.sel_txt_uc_alpha.Disable()
    self.sel_ctr_uc_alpha.Disable()
    self.sel_txt_uc_beta.Disable()
    self.sel_ctr_uc_beta.Disable()
    self.sel_txt_uc_gamma.Disable()
    self.sel_ctr_uc_gamma.Disable()
    sel_box_4a.Add(self.sel_txt_uc_alpha, flag=wx.LEFT, border=uc_spacer_2)
    sel_box_4a.Add(self.sel_ctr_uc_alpha, flag=wx.LEFT, border=5)
    sel_box_4a.Add(self.sel_txt_uc_beta, flag=wx.LEFT, border=5)
    sel_box_4a.Add(self.sel_ctr_uc_beta, flag=wx.LEFT, border=5)
    sel_box_4a.Add(self.sel_txt_uc_gamma, flag=wx.LEFT, border=5)
    sel_box_4a.Add(self.sel_ctr_uc_gamma, flag=wx.LEFT, border=5)
    cctbx_box.Add((-1, 10))
    cctbx_box.Add(sel_box_4a, flag=wx.LEFT | wx.RIGHT, border=5)

    sel_box_4b = wx.BoxSizer(wx.HORIZONTAL)
    self.sel_txt_uc_tol = wx.StaticText(self, label='tolerance:')
    self.sel_ctr_uc_tol = wx.TextCtrl(self, size=(60, -1))
    self.sel_txt_uc_tol.Disable()
    self.sel_ctr_uc_tol.Disable()
    sel_box_4b.Add(self.sel_txt_uc_tol, flag=wx.LEFT, border=uc_spacer_2)
    sel_box_4b.Add(self.sel_ctr_uc_tol, flag=wx.LEFT, border=5)
    cctbx_box.Add((-1, 10))
    cctbx_box.Add(sel_box_4b, flag=wx.LEFT | wx.RIGHT, border=5)

    sel_box_5 = wx.BoxSizer(wx.HORIZONTAL)
    self.sel_chk_minref = wx.CheckBox(self, label='Number of reflections:')
    self.sel_ctr_minref = wx.TextCtrl(self, size=(60, -1))
    self.sel_chk_res = wx.CheckBox(self, label='Resolution')
    self.sel_ctr_res = wx.TextCtrl(self, size=(60, -1))
    self.sel_txt_res = wx.StaticText(self, label=u'\N{ANGSTROM SIGN}')
    sel_box_5.Add(self.sel_chk_minref, flag=wx.LEFT, border=uc_spacer_1)
    sel_box_5.Add(self.sel_ctr_minref, flag=wx.LEFT, border=5)
    sel_box_5.Add(self.sel_chk_res, flag=wx.LEFT, border=15)
    sel_box_5.Add(self.sel_ctr_res, flag=wx.LEFT, border=5)
    sel_box_5.Add(self.sel_txt_res)
    self.sel_ctr_minref.Disable()
    self.sel_ctr_res.Disable()
    cctbx_box.Add((-1, 10))
    cctbx_box.Add(sel_box_5, flag=wx.LEFT | wx.RIGHT, border=5)
    cctbx_box.Add((-1, 10))

    # Button bindings
    self.cctbx_btn_target.Bind(wx.EVT_BUTTON, self.onTargetBrowse)
    self.sel_btn_objpath.Bind(wx.EVT_BUTTON, self.onSelOnlyBrowse)
    self.gs_rb1_type.Bind(wx.EVT_RADIOBUTTON, self.onGSType)
    self.gs_rb2_type.Bind(wx.EVT_RADIOBUTTON, self.onGSType)
    self.gs_rb3_type.Bind(wx.EVT_RADIOBUTTON, self.onGSType)
    self.sel_chk_only.Bind(wx.EVT_CHECKBOX, self.onSelOnly)
    self.sel_chk_lattice.Bind(wx.EVT_CHECKBOX, self.onFilterCheck)
    self.sel_chk_unitcell.Bind(wx.EVT_CHECKBOX, self.onFilterCheck)
    self.sel_chk_minref.Bind(wx.EVT_CHECKBOX, self.onFilterCheck)
    self.sel_chk_res.Bind(wx.EVT_CHECKBOX, self.onFilterCheck)

    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
    main_sizer.Add(cctbx_box, flag=wx.ALL, border=15)
    main_sizer.Add(dialog_box, flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.LEFT | wx.RIGHT | wx.TOP | wx.BOTTOM, border=10)

    self.SetSizer(main_sizer)

  def onFilterCheck(self, e):
    # Controls prefilter options

    if e.GetId() == self.sel_chk_lattice.GetId():
      if self.sel_chk_lattice.GetValue() == True:
        self.sel_ctr_lattice.Enable()
        self.sel_ctr_lattice.SetValue('P4')
      else:
        self.sel_ctr_lattice.Disable()
        self.sel_ctr_lattice.SetValue('')

    if e.GetId() == self.sel_chk_unitcell.GetId():
      if self.sel_chk_unitcell.GetValue() == True:
        self.sel_txt_uc_a.Enable()
        self.sel_ctr_uc_a.Enable()
        self.sel_ctr_uc_a.SetValue('79.4')
        self.sel_txt_uc_b.Enable()
        self.sel_ctr_uc_b.Enable()
        self.sel_ctr_uc_b.SetValue('79.4')
        self.sel_txt_uc_c.Enable()
        self.sel_ctr_uc_c.Enable()
        self.sel_ctr_uc_c.SetValue('38.1')
        self.sel_txt_uc_alpha.Enable()
        self.sel_ctr_uc_alpha.Enable()
        self.sel_ctr_uc_alpha.SetValue('90.0')
        self.sel_txt_uc_beta.Enable()
        self.sel_ctr_uc_beta.Enable()
        self.sel_ctr_uc_beta.SetValue('90.0')
        self.sel_txt_uc_gamma.Enable()
        self.sel_ctr_uc_gamma.Enable()
        self.sel_ctr_uc_gamma.SetValue('90.0')
        self.sel_txt_uc_tol.Enable()
        self.sel_ctr_uc_tol.Enable()
        self.sel_ctr_uc_tol.SetValue('0.05')
      else:
        self.sel_txt_uc_a.Disable()
        self.sel_ctr_uc_a.Disable()
        self.sel_ctr_uc_a.SetValue('')
        self.sel_txt_uc_b.Disable()
        self.sel_ctr_uc_b.Disable()
        self.sel_ctr_uc_b.SetValue('')
        self.sel_txt_uc_c.Disable()
        self.sel_ctr_uc_c.Disable()
        self.sel_ctr_uc_c.SetValue('')
        self.sel_txt_uc_alpha.Disable()
        self.sel_ctr_uc_alpha.Disable()
        self.sel_ctr_uc_alpha.SetValue('')
        self.sel_txt_uc_beta.Disable()
        self.sel_ctr_uc_beta.Disable()
        self.sel_ctr_uc_beta.SetValue('')
        self.sel_txt_uc_gamma.Disable()
        self.sel_ctr_uc_gamma.Disable()
        self.sel_ctr_uc_gamma.SetValue('')
        self.sel_txt_uc_tol.Disable()
        self.sel_ctr_uc_tol.Disable()
        self.sel_ctr_uc_tol.SetValue('')

    if e.GetId() == self.sel_chk_minref.GetId():
      if self.sel_chk_minref.GetValue() == True:
        self.sel_ctr_minref.Enable()
        self.sel_ctr_minref.SetValue('100')
      else:
        self.sel_ctr_minref.Disable()
        self.sel_ctr_minref.SetValue('')

    if e.GetId() == self.sel_chk_res.GetId():
      if self.sel_chk_res.GetValue() == True:
        self.sel_ctr_res.Enable()
        self.sel_ctr_res.SetValue('2.5')
      else:
        self.sel_ctr_res.Disable()
        self.sel_ctr_res.SetValue('')

  def onSelOnly(self, e):
    # Controls selection-only option

    if self.sel_chk_only.GetValue():
      self.sel_txt_objpath.Enable()
      self.sel_ctr_objpath.Enable()
      self.sel_btn_objpath.Enable()
    else:
      self.sel_txt_objpath.Disable()
      self.sel_ctr_objpath.Disable()
      self.sel_btn_objpath.Disable()

  def onSelOnlyBrowse(self, event):
    """ On clincking the Browse button: show the DirDialog and populate 'Input' box w/ selection """
    dlg = wx.DirDialog(self, "Choose directory with image objects:", style=wx.DD_DEFAULT_STYLE)
    if dlg.ShowModal() == wx.ID_OK:
      self.sel_ctr_objpath.SetValue(dlg.GetPath())
    dlg.Destroy()

  def onGSType(self, e):
    # Determines which grid search options are available

    if self.gs_rb1_type.GetValue():
      self.gs_ctr_hrange.SetValue('0')
      self.gs_ctr_hrange.Disable()
      self.gs_ctr_arange.SetValue('0')
      self.gs_ctr_arange.Disable()
      self.gs_chk_sih.Disable()

    if self.gs_rb2_type.GetValue():
      self.gs_ctr_hrange.Enable()
      self.gs_ctr_arange.Enable()
      self.gs_chk_sih.Enable()

    if self.gs_rb3_type.GetValue():
      self.gs_ctr_hrange.SetValue('1')
      self.gs_ctr_hrange.Disable()
      self.gs_ctr_arange.SetValue('1')
      self.gs_ctr_arange.Disable()
      self.gs_chk_sih.Enable()

  def onTargetBrowse(self, e):
    # Opens file dialog for target file browsing

    dlg = wx.FileDialog(
      self, message="Select CCTBX.XFEL target file",
      defaultDir=os.curdir,
      defaultFile="*.phil",
      wildcard="*",
      style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR
    )
    if dlg.ShowModal() == wx.ID_OK:
      filepath = dlg.GetPaths()[0]
      self.cctbx_ctr_target.SetValue(filepath)


class DIALSOptions(wx.Dialog):
  # DIALS options

  def __init__(self, *args, **kwargs):
    super(DIALSOptions, self).__init__(*args, **kwargs)

    main_sizer = wx.BoxSizer(wx.VERTICAL)

    main_box = wx.StaticBox(self, label='DIALS Options')
    dials_box = wx.StaticBoxSizer(main_box, wx.VERTICAL)
    target_box = wx.BoxSizer(wx.HORIZONTAL)

    self.dials_txt_target = wx.StaticText(self, label='DIALS target file:')
    self.dials_btn_target = wx.Button(self, label='Browse...')
    ctr_length = 540 - self.dials_btn_target.GetSize()[0] - self.dials_txt_target.GetSize()[0]
    self.dials_ctr_target = wx.TextCtrl(self, size=(ctr_length, -1))

    target_box.Add(self.dials_txt_target)
    target_box.Add(self.dials_ctr_target, flag=wx.LEFT, border=5)
    target_box.Add(self.dials_btn_target, flag=wx.LEFT, border=5)
    dials_box.Add((-1, 25))
    dials_box.Add(target_box, flag=wx.LEFT | wx.RIGHT, border=5)

    # Options separator
    # opt_sep_box = wx.BoxSizer(wx.HORIZONTAL)
    # opt_stl = wx.StaticLine(self, size=(550, -1))
    # opt_sep_box.Add(opt_stl, flag=wx.ALIGN_BOTTOM)
    # dials_box.Add((-1, 10))
    # dials_box.Add(opt_sep_box, flag=wx.LEFT|wx.TOP|wx.BOTTOM, border=5)

    # DIALS options
    d_opt_box = wx.BoxSizer(wx.HORIZONTAL)
    self.d_txt_smin = wx.StaticText(self, label='Minimum spot size: ')
    self.d_ctr_smin = wx.TextCtrl(self, size=(60, -1))
    self.d_txt_gthr = wx.StaticText(self, label='Global threshold:')
    self.d_ctr_gthr = wx.TextCtrl(self, size=(60, -1))
    d_opt_box.Add(self.d_txt_smin)
    d_opt_box.Add(self.d_ctr_smin, flag=wx.LEFT | wx.RIGHT, border=5)
    d_opt_box.Add(self.d_txt_gthr, flag=wx.LEFT | wx.RIGHT, border=5)
    d_opt_box.Add(self.d_ctr_gthr, flag=wx.LEFT | wx.RIGHT, border=5)
    dials_box.Add((-1, 10))
    dials_box.Add(d_opt_box, flag=wx.LEFT | wx.RIGHT, border=5)
    dials_box.Add((-1, 10))

    # Button bindings
    self.dials_btn_target.Bind(wx.EVT_BUTTON, self.onTargetBrowse)

    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
    main_sizer.Add(dials_box, flag=wx.ALL, border=15)
    main_sizer.Add(dialog_box, flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.LEFT | wx.RIGHT | wx.TOP | wx.BOTTOM, border=10)

    self.SetSizer(main_sizer)

  def onTargetBrowse(self, e):
    # Opens file dialog for target file browsing

    dlg = wx.FileDialog(
      self, message="Select DIALS target file",
      defaultDir=os.curdir,
      defaultFile="*.phil",
      wildcard="*",
      style=wx.OPEN | wx.CHANGE_DIR
    )
    if dlg.ShowModal() == wx.ID_OK:
      filepath = dlg.GetPaths()[0]
      self.dials_ctr_target.SetLabel(filepath)


class AnalysisWindow(wx.Dialog):
  # Import window - image import, modification and triage

  def __init__(self, *args, **kwargs):
    super(AnalysisWindow, self).__init__(*args, **kwargs)

    main_sizer = wx.BoxSizer(wx.VERTICAL)

    main_box = wx.StaticBox(self, label='Analysis Options')
    vbox = wx.StaticBoxSizer(main_box, wx.VERTICAL)

    # Clustering options - DISABLED FOR NOW
    an_box_1 = wx.BoxSizer(wx.HORIZONTAL)
    self.an_chk_cluster = wx.CheckBox(self, label='Unit cell clustering', style=wx.ALIGN_TOP)
    self.an_chk_cluster.Disable()
    self.an_txt_cluster = wx.StaticText(self, label='Threshold:')
    self.an_ctr_cluster = wx.TextCtrl(self, size=(100, -1))
    self.an_txt_cluster.Disable()
    self.an_ctr_cluster.Disable()
    self.an_ctr_cluster.SetValue('5000')
    an_box_1.Add(self.an_chk_cluster, flag=wx.ALIGN_CENTER)
    an_box_1.Add(self.an_txt_cluster, flag=wx.ALIGN_CENTER|wx.LEFT, border=25)
    an_box_1.Add(self.an_ctr_cluster, flag=wx.LEFT, border=5)

    # Visualization options
    an_box_2 = wx.BoxSizer(wx.HORIZONTAL)
    self.an_txt_viz = wx.StaticText(self, label='Visualization: ')
    viz = ['None', 'Integration', 'CV vectors']
    self.an_cbx_viz = wx.Choice(self, choices=viz)
    self.an_chk_charts = wx.CheckBox(self, label='Output processing charts')
    an_box_2.Add(self.an_txt_viz)
    an_box_2.Add(self.an_cbx_viz, flag=wx.ALIGN_CENTER|wx.LEFT, border=5)
    an_box_2.Add(self.an_chk_charts, flag=wx.ALIGN_CENTER|wx.LEFT, border=25)

    vbox.Add(an_box_1, flag=wx.LEFT | wx.RIGHT | wx.TOP | wx.BOTTOM, border=5)
    vbox.Add(an_box_2, flag=wx.LEFT | wx.RIGHT | wx.TOP | wx.BOTTOM, border=5)
    vbox.Add((-1, 10))

    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
    main_sizer.Add(vbox, flag=wx.ALL, border=15)
    main_sizer.Add(dialog_box, flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.LEFT | wx.RIGHT | wx.TOP | wx.BOTTOM, border=10)

    self.SetSizer(main_sizer)

    # Button bindings
    self.an_chk_cluster.Bind(wx.EVT_CHECKBOX, self.onCluster)

  def onCluster(self, e):
    if self.an_chk_cluster.GetValue():
      self.an_ctr_cluster.Enable()
      self.an_txt_cluster.Enable()
    else:
      self.an_ctr_cluster.Disable()
      self.an_txt_cluster.Disable()

class ProgressWindow(wx.Panel):
  # Output window - will show gauges for every process

  def __init__(self, parent):
    super(ProgressWindow, self).__init__(parent)

    progress_box = wx.StaticBox(self, label='Progress', pos=(5, 5), size=(600, 300))
    vbox = wx.StaticBoxSizer(progress_box, wx.VERTICAL)

    self.SetSizer(vbox)


class LogWindow(wx.Panel):
  # Output text log - not yet sure if this is necessary

  def __init__(self, parent):
    super(LogWindow, self).__init__(parent)

    log_box = wx.StaticBox(self, label='Log', pos=(5, 5), size=(600, 600))
    vbox = wx.StaticBoxSizer(log_box, wx.VERTICAL)

    self.SetSizer(vbox)


# ------------------------------------- Main Frame ------------------------------------- #

class MainWindow(wx.Frame):
  def __init__(self, parent, id, title):
    wx.Frame.__init__(self, parent, id, title, pos=wx.DefaultPosition, size=(800, 500))

    # font = wx.Font(12, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL, True)
    # self.SetFont(font)
    self.InitUI()
    self.Show()

  def InitUI(self):
    # Menu bar
    menubar = wx.MenuBar()

    # Status bar
    self.sb = self.CreateStatusBar()

    # Help menu item with the about dialog
    help = wx.Menu()
    help.Append(100, '&About')
    self.Bind(wx.EVT_MENU, self.OnAboutBox, id=100)
    menubar.Append(help, '&Help')
    self.SetMenuBar(menubar)

    # main_grid = wx.GridBagSizer(15, 15)
    # spacer_left = wx.BoxSizer(wx.HORIZONTAL)
    # spacer_left.Add((600, -1))
    # spacer_right = wx.BoxSizer(wx.HORIZONTAL)
    # spacer_right.Add((600, -1))

    main_box = wx.BoxSizer(wx.VERTICAL)

    # Toolbar
    self.toolbar = self.CreateToolBar(wx.TB_TEXT)
    self.texit = self.toolbar.AddLabelTool(wx.ID_EXIT, label='Quit',
                                           bitmap=wx.Bitmap('{}/32x32/exit.png'.format(icons)),
                                           shortHelp='Quit',
                                           longHelp='Quit IOTA')

    self.toolbar.AddSeparator()
    self.tb_btn_run = self.toolbar.AddLabelTool(wx.ID_ANY, label='Run',
                                                bitmap=wx.Bitmap('{}/32x32/run.png'.format(icons)),
                                                shortHelp='Run', longHelp='Run all stages of refinement')
    self.tb_btn_stop = self.toolbar.AddLabelTool(wx.ID_ANY, label='Abort',
                                                 bitmap=wx.Bitmap('{}/32x32/stop.png'.format(icons)),
                                                 shortHelp='Abort', longHelp='Abort all running processes')
    self.toolbar.EnableTool(self.tb_btn_run.GetId(), False)
    self.toolbar.EnableTool(self.tb_btn_stop.GetId(), False)
    self.toolbar.AddSeparator()
    self.tb_btn_script = self.toolbar.AddLabelTool(wx.ID_ANY, label='Output Script',
                                                   bitmap=wx.Bitmap('{}/32x32/write.png'.format(icons)),
                                                   shortHelp='Output Script',
                                                   longHelp='Output current settings to script file')
    self.tb_btn_imglist = self.toolbar.AddLabelTool(wx.ID_ANY,
                                                    label='Write Image List',
                                                    bitmap=wx.Bitmap('{}/32x32/list.png'.format(icons)),
                                                    shortHelp='Write Image List',
                                                    longHelp='Collect list of raw images and output to file')
    self.tb_btn_convert = self.toolbar.AddLabelTool(wx.ID_ANY,
                                                    label='Convert Images',
                                                    bitmap=wx.Bitmap('{}/32x32/convert.png'.format(icons)),
                                                    shortHelp='Convert Image Files',
                                                    longHelp='Collect list of raw images and output to file')
    self.toolbar.Realize()

    # Instantiate windows
    self.input_window = InputWindow(self)

    # self.import_window = ImportWindow(self)
    # self.process_window = ProcessWindow(self)
    # self.analysis_window = AnalysisWindow(self)
    # self.progress_window = ProgressWindow(self)
    # self.log_window = LogWindow(self)

    # Single input window
    main_box.Add(self.input_window, flag=wx.LEFT | wx.RIGHT | wx.TOP | wx.BOTTOM, border=10)
    main_box.Add((-1, 20))

    # Input window - data input, description of project
    # main_grid.Add(spacer_left, pos=(0, 0), flag=wx.EXPAND | wx.LEFT, border=15)
    # main_grid.Add(self.input_window, pos=(1, 0), flag=wx.EXPAND | wx.LEFT, border=15)
    # main_grid.Add(self.import_window, pos=(2, 0), flag=wx.EXPAND | wx.LEFT, border=15)
    # main_grid.Add(self.process_window, pos=(3, 0), flag=wx.EXPAND | wx.LEFT, border=15)
    # main_grid.Add(self.analysis_window, pos=(4, 0), flag=wx.EXPAND | wx.LEFT, border=15)
    #
    # main_grid.Add(spacer_right, pos=(0, 1), flag=wx.EXPAND | wx.LEFT, border=15)
    # main_grid.Add(self.progress_window, pos=(1, 1), span=(2, 1), flag=wx.EXPAND | wx.LEFT, border=15)
    # main_grid.Add(self.log_window, pos=(3, 1), span=(2, 1), flag=wx.EXPAND | wx.LEFT, border=15)
    # self.SetSizer(main_grid)

    self.SetSizer(main_box)

    # Toolbar button bindings
    self.Bind(wx.EVT_TOOL, self.onQuit, self.texit)
    self.Bind(wx.EVT_BUTTON, self.onInput, self.input_window.inp_btn_browse)

  def OnAboutBox(self, e):
    ''' About dialog '''
    info = wx.AboutDialogInfo()
    # info.SetIcon(wx.Icon('icons/iota_100.png', wx.BITMAP_TYPE_PNG))
    info.SetName('IOTA')
    info.SetVersion(iota_version)
    info.SetDescription(description)
    info.SetCopyright('(C) 2015 - 2016 Art Lyubimov')
    info.SetWebSite('http://cci.lbl.gov/xfel')
    info.SetLicense(license)
    info.AddDeveloper('Art Lyubimov')
    info.AddDocWriter('Art Lyubimov')
    info.AddArtist('The Internet')
    info.AddTranslator('Art Lyubimov')
    wx.AboutBox(info)

  def onInput(self, e):
    if self.input_window.inp_ctr.GetValue() != '':
      self.toolbar.EnableTool(self.tb_btn_run.GetId(), True)
    else:
      self.toolbar.EnableTool(self.tb_btn_run.GetId(), False)

  def onQuit(self, e):
    self.Close()
