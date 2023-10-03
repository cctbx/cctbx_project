# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# $Id$

from __future__ import absolute_import, division, print_function
from six.moves import range

import wx
from scitbx.matrix import col


class SBSettingsFrame(wx.MiniFrame):
  def __init__(self, *args, **kwds):
    super(SBSettingsFrame, self).__init__(*args, **kwds)
    szr = wx.BoxSizer(wx.VERTICAL)
    panel = SBSettingsPanel(self)
    self.SetSizer(szr)
    szr.Add(panel, 1, wx.EXPAND)
    szr.Fit(panel)
    self.panel = panel
    self.sizer = szr
    self.Fit()
    self.Bind(wx.EVT_CLOSE, lambda evt : self.Destroy(), self)

  # XXX Could have a set_image() function instead of referring back to
  # the frame all the time?


class SBSettingsPanel(wx.Panel):
  # XXX Names: they're not really settings.  XXX Allow for setting
  # rotation, and provide a hierarchical drop-down menu to play with
  # detector, panel, sensor and ASIC.

  def __init__(self, *args, **kwds):
    super(SBSettingsPanel, self).__init__(*args, **kwds)
    sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(sizer)

    # Number of decimal digits for distances.
    self.digits = 2

    # Quad translation controls
    from wx.lib.agw.floatspin import EVT_FLOATSPIN, FloatSpin

    img = self.GetParent().GetParent().pyslip.tiles.raw_image
    d = img.get_detector()
    self._quad_spinners = []
    for serial in range(4):
      fast, slow = d.hierarchy()[serial].get_origin()[0:2]
      name_quadrant = ["Q0", "Q1", "Q2", "Q3"][serial]
      box = wx.BoxSizer(wx.HORIZONTAL)

      for (name_direction, value) in [("fast", fast), ("slow", slow)]:
        name_ctrl = name_quadrant + "_" + name_direction + "_ctrl"

        spinner = FloatSpin(
          self, digits=self.digits, name=name_ctrl, value=value)
        self.Bind(EVT_FLOATSPIN, self.OnUpdateQuad, spinner)

        box.Add(spinner,
                0, wx.RIGHT|wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_VERTICAL, 5)
        box.Add(wx.StaticText(self, label=name_quadrant + " " + name_direction),
                0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

        setattr(self, "_" + name_ctrl, spinner)
        self._quad_spinners.append(spinner)

      sizer.Add(box)

    #Spinner amount control
    box = wx.BoxSizer(wx.HORIZONTAL)
    self._spinner_amt_control = FloatSpin(
      self, digits=self.digits, name="spin_amount", value=1, min_val= 0.1, increment=0.1)
    self.Bind(EVT_FLOATSPIN, self.OnSpinAmount, self._spinner_amt_control)
    box.Add(self._spinner_amt_control,
            0, wx.RIGHT|wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(wx.StaticText(self, label="Spinner increment (mm)"),
            0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    sizer.Add(box)

    box = wx.BoxSizer(wx.HORIZONTAL)

    btn = wx.Button(self, label="Restore metrology")
    box.Add(btn, flag=wx.ALL, border=5)
    self.Bind(wx.EVT_BUTTON, self.OnRestoreMetrology, btn)

    btn = wx.Button(self, label="Save current metrology")
    box.Add(btn, flag=wx.ALL, border=5)
    self.Bind(wx.EVT_BUTTON, self.OnSaveMetrology, btn)

    sizer.Add(box, flag=wx.ALIGN_CENTER)

    # XXX Rename to metrology tool?


  def OnRestoreMetrology(self, event):
    print("Not implemented")
    return

    dialog = wx.FileDialog(
      self,
      defaultDir="",
      message="Restore metrology file",
      style=wx.FD_OPEN,
      wildcard="Phil files (*.eff; *.def)|*.eff;*.def")
    if dialog.ShowModal() == wx.ID_OK:
      path = dialog.GetPath()
      if (path != ""):
        from serialtbx.detector.legacy_metrology.metrology import \
          master_phil, metrology_as_transformation_matrices
        from libtbx import phil

        frame = self.GetParent().GetParent()
        stream = open(path)
        metrology_phil = master_phil.fetch(sources=[phil.parse(stream.read())])
        stream.close()

        # Merge restored metrology into the raw image
        from libtbx.phil import experimental
        experimental.merge_params_by_key(
          frame.pyslip.tiles.raw_image._metrology_params,
          metrology_phil.extract(),
          'serial')

        img = frame.pyslip.tiles.raw_image
        img.apply_metrology_from_matrices(metrology_as_transformation_matrices(
          metrology_phil.extract()))

        # Update the view, trigger redraw.  XXX Duplication
        # w.r.t. OnUpdateQuad().
        tiles = frame.pyslip.tiles
        tiles.flex_image = frame.pyslip.tiles.raw_image.get_flex_image(
          brightness=tiles.current_brightness / 100)
        tiles.flex_image.adjust(color_scheme=tiles.current_color_scheme)

        tiles.reset_the_cache()
        tiles.tile_cache = tiles.cache[tiles.zoom_level]
        tiles.tile_list = tiles.lru[tiles.zoom_level]
        frame.pyslip.Update()

        # Update the controls, remember to reset the default values
        # for the spinners.
        for serial in range(4):
          fast, slow = img.get_panel_fast_slow(serial)
          name_quadrant = ["Q0", "Q1", "Q2", "Q3"][serial]

          spinner = getattr(self, "_" + name_quadrant + "_fast_ctrl")
          spinner.SetDefaultValue(fast)
          spinner.SetValue(fast)

          spinner = getattr(self, "_" + name_quadrant + "_slow_ctrl")
          spinner.SetDefaultValue(slow)
          spinner.SetValue(slow)


  def OnSaveMetrology(self, event):
    import pycbf, os

    dialog = wx.FileDialog(
      self,
      defaultDir=os.curdir,
      defaultFile="quadrants.def",
      message="Save metrology file",
      style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT,
      wildcard="Phil files (*.def)|*.def")
    if dialog.ShowModal() == wx.ID_OK:
      path = str(dialog.GetPath())
      if (path != ""):
        # The detector object of the format instance is adjusted when the quadrant calibration
        # arrows are clicked.  Sync those adjustments to a new cbf handle, drop uneeded categories
        # (categories frame specific but not metrology specific) and write the file.
        frame = self.GetParent().GetParent()
        img = frame.pyslip.tiles.raw_image
        header=img.image_set.get_format_class()(img.full_path)
        header.sync_detector_to_cbf(img.get_detector())
        cbf = header._cbf_handle
        cbf.find_category("array_data")                  ; cbf.remove_category()
        cbf.find_category("array_structure")             ; cbf.remove_category()
        cbf.find_category("array_intensities")           ; cbf.remove_category()
        cbf.find_category("diffrn_radiation")            ; cbf.remove_category()
        cbf.find_category("diffrn_radiation_wavelength") ; cbf.remove_category()
        cbf.find_category("diffrn_measurement")          ; cbf.remove_category()
        cbf.find_category("diffrn_scan")                 ; cbf.remove_category()
        cbf.find_category("diffrn_scan_frame")           ; cbf.remove_category()

        cbf.write_widefile(path,pycbf.CBF,\
              pycbf.MIME_HEADERS|pycbf.MSG_DIGEST|pycbf.PAD_4K,0)

        print("Saved cbf header to", path)


  def OnUpdateQuad(self, event):
    # Get the name of the spinner and its delta, the deviation from
    # the default value.  Update the default for the next event.
    obj = event.EventObject
    name = obj.GetName()
    value = obj.GetValue()
    delta = float(value - obj.GetDefaultValue())
    obj.SetDefaultValue(value)

    # Update the frame's effective metrology parameters.
    frame = self.GetParent().GetParent()
    img = frame.pyslip.tiles.raw_image

    quads = img.get_detector().hierarchy()

    if   (name == "Q0_fast_ctrl"): quad, delta = (quads[0], col((delta,0,0)))
    elif (name == "Q0_slow_ctrl"): quad, delta = (quads[0], col((0,delta,0)))
    elif (name == "Q1_fast_ctrl"): quad, delta = (quads[1], col((delta,0,0)))
    elif (name == "Q1_slow_ctrl"): quad, delta = (quads[1], col((0,delta,0)))
    elif (name == "Q2_fast_ctrl"): quad, delta = (quads[2], col((delta,0,0)))
    elif (name == "Q2_slow_ctrl"): quad, delta = (quads[2], col((0,delta,0)))
    elif (name == "Q3_fast_ctrl"): quad, delta = (quads[3], col((delta,0,0)))
    elif (name == "Q3_slow_ctrl"): quad, delta = (quads[3], col((0,delta,0)))
    else:
      raise RuntimeError("Unknown control name " + name)

    ldm = quad.get_local_d_matrix()
    fast = (ldm[0],ldm[3],ldm[6])
    slow = (ldm[1],ldm[4],ldm[7])

    orig = col((ldm[2],ldm[5],ldm[8])) + delta

    quad.set_local_frame(fast,slow,orig)

    # Update the view, trigger redraw.
    tiles = frame.pyslip.tiles
    tiles.set_image(tiles.raw_image)
    tiles.flex_image.adjust(color_scheme=tiles.current_color_scheme)

    tiles.reset_the_cache()
    tiles.tile_cache = tiles.cache[tiles.zoom_level]
    tiles.tile_list = tiles.lru[tiles.zoom_level]
    frame.pyslip.Update()

  def OnSpinAmount(self, event):
    obj = event.EventObject
    for spinner in self._quad_spinners:
      spinner.SetIncrement(obj.GetValue())
