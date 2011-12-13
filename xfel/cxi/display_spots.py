# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8; -*-

import os, sys, copy
from iotbx.detectors.npy import NpyImage

class empty:
  pass

parameters=empty()

def ImageFactory(filename):
  global parameters
  if type(filename)==type("") and os.path.isfile(filename):
    I = NpyImage(filename)
  else:
    print "This is not a file; assume the data are in the defined dictionary format"
    I = NpyImage(filename, source_data=parameters.horizons_phil.indexing.data)
  I.readHeader(parameters.horizons_phil)

  I.translate_tiles(parameters.horizons_phil)

  #from xfel.cxi.display_powder_arcs import apply_gaussian_noise
  #apply_gaussian_noise(I,parameters.horizons_phil)
  if (parameters.horizons_phil.viewer.powder_arcs.show):
    from xfel.cxi.display_powder_arcs import superimpose_powder_arcs
    superimpose_powder_arcs(I,parameters.horizons_phil)
  return I


from iotbx import detectors
detectors.ImageFactory = ImageFactory

from spotfinder.diffraction.imagefiles import FileName
FileName.exts.append("pickle")

from spotfinder.applications.xfel import cxi_phil

show_bg_boxes=True
class wrapper_of_callback(object):
  def __init__(self, info=None):
    if info != None:
      self.spotfinder = info.S
      self.frames     = info.frames
      self.Files      = info.Files
      self.phil_params= info.phil_params

  def display(self,path):
    import wx
    from rstbx.viewer.frame import XrayFrame

    app   = wx.App(0)
    frame = XrayFrame(None, -1, "X-ray image display", size=(1200,1080))
    frame.SetSize((1024,780))
    frame.load_image(path)
    frame.Show()
    app.MainLoop()

  def display_with_callback(self,path):
    from rstbx.viewer       import display
    display.user_callback = self.user_callback
    self.display(path)

  def user_callback(self,dc,wxpanel,wx):
    # arguments are a wx Device Context, an Xray Frame, and the wx Module itself
    """
    for spot in self.spotfinder.images[self.frames[0]]["spots_total"]:
      for pxl in spot.bodypixels:
        x,y = wxpanel._img.image_coords_as_screen_coords(
          pxl.y,
          pxl.x)
        dc.SetPen(wx.Pen('red'))
        dc.SetBrush(wx.RED_BRUSH)
        dc.DrawCircle(x,y,1)

      x,y = wxpanel._img.image_coords_as_screen_coords(
        spot.ctr_mass_y(),
        spot.ctr_mass_x())
      dc.SetPen(wx.Pen('green'))
      dc.SetBrush(wx.GREEN_BRUSH)
      dc.DrawCircle(x,y,1)
    """
    """
    for spot in self.spotfinder.images[self.frames[0]]["lo_pass_resolution_spots"]:
      for pxl in spot.bodypixels:
        x,y = wxpanel._img.image_coords_as_screen_coords(
          pxl.y,
          pxl.x)
        dc.SetPen(wx.Pen('yellow'))
        dc.DrawCircle(x,y,1)
    """
    for spot in self.spotfinder.images[self.frames[0]]["spots_inlier"]:
      for pxl in spot.bodypixels:
        x,y = wxpanel._img.image_coords_as_screen_coords(
          pxl.y,
          pxl.x)
        dc.SetPen(wx.Pen('cyan'))
        dc.DrawCircle(x,y,1)

    if  show_bg_boxes:
      imgobj = self.Files.imageindex(self.frames[0])
      from iotbx.detectors.npy import tile_manager
      aa = tile_manager(self.phil_params).effective_tiling_as_flex_int(
        reapply_peripheral_margin = True, beam = (
        imgobj.beamx/imgobj.pixel_size,imgobj.beamy/imgobj.pixel_size))

      dc.SetPen(wx.Pen('orange'))
      dc.SetBrush(wx.Brush('red', wx.TRANSPARENT))
      for i in xrange(0, len(aa), 4):
        p = wxpanel._img.image_coords_as_screen_coords(aa[i + 1], aa[i + 0])
        p2 = wxpanel._img.image_coords_as_screen_coords(aa[i + 3]-1, aa[i + 2]-1)
        dc.DrawRectangle(x=p[0], y=p[1], width=p2[0]-p[0],height=p2[1]-p[1])

    if False:
      x, y = wxpanel._img.image_coords_as_screen_coords(850, 850)

      dc.SetPen(wx.Pen('red'))
      dc.SetBrush(wx.Brush('red', wx.TRANSPARENT))
      for i in range(24):
        r, s = wxpanel._img.image_coords_as_screen_coords(850 + (i + 1) * 25, 850)
        dc.DrawCircle(x, y, r  - x)

def view_raw_image(path, *command_line, **kwargs):
  args = [path,
          "distl.detector_format_version=CXI 5.1",
          "viewer.powder_arcs.show=False",
          "viewer.powder_arcs.code=3n9c",
         ]

  horizons_phil = cxi_phil.cxi_versioned_extract(
                    copy.deepcopy(args),list(command_line))

  global parameters
  parameters.horizons_phil = horizons_phil

  if horizons_phil.viewer.calibrate_silver==True:
    from rstbx.viewer.calibration import sb_wrapper
    sb_wrapper(horizons_phil).display(path)
    return

  wrapper_of_callback().display(path)

def run_one(path, *command_line, **kwargs):
  args = ["distl.image=%s"%path,
          "distl.res.outer=2.1",
          "distl.detector_format_version=CXI 5.1",
          ]

  horizons_phil = cxi_phil.cxi_versioned_extract(
                    copy.deepcopy(args),list(command_line))

  global parameters
  parameters.horizons_phil = horizons_phil

  from spotfinder.applications import signal_strength
  info = signal_strength.run_signal_strength(horizons_phil)

  if kwargs.get("display",False):

    work = wrapper_of_callback(info)
    work.display_with_callback(path)

def run_one_index_core(horizons_phil):
  global parameters
  parameters.horizons_phil = horizons_phil

  from rstbx.new_horizons.index import pre_indexing_validation,pack_names,new_horizons_state
  pre_indexing_validation(horizons_phil)
  imagefile_arguments = pack_names(horizons_phil)
  info = new_horizons_state(horizons_phil,imagefile_arguments)

  info.process()

  info.S = info.spotfinder_results
  return info

def run_one_index(path, *arguments, **kwargs):

  assert arguments[0].find("target=")==0
  target = arguments[0].split("=")[1]
  import xfel_targets

  args = ["indexing.data=%s"%path,
          "distl.detector_format_version=CXI 5.1",
          "beam_search_scope=0.5",
          "lepage_max_delta = 3.0",
          "spots_pickle = None",
          "subgroups_pickle = None",
          "refinements_pickle = None",
          "rmsd_tolerance = 5.0",
          "mosflm_rmsd_tolerance = 5.0",
          "difflimit_sigma_cutoff=2.0",
          #"indexing.verbose_cv=True",
          "indexing.open_wx_viewer=True"
          ] + targets[target] + list(arguments[1:])

  horizons_phil = cxi_phil.cxi_versioned_extract(
                    copy.deepcopy(args))

  info = run_one_index_core(horizons_phil)
  info.Files = info.organizer.Files
  info.phil_params = info.horizons_phil
  work = wrapper_of_callback(info)

  if kwargs.get("display",False):
      import wx
      from rstbx.viewer       import display
      from rstbx.viewer.frame import XrayFrame
      display.user_callback = work.user_callback

      app   = wx.App(0)
      frame = XrayFrame(None, -1, "X-ray image display", size=(1200,1080))
      frame.SetSize((1024,780))
      frame.load_image(path)
      frame.Show()
      app.MainLoop()

if __name__ == "__main__":
  function = sys.argv[1] # either view, spots, or index
  selector = {"view":view_raw_image, "spots":run_one, "index":run_one_index}
  files = [arg for arg in sys.argv[2:] if os.path.isfile(arg)]
  arguments = [arg for arg in sys.argv[2:] if not os.path.isfile(arg)]
  for file in files:
    selector[function](file, *arguments, **({'display':True}))
