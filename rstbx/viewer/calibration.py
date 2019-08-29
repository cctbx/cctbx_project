from __future__ import absolute_import, division, print_function
from six.moves import range
from scitbx.array_family import flex
import math
from six.moves import cPickle as pickle

class sb_wrapper:
  def __init__(self,working_phil):
      self.working_phil = working_phil
      self.two_theta_experimental = flex.double([
        1.5114,3.0160,4.5289,6.0422,7.5575,9.0733,10.5922,12.1106,13.6320,15.1563,16.6836,18.2107])
        # first 12 orders of silver behenate rings reported by Huang (1993) J Appl Cryst 26, 180-184.
      copper_Kalpha = 1.5418 # Angstroms
      d = flex.double(len(self.two_theta_experimental), copper_Kalpha/2.)
      self.experimental_d = d/flex.sin((math.pi/360.)*self.two_theta_experimental)
      # should correspond to repeat spacing of 58.7 Angstrom according to Huang.

  def display(self,path):
    import wx
    from rstbx.viewer.calibration_frame import SBFrame
    from rstbx.viewer import display
    display.user_callback  = self.user_callback

    app   = wx.App(0)
    frame = SBFrame(None, -1, "X-ray image display", size=(1200,1080),
      pos=(100,100), horizons_phil=self.working_phil)
    frame.SetSize((1024,780))
    frame.load_image(path)
    frame.path = path
    self.path = path
    frame.settings_frame.panel.distance_ctrl.SetFloat(frame._img._raw.distance)
    frame.settings_frame.panel.distance_ctrl.DoSendEvent()
    frame.Show()
    app.MainLoop()

  def user_callback(self,dc,panel,wx):
    center_x, center_y = panel._img.get_beam_center()
    xc, yc = panel._img.image_coords_as_screen_coords(center_x, center_y)
    dc.SetPen(wx.Pen('red'))
    dc.SetBrush(wx.TRANSPARENT_BRUSH)

    wavelength = panel._img._raw.wavelength #should be this
    wavelength_from_avg_file = False
    if wavelength_from_avg_file:
      # go through hoops to get the proper wavelength corresponding to this run
      avepath = self.path.replace("stddev","avg")
      info = pickle.load(open(avepath,"rb"))
      wavelength = info["WAVELENGTH"]

    twotheta = 2.* flex.asin(
                 flex.double(len(self.two_theta_experimental), wavelength/2.)/
                 self.experimental_d)
    L_mm = panel.settings.distance * flex.tan(twotheta)
    L_pixels = L_mm / panel._img._raw.pixel_size

    [ dc.DrawCircle(xc, yc, panel._img.get_scale() * pxl) for pxl in L_pixels ]

class pdb_code_wrapper(sb_wrapper):

  def __init__(self,working_phil):
      self.working_phil = working_phil
      # should corresponds to low angle scattering of crystal structure---up to 20 Angstroms
      from xfel.cxi.display_powder_arcs import get_mmtbx_icalc
      intensities = get_mmtbx_icalc(
        code = working_phil.viewer.calibrate_pdb.code,
        d_min = working_phil.viewer.calibrate_pdb.d_min,
        anomalous_flag=False)
      self.hkl_list = intensities.indices()
      self.uc = intensities.unit_cell()
      spacings = self.uc.d(self.hkl_list)
      rev_order = flex.sort_permutation(spacings,reverse = True)
      for x in range(len(rev_order)):
        print(self.hkl_list[rev_order[x]], spacings[rev_order[x]])
      self.experimental_d =  spacings.select(rev_order)

  def user_callback(self,dc,panel,wx):
    center_x, center_y = panel._img.get_beam_center()
    xc, yc = panel._img.image_coords_as_screen_coords(center_x, center_y)
    dc.SetPen(wx.Pen('red'))
    dc.SetBrush(wx.TRANSPARENT_BRUSH)

    wavelength = panel._img._raw.wavelength #should be this
    wavelength_from_avg_file = False
    if wavelength_from_avg_file:
      # go through hoops to get the proper wavelength corresponding to this run
      avepath = self.path.replace("stddev","avg")
      info = pickle.load(open(avepath,"rb"))
      wavelength = info["WAVELENGTH"]

    twotheta = self.uc.two_theta(miller_indices = self.hkl_list, wavelength = wavelength)
    L_mm = panel.settings.distance * flex.tan(twotheta)
    L_pixels = L_mm / panel._img._raw.pixel_size

    [ dc.DrawCircle(xc, yc, panel._img.get_scale() * pxl) for pxl in L_pixels ]

class unit_cell_wrapper(sb_wrapper):

  def __init__(self,working_phil):
    self.working_phil = working_phil

    from cctbx.crystal import symmetry
    import cctbx.miller
    self.uc = symmetry(unit_cell=self.working_phil.viewer.calibrate_unitcell.unitcell,
                       space_group_symbol=self.working_phil.viewer.calibrate_unitcell.spacegroup)
    self.hkl_list = cctbx.miller.build_set(self.uc, False, d_min=working_phil.viewer.calibrate_unitcell.d_min)

    spacings = list(self.hkl_list.d_spacings())
    print("Printing spacings, len: %s"%len(spacings))

    def cmp(a,b):
      if a[1] > b[1]: return 1
      elif a[1] < b[1]: return -1
      return 0

    spacings = sorted(spacings, cmp=cmp, reverse=True)

    for d in spacings:
      print(d)

  def user_callback(self,dc,panel,wx):
    if not hasattr(panel.settings, "distance"): return # fixes a crash on exit

    center_x, center_y = panel._img.get_beam_center()
    xc, yc = panel._img.image_coords_as_screen_coords(center_x, center_y)
    dc.SetPen(wx.Pen('red'))
    dc.SetBrush(wx.TRANSPARENT_BRUSH)

    wavelength = panel._img._raw.wavelength #should be this
    wavelength_from_avg_file = False
    if wavelength_from_avg_file:
      # go through hoops to get the proper wavelength corresponding to this run
      avepath = self.path.replace("stddev","avg")
      info = pickle.load(open(avepath,"rb"))
      wavelength = info["WAVELENGTH"]

    twotheta = self.hkl_list.two_theta(wavelength = wavelength)
    L_mm = []
    L_pixels = []
    for tt in twotheta: L_mm.append(panel.settings.distance * math.tan(tt[1]))
    for lmm in L_mm: L_pixels.append(lmm/panel._img._raw.pixel_size)

    [ dc.DrawCircle(xc, yc, panel._img.get_scale() * pxl) for pxl in L_pixels ]
