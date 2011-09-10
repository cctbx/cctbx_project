from scitbx.array_family import flex
import math
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
      horizons_phil=self.working_phil)
    frame.SetSize((1024,780))
    frame.load_image(path)
    frame.path = path
    self.path = path
    frame.Show()
    app.MainLoop()

  def user_callback(self,dc,panel,wx):
    center_x, center_y = panel._img.get_beam_center()
    xc, yc = panel._img.image_coords_as_screen_coords(center_x, center_y)
    dc.SetPen(wx.Pen('red'))
    #dc.SetBrush(wx.TRANSPARENT_BRUSH)

    wavelength = panel._img._raw.wavelength #should be this
    wavelength_from_avg_file = True
    if wavelength_from_avg_file:
      # go through hoops to get the proper wavelength corresponding to this run
      avepath = self.path.replace("stddev","avg")
      import pickle
      info = pickle.load(open(avepath,"rb"))
      wavelength = info["WAVELENGTH"]

    twotheta = 2.* flex.asin(
                 flex.double(len(self.two_theta_experimental), wavelength/2.)/
                 self.experimental_d)
    L_mm = panel.settings.distance * flex.atan(twotheta)
    L_pixels = L_mm / panel._img._raw.pixel_size

    [ dc.DrawCircle(xc, yc, panel._img.zoom * pxl) for pxl in L_pixels ]
