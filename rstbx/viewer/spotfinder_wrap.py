from __future__ import absolute_import, division, print_function
class spot_wrapper:
  def __init__(self,working_phil):
      self.working_phil = working_phil

  def display(self,path,organizer):
    import wx
    from rstbx.viewer.spotfinder_frame import SpotFrame
    from rstbx.viewer import display

    app   = wx.App(0)
    frame = SpotFrame(None, -1, "X-ray image display", size=(1200,1080),
      pos=(100,100), horizons_phil=self.working_phil,
      spot_organizer = organizer)
    frame.SetSize((1024,780))
    frame.load_image(path)
    frame.path = path
    self.path = path
    frame.Show()
    app.MainLoop()
