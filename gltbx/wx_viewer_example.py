from gltbx import wx_viewer
import wx
from gltbx.gl import *
from scitbx.math import minimum_covering_sphere
from scitbx.array_family import flex

class MyGLWindow(wx_viewer.show_points_and_lines_mixin):

  def __init__(self, *args, **kwds):
    super(MyGLWindow, self).__init__(*args, **kwds)
    self.points = flex.vec3_double([ (-5,-5,-5), (-4,0,0), (0,-8,0), (0,0,-11) ])
    self.line_i_seqs = [ (0,1), (0,2), (0,3), (1,2), (1,3), (2,3) ]
    self.spheres = [ ((0,0,0), 1) ]
    self.flag_show_minimum_covering_sphere = False
    self.minimum_covering_sphere = minimum_covering_sphere(self.points)

class MyApp(wx_viewer.App):

  def init_view_objects(self):
    box = wx.BoxSizer(wx.VERTICAL)
    self.view_objects = MyGLWindow(self.frame, size=(600,600))
    box.Add(self.view_objects, wx.EXPAND, wx.EXPAND)
    self.frame.SetSizer(box)
    box.SetSizeHints(self.frame)

if __name__ == '__main__':
  a = MyApp(title="An example of using gltbx.wx_viewer")
  a.MainLoop()
