"""\
Example:
  - use results from a function with callback (motion() in this example)
  - return from callback only after an event in the GUI (Tab-key)
"""

from gltbx import wx_viewer
from libtbx.thread_utils import thread_with_callback_and_wait
from scitbx.rigid_body.proto.free_motion_reference_impl import \
  create_triangle_with_center_of_mass_at_origin
from scitbx.math import minimum_covering_sphere, sphere_3d
from scitbx.array_family import flex
from scitbx import matrix
import wx
import sys

class viewer(wx_viewer.show_points_and_lines_mixin):

  def __init__(self, *args, **kwds):
    super(viewer, self).__init__(*args, **kwds)

  def set_points_and_lines(self):
    self.child_thread = thread_with_callback_and_wait(
      run = motion,
      callback = self.motion_callback,
      first_callback = self.first_motion_callback)
    self.child_thread.start_and_wait_for_first_callback()

  def first_motion_callback(self, points):
    self.points = flex.vec3_double(points)
    self.labels = ["A", "B", "C"]
    def add_line(i, j, color):
      line = (i,j)
      self.line_i_seqs.append(line)
      self.line_colors[line] = color
    add_line(0, 1, (1,0,0))
    add_line(1, 2, (0,1,0))
    add_line(2, 0, (0,0,1))
    mcs = minimum_covering_sphere(self.points, epsilon=1.e-2)
    self.minimum_covering_sphere = sphere_3d(
      center=matrix.col(mcs.center())+matrix.col((0.5,0.5,0.5)),
      radius=mcs.radius()*2.0)
    self.flag_show_minimum_covering_sphere = False
    self.flag_show_rotation_center = False

  def motion_callback(self, points):
    wx.PostEvent(self, wx_viewer.ViewerUpdateEvent(points))

  def OnUpdate (self, event) :
    for i in xrange(len(event.data)):
      self.points[i] = event.data[i]
    self.labels_display_list = None
    self.lines_display_list = None
    self.points_display_list = None
    self.OnRedraw()

  def tab_callback(self, shift_down=False, control_down=False):
    self.child_thread.resume()

  def CleanupBeforeFrameClose(self):
    self.child_thread.resume(last_iteration=True)

class App(wx_viewer.App):

  def __init__(self, args):
    assert len(args) == 0
    super(App, self).__init__(title="Motion Viewer")

  def init_view_objects(self):
    box = wx.BoxSizer(wx.VERTICAL)
    self.view_objects = viewer(self.frame, size=(600,600))
    self.view_objects.set_points_and_lines()
    box.Add(self.view_objects, wx.EXPAND, wx.EXPAND)
    self.frame.SetSizer(box)
    box.SetSizeHints(self.frame)

def run(args):
  App(args).MainLoop()

def motion(n_steps=20, callback=None):
  points = create_triangle_with_center_of_mass_at_origin()
  shift = matrix.col((0.1,0.1,0.1))
  i_step = 0
  while (callback is not None or i_step != n_steps):
    if (callback is None):
      print [point.elems for point in points]
    else:
      status = callback(points)
      if status == False :
        break
    for i in xrange(len(points)):
      points[i] += shift
    i_step += 1
    if (i_step % 10 == 0):
      shift = -shift

if (__name__ == "__main__"):
  run(sys.argv[1:])
