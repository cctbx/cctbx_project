from mmtbx.refinement import tst_tardy_pdb
from gltbx import wx_viewer
from libtbx.thread_utils import thread_with_callback_and_wait
from scitbx.math import minimum_covering_sphere, sphere_3d
from scitbx.array_family import flex
import wx
import sys

class viewer(wx_viewer.show_points_and_lines_mixin):

  def __init__(self, *args, **kwds):
    super(viewer, self).__init__(*args, **kwds)

  def set_points_and_lines(self, args):
    self.child_thread = thread_with_callback_and_wait(
      run=tst_tardy_pdb.run,
      run_kwds={"args": args},
      callback=self.action_callback,
      first_callback=self.first_action_callback)
    self.child_thread.start_and_wait_for_first_callback()

  def first_action_callback(self, sim):
    self.points = flex.vec3_double(sim.sites_moved())
    self.labels = sim.labels
    for line,color in sim.tardy_tree.viewer_lines_with_colors(
          include_loop_edge_bendings=False):
      self.line_i_seqs.append(line)
      self.line_colors[line] = color
    print "\n".join(sim.tardy_tree.viewer_lines_with_colors_legend(
      include_loop_edge_bendings=False))
    mcs = minimum_covering_sphere(self.points, epsilon=1.e-2)
    self.minimum_covering_sphere = sphere_3d(
      center=mcs.center(),
      radius=mcs.radius()*2.0)
    self.flag_show_minimum_covering_sphere = False
    self.flag_show_rotation_center = False

  def action_callback(self, sim):
    wx.PostEvent(self, wx_viewer.ViewerUpdateEvent(sim))

  def OnUpdate (self, event) :
    sim = event.data
    self.points = flex.vec3_double(sim.sites_moved())
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
    self.args = args
    super(App, self).__init__(title="wx_tardy_pdb")

  def init_view_objects(self):
    box = wx.BoxSizer(wx.VERTICAL)
    self.view_objects = viewer(self.frame, size=(600,600))
    self.view_objects.set_points_and_lines(args=self.args)
    box.Add(self.view_objects, wx.EXPAND, wx.EXPAND)
    self.frame.SetSizer(box)
    box.SetSizeHints(self.frame)

def run(args):
  App(args=args).MainLoop()

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
