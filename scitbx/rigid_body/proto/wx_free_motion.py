from __future__ import absolute_import, division, print_function
from scitbx.rigid_body.proto.free_motion_reference_impl import simulation
from gltbx import wx_viewer
from scitbx.math import minimum_covering_sphere, sphere_3d
from scitbx.array_family import flex
import wx
import sys
from six.moves import range

class viewer(wx_viewer.show_points_and_lines_mixin):

  def __init__(self, *args, **kwds):
    super(viewer, self).__init__(*args, **kwds)

  def set_points_and_lines(self):
    self.sim_as = simulation()
    self.sim_ac = simulation()
    self.points = flex.vec3_double(self.sim_as.sites_cart_moved_F01)
    self.points.extend(flex.vec3_double(self.sim_ac.sites_cart_moved_F01))
    self.points.extend(flex.vec3_double(self.sim_as.sites_cart_wells_F01))
    def add_line(i, j, color):
      line = (i,j)
      self.line_i_seqs.append(line)
      self.line_colors[line] = color
    self.labels = []
    n = len(self.sim_as.sites_cart_F1)
    offs = 0
    for prefix,color in [("S",(1,0,0)),("C",(0,0,1)),("W",(0,1,0))]:
      for i in range(n):
        add_line(offs+i, offs+(i+1)%n, color)
        self.labels.append(prefix+str(i))
      offs += n
    mcs = minimum_covering_sphere(self.points, epsilon=1.e-2)
    self.minimum_covering_sphere = sphere_3d(
      center=mcs.center(), radius=mcs.radius()*1.3)
    self.flag_show_minimum_covering_sphere = False
    self.flag_show_rotation_center = False
    self.steps_per_tab = 8
    print("Press and hold Tab key to run the simulation.")
    print("Press Shift-Tab to increase speed.")
    print("Press Ctrl-Tab  to decrease speed.")

  def tab_callback(self, shift_down=False, control_down=False):
    if (shift_down or control_down):
      if (shift_down):
        self.steps_per_tab = min(256, self.steps_per_tab * 2)
      else:
        self.steps_per_tab = max(1, self.steps_per_tab // 2)
      print("Steps per Tab:", self.steps_per_tab)
      return
    ip = 0
    for sim in [self.sim_as, self.sim_ac]:
      use_classical_accel = (sim is self.sim_ac)
      for ids in range(self.steps_per_tab):
        sim.dynamics_step(
          delta_t=0.01,
          use_classical_accel=use_classical_accel)
      for site in sim.sites_cart_moved_F01:
        self.points[ip] = site
        ip += 1
    self.labels_display_list = None
    self.lines_display_list = None
    self.points_display_list = None
    self.OnRedraw()

class App(wx_viewer.App):

  def __init__(self, args):
    assert len(args) == 0
    super(App, self).__init__(title="Free Motion Viewer")

  def init_view_objects(self):
    box = wx.BoxSizer(wx.VERTICAL)
    self.view_objects = viewer(self.frame, size=(600,600))
    self.view_objects.set_points_and_lines()
    box.Add(self.view_objects, wx.EXPAND, wx.EXPAND)
    self.frame.SetSizer(box)
    box.SetSizeHints(self.frame)

def run(args):
  App(args).MainLoop()

if (__name__ == "__main__"):
  run(sys.argv[1:])
