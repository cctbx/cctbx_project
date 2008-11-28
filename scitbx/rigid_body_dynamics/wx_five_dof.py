from scitbx.rigid_body_dynamics import tst_free_motion_hard
from scitbx.math import minimum_covering_sphere, sphere_3d
from scitbx.array_family import flex
from scitbx import matrix
from gltbx import wx_viewer
import wx
import sys

class viewer(wx_viewer.show_points_and_lines_mixin):

  def __init__(self, *args, **kwds):
    super(viewer, self).__init__(*args, **kwds)

  def set_points(self):
    self.points.clear()
    for sim in [self.sim5, self.sim56]:
      self.points.extend(flex.vec3_double(sim.sites_moved()))
    assert self.sim56.wells is self.sim5.wells
    self.points.extend(flex.vec3_double(self.sim5.wells))

  def set_points_and_lines(self):
    self.sim5 = tst_free_motion_hard.five_dof_simulation(
      r_is_qr=True,
      mersenne_twister=flex.mersenne_twister(seed=0))
    self.sim56 = tst_free_motion_hard.five_six_dof_simulation(
      six_dof_type="euler_angles_xyz",
      sim5=self.sim5)
    self.sim5.qd = matrix.zeros(n=5)
    self.sim5.energies_and_accelerations_update()
    self.sim56.qd = matrix.zeros(n=6)
    self.sim56.energies_and_accelerations_update()
    self.points = flex.vec3_double()
    self.set_points()
    assert self.points.size() == 6
    def add_line(i, j, color):
      line = (i,j)
      self.line_i_seqs.append(line)
      self.line_colors[line] = color
    self.labels = ["5a", "5b", "6a", "6b", "wa", "wb"]
    add_line(0, 1, (1,0,0))
    add_line(2, 3, (0,0,1))
    add_line(4, 5, (0,1,0))
    mcs = minimum_covering_sphere(self.points, epsilon=1.e-2)
    self.minimum_covering_sphere = sphere_3d(
      center=mcs.center(), radius=mcs.radius()*2.5)
    self.flag_show_minimum_covering_sphere = False
    self.flag_show_rotation_center = False
    self.steps_per_tab = 8
    print "Press and hold Tab key to run the simulation."
    print "Press Shift-Tab to increase speed."
    print "Press Ctrl-Tab  to decrease speed."

  def tab_callback(self, shift_down=False, control_down=False):
    if (shift_down or control_down):
      if (shift_down):
        self.steps_per_tab = min(256, self.steps_per_tab * 2)
      else:
        self.steps_per_tab = max(1, self.steps_per_tab // 2)
      print "Steps per Tab:", self.steps_per_tab
      return
    for sim in [self.sim5, self.sim56]:
      sim.dynamics_step(delta_t=0.05)
    self.set_points()
    self.labels_display_list = None
    self.lines_display_list = None
    self.points_display_list = None
    self.OnRedraw()

class App(wx_viewer.App):

  def __init__(self, args):
    assert len(args) == 0
    super(App, self).__init__(title="five_dof")

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
