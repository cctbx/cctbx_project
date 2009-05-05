from scitbx.rigid_body.essence import tst_molecules
from scitbx.math import minimum_covering_sphere, sphere_3d
from scitbx.array_family import flex
from gltbx import wx_viewer
from libtbx.option_parser import libtbx_option_parser
from libtbx.utils import Usage
import wx
import sys

class viewer(wx_viewer.show_points_and_lines_mixin):

  def __init__(self, *args, **kwds):
    super(viewer, self).__init__(*args, **kwds)

  def set_points(self):
    self.points = flex.vec3_double(self.sim.sites_moved())
    self.labels_display_list = None
    self.lines_display_list = None
    self.points_display_list = None

  def set_points_and_lines(self,
        sim,
        velocity_scaling=False,
        e_kin_per_dof=1,
        minimum_covering_sphere_view_scale=1.3,
        show_loop_edge_bendings=True):
    self.sim = sim
    if (e_kin_per_dof is None):
      self.e_kin_target = sim.e_kin() / max(1, sim.degrees_of_freedom)
    else:
      self.e_kin_target = e_kin_per_dof * sim.degrees_of_freedom
      sim.assign_random_velocities(e_kin_target=self.e_kin_target)
    self.velocity_scaling = velocity_scaling
    self.labels = self.sim.labels
    self.set_points()
    for line,color in sim.tardy_tree.viewer_lines_with_colors(
          include_loop_edge_bendings=show_loop_edge_bendings):
      self.line_i_seqs.append(line)
      self.line_colors[line] = color
    print "\n".join(sim.tardy_tree.viewer_lines_with_colors_legend(
      include_loop_edge_bendings=show_loop_edge_bendings))
    mcs = minimum_covering_sphere(self.points, epsilon=1.e-2)
    self.minimum_covering_sphere = sphere_3d(
      center=mcs.center(),
      radius=mcs.radius()*minimum_covering_sphere_view_scale)
    self.flag_show_minimum_covering_sphere = False
    self.flag_show_rotation_center = False
    self.steps_per_tab = 1
    self.show_key_stroke_help()

  def show_key_stroke_help(self):
    print "Press and hold Tab key to run the simulation."
    print "Press Shift-Tab to increase speed."
    print "Press Ctrl-Tab  to decrease speed."
    print "Press M for minimization."

  def process_key_stroke(self, key):
    if (key == ord("M")):
      return self.minimization()
    print "No action for this key stroke."
    self.show_key_stroke_help()

  def tab_callback(self, shift_down=False, control_down=False):
    if (shift_down or control_down):
      if (shift_down):
        self.steps_per_tab = min(256, self.steps_per_tab * 2)
      else:
        self.steps_per_tab = max(1, self.steps_per_tab // 2)
      print "Steps per Tab:", self.steps_per_tab
      return
    self.sim.dynamics_step(delta_t=0.05)
    if (self.velocity_scaling):
      self.sim.reset_e_kin(e_kin_target=self.e_kin_target)
    print "e_kin+e_pot: %12.6g + %12.6g = %12.6g" % (
      self.sim.e_kin(), self.sim.e_pot(), self.sim.e_tot())
    self.set_points()
    self.OnRedraw()

  def minimization(self):
    print "Minimization:"
    print "  start e_pot:", self.sim.e_pot()
    self.sim.minimization(
      max_iterations=10,
      callback_after_step=self.minimization_callback)
    print "  final e_pot:", self.sim.e_pot()

  def minimization_callback(self, minimizer):
    print "        e_pot:", self.sim.e_pot()
    self.set_points()
    self.OnRedraw()

class App(wx_viewer.App):

  def __init__(self, args):
    n = tst_molecules.n_test_simulations
    command_line = (libtbx_option_parser(
      usage="""\
scitbx.python wx_molecules.py [options] sim_index
  sim_index range: 0 ... %d\
""" % (n-1))
      .option(None, "--i_seq_labels",
        action="store_true",
        default=False)
      .option(None, "--velocity_scaling",
        action="store_true",
        default=False)
      .option(None, "--e_kin_per_dof",
        type="float",
        default=1.0,
        metavar="FLOAT")
      .option(None, "--view_scale",
        type="float",
        default=1.3,
        metavar="FLOAT")
    ).process(args=args, nargs=1)
    co = command_line.options
    self.i_seq_labels = co.i_seq_labels
    self.velocity_scaling = co.velocity_scaling
    self.e_kin_per_dof = co.e_kin_per_dof
    self.view_scale = co.view_scale
    self.simulation_index = int(command_line.args[0])
    assert 0 <= self.simulation_index < n
    super(App, self).__init__(title="wx_molecules")

  def init_view_objects(self):
    box = wx.BoxSizer(wx.VERTICAL)
    self.view_objects = viewer(self.frame, size=(600,600))
    sim = tst_molecules.get_test_simulation_by_index(i=self.simulation_index)
    if (self.i_seq_labels):
      sim.labels = [str(i) for i in xrange(len(sim.labels))]
    self.view_objects.set_points_and_lines(
      sim=sim,
      velocity_scaling=self.velocity_scaling,
      e_kin_per_dof=self.e_kin_per_dof,
      minimum_covering_sphere_view_scale=self.view_scale)
    box.Add(self.view_objects, wx.EXPAND, wx.EXPAND)
    self.frame.SetSizer(box)
    box.SetSizeHints(self.frame)

def run(args):
  App(args).MainLoop()

if (__name__ == "__main__"):
  run(sys.argv[1:])
