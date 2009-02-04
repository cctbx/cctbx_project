from scitbx.rigid_body.essence import tst_molecules
from scitbx.math import minimum_covering_sphere, sphere_3d
from scitbx.array_family import flex
from gltbx import wx_viewer
from libtbx.utils import Usage
import wx
import sys

class viewer(wx_viewer.show_points_and_lines_mixin):

  def __init__(self, *args, **kwds):
    super(viewer, self).__init__(*args, **kwds)

  def set_points(self):
    self.points = flex.vec3_double(self.sim.sites_moved)
    self.labels_display_list = None
    self.lines_display_list = None
    self.points_display_list = None

  def set_points_and_lines(self, sim):
    self.sim = sim
    self.labels = self.sim.labels
    self.set_points()
    cm = self.sim.cluster_manager
    he = set()
    for e in cm.hinge_edges:
      if (e[0] == -1): continue
      he.add(tuple(sorted(e)))
    le = set([tuple(sorted(e)) for e in cm.loop_edges])
    assert len(he.intersection(le)) == 0
    for line in self.sim.bonds:
      self.line_i_seqs.append(line)
      if (line in he):
        color = (0,1,0)
      elif (line in le):
        color = (1,0,0)
      else:
        cii, cij = [cm.cluster_indices[i] for i in line]
        if (cii == cij and cm.hinge_edges[cii][0] == -1):
          color = (0,1,1)
        else:
          color = (0,0,1)
      self.line_colors[line] = color
    for line in cm.loop_edge_bendings:
      self.line_i_seqs.append(line)
      self.line_colors[line] = (0.5,0,0.5)
    mcs = minimum_covering_sphere(self.points, epsilon=1.e-2)
    self.minimum_covering_sphere = sphere_3d(
      center=mcs.center(), radius=mcs.radius()*1.3)
    self.flag_show_minimum_covering_sphere = False
    self.flag_show_rotation_center = False
    self.steps_per_tab = 1
    print "Edge colors:"
    print "  turquoise: intra-base-cluster, six degrees of freedom"
    print "  green:     rotatable bond, one degree of freedom"
    print "  blue:      intra-cluster"
    print "  red:       loop edge (restrained only)"
    print "  purple:    loop bending edge (restrained only)"
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
    self.set_points()
    self.OnRedraw()

  def minimization(self):
    print "Minimization:"
    print "  start e_pot:", self.sim.e_pot
    self.sim.minimization(
      max_iterations=10,
      callback_after_step=self.minimization_callback)
    print "  final e_pot:", self.sim.e_pot

  def minimization_callback(self, minimizer):
    print "        e_pot:", self.sim.e_pot
    self.set_points()
    self.OnRedraw()

class App(wx_viewer.App):

  def __init__(self, args):
    n = tst_molecules.n_test_simulations
    if (len(args) != 1):
      raise Usage("""\
scitbx.python wx_molecules.py sim_index
  sim_index range: 0 ... %d
""" % (n-1))
    self.simulation_index = int(args[0])
    assert 0 <= self.simulation_index < n
    super(App, self).__init__(title="wx_molecules")

  def init_view_objects(self):
    box = wx.BoxSizer(wx.VERTICAL)
    self.view_objects = viewer(self.frame, size=(600,600))
    self.view_objects.set_points_and_lines(
      sim=tst_molecules.get_test_simulation_by_index(self.simulation_index))
    box.Add(self.view_objects, wx.EXPAND, wx.EXPAND)
    self.frame.SetSizer(box)
    box.SetSizeHints(self.frame)

def run(args):
  App(args).MainLoop()

if (__name__ == "__main__"):
  run(sys.argv[1:])
