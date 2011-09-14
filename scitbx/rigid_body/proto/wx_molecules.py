from scitbx.rigid_body.proto import tst_molecules
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
    self.points.clear()
    for B,AJA in zip(self.sim.bodies, self.sim.AJA_accu):
      for s in B.sites:
        self.points.append(AJA * s)
    self.labels_display_list = None
    self.lines_display_list = None
    self.points_display_list = None

  def set_points_and_lines(self, simulation_factory_index, n_zigzag):
    if (simulation_factory_index == 0):
      self.sim = tst_molecules.simulation_zigzag(NB=n_zigzag)
    else:
      self.sim = tst_molecules.simulation_factories[simulation_factory_index]()
    self.points = flex.vec3_double()
    self.set_points()
    def add_line(i, j, color):
      line = (i,j)
      self.line_i_seqs.append(line)
      self.line_colors[line] = color
    self.labels = []
    B_off = []
    for B in self.sim.bodies:
      B_off.append(len(self.labels))
      self.labels.extend(B.labels)
      def add_off(i):
        if (i < 0): return B_off[B.parent+1] + i
        else:       return B_off[-1] + i
      for bond in B.bonds:
        i,j = [add_off(b) for b in bond]
        add_line(i, j, (1,0,0))
    mcs = minimum_covering_sphere(self.points, epsilon=1.e-2)
    self.minimum_covering_sphere = sphere_3d(
      center=mcs.center(), radius=mcs.radius()*1.3)
    self.flag_show_minimum_covering_sphere = False
    self.flag_show_rotation_center = False
    self.steps_per_tab = 1
    self.show_key_stroke_help()

  def show_key_stroke_help(self):
    print "Press and hold Tab key to run the simulation."
    print "Press Shift-Tab to increase speed."
    print "Press Ctrl-Tab  to decrease speed."
    print "Press [0-9A-F] for sensitivity test with this many significant" \
          " digits."
    print "Press M for minimization."

  def process_key_stroke(self, key):
    if (key == ord("M")):
      return self.minimization()
    for n,digit in enumerate("0123456789ABCDEF"):
      if (key == ord(digit)):
        return self.sensitivity_test(n_significant_digits=n)
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

  def sensitivity_test(self, n_significant_digits):
    if (n_significant_digits == 0):
      print "Sensitivity test (full precision):"
      n_significant_digits = None
    else:
      print "Sensitivity test (%d significant digits):" % n_significant_digits
    qdd = self.sim.sensitivity_test(n_significant_digits=n_significant_digits)
    flex.double(qdd).min_max_mean().show(prefix=" ")

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
    n = len(tst_molecules.simulation_factories)
    if (len(args) not in [1,2]):
      raise Usage("scitbx.python wx_molecules.py sim_index [n_zigzag]")
    self.simulation_factory_index = int(args[0])
    if (len(args) > 1):
      self.n_zigzag = int(args[1])
    else:
      self.n_zigzag = 10
    assert 0 <= self.simulation_factory_index < n
    super(App, self).__init__(title="wx_molecules")

  def init_view_objects(self):
    box = wx.BoxSizer(wx.VERTICAL)
    self.view_objects = viewer(self.frame, size=(600,600))
    self.view_objects.set_points_and_lines(
      simulation_factory_index=self.simulation_factory_index,
      n_zigzag=self.n_zigzag)
    box.Add(self.view_objects, wx.EXPAND, wx.EXPAND)
    self.frame.SetSizer(box)
    box.SetSizeHints(self.frame)

def run(args):
  App(args).MainLoop()

if (__name__ == "__main__"):
  run(sys.argv[1:])
