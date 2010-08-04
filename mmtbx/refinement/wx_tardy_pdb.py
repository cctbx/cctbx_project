from mmtbx.refinement import tst_tardy_pdb
from gltbx import wx_viewer
from gltbx import gl_managed
from gltbx.gl import *
import gltbx.util
from libtbx.thread_utils import thread_with_callback_and_wait
import scitbx.iso_surface
from scitbx.math import minimum_covering_sphere, sphere_3d
from scitbx.array_family import flex
import wx
import sys

class draw_map(object):

  def __init__(O, iso_level=2):
    O.iso_level = iso_level
    O.set_unit_cell_and_density_map(unit_cell=None, density_map=None)

  def set_unit_cell_and_density_map(O, unit_cell, density_map):
    O.unit_cell = unit_cell
    O.density_map = density_map
    O.display_list = None

  def change_iso_level(O, direction):
    assert direction in ["up", "down"]
    old = O.iso_level
    if (old < 2): shift = 0.5
    else:         shift = 1
    if (direction == "down"):
      shift *= -1
    O.iso_level = min(max(0.5, old + shift), 10)
    print "New iso-level: %.2f" % O.iso_level
    O.display_list = None

  def draw_unit_cell(O):
    glLineWidth(1)
    glColor3f(1,0,0)
    glBegin(GL_LINES)
    glVertex3f(0,0,0)
    glVertex3f(1,0,0)
    glEnd()
    glColor3f(0,1,0)
    glBegin(GL_LINES)
    glVertex3f(0,0,0)
    glVertex3f(0,1,0)
    glEnd()
    glColor3f(0,0,1)
    glBegin(GL_LINES)
    glVertex3f(0,0,0)
    glVertex3f(0,0,1)
    glEnd()
    glColor3f(1,1,1)
    glBegin(GL_LINE_LOOP)
    glVertex3f(1,0,0)
    glVertex3f(1,0,1)
    glVertex3f(1,1,1)
    glVertex3f(1,1,0)
    glEnd()
    glBegin(GL_LINES)
    glVertex3f(0,0,1)
    glVertex3f(1,0,1)
    glVertex3f(0,0,1)
    glVertex3f(0,1,1)
    glVertex3f(0,1,0)
    glVertex3f(0,1,1)
    glVertex3f(0,1,1)
    glVertex3f(1,1,1)
    glVertex3f(0,1,0)
    glVertex3f(1,1,0)
    glEnd()

  def draw_triangulation(O):
    triangulation = scitbx.iso_surface.triangulation(
      map=O.density_map,
      iso_level=O.iso_level,
      map_extent=(1,1,1),
      from_here=(0,0,0),
      to_there=(1,1,1),
      periodic=True,
      ascending_normal_direction=False)
    glLineWidth(1)
    glColor3f(*[0.2]*3)
    vertices = triangulation.vertices
    for triangle in triangulation.triangles:
      glBegin(GL_LINE_LOOP)
      for i in triangle:
        glVertex3f(*vertices[i])
      glEnd()

  def __call__(O):
    if (O.unit_cell is None): return
    if (O.display_list is None):
      O.display_list = gltbx.gl_managed.display_list()
      O.display_list.compile()
      glMatrixMode(GL_MODELVIEW)
      glPushMatrix()
      glMultTransposeMatrixd(
        gltbx.util.augment_3x3(m=O.unit_cell.orthogonalization_matrix()))
      O.draw_unit_cell()
      O.draw_triangulation()
      glPopMatrix()
      O.display_list.end()
      gltbx.util.handle_error()
    O.display_list.call()

class viewer(wx_viewer.show_points_and_lines_mixin):

  def __init__(O, *args, **kwds):
    super(viewer, O).__init__(*args, **kwds)
    O.draw_map = draw_map()
    O.first_first = True

  def DrawGL(O):
    super(viewer, O).DrawGL()
    O.draw_map()

  def set_points_and_lines(O, args=None):
    if (args is not None):
      O.args = args
    O.draw_map.set_unit_cell_and_density_map(unit_cell=None, density_map=None)
    O.child_thread = thread_with_callback_and_wait(
      run=tst_tardy_pdb.run,
      run_kwds={"args": O.args},
      callback=O.action_callback,
      first_callback=O.first_action_callback)
    O.child_thread.start_and_wait_for_first_callback()

  def first_action_callback(O, tardy_model, rmsd_calculator):
    O.tardy_model = tardy_model
    tpo = tardy_model.potential_obj
    O.draw_map.set_unit_cell_and_density_map(
      unit_cell=tpo.geo_manager.crystal_symmetry.unit_cell(),
      density_map=tpo.density_map)
    O.points = tardy_model.sites_moved().deep_copy()
    if (tpo.ideal_sites_cart is not None):
      O.points.extend(tpo.ideal_sites_cart)
    if (O.points.size() < 20):
      if (tpo.ideal_sites_cart is None):
        O.labels = tardy_model.labels
      else:
        O.labels = tardy_model.labels + [""] * len(tardy_model.labels)
    def draw_ideal_line():
      if (tpo.ideal_sites_cart is None): return
      n = tardy_model.tardy_tree.n_vertices
      ideal_line = tuple([i+n for i in line])
      O.line_i_seqs.append(ideal_line)
      O.line_colors[ideal_line] = [0.6]*3
    if (tardy_model.potential_obj.reduced_geo_manager is not None):
      for line,color in tardy_model.tardy_tree.viewer_lines_with_colors(
            include_loop_edge_bendings=False):
        draw_ideal_line()
        O.line_i_seqs.append(line)
        O.line_colors[line] = color
      print "\n".join(tardy_model.tardy_tree.viewer_lines_with_colors_legend(
        include_loop_edge_bendings=False))
    else:
      for line in tardy_model.potential_obj.geo_manager.simple_edge_list():
        draw_ideal_line()
        O.line_i_seqs.append(line)
        O.line_colors[line] = (1,0,0)
    mcs = minimum_covering_sphere(O.points, epsilon=1.e-2)
    O.minimum_covering_sphere = sphere_3d(
      center=mcs.center(),
      radius=mcs.radius()*2.0)
    O.flag_show_minimum_covering_sphere = False
    O.flag_show_rotation_center = False
    O.show_key_stroke_help()
    if (O.first_first): O.first_first = False
    else:               O.action_callback()

  def show_key_stroke_help(O):
    print "Press and hold Tab key to continue execution."
    print "Press i to change iso-surface level down."
    print "Press I to change iso-surface level up."
    print "Press R to re-run from the start."

  def process_key_stroke(O, key):
    if   (key == ord("i")):
      O.draw_map.change_iso_level(direction="down")
    elif (key == ord("I")):
      O.draw_map.change_iso_level(direction="up")
    elif (key == ord("R")):
      O.child_thread.resume(last_iteration=True)
      O.set_points_and_lines()
    else:
      print "No action for this key stroke."
      O.show_key_stroke_help()
    O.OnRedraw()

  def action_callback(O):
    wx.PostEvent(O, wx_viewer.ViewerUpdateEvent(
      data=O.tardy_model.sites_moved()))

  def OnUpdate(O, event) :
    sites_moved = event.data
    for i in xrange(len(sites_moved)):
      O.points[i] = sites_moved[i]
    O.labels_display_list = None
    O.lines_display_list = None
    O.points_display_list = None
    O.OnRedraw()

  def tab_callback(O, shift_down=False, control_down=False):
    O.child_thread.resume()

  def CleanupBeforeFrameClose(O):
    O.child_thread.resume(last_iteration=True)

class App(wx_viewer.App):

  def __init__(O, args):
    O.args = args
    super(App, O).__init__(title="wx_tardy_pdb")

  def init_view_objects(O):
    box = wx.BoxSizer(wx.VERTICAL)
    O.view_objects = viewer(O.frame, size=(600,600))
    O.view_objects.set_points_and_lines(args=O.args)
    box.Add(O.view_objects, wx.EXPAND, wx.EXPAND)
    O.frame.SetSizer(box)
    box.SetSizeHints(O.frame)

def run(args):
  App(args=args).MainLoop()

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
