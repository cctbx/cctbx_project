from PyQt4 import QtGui, QtCore, QtOpenGL
from PyQt4.QtCore import Qt
import gltbx.util
from crys3d import qttbx
import crys3d.qttbx.map_viewer_controls
from gltbx import quadrics, gl_managed
from gltbx.gl import *
from gltbx.glu import *
from scitbx import matrix as mat
from scitbx import iso_surface
from scitbx.array_family import flex
from cctbx import maptbx, uctbx
import math
import sys

def display(window_title="Map Viewer", **kwds):
  app = QtGui.QApplication([])
  view = map_viewer(**kwds)
  view.setWindowTitle(window_title)
  view.resize(800, 800)
  ctrls = map_viewer_controls(view)
  ctrls.show()
  view.show()
  app.exec_()


class map_viewer_controls(qttbx.widget_control_mixin,
                          qttbx.map_viewer_controls.Ui_Form):

  def __init__(self, map_viewer):
    qttbx.widget_control_mixin.__init__(self, map_viewer, None, Qt.Tool)
    self.move(0,50)
    self.wiresBox.setChecked(self.view.wires)
    self.posIsoLevelSlider.setRange(0, 100)

    self.wiresBox.stateChanged[int].connect(self.view.set_wires)
    self.posIsoLevelSlider.valueChanged[int].connect(
      self.posIsoLevelLCD.display)
    self.posIsoLevelSlider.setValue(
      int(self.view.positive_iso_level/self.view.max_density*100))
    self.posIsoLevelSlider.valueChanged[int].connect(
      lambda x: self.view.set_positive_iso_level(self.view.max_density*x/100))


class map_viewer(qttbx.widget):

  def __init__(self,
               fft_map=None,
               unit_cell=None, raw_map=None,
               from_here=None,
               to_there=None,
               periodic=False,
               positive_iso_level=None,
               iso_level_positive_range_fraction=None,
               negative_iso_level=None,
               iso_level_negative_range_fraction=None,
               wires=True,
               **kwds):
    if fft_map is not None:
      unit_cell = fft_map.unit_cell()
    super(map_viewer, self).__init__(unit_cell=unit_cell,
                                     from_here=from_here,
                                     to_there=to_there,
                                     light_position=(1, 1, 1, 0),
                                     **kwds)
    assert (fft_map is not None
            or (unit_cell is not None and raw_map is not None))
    assert (positive_iso_level is not None
            or iso_level_positive_range_fraction is not None)
    assert (bool(negative_iso_level is not None)
            ^ bool(iso_level_negative_range_fraction is None))
    if fft_map is not None:
      self.rho = fft_map.real_map()
    else:
      self.rho = raw_map
    density_stats = maptbx.statistics(self.rho)
    self.min_density = density_stats.min()
    self.max_density = density_stats.max()
    if positive_iso_level is not None:
      self.positive_iso_level = positive_iso_level
    else:
      self.positive_iso_level = (
        iso_level_positive_range_fraction*self.max_density)
    if (negative_iso_level is None
        and iso_level_negative_range_fraction is not None):
      negative_iso_level = (
        iso_level_negative_range_fraction*self.min_density)
    self.negative_iso_level = negative_iso_level

    self.periodic = periodic
    self.wires = wires
    self.compute_triangulation()

  def set_positive_iso_level(self, level):
    if self.positive_iso_level != level:
      self.positive_iso_level = level
      self.compute_triangulation()
      self.updateGL()
    return self

  def set_perspective(self, flag):
    if self.orthographic != (not flag):
      self.orthographic = not flag
      self.resizeGL(self.width(), self.height())
      self.updateGL()
    return self

  def set_wires(self, flag):
    if self.wires != flag:
      self.wires = flag
      self.updateGL()
    return self

  def compute_triangulation(self):
    self.triangulation = iso_surface.triangulation(
      self.rho, self.positive_iso_level,
      map_extent=(1,1,1),
      from_here=self.from_here,
      to_there=self.to_there,
      periodic=self.periodic,
      ascending_normal_direction=False)

  def initialise_opengl(self):
    glEnableClientState(GL_VERTEX_ARRAY)
    glEnableClientState(GL_NORMAL_ARRAY)
    self.material = gl_managed.material_model()

  def draw_object(self):
    self.draw_triangulation()

  def draw_triangulation(self):
    self.material.execute(specular=not self.wires)
    if self.wires:
      self.light_with_only_ambient(1., 1., 1.)
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
      glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE)
    else:
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
      glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE)
    va = gltbx.util.vertex_array(self.triangulation.vertices,
                                 self.triangulation.normals)
    va.draw_triangles(self.triangulation.triangles)
    gltbx.util.handle_error()

if __name__ == '__main__':
  from scitbx.iso_surface import tst_iso_surface
  uc = uctbx.unit_cell((1,1,1,60,120,90))
  case = tst_iso_surface.triangulation_test_case(
    func=tst_iso_surface.elliptic(),
    grid_size=(50, 40, 30),
    periodic=False,
    lazy_normals=False,
    descending_normals=True)
  display(unit_cell=uc,
          raw_map=case.map,
          iso_level_positive_range_fraction=0.3,
          wires=False,
          orthographic=True)
