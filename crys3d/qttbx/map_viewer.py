from __future__ import absolute_import, division, print_function
from PyQt4 import QtGui
from PyQt4.QtCore import Qt
import gltbx.util
from crys3d import qttbx
import crys3d.qttbx.map_viewer_controls
from gltbx import gl_managed
from gltbx.gl import *
from gltbx.glu import *
from scitbx import iso_surface
from cctbx import maptbx, uctbx

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

  positive_iso_level_fmt = "%.3g"

  def __init__(self, map_viewer):
    qttbx.widget_control_mixin.__init__(self, map_viewer, None, Qt.Tool)
    self.move(0,50)
    self.wiresBox.setChecked(self.view.wires)
    self.posIsoLevelSlider.setRange(0, 100)
    iso_level_validator = QtGui.QDoubleValidator(self)
    iso_level_validator.setRange(0, self.view.max_density, decimals=3)
    iso_level_validator.setNotation(QtGui.QDoubleValidator.ScientificNotation)
    self.posIsoLevel.setValidator(iso_level_validator)
    self.posIsoLevel.setText(self.positive_iso_level_fmt
                             % self.view.positive_iso_level)

    self.wiresBox.stateChanged[int].connect(self.view.set_wires)
    self.posIsoLevelSlider.valueChanged[int].connect(
      self.posIsoLevelLCD.display)
    self.posIsoLevelSlider.setValue(
      int(self.view.positive_iso_level/self.view.max_density*100))
    self.posIsoLevelSlider.sliderMoved[int].connect(
      self.on_positive_iso_level_slider_moved)
    self.posIsoLevelSlider.valueChanged[int].connect(
      self.posIsoLevelLCD.display)
    self.posIsoLevel.editingFinished.connect(
      self.on_iso_level_field_edited)

  def on_positive_iso_level_slider_moved(self, x):
    iso_level = self.view.max_density*x/100
    self.view.set_positive_iso_level(iso_level)
    self.posIsoLevel.setText(self.positive_iso_level_fmt % iso_level)

  def on_iso_level_field_edited(self):
    iso_level = float(self.posIsoLevel.text())
    self.view.set_positive_iso_level(iso_level)
    self.posIsoLevelSlider.setValue(iso_level/self.view.max_density*100)



class map_viewer(qttbx.widget):

  def __init__(self,
               fft_map=None,
               unit_cell=None, raw_map=None,
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
                                     light_position=(-1, 1, 1, 0),
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
    self.material = gl_managed.material_model(
      ambient_front_colour=(0, 0.25, 0, 1),
      ambient_back_colour=(0., 0.25, 0, 1),
      diffuse_front_colour=(0, 1, 0, 1),
      diffuse_back_colour=(0, 0.75, 0, 1),
      specular_front_colour=(0.5, 0, 0.5, 1),
      specular_focus=100)
    self.wire_colour = (0, 1, 0, 1)

  def draw_object_in_fractional_coordinates(self):
    self.draw_triangulation()

  def draw_triangulation(self):
    if self.wires:
      glPushAttrib(GL_LIGHTING_BIT)
      glDisable(GL_LIGHTING)
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
      glColor4fv(self.wire_colour)
    else:
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
      self.material.execute(specular=not self.wires)
    va = gltbx.util.vertex_array(self.triangulation.vertices,
                                 self.triangulation.normals)
    va.draw_triangles(self.triangulation.triangles)
    if self.wires:
      glPopAttrib()
    gltbx.util.handle_error()

if __name__ == '__main__':
  from scitbx.iso_surface import tst_iso_surface
  uc = uctbx.unit_cell((1,1,1,60,120,90))
  case = tst_iso_surface.triangulation_test_case(
    func=tst_iso_surface.elliptic(),
    grid_size=(5, 3, 2),
    periodic=False,
    lazy_normals=False,
    descending_normals=True)
  display(unit_cell=uc,
          raw_map=case.map,
          positive_iso_level=2.26,
          wires=True,
          is_unit_cell_shown=False,
          orthographic=True)
