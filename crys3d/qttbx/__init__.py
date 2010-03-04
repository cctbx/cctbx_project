from __future__ import division
import gltbx
import gltbx.util
from gltbx.gl import *
from gltbx.glu import *
from PyQt4.QtOpenGL import *
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from scitbx import matrix as mat
import math

f = QGLFormat(QGL.SampleBuffers)
QGLFormat.setDefaultFormat(f)
del f

class widget(QGLWidget):

  def __init__(self,
               unit_cell,
               light_position,
               from_here=None,
               to_there=None,
               clear_colour=(0, 0, 0, 1),
               fovy=30, orthographic=False,
               mouse_rotation_scale=0.6, mouse_wheel_scale=0.1,
               mouse_translation_scale=0.01,
               unit_cell_axis_label_font=None,
               *args, **kwds):
    super(widget, self).__init__(*args, **kwds)
    self.unit_cell = unit_cell
    self.fovy = fovy
    if from_here is None: from_here = (0,0,0)
    if to_there is None: to_there = (1,1,1)
    self.set_extent(from_here, to_there)
    self.orthogonaliser = gltbx.util.matrix(
      unit_cell.orthogonalization_matrix())
    self.orthographic = orthographic
    self.light_position = light_position
    self.clear_colour = clear_colour
    self.mouse_rotation_scale = mouse_rotation_scale
    self.mouse_wheel_scale = mouse_wheel_scale
    self.mouse_translation_scale = mouse_translation_scale
    self.orbiting = None
    self.dolly = None
    self.mouse_position = None
    if unit_cell_axis_label_font is None:
      unit_cell_axis_label_font = QFont("Helvetica", pointSize=16)
    self.unit_cell_axis_label_font = unit_cell_axis_label_font

  def set_extent(self, from_here, to_there):
    self.from_here = mat.col(from_here)
    self.to_there = mat.col(to_there)
    self.object_centre_wrt_frac = (self.from_here + self.to_there)/2
    x0, y0, z0 = self.from_here
    x1, y1, z1 = self.to_there
    dx, dy, dz = x1 - x0, y1 - y0, z1 - z0
    diameter = max([
      self.unit_cell.length(u)
      for u in [(dx, dy, dz), (dx, dy, -dz), (dx, -dy, dz), (-dx, dy, dz) ] ])
    self.object_radius = diameter/2

  def eye_distance(self):
    return self.object_radius / math.tan(self.fovy/2*math.pi/180)
  eye_distance = property(eye_distance)

  def initializeGL(self):
    glEnable(GL_DEPTH_TEST)
    glEnable(GL_LIGHTING)
    glEnable(GL_LIGHT0)
    glEnable(GL_NORMALIZE)
    glClearColor(*self.clear_colour)
    glMatrixMode(GL_MODELVIEW)
    glPushMatrix()
    glLoadIdentity()
    self.orbiting = gltbx.util.matrix().get()
    self.dolly = gltbx.util.matrix().get()
    glPopMatrix()
    gltbx.util.handle_error()
    self.initialise_opengl()

  def resizeGL(self, w, h):
    w = max(w, 1)
    h = max(h, 1)
    aspect = w/h
    glViewport(0, 0, w, h)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    if self.orthographic:
      left = bottom = -self.object_radius
      right = top = self.object_radius
      if aspect < 1:
        bottom /= aspect
        top /= aspect
      else:
        left *= aspect
        right *= aspect
      glOrtho(left, right, bottom, top, 0, 100)
    else:
      gluPerspective(aspect=w/h, fovy=self.fovy, zNear=0.1, zFar=100)
    gltbx.util.handle_error()

  def paintGL(self):
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
    glTranslatef(0, 0, -self.eye_distance)
    self.dolly.multiply()
    self.orbiting.multiply()
    self.orthogonaliser.multiply()
    glTranslatef(*-self.object_centre_wrt_frac)
    self.draw_unit_cell()
    self.draw_object()
    gltbx.util.handle_error()

  def draw_unit_cell(self):
    glLightfv(GL_LIGHT0, GL_AMBIENT, [0.5, 0.5, 0.5, 1.])
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, (0., 0., 0., 1.))
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, (0, 0, 0, 1))
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 0)

    (x0, y0, z0), (x1, y1, z1) = self.triangulation.bounds

    r,g,b = 0.9, 0.4, 0.3
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, (r, g, b, 1.))
    e = -0.02
    for pos, label in zip([((x0 + x1)/2, y0 + e, z0 + e),
                           (x0 + e, (y0 + y1)/2, z0 + e),
                           (x0 + e, y0 + e, (z0 + z1)/2)],
                          ['a','b','c']):
      self.renderText(pos[0], pos[1], pos[2], label,
                      self.unit_cell_axis_label_font)

    lw = [0.]
    glGetFloatv(GL_LINE_WIDTH, lw)
    glLineWidth(3)

    r,g,b = (0.6,)*3
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, (r, g, b, 1.))
    glBegin(GL_LINE_LOOP)
    glVertex3f(x0,y0,z0)
    glVertex3f(x1,y0,z0)
    glVertex3f(x1,y1,z0)
    glVertex3f(x0,y1,z0)
    glEnd()
    glBegin(GL_LINE_LOOP)
    glVertex3f(x0,y0,z1)
    glVertex3f(x1,y0,z1)
    glVertex3f(x1,y1,z1)
    glVertex3f(x0,y1,z1)
    glEnd()
    glBegin(GL_LINES)
    glVertex3f(x0,y0,z0)
    glVertex3f(x0,y0,z1)
    glVertex3f(x1,y0,z0)
    glVertex3f(x1,y0,z1)
    glVertex3f(x1,y1,z0)
    glVertex3f(x1,y1,z1)
    glVertex3f(x0,y1,z0)
    glVertex3f(x0,y1,z1)
    glEnd()
    glLineWidth(lw[0])

  def mousePressEvent(self, event):
    if event.button() == Qt.LeftButton:
      self.mouse_position = event.pos()

  def mouseReleaseEvent(self, event):
    if event.button() == Qt.LeftButton:
      self.mouse_position = None

  def mouseMoveEvent(self, event):
    if self.mouse_position is None: return
    p, p0 = event.pos(), self.mouse_position
    self.mouse_position = p
    delta_x, delta_y = p.x() - p0.x(), p.y() - p0.y()
    if event.modifiers() == Qt.NoModifier:
      if not self.orbiting: return
      s = self.mouse_rotation_scale
      glMatrixMode(GL_MODELVIEW)
      glPushMatrix()
      glLoadIdentity()
      glRotatef(s*delta_x, 0, 1, 0)
      glRotatef(s*delta_y, 1, 0, 0)
      self.orbiting.multiply()
      self.orbiting.get()
      glPopMatrix()
    elif event.modifiers() == Qt.AltModifier:
      if not self.dolly: return
      s = self.mouse_translation_scale
      glMatrixMode(GL_MODELVIEW)
      glPushMatrix()
      glLoadIdentity()
      glTranslatef(s*delta_x, -s*delta_y, 0)
      self.dolly.multiply()
      self.dolly.get()
      glPopMatrix()
    self.updateGL()

  def wheelEvent(self, event):
    if not self.dolly: return
    wheel_rotation = event.delta()/8 #degrees
    wheel_rotation *= self.mouse_wheel_scale
    glMatrixMode(GL_MODELVIEW)
    glPushMatrix()
    glLoadIdentity()
    if self.orthographic:
      s = 1 + math.tanh(wheel_rotation)
      glScalef(s, s, s)
    else:
      glTranslatef(0, 0, wheel_rotation)
    self.dolly.multiply()
    self.dolly.get()
    glPopMatrix()
    self.updateGL()
