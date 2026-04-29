from __future__ import absolute_import, division, print_function
import wx
import gltbx.util
from gltbx.gl import *
from gltbx.glu import *
from gltbx import wx_viewer
from gltbx import quadrics
import unicodedata
from six.moves import zip

class MyGLWindow(wx_viewer.wxGLWindow):

  def __init__(self, *args, **kwds):
    super(MyGLWindow, self).__init__(*args, **kwds)

    # location of the ellipsoids
    self.locations = [ (1,0,0), (-1,0,0),
                       (0,1,0), (0,-1,0),
                       (0,0,1), (0,0,-1) ]
    self.set_minimum_covering_sphere(self.locations)
    self.minimum_covering_sphere = (
      self.minimum_covering_sphere.expand_relative(0.2))

    # each element is a list of 4 symmetric matrices,
    # each defining an ellipsoid
    # one can move from tests to tests by using the left and right arrow keys
    a, b, c = (0.05, 0.1, 0.2)
    self.tests = [
      [ (a, b, c, 0, 0, 0),
        (a, c, b, 0, 0, 0),
        (b, a, c, 0, 0, 0),
        (b, c, a, 0, 0, 0),
        (c, a, b, 0, 0, 0),
        (c, b, a, 0, 0, 0) ],
    ]
    self.test_index = 0

  def InitGL(self):
    gltbx.util.handle_error()

    self.initialize_modelview(angle=45) # that does glMatrixMode(GL_MODELVIEW)
                                        # so we can safely position the light
                                        # later in this method
    self.buffer_factor = 2.

    glEnable(GL_NORMALIZE) # mighty important because we use non-uniform
                           # scaling to render cylinders and ellipsoids,
                           # which also means thatGL_RESCALE_NORMALS
                           # would not do here

    # Colours and lighting
    glClearColor(0., 0., 0., 0.)
    glShadeModel(GL_SMOOTH)
    glEnable(GL_DEPTH_TEST)
    glEnable(GL_LIGHTING)
    glEnable(GL_LIGHT0)
    glLightfv(GL_LIGHT0, GL_POSITION, [0, 0, 1, 0])
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE)
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, [0.5, 0.5, 0.5, 1.0])
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, (0.2, 0.2, 0.2, 1.))
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, (1, 1, 1, 1.))
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 50.)

    # We build our quadrics
    self.proto_cylinder = quadrics.proto_cylinder(slices=16)
    self.proto_ellipsoid = quadrics.proto_ellipsoid(slices=32, stacks=32)

    # We build the texture to paint the principal ellipses on ellipsoids
    self.principal_ellipses_tex = \
      quadrics.ellipsoid_principal_sections_texture(darkening=0.75,
                                                    n_s=64, n_t=64)

    # Enable texturing and specify how to lay the texture on the ellipsoids
    glEnable(GL_TEXTURE_2D)
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE)

    gltbx.util.handle_error()

  def DrawGL(self):
    glMatrixMode(GL_MODELVIEW) # don't forget!

    radius=0.08
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, (1., 0., 0., 1.))
    self.proto_cylinder.draw(self.locations[0], self.locations[1],
                             base_radius=radius)
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, (0., 1., 0., 1.))
    self.proto_cylinder.draw(self.locations[2], self.locations[3],
                             base_radius=radius)
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, (0., 0., 1., 1.))
    self.proto_cylinder.draw(self.locations[4], self.locations[5],
                             base_radius=radius)

    # Let's draw the following with our texture
    self.principal_ellipses_tex.bind()
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, (0.2, 0.4, 0.6, 1.))
    for x, m in zip(self.locations, self.tests[self.test_index]):
      self.proto_ellipsoid.draw(x, m)
    # End of drawing with our texture
    self.principal_ellipses_tex.unbind()

  def OnChar(self, event):
    key = event.GetKeyCode()
    if key == wx.WXK_LEFT:
      self.test_index += 1
    elif key == wx.WXK_RIGHT:
      self.test_index -= 1
    self.test_index %= len(self.tests)
    super(MyGLWindow, self).OnChar(event)


class MyApp(wx_viewer.App):

  def init_view_objects(self):
    box = wx.BoxSizer(wx.VERTICAL)
    instructions = wx.StaticText(
      self.frame,
      label="Use %s and %s to move up and down the test cases" % (
        unicodedata.lookup('LEFTWARDS ARROW'),
        unicodedata.lookup('RIGHTWARDS ARROW')),
      style=wx.ALIGN_CENTER,
    )
    box.Add(instructions, 0, wx.EXPAND)
    self.view_objects = MyGLWindow(self.frame, size=(1280, 800))
    self.view_objects.SetFocus()
    box.Add(self.view_objects, wx.EXPAND, wx.EXPAND)
    self.frame.SetSizer(box)
    box.SetSizeHints(self.frame)

if __name__ == '__main__':
  a = MyApp(title="Ellipsoids")
  a.MainLoop()
