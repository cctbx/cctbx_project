# This code is based on:
#   http://lists.wxwidgets.org/archive/wxPython-users/msg11078.html

import gltbx.util
from gltbx.gl import *
from gltbx.glu import *
import gltbx.gl_managed
import gltbx.fonts
from scitbx import matrix
import wx
import wx.glcanvas

def v3distsq(a, b):
  result = 0
  for x,y in zip(a,b): result += (x-y)**2
  return result

class line_given_points:
  "http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html"

  def __init__(self, points):
    self.points = [matrix.col(point) for point in points]
    self.delta = self.points[1] - self.points[0]
    self.delta_norm_sq = self.delta.norm_sq()
    assert self.delta_norm_sq != 0

  def distance_sq(self, point):
    return self.delta.cross(point - self.points[0]).norm_sq() \
         / self.delta_norm_sq

class wxGLWindow(wx.glcanvas.GLCanvas):

  def InitGL(self):
    raise NotImplemented

  def DrawGL(self):
    raise NotImplemented

  def __init__(self, parent, *args, **kw):
    if (kw.has_key("autospin_allowed")):
      self.autospin_allowed = kw["autospin_allowed"]
      del kw["autospin_allowed"]
    else:
      self.autospin_allowed = 0
    self.GL_uninitialised = 1
    wx.glcanvas.GLCanvas.__init__(*((self, parent)+args), **kw)

    self.Bind(wx.EVT_ERASE_BACKGROUND, self.OnEraseBackground)
    self.Bind(wx.EVT_SIZE, self.OnSize)
    self.Bind(wx.EVT_PAINT, self.OnPaint)
    self.Bind(wx.EVT_CHAR, self.OnChar)
    self.Bind(wx.EVT_LEFT_DOWN, self.OnLeftClick)
    self.Bind(wx.EVT_LEFT_DCLICK, self.OnLeftDClick)
    self.Bind(wx.EVT_LEFT_UP, self.OnLeftUp)
    self.Bind(wx.EVT_MIDDLE_DOWN, self.OnMiddleClick)
    self.Bind(wx.EVT_RIGHT_DOWN, self.OnRightClick)
    self.Bind(wx.EVT_RIGHT_DCLICK, self.OnRightDClick)
    self.Bind(wx.EVT_MOTION, self.OnMouseMotion)
    self.Bind(wx.EVT_IDLE, self.OnIdle)

    self.w, self.h = self.GetClientSizeTuple()

    # The _back color
    self.r_back = 0
    self.g_back = 0
    self.b_back = 0

    # Where the eye is
    self.base_distance = self.distance = 10.0

    # Field of view in y direction
    self.fovy = 30.0

    # Position of clipping planes.
    self.near = 0.1
    self.far = 1000.0

    # Where we are centering.
    self.xcenter = 0.0
    self.ycenter = 0.0
    self.zcenter = 0.0

    self.translation = [0,0,0]

    self.parent = parent
    # Current coordinates of the mouse.
    self.xmouse = 0
    self.ymouse = 0

    self.xspin = 0
    self.yspin = 0

    # Is the widget currently autospinning?
    self.autospin = 0

    self.initLeft = (0,0)
    self.was_dragged = False
    self.pick_points = None

  def OnEraseBackground(self, event):
    pass # Do nothing, to avoid flashing on MSW.

  def OnSize(self, event=None):
    self.w, self.h = self.GetClientSizeTuple()
    if self.GetContext():
      self.SetCurrent()
      glViewport(0, 0, self.w, self.h)

  def OnPaint(self, event=None):
    dc = wx.PaintDC(self)
    self.OnRedrawGL(event)

  def OnIdle(self,event):
    if self.autospin:
      wx.WakeUpIdle()
      self.do_AutoSpin(event)
      event.Skip(1)

  def OnChar(self,event):
    key = event.GetKeyCode()
    if key == ord('a'):
      self.autospin_allowed = not self.autospin_allowed
    if self.autospin:
      self.autospin = 0
    elif key == ord('q'):
      self.parent.Destroy()

  def OnLeftClick(self,event):
    self.OnRecordMouse(event)
    self.initLeft = event.GetX(),event.GetY()
    self.was_dragged = False

  def OnLeftDClick(self,event):
    self.OnRecordMouse(event)
    self.reset()

  def OnLeftUp(self,event):
    if (not self.was_dragged):
      self.get_pick_points(self.initLeft)
      self.process_pick_points()
      self.OnRedraw()
    else:
      self.was_dragged = False
      if not event.m_shiftDown:
        self.OnAutoSpin(event)

  def OnMiddleClick(self,event):
    self.OnRecordMouse(event)

  def OnRightClick(self,event):
    self.OnRecordMouse(event)

  def OnRightDClick(self,event):
    self.OnRecordMouse(event)
    self.distance=self.base_distance
    self.OnRedraw()

  def OnLeftDrag(self,event):
    self.was_dragged = True
    self.OnRotate(event)

  def OnMiddleDrag(self,event):
    self.OnTranslate(event)

  def OnRightDrag(self,event):
    self.OnScale(event)

  def OnMouseMotion(self,event):
    if not event.Dragging():
      return
    if event.LeftIsDown():
      self.OnLeftDrag(event)
    elif event.MiddleIsDown():
      self.OnMiddleDrag(event)
    elif event.RightIsDown():
      self.OnRightDrag(event)

  def SetBgColour(self, r, g, b):
    self.r_back = r
    self.g_back = g
    self.b_back = b

  def SetCenterpoint(self, x, y, z):
    self.xcenter = x
    self.ycenter = y
    self.zcenter = z

  def set_base_distance(self, distance):
    self.base_distance = distance

  def set_distance(self, distance):
    self.distance = distance

  def reset(self):
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
    self.translation = [0,0,0]
    self.OnRedraw()

  def OnRecordMouse(self, event):
    self.xmouse = event.GetX()
    self.ymouse = event.GetY()

  def OnStartRotate(self, event):
    self.autospin = 0
    self.OnRecordMouse(event)

  def OnScale(self, event):
    scale = 1 + 0.01 * (event.GetY() - self.ymouse)
    self.distance = self.distance * scale
    self.OnRedraw()
    self.OnRecordMouse(event)

  def do_AutoSpin(self):
    spin_factor = 0.05
    gltbx.util.RotateScene(
      spin_factor, self.xcenter, self.ycenter, self.zcenter,
      self.yspin, self.xspin, 0, 0)
    self.OnRedraw()

  def OnAutoSpin(self, event):
    if self.autospin_allowed:
      self.autospin = 1
      self.yspin = .1 * (event.GetX()-self.initLeft[0])
      self.xspin = .1 * (event.GetY()-self.initLeft[1])
      if self.xspin == 0 and self.yspin == 0:
        self.autospin = 0
      else:
        self.do_AutoSpin()

  def OnRotate(self, event):
    xp = event.GetX()
    yp = event.GetY()
    if not event.m_shiftDown:
      gltbx.util.RotateScene(0.5,
                    self.xcenter, self.ycenter, self.zcenter,
                    xp, yp, self.xmouse, self.ymouse)
    else:
      # rotate about z
      sz = self.GetClientSizeTuple()
      sz = (sz[0]/2, sz[1]/2)
      dy = (self.ymouse-yp)
      dx = (self.xmouse-xp)
      if yp > sz[1]:
        dx = dx * -1
      if xp < sz[0]:
        dy = dy * -1
      d = dx + dy
      gltbx.util.RotateAboutVector(
        xcenter=self.xcenter,
        ycenter=self.ycenter,
        zcenter=self.zcenter,
        xvector=0,
        yvector=0,
        zvector=1,
        angle=.5*d)
    self.OnRedraw()
    self.OnRecordMouse(event)

  def OnTranslate(self, event):
    win_height = max(1, self.w)
    obj_center = (self.xcenter, self.ycenter, self.zcenter)
    model = gltbx.util.get_gl_modelview_matrix()
    proj = gltbx.util.get_gl_projection_matrix()
    view = gltbx.util.get_gl_viewport()
    winx = []
    winy = []
    winz = []
    assert gluProject(
      obj_center[0], obj_center[1], obj_center[2],
      model, proj, view,
      winx, winy, winz)
    objx = []
    objy = []
    objz = []
    assert gluUnProject(
      winx[0], winy[0]+0.5*win_height, winz[0],
      model, proj, view,
      objx, objy, objz)
    dist = v3distsq( (objx[0],objy[0],objz[0]), obj_center )**0.5
    scale = abs( dist / ( 0.5 * win_height ) )
    x,y = event.GetX(), event.GetY()
    translation = (
      (scale * (x - self.xmouse)), (scale * (self.ymouse - y)), 0.0)
    self.translation = [old-new
      for old,new in zip(self.translation, translation)]
    self.OnRedraw()
    self.OnRecordMouse(event)

  def get_pick_points(self, mouse_xy):
    model = gltbx.util.get_gl_modelview_matrix()
    proj = gltbx.util.get_gl_projection_matrix()
    view = gltbx.util.get_gl_viewport()
    self.pick_points = []
    for winz in [0.0, 1.0]:
      objx = []
      objy = []
      objz = []
      ok = gluUnProject(
        mouse_xy[0], self.h-mouse_xy[1], winz,
        model, proj, view,
        objx, objy, objz)
      if (not ok):
        self.pick_points = None
        break
      self.pick_points.append((objx[0], objy[0], objz[0]))

  def OnRedraw(self, event=None):
    dc = wx.ClientDC(self)
    self.OnRedrawGL(event)

  def OnRedrawGL(self, event=None):
    self.SetCurrent()
    if self.GL_uninitialised:
      glViewport(0, 0, self.w, self.h)
      self.InitGL()
      self.GL_uninitialised = 0
    glClearColor(self.r_back, self.g_back, self.b_back, 0.0)
    glClear(GL_COLOR_BUFFER_BIT)
    glClear(GL_DEPTH_BUFFER_BIT)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(self.fovy, float(self.w)/float(self.h), self.near, self.far)
    t = self.translation
    gluLookAt(self.xcenter, self.ycenter, self.zcenter + self.distance,
              self.xcenter+t[0], self.ycenter+t[1], self.zcenter+t[2],
              0., 1., 0.)
    glMatrixMode(GL_MODELVIEW)
    glPushMatrix()
    self.DrawGL()
    glPopMatrix()
    glFlush()
    self.SwapBuffers()
    if (event is not None): event.Skip()

class show_cube(wxGLWindow):

  def InitGL(self):
    self.cube_display_list = None
    self.SetCenterpoint(0.5, 0.5, 0.5)

  def DrawGL(self):
    self.draw_cube()
    self.draw_center_of_rotation()

  def draw_cube(self, f=1):
    if (self.cube_display_list is None):
      self.cube_display_list = gltbx.gl_managed.display_list()
      self.cube_display_list.compile()
      glBegin(GL_LINES)
      glColor3f(0,f,0)
      glVertex3f(0,0,0)
      glVertex3f(0,f,0)
      glColor3f(f,f,f)
      glVertex3f(0,f,0)
      glVertex3f(f,f,0)
      glVertex3f(f,f,0)
      glVertex3f(f,0,0)
      glColor3f(f,0,0)
      glVertex3f(f,0,0)
      glVertex3f(0,0,0)
      glColor3f(f,f,f)
      glVertex3f(0,0,f)
      glVertex3f(0,f,f)
      glVertex3f(0,f,f)
      glVertex3f(f,f,f)
      glVertex3f(f,f,f)
      glVertex3f(f,0,f)
      glVertex3f(f,0,f)
      glVertex3f(0,0,f)
      glColor3f(0,0,f)
      glVertex3f(0,0,0)
      glVertex3f(0,0,f)
      glColor3f(f,f,f)
      glVertex3f(f,0,0)
      glVertex3f(f,0,f)
      glVertex3f(0,f,0)
      glVertex3f(0,f,f)
      glVertex3f(f,f,0)
      glVertex3f(f,f,f)
      glEnd()
      self.cube_display_list.end()
    self.cube_display_list.call()
    font = gltbx.fonts.ucs_bitmap_10x20
    glRasterPos3f(0, 0, 0)
    font.render_string("O")
    glRasterPos3f(1, 0, 0)
    font.render_string("x")
    glRasterPos3f(0, 1, 0)
    font.render_string("y")
    glRasterPos3f(0, 0, 1)
    font.render_string("z")

  def draw_cross_at(self, (x,y,z), f=0.1):
    glBegin(GL_LINES)
    glColor3f(1,1,1)
    glVertex3f(x-f,y,z)
    glVertex3f(x+f,y,z)
    glVertex3f(x,y-f,z)
    glVertex3f(x,y+f,z)
    glVertex3f(x,y,z-f)
    glVertex3f(x,y,z+f)
    glEnd()

  def draw_center_of_rotation(self):
    self.draw_cross_at((self.xcenter, self.ycenter, self.zcenter))

  def process_pick_points(self):
    line = line_given_points(self.pick_points)
    min_dist_sq = 1**2
    closest_point = None
    for x in [0,1]:
      for y in [0,1]:
        for z in [0,1]:
          point = matrix.col((x,y,z))
          dist_sq = line.distance_sq(point=point)
          if (min_dist_sq > dist_sq):
            min_dist_sq = dist_sq
            closest_point = point
    if (closest_point is not None):
      old_center = matrix.col((self.xcenter, self.ycenter, self.zcenter))
      old_transl = matrix.col(self.translation)
      new_center = closest_point
      new_transl = old_transl - (new_center - old_center)
      self.xcenter, self.ycenter, self.zcenter = new_center
      self.translation = list(new_transl)
      # XXX TODO update self.distance

class App(wx.App):

  def OnInit(self):
    frame = wx.Frame(
      None, -1,
      "Simple Cube",
      wx.DefaultPosition,
      wx.Size(400,400))
    show_cube(frame, -1, wx.Point(0,0), wx.Size(400,400))
    frame.Show(True)
    self.SetTopWindow(frame)
    return True

def run():
  App(0).MainLoop()

if (__name__ == "__main__"):
  run()
