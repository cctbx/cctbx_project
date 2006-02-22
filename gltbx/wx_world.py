from __future__ import division
from __future__ import generators

# This code is based on:
#   http://lists.wxwidgets.org/archive/wxPython-users/msg11078.html

import gltbx.util
from gltbx.gl import *
from gltbx.glu import *
import gltbx.gl_managed
import gltbx.fonts
import gltbx.images
from scitbx.array_family import flex
import scitbx.math
from scitbx import matrix
import wx
import wx.glcanvas
import math
import time
import sys, os

def animation_stepper(time_move=1.0, move_factor=1, frames_per_second=100):
  time_move *= move_factor
  n_steps = max(1, int(time_move * frames_per_second + 0.5))
  time_per_frame = time_move / n_steps
  i_step = 1
  t0 = time.time()
  while True:
    f = min(1, i_step/n_steps)
    yield f
    if (f == 1): break
    ideal_t = i_step * time_per_frame
    i_step += 1
    delta_t = time.time() - t0
    if (delta_t < ideal_t):
      time.sleep(ideal_t - delta_t)
    else:
      i_step = max(i_step, int(delta_t/time_per_frame + 0.5))

def v3distsq(a, b):
  result = 0
  for x,y in zip(a,b): result += (x-y)**2
  return result

class line_given_points:

  def __init__(self, points):
    self.points = [matrix.col(point) for point in points]
    self.delta = self.points[1] - self.points[0]
    self.delta_norm_sq = self.delta.norm_sq()
    assert self.delta_norm_sq != 0

  def distance_sq(self, point):
    "http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html"
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
    self.Bind(wx.EVT_LEFT_UP, self.OnLeftUp)
    self.Bind(wx.EVT_MIDDLE_DOWN, self.OnMiddleClick)
    self.Bind(wx.EVT_RIGHT_DOWN, self.OnRightClick)
    self.Bind(wx.EVT_MOTION, self.OnMouseMotion)
    self.Bind(wx.EVT_IDLE, self.OnIdle)

    self.w, self.h = self.GetClientSizeTuple()

    # The _back color
    self.r_back = 0
    self.g_back = 0
    self.b_back = 0

    self.field_of_view_y = 10.0
    self.min_near = 1
    self.min_viewport_use_fraction = 0.01

    self.rotation_center = (0,0,0)

    self.parent = parent
    # Current coordinates of the mouse.
    self.xmouse = 0
    self.ymouse = 0

    self.xspin = 0
    self.yspin = 0

    # Is the widget currently autospinning?
    self.autospin = False

    self.initLeft = (0,0)
    self.was_dragged = False
    self.pick_points = None

  def OnEraseBackground(self, event):
    pass # Do nothing, to avoid flashing on MSW.

  def OnSize(self, event=None):
    self.w, self.h = self.GetClientSizeTuple()
    if (self.GetContext()):
      self.SetCurrent()
      glViewport(0, 0, self.w, self.h)

  def OnIdle(self,event):
    if (self.autospin):
      wx.WakeUpIdle()
      self.do_AutoSpin()
      event.Skip(1)

  def OnChar(self,event):
    key = event.GetKeyCode()
    if   (key == ord('c')):
      self.move_to_center_of_viewport(self.rotation_center)
    elif (key == ord('f')):
      self.fit_into_viewport()
    elif (key == ord('a')):
      self.reset_rotation()
    elif (key == ord('s')):
      self.autospin_allowed = not self.autospin_allowed
    elif (key == ord('V')):
      gltbx.util.show_versions()
    self.autospin = False

  def OnLeftClick(self,event):
    self.OnRecordMouse(event)
    self.initLeft = event.GetX(),event.GetY()
    self.was_dragged = False
    self.autospin = False

  def OnLeftUp(self,event):
    if (not self.was_dragged):
      self.get_pick_points(self.initLeft)
      self.process_pick_points()
      self.OnRedraw()
    else:
      self.was_dragged = False
      if (not event.m_shiftDown):
        self.OnAutoSpin(event)

  def OnMiddleClick(self, event):
    self.OnRecordMouse(event)

  def OnRightClick(self, event):
    self.OnRecordMouse(event)

  def OnLeftDrag(self,event):
    self.was_dragged = True
    self.OnRotate(event)

  def OnMiddleDrag(self,event):
    self.OnTranslate(event)

  def OnRightDrag(self,event):
    self.OnScale(event)

  def OnMouseMotion(self,event):
    if (not event.Dragging()):
      return
    if (event.LeftIsDown()):
      self.OnLeftDrag(event)
    elif (event.MiddleIsDown()):
      self.OnMiddleDrag(event)
    elif (event.RightIsDown()):
      self.OnRightDrag(event)

  def setup_viewing_volume(self, orthographic=False):
    aspect = self.w / max(1,self.h)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    if (orthographic):
      s = self.minimum_covering_sphere
      c = s.center()
      r = s.radius()
      rf = self.buffer_factor * r
      left = c[0] - rf
      right = c[0] + rf
      bottom = c[1] - rf
      top = c[1] + rf
      if (aspect < 1):
        bottom /= aspect
        top /= aspect
      else:
        left *= aspect
        right *= aspect
      glOrtho(left, right, bottom, top, self.near, self.far)
    else:
      gluPerspective(self.field_of_view_y, aspect, self.near, self.far)

  def compute_home_translation(self):
    s = self.minimum_covering_sphere
    x,y,z = [-v for v in gltbx.util.object_as_eye_coordinates(s.center())]
    h = s.radius()
    if (self.w < self.h):
      h *= self.h / max(1,self.w)
    assert 0 < self.field_of_view_y < 90
    z -= h / math.tan(self.field_of_view_y*math.pi/180/2)
    return x,y,z

  def initialize_modelview(self):
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
    gluLookAt(0,0,0, 0,0,-1, 0,1,0)
    glTranslated(*self.compute_home_translation())
    self.rotation_center = self.minimum_covering_sphere.center()

  def rotation_move_factor(self, rotation_angle):
    return abs(rotation_angle)/180

  def translation_move_factor(self, translation_vector):
    return min(2,  abs(matrix.col(translation_vector))
                 / max(1,self.minimum_covering_sphere.radius()))

  def fit_into_viewport(self):
    dx,dy,dz = self.compute_home_translation()
    mvm = gltbx.util.get_gl_modelview_matrix()
    for f in animation_stepper(move_factor=self.translation_move_factor(
                                translation_vector=(dx,dy,dz))):
      glMatrixMode(GL_MODELVIEW)
      glLoadIdentity()
      glTranslated(f*dx, f*dy, f*dz)
      glMultMatrixd(mvm)
      self.OnRedraw()

  def reset_rotation(self):
    co = gltbx.util.object_as_eye_coordinates(self.rotation_center)
    rc = self.rotation_center
    aa = scitbx.math.r3_rotation_axis_and_angle_from_matrix(
      r=gltbx.util.extract_rotation_from_gl_modelview_matrix())
    u,v,w = aa.axis
    angle = -aa.angle(deg=True)
    mvm = gltbx.util.get_gl_modelview_matrix()
    for f in animation_stepper(move_factor=self.rotation_move_factor(
                                 rotation_angle=angle)):
      glMatrixMode(GL_MODELVIEW)
      glLoadMatrixd(mvm)
      gltbx.util.rotate_object_about_eye_vector(
        xcenter=rc[0], ycenter=rc[1], zcenter=rc[2],
        xvector=u, yvector=v, zvector=w,
        angle=f*angle)
      self.OnRedraw()
    # just to eliminate round-off errors
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
    cn = gltbx.util.object_as_eye_coordinates(self.rotation_center)
    mvm = gltbx.util.get_gl_modelview_matrix()
    glLoadIdentity()
    glTranslated(*[o-n for o,n in zip(co,cn)])
    glMultMatrixd(mvm)
    self.OnRedraw()

  def move_to_center_of_viewport(self, obj_coor):
    dx,dy = [-x for x in gltbx.util.object_as_eye_coordinates(obj_coor)[:2]]
    move_factor = min(1,  abs(matrix.col((dx,dy,0)))
                        / max(1,self.minimum_covering_sphere.radius()))
    mvm = gltbx.util.get_gl_modelview_matrix()
    for f in animation_stepper(move_factor=self.translation_move_factor(
                                translation_vector=(dx,dy,0))):
      glMatrixMode(GL_MODELVIEW)
      glLoadIdentity()
      glTranslated(f*dx, f*dy, 0)
      glMultMatrixd(mvm)
      self.OnRedraw()

  def OnRecordMouse(self, event):
    self.xmouse = event.GetX()
    self.ymouse = event.GetY()

  def OnStartRotate(self, event):
    self.autospin = False
    self.OnRecordMouse(event)

  def OnScale(self, event):
    s = self.minimum_covering_sphere
    r = (1+1.e-6)*s.radius()
    d = -gltbx.util.object_as_eye_coordinates(self.rotation_center)[2]
    dr = d + r
    scale = 0.02 * (event.GetY() - self.ymouse)
    if (scale > 0 and dr <= self.min_near):
      pass # near limit
    elif (scale < 0 and r < d * math.sin(self.field_of_view_y*math.pi/180/2)
                              * self.min_viewport_use_fraction):
      pass # far limit
    else:
      gltbx.util.translate_object(0,0,dr*scale)
      self.OnRedraw()
    self.OnRecordMouse(event)

  def OnAutoSpin(self, event):
    if (self.autospin_allowed):
      self.autospin = True
      self.yspin = 0.1 * (event.GetX()-self.initLeft[0])
      self.xspin = 0.1 * (event.GetY()-self.initLeft[1])
      if (self.xspin == 0 and self.yspin == 0):
        self.autospin = False
      else:
        self.do_AutoSpin()

  def do_AutoSpin(self):
    spin_factor = 0.05
    rc = self.rotation_center
    glMatrixMode(GL_MODELVIEW)
    gltbx.util.rotate_object_about_eye_x_and_y(
      spin_factor, rc[0], rc[1], rc[2],
      self.yspin, self.xspin, 0, 0)
    self.OnRedraw()

  def OnRotate(self, event):
    xp = event.GetX()
    yp = event.GetY()
    rc = self.rotation_center
    if (not event.m_shiftDown):
      gltbx.util.rotate_object_about_eye_x_and_y(
        0.5, rc[0], rc[1], rc[2],
        xp, yp, self.xmouse, self.ymouse)
    else:
      sz = self.GetClientSizeTuple()
      sz = (sz[0]/2, sz[1]/2)
      dy = (self.ymouse-yp)
      dx = (self.xmouse-xp)
      if (yp > sz[1]): dx *= -1
      if (xp < sz[0]): dy *= -1
      angle = (dx + dy)/2
      gltbx.util.rotate_object_about_eye_vector(
        xcenter=rc[0], ycenter=rc[1], zcenter=rc[2],
        xvector=0, yvector=0, zvector=1,
        angle=angle)
    self.OnRedraw()
    self.OnRecordMouse(event)

  def OnTranslate(self, event):
    model = gltbx.util.get_gl_modelview_matrix()
    proj = gltbx.util.get_gl_projection_matrix()
    view = gltbx.util.get_gl_viewport()
    winx = []
    winy = []
    winz = []
    rc = self.rotation_center
    assert gluProject(
      rc[0], rc[1], rc[2],
      model, proj, view,
      winx, winy, winz)
    objx = []
    objy = []
    objz = []
    win_height = max(1, self.w)
    assert gluUnProject(
      winx[0], winy[0]+0.5*win_height, winz[0],
      model, proj, view,
      objx, objy, objz)
    dist = v3distsq((objx[0],objy[0],objz[0]), rc)**0.5
    scale = abs(dist / (0.5 * win_height))
    x,y = event.GetX(), event.GetY()
    gltbx.util.translate_object(scale, x, y, self.xmouse, self.ymouse)
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

  def OnPaint(self, event=None):
    wx.PaintDC(self)
    self.SetCurrent()
    if (self.GL_uninitialised):
      glViewport(0, 0, self.w, self.h)
      self.InitGL()
      self.GL_uninitialised = 0
    self.OnRedrawGL(event)

  def OnRedraw(self, event=None):
    wx.ClientDC(self)
    self.OnRedrawGL(event)

  def OnRedrawGL(self, event=None):
    gltbx.util.handle_error()
    s = self.minimum_covering_sphere
    r = (1+1.e-6)*s.radius()
    z = -gltbx.util.object_as_eye_coordinates(s.center())[2]
    self.near = max(self.min_near, z-r)
    self.far = max(self.near*(1.e-6), z+r)
    self.setup_viewing_volume()
    gltbx.util.handle_error()
    glClear(GL_COLOR_BUFFER_BIT)
    glClear(GL_DEPTH_BUFFER_BIT)
    mvm = gltbx.util.get_gl_modelview_matrix()
    self.DrawGL()
    gltbx.util.handle_error()
    glFlush()
    self.SwapBuffers()
    if (event is not None): event.Skip()

class show_cube(wxGLWindow):

  def __init__(self, *args, **keyword_args):
    wxGLWindow.__init__(self, *args, **keyword_args)
    self.labels = None
    self.points = None
    self.lines = None

  def set_points(self, atom_attributes_list):
    self.labels = [atom.pdb_format() for atom in atom_attributes_list]
    self.points = flex.vec3_double([atom.coordinates
      for atom in atom_attributes_list])
    s = scitbx.math.minimum_covering_sphere_3d(points=self.points)
    self.minimum_covering_sphere = s
    self.buffer_factor = 1.20
    self.labels_display_list = None
    self.points_display_list = None

  def set_lines(self, bond_proxies):
    self.lines_display_list = None
    self.bond_proxies = bond_proxies

  def InitGL(self):
    gltbx.util.handle_error()
    glClearColor(self.r_back, self.g_back, self.b_back, 0.0)
    self.cube_display_list = None
    self.labels_display_list = None
    self.points_display_list = None
    self.lines_display_list = None
    self.minimum_covering_sphere_display_list = None
    self.initialize_modelview()
    gltbx.util.handle_error()

  def DrawGL(self):
    #self.draw_cube()
    self.draw_minimum_covering_sphere()
    self.draw_points()
    self.draw_lines()
    self.draw_rotation_center()
    #self.draw_labels()

  def draw_minimum_covering_sphere(self):
    if (self.minimum_covering_sphere_display_list is None):
      self.minimum_covering_sphere_display_list=gltbx.gl_managed.display_list()
      self.minimum_covering_sphere_display_list.compile()
      s = self.minimum_covering_sphere
      c = s.center()
      r = s.radius()
      gray = 0.3
      glColor3f(gray,gray,gray)
      glBegin(GL_POLYGON)
      for i in xrange(360):
        a = i * math.pi / 180
        rs = r * math.sin(a)
        rc = r * math.cos(a)
        glVertex3f(c[0]+rs,c[1]+rc,c[2])
      glEnd()
      self.draw_cross_at(c, color=(1,0,0))
      self.minimum_covering_sphere_display_list.end()
    self.minimum_covering_sphere_display_list.call()

  def draw_points(self):
    if (self.points_display_list is None):
      self.points_display_list = gltbx.gl_managed.display_list()
      self.points_display_list.compile()
      for point in self.points:
        self.draw_cross_at(point)
      self.points_display_list.end()
    self.points_display_list.call()

  def draw_lines(self):
    if (self.lines_display_list is None):
      self.lines_display_list = gltbx.gl_managed.display_list()
      self.lines_display_list.compile()
      glColor3f(1,0,1)
      glBegin(GL_LINES)
      for proxy in self.bond_proxies.simple:
        i_seq, j_seq = proxy.i_seqs
        glVertex3f(*self.points[i_seq])
        glVertex3f(*self.points[j_seq])
      glEnd()
      self.lines_display_list.end()
    self.lines_display_list.call()

  def draw_labels(self):
    if (self.labels_display_list is None):
      font = gltbx.fonts.ucs_bitmap_8x13
      font.setup_call_lists()
      self.labels_display_list = gltbx.gl_managed.display_list()
      self.labels_display_list.compile()
      for label,point in zip(self.labels, self.points):
        glRasterPos3f(*point)
        font.render_string(label)
      self.labels_display_list.end()
    self.labels_display_list.call()

  def draw_cube(self, f=1):
    if (self.cube_display_list is None):
      font = gltbx.fonts.ucs_bitmap_10x20
      font.setup_call_lists()
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
      glRasterPos3f(0, 0, 0)
      font.render_string("O")
      glRasterPos3f(1, 0, 0)
      font.render_string("x")
      glRasterPos3f(0, 1, 0)
      font.render_string("y")
      glRasterPos3f(0, 0, 1)
      font.render_string("z")
      self.cube_display_list.end()
    self.cube_display_list.call()

  def draw_cross_at(self, (x,y,z), color=(1,1,1), f=0.1):
    glBegin(GL_LINES)
    glColor3f(*color)
    glVertex3f(x-f,y,z)
    glVertex3f(x+f,y,z)
    glVertex3f(x,y-f,z)
    glVertex3f(x,y+f,z)
    glVertex3f(x,y,z-f)
    glVertex3f(x,y,z+f)
    glEnd()

  def draw_rotation_center(self):
    self.draw_cross_at(self.rotation_center, color=(0,1,0))

  def process_pick_points(self):
    line = line_given_points(self.pick_points)
    min_dist_sq = 1**2
    closest_point = None
    for point in self.points:
      dist_sq = line.distance_sq(point=matrix.col(point))
      if (min_dist_sq > dist_sq):
        min_dist_sq = dist_sq
        closest_point = point
    if (closest_point is not None):
      self.rotation_center = closest_point

def pdb_interpretation(file_name):
  from mmtbx import monomer_library
  import mmtbx.monomer_library.server
  import mmtbx.monomer_library.pdb_interpretation
  #
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  processed_pdb = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=file_name,
    strict_conflict_handling=False,
    keep_monomer_mappings=False,
    log=sys.stdout)
  processed_pdb.geometry_restraints_manager()
  return processed_pdb

class App(wx.App):

  def __init__(self, args):
    assert len(args) == 1
    self.args = args
    self.processed_pdb = pdb_interpretation(file_name=args[0])
    wx.App.__init__(self, 0)

  def OnInit(self):
    # XXX gltbx.images.create_encoded("align.bmp")
    self.frame = wx.Frame(
      None, -1,
      "Wire World",
      wx.DefaultPosition,
      wx.Size(400,400))

    self.frame.CreateStatusBar()

    tb = self.frame.CreateToolBar(
      style = wx.TB_HORIZONTAL | wx.NO_BORDER | wx.TB_FLAT | wx.TB_TEXT)
    tb.SetToolBitmapSize(wx.Size(20,20))

    center_bmp = gltbx.images.center_img.as_wx_Bitmap()
    fit_bmp = gltbx.images.fit_img.as_wx_Bitmap()
    align_bmp = gltbx.images.align_img.as_wx_Bitmap()
    spiral_bmp = gltbx.images.spiral_img.as_wx_Bitmap()

    tb.AddSimpleTool(10, center_bmp, "Center",
      "Moves center of rotation to center of window. Keyboard shortcut: c")
    self.Bind(wx.EVT_TOOL, self.OnToolClick, id=10)

    tb.AddSimpleTool(20, fit_bmp, "Fit size",
      "Resize object to fit into window. Keyboard shortcut: f")
    self.Bind(wx.EVT_TOOL, self.OnToolClick, id=20)

    tb.AddSimpleTool(30, align_bmp, "Align",
      "Aligns object and eye coordinate systems. Keyboard shortcut: a")
    self.Bind(wx.EVT_TOOL, self.OnToolClick, id=30)

    tb.AddSimpleTool(40, spiral_bmp, "Spin on/off",
      "Turns auto-spin on/off. Keyboard shortcut: s")
    self.Bind(wx.EVT_TOOL, self.OnToolClick, id=40)

    tb.AddSeparator()
    tb.AddControl(wx.StaticText(tb, -1, "Pick:"))
    pick_action_choices = [
      "Center of Rotation",
      "Atom",
      "Residue",
      "Chain",
      "Conformer",
      "Model"]
    pick_action_combobox = wx.ComboBox(
      tb, -1,
      pick_action_choices[0],
      choices=pick_action_choices,
      size=(-1,-1), style=wx.CB_DROPDOWN | wx.CB_READONLY)
    tb.AddControl(pick_action_combobox)
    self.Bind(wx.EVT_COMBOBOX, self.OnPickActionSelect, pick_action_combobox)

    tb.Realize()

    menuBar = wx.MenuBar()
    file_menu = wx.Menu()
    item = file_menu.Append(-1, "E&xit\tAlt-X", "Exit demo")
    self.Bind(wx.EVT_MENU, self.OnExitApp, item)
    menuBar.Append(file_menu, "&File")

    self.frame.SetMenuBar(menuBar)
    self.frame.Show(True)
    self.cube = show_cube(self.frame, -1, wx.Point(0,0), wx.Size(400,400))
    self.cube.set_points(
      self.processed_pdb.all_chain_proxies.stage_1.atom_attributes_list)
    self.cube.set_lines(
      self.processed_pdb.geometry_restraints_manager()
        .pair_proxies().bond_proxies)
    self.cube.SetFocus()
    self.SetTopWindow(self.frame)
    return True

  def OnExitApp(self, event):
    self.frame.Close(True)

  def OnToolClick(self, event):
    id = event.GetId()
    if (id == 10):
      self.cube.move_to_center_of_viewport(self.cube.rotation_center)
    elif (id == 20):
      self.cube.fit_into_viewport()
    if (id == 30):
      self.cube.reset_rotation()
    elif (id == 40):
      self.cube.autospin_allowed = not self.cube.autospin_allowed
      self.cube.autospin = False
      self.frame.SetStatusText(
        "Auto Spin %s" % ["Off", "On"][int(self.cube.autospin_allowed)])

  def OnPickActionSelect(self, event):
    self.frame.SetStatusText(event.GetString())

def run(args):
  App(args).MainLoop()

if (__name__ == "__main__"):
  run(sys.argv[1:])
