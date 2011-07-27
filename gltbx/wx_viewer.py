from __future__ import division

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
import os

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

VIEWER_UPDATE_ID = wx.NewId()
class ViewerUpdateEvent (wx.PyEvent) :
  def __init__ (self, data=None, recenter=False) :
    wx.PyEvent.__init__(self)
    self.data = data
    self.recenter = recenter
    self.SetEventType(VIEWER_UPDATE_ID)

class wxGLWindow(wx.glcanvas.GLCanvas):

  def InitGL(self):
    raise NotImplemented

  def DrawGL(self):
    raise NotImplemented

  def process_keyword_arguments(self,
                                auto_spin_allowed=False,
                                orthographic=False,
                                animation_time=1, #seconds
                                background_rgb=(0,0,0),
                                **kw):
    self.autospin_allowed = auto_spin_allowed
    self.orthographic = orthographic
    self.animation_time = animation_time
    self.background_rgb = background_rgb
    return kw

  def __init__(self, parent, *args, **kw):
    kw = self.process_keyword_arguments(**kw)
    self.GL_uninitialised = 1
    wx.glcanvas.GLCanvas.__init__(*((self, parent)+args), **kw)
    self.context = None
    if (wx.VERSION[1] >= 9) : # wxPython 2.9.*
      self.context = wx.glcanvas.GLContext(self)

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
    self.Bind(wx.EVT_MOUSEWHEEL, self.OnMouseWheel)
    self.Connect(-1, -1, VIEWER_UPDATE_ID, self.OnUpdate)

    self.w, self.h = self.GetClientSizeTuple()

    self.field_of_view_y = 10.0
    self.min_near = 1
    self.min_dist = -100
    self.min_viewport_use_fraction = 0.01
    self.slab_scale = 1.0
    self.scale_max = 0.1
    self.z_min = 4.0
    self.fog_scale_factor = 0.5
    self.flag_show_fog = False # leave off by default
    self._settings_widget = None

    self.flag_use_lights = False

    self.rotation_center = (0,0,0)
    self.marked_rotation = None

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
    self.light0_position = [0, 0, 100, 0]

  def OnEraseBackground(self, event):
    pass # Do nothing, to avoid flashing on MSW.

  def OnSize(self, event=None):
    self.w, self.h = self.GetClientSizeTuple()
    if (self.GetParent().IsShown()) :
      if (self.GetContext()) or (wx.VERSION[1] >= 9) :
        if (self.context is not None) :
          self.SetCurrent(self.context)
        else :
          self.SetCurrent()
        glViewport(0, 0, self.w, self.h)

  def OnIdle(self,event):
    if (self.autospin):
      wx.WakeUpIdle()
      self.do_AutoSpin()
      event.Skip(1)

  def OnChar(self,event):
    key = event.GetKeyCode()
    if   (key == ord('m')):
      self.move_rotation_center_to_mcs_center()
    elif (key == ord('c')):
      self.move_to_center_of_viewport(self.rotation_center)
    elif (key == ord('f')):
      self.fit_into_viewport()
    elif key == ord('k'):
      self.mark_rotation()
    elif (key == ord('a')):
      self.snap_back_rotation()
    elif (key == ord('s')):
      self.autospin_allowed = not self.autospin_allowed
    elif (key == ord('l')) :
      self.flag_show_labels = not self.flag_show_labels
      self.OnRedraw()
    elif (key == ord('S')) :
      self.save_screen_shot()
    elif (key == ord('V')):
      gltbx.util.show_versions()
    elif (key == ord('O')):
      self.edit_opengl_settings()
    elif (key == ord('\t')):
      callback = getattr(self, "tab_callback", None)
      if (callback is None):
        print "Tab callback not available."
      else:
        kwargs = {"shift_down": event.ShiftDown() }
        if (event.ControlDown()): kwargs["control_down"] = True
        callback(**kwargs)
    else:
      callback = getattr(self, "process_key_stroke", None)
      if (callback is None):
        print "No action for this key stroke."
      else:
        if (callback(key=key) == False):
          print "No action for this key stroke."
    self.autospin = False

  def OnMouseWheel(self, event):
    scale = 0.002*event.GetWheelRotation()
    self.OnScale(scale)

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
      if (not event.ShiftDown()):
        self.OnAutoSpin(event)

  def OnMiddleClick(self, event):
    self.OnRecordMouse(event)

  def OnRightClick(self, event):
    self.OnRecordMouse(event)

  def OnLeftDrag(self,event):
    self.was_dragged = True
    if event.AltDown():
      self.OnTranslate(event)
    else:
      self.OnRotate(event)

  def OnMiddleDrag(self,event):
    self.OnTranslate(event)

  def OnRightDrag(self,event):
    scale = 0.02 * (event.GetY() - self.ymouse)
    self.OnScale(scale)
    self.OnRecordMouse(event)

  def OnMouseMotion(self,event):
    if (not event.Dragging()):
      return
    if (event.LeftIsDown()):
      self.OnLeftDrag(event)
    elif (event.MiddleIsDown()):
      self.OnMiddleDrag(event)
    elif (event.RightIsDown()):
      self.OnRightDrag(event)

  def setup_viewing_volume(self):
    aspect = self.w / max(1,self.h)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    near, far = self.get_clipping_distances()
    if self.orthographic:
      s = self.minimum_covering_sphere
      #c = s.center()
      c = self.rotation_center
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
      glOrtho(left, right, bottom, top, near, far)
    else:
      gluPerspective(self.field_of_view_y, aspect, near, far)
    self.set_lights()
    self.setup_fog()

  def get_clipping_distances (self) :
    slab = self.far - self.near
    clip = (1.0 - self.slab_scale) * (slab / 2.0)
    near = self.near + clip
    far = self.far - clip
    if near > far :
      near = far - 1
    if near < self.min_near :
      near = self.min_near
    return (near, far)

  def setup_lighting (self) :
    if self.flag_use_lights :
      glEnable(GL_LIGHTING)
      glEnable(GL_LIGHT0)
      glEnable(GL_BLEND)
      #glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
      glBlendFunc(GL_SRC_ALPHA,GL_ZERO)
      glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE)
      glLightfv(GL_LIGHT0, GL_AMBIENT, [0.0, 0.0, 0.0, 1.0])
      glLightfv(GL_LIGHT0, GL_DIFFUSE, [1, 1, 1, 1])
      glLightfv(GL_LIGHT0, GL_SPECULAR, [0.5, 0.5, 0.5, 1.0])

  def set_lights (self) :
    if self.flag_use_lights :
      glMatrixMode(GL_MODELVIEW)
      glPushMatrix()
      glLoadIdentity()
      glEnable(GL_LIGHTING)
      glEnable(GL_LIGHT0)
      glLightfv(GL_LIGHT0, GL_POSITION, self.light0_position)
      glPopMatrix()

  def setup_fog (self) :
    if self.flag_show_fog :
      near, far = self.get_clipping_distances()
      fog_start = near + self.fog_scale_factor*(far - near)
      fog_end = max(fog_start + 5, far)
      glMatrixMode(GL_MODELVIEW)
      glEnable(GL_FOG)
      glFogi(GL_FOG_MODE, GL_LINEAR)
      glFogf(GL_FOG_START, fog_start)
      glFogf(GL_FOG_END, fog_end)
      b = self.background_rgb
      glFogfv(GL_FOG_COLOR, [b[0], b[1], b[2], 1.0])
    else :
      glDisable(GL_FOG)

  def set_minimum_covering_sphere(self, atoms=None):
    if (atoms is None):
      points = self.points
    else:
      points = flex.vec3_double()
      for atom in atoms:
        points.append(atom)
    if (points is not None and len(points) > 1):
      s = scitbx.math.minimum_covering_sphere_3d(points=points)
    else:
      if (points is None or len(points) == 0):
        center = (0,0,0)
      else:
        center = points[0]
      s = scitbx.math.sphere_3d(center=center, radius=1)
    self.minimum_covering_sphere = s

  def compute_home_translation(self):
    s = self.minimum_covering_sphere
    x,y,z = [-v for v in gltbx.util.object_as_eye_coordinates(s.center())]
    h = s.radius()
    if (self.w < self.h):
      h *= self.h / max(1,self.w)
    assert 0 < self.field_of_view_y < 90
    z -= h / math.tan(self.field_of_view_y*math.pi/180/2)
    return x,y,z

  def initialize_modelview(self, eye_vector=None, angle=None):
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
    self.setup_lighting()
    gluLookAt(0,0,0, 0,0,-1, 0,1,0)
    translation = self.compute_home_translation()
    glTranslated(*translation)
    rc = self.minimum_covering_sphere.center()
    self.rotation_center = rc
    if eye_vector is None: eye_vector = (1,1,1)
    if angle is None: angle = -120
    gltbx.util.rotate_object_about_eye_vector(
      xcenter=rc[0], ycenter=rc[1], zcenter=rc[2],
      xvector=eye_vector[0], yvector=eye_vector[1], zvector=eye_vector[2],
      angle=angle)

  def rotation_move_factor(self, rotation_angle):
    return abs(rotation_angle)/180

  def translation_move_factor(self, translation_vector):
    return min(2,  abs(matrix.col(translation_vector))
                 / max(1,self.minimum_covering_sphere.radius()))

  def fit_into_viewport(self):
    dx,dy,dz = self.compute_home_translation()
    move_factor=self.translation_move_factor((dx,dy,dz))
    mvm = gltbx.util.get_gl_modelview_matrix()
    for f in animation_stepper(time_move=self.animation_time,
                               move_factor=move_factor):
      glMatrixMode(GL_MODELVIEW)
      glLoadIdentity()
      glTranslated(f*dx, f*dy, f*dz)
      glMultMatrixd(mvm)
      self.OnRedraw()

  def mark_rotation(self):
    self.marked_rotation = matrix.sqr(
      gltbx.util.extract_rotation_from_gl_modelview_matrix())

  def snap_back_rotation(self):
    rc = self.rotation_center
    rotation_to_undo = matrix.sqr(
      gltbx.util.extract_rotation_from_gl_modelview_matrix())
    if self.marked_rotation is not None:
      rotation_to_undo *= self.marked_rotation.inverse()
    aa = scitbx.math.r3_rotation_axis_and_angle_from_matrix(
      r=rotation_to_undo.as_mat3())
    u,v,w = aa.axis
    angle = -aa.angle(deg=True)
    mvm = gltbx.util.get_gl_modelview_matrix()
    for f in animation_stepper(time_move=self.animation_time,
                               move_factor=self.rotation_move_factor(angle)):
      glMatrixMode(GL_MODELVIEW)
      glLoadMatrixd(mvm)
      gltbx.util.rotate_object_about_eye_vector(
        xcenter=rc[0], ycenter=rc[1], zcenter=rc[2],
        xvector=u, yvector=v, zvector=w,
        angle=f*angle)
      self.OnRedraw()

  def move_rotation_center_to_mcs_center(self):
    self.rotation_center = self.minimum_covering_sphere.center()

  def move_to_center_of_viewport(self, obj_coor):
    dx,dy = [-x for x in gltbx.util.object_as_eye_coordinates(obj_coor)[:2]]
    move_factor = self.translation_move_factor((dx,dy,0))
    mvm = gltbx.util.get_gl_modelview_matrix()
    for f in animation_stepper(time_move=self.animation_time,
                               move_factor=move_factor):
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

  def OnScale(self, scale):
    if (abs(scale) > self.scale_max) :
      if (scale < 0) :
        scale = -self.scale_max
      else :
        scale = self.scale_max
    s = self.minimum_covering_sphere
    r = (1+1.e-6)*s.radius()
    d = -gltbx.util.object_as_eye_coordinates(self.rotation_center)[2]
    dr = d + r
    if (scale > 0 and (dr <= self.min_near or d <= self.min_dist)):
      pass # near limit
    elif (scale < 0 and r < d * math.sin(self.field_of_view_y*math.pi/180/2)
                              * self.min_viewport_use_fraction):
      pass # far limit
    else:
      gltbx.util.translate_object(0,0,dr*scale)
      self.OnRedraw()

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
    if (not event.ShiftDown()):
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
    rc_eye = gltbx.util.object_as_eye_coordinates(rc)
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
    self.rotation_center = tuple(
      gltbx.util.modelview_matrix_as_rt().inverse() * matrix.col(rc_eye))
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
    if (self.context is not None) :
      self.SetCurrent(self.context)
    else :
      self.SetCurrent()
    if (self.GL_uninitialised):
      glViewport(0, 0, self.w, self.h)
      self.InitGL()
      self.GL_uninitialised = 0
    self.OnRedrawGL(event)

  def OnRedraw(self, event=None):
    wx.ClientDC(self)
    self.OnRedrawGL(event)

  def setup_distances (self) :
    s = self.minimum_covering_sphere
    r = self.buffer_factor*s.radius()
    #z = -gltbx.util.object_as_eye_coordinates(s.center())[2]
    z = -gltbx.util.object_as_eye_coordinates(self.rotation_center)[2]
    if (z < self.z_min) :
      z = self.z_min
    self.near = max(self.min_near, z-r)
    self.far = max(self.near+1, z+r)

  def OnRedrawGL(self, event=None):
    gltbx.util.handle_error()
    self.setup_distances()
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

  def show_stack_sizes (self) :
    mv_depth = [0]
    pr_depth = [0]
    tx_depth = [0]
    glGetIntegerv(GL_MODELVIEW_STACK_DEPTH, mv_depth)
    glGetIntegerv(GL_PROJECTION_STACK_DEPTH, pr_depth)
    glGetIntegerv(GL_TEXTURE_STACK_DEPTH, tx_depth)
    print "Modelview: %d  Projection: %d  Texture: %d" % (mv_depth[0],
      pr_depth[0], tx_depth[0])

  def OnUpdate (self, event=None) :
    pass

  def force_update (self, recenter=False) :
    wx.PostEvent(self, ViewerUpdateEvent(data=None, recenter=recenter))

  def edit_opengl_settings (self, event=None) :
    if self._settings_widget is None :
      self._settings_widget = OpenGLSettingsToolbox(self)
      self._settings_widget.Show()

  def save_screen_shot_via_pil(self,
        file_name="wx_viewer",
        extensions=["png", "jpg", "tiff", "eps", "pdf"]):
    import gltbx.viewer_utils
    from libtbx.utils import Sorry
    from libtbx.str_utils import show_string
    pil_image = gltbx.viewer_utils.read_pixels_to_pil_image(
      x=0, y=0, width=self.w, height=self.h)
    if (pil_image is None):
      print \
        "Cannot save screen shot to file:" \
        " Python Imaging Library (PIL) not available."
      return 0
    print "Screen shot width x height: %d x %d" % (self.w, self.h)
    save = pil_image.save
    def try_save(file_name_ext):
      try: save(file_name_ext)
      except KeyboardInterrupt: raise
      except Exception: return False
      return True
    for ext in extensions:
      if (file_name.endswith("."+ext)):
        print "Writing file: %s" % show_string(os.path.abspath(file_name))
        if (not try_save(file_name_ext=file_name)):
          print "Failure saving screen shot as %s file." % ext.upper()
        return 1
    n_written = 0
    for ext in extensions:
      file_name_ext = file_name + "."+ext
      if (not try_save(file_name_ext=file_name_ext)):
        print "Image output format not available: %s" % ext.upper()
      else:
        print "Wrote file: %s" % show_string(os.path.abspath(file_name_ext))
        n_written += 1
    return n_written

  def save_screen_shot_via_gl2ps(self, file_name):
    from libtbx.str_utils import show_string
    gl2ps = gltbx.util.gl2ps_interface
    if (not gl2ps(file_name=None, draw_background=False, callback=None)):
      print "PDF output via gl2ps not available: cannot write file %s" \
        % file_name
      return 0
    try:
      # preempt potential error in C++, for better reporting here
      open(file_name, "wb")
    except KeyboardInterrupt: raise
    except Exception:
      print "Error opening file for writing: %s" % \
        show_string(os.path.abspath(file_name))
      return 0
    gl2ps(file_name=file_name, draw_background=False, callback=self.OnRedraw)
    print "Wrote file: %s" % show_string(os.path.abspath(file_name))
    return 1

  def save_screen_shot(self,
        file_name="wx_viewer",
        extensions=["png", "jpg", "tiff", "eps", "pdf"]):
    extensions_pil = []
    save_pdf = file_name.endswith(".pdf")
    gl2ps_file_name = file_name
    if (not save_pdf):
      for ext in extensions:
        if (ext == "pdf"):
          save_pdf = True
          gl2ps_file_name += "."+ext
        else:
          extensions_pil.append(ext)
    n_written = 0
    if (len(extensions_pil) != 0):
      n_written += self.save_screen_shot_via_pil(
        file_name=file_name, extensions=extensions_pil)
    if (save_pdf):
      n_written += self.save_screen_shot_via_gl2ps(file_name=gl2ps_file_name)
    if (n_written == 0):
      raise Sorry(
        "Cannot save screen shot in any of the formats specified.")
    return n_written

class show_points_and_lines_mixin(wxGLWindow):

  def __init__(self, *args, **keyword_args):
    wxGLWindow.__init__(self, *args, **keyword_args)
    self.buffer_factor = 2.0
    self.flag_show_minimum_covering_sphere = True
    self.flag_show_labels = True
    self.flag_show_points = True
    self.flag_show_lines = True
    self.flag_show_rotation_center = True
    self.flag_show_spheres = True
    self.long_labels = None
    self.labels = []
    self.points = flex.vec3_double()
    self.line_i_seqs = []
    self.line_colors = {}
    self.spheres = []
    self.line_width = 1
    self.labels_display_list = None
    self.points_display_list = None
    self.lines_display_list = None

  def InitGL(self):
    gltbx.util.handle_error()
    b = self.background_rgb
    glClearColor(b[0], b[1], b[2], 0.0)
    self.minimum_covering_sphere_display_list = None
    self.initialize_modelview()
    gltbx.util.handle_error()

  def DrawGL(self):
    if (self.flag_show_minimum_covering_sphere):
      self.draw_minimum_covering_sphere()
    if (self.flag_show_points):
      self.draw_points()
    if (self.flag_show_lines):
      self.draw_lines()
    if (self.flag_show_rotation_center):
      self.draw_rotation_center()
    if (self.flag_show_labels):
      self.draw_labels()
    if (self.flag_show_spheres):
      self.draw_spheres()

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
        glVertex3f(c[0],c[1]+rs,c[2]+rc)
      glEnd()
      self.draw_cross_at(c, color=(1,0,0))
      self.minimum_covering_sphere_display_list.end()
    self.minimum_covering_sphere_display_list.call()

  def draw_points(self):
    if (self.points_display_list is None):
      self.points_display_list = gltbx.gl_managed.display_list()
      self.points_display_list.compile()
      glLineWidth(1)
      for point in self.points:
        self.draw_cross_at(point)
      self.points_display_list.end()
    self.points_display_list.call()

  def draw_lines(self):
    if (self.lines_display_list is None):
      self.lines_display_list = gltbx.gl_managed.display_list()
      self.lines_display_list.compile()
      assert self.line_width > 0
      for i_seqs in self.line_i_seqs:
        color = self.line_colors.get(tuple(i_seqs))
        if (color is None):
          color = self.line_colors.get(tuple(reversed(i_seqs)))
          if (color is None):
            color = (1,0,1)
        glColor3f(*color)
        glLineWidth(self.line_width)
        glBegin(GL_LINES)
        glVertex3f(*self.points[i_seqs[0]])
        glVertex3f(*self.points[i_seqs[1]])
        glEnd()
      self.lines_display_list.end()
    self.lines_display_list.call()

  def draw_labels(self, color=(1,1,1)):
    if (self.labels_display_list is None):
      font = gltbx.fonts.ucs_bitmap_8x13
      font.setup_call_lists()
      self.labels_display_list = gltbx.gl_managed.display_list()
      self.labels_display_list.compile()
      glColor3f(*color)
      for label,point in zip(self.labels, self.points):
        glRasterPos3f(*point)
        font.render_string(label)
      self.labels_display_list.end()
    self.labels_display_list.call()

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

  def draw_spheres(self, solid=False):
    glMatrixMode(GL_MODELVIEW)
    gray = 0.3
    glColor3f(gray,gray,gray)
    if (solid):
      glEnable(GL_LIGHTING)
      glEnable(GL_LIGHT0)
      glLightfv(GL_LIGHT0, GL_AMBIENT, [1, 1, 1, 1])
      glLightfv(GL_LIGHT0, GL_POSITION, [0, 0, 1, 0])
      glEnable(GL_BLEND)
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
      glMaterialfv(GL_FRONT, GL_DIFFUSE, [1,1,1,0.5])
      sphere = gltbx.util.SolidSphere
      grid = 50
    else:
      sphere = gltbx.util.WireSphere
      grid = 20
    for i, (x,r) in enumerate(self.spheres):
      glPushMatrix()
      glTranslated(*(x))
      sphere(radius=r, slices=grid, stacks=grid)
      glPopMatrix()
    if (solid):
      glDisable(GL_LIGHTING)
      glDisable(GL_LIGHT0)
      glDisable(GL_BLEND)

  def process_pick_points(self):
    line = scitbx.math.line_given_points(self.pick_points)
    min_dist_sq = 1**2
    i_point_closest = None
    for i_point,point in enumerate(self.points):
      dist_sq = line.distance_sq(point=matrix.col(point))
      if (min_dist_sq > dist_sq):
        min_dist_sq = dist_sq
        i_point_closest = i_point
    if (i_point_closest is not None):
      self.rotation_center = self.points[i_point_closest]
      if (self.long_labels is not None):
        assert len(self.long_labels) == len(self.points)
        lbl = self.long_labels[i_point_closest]
      elif (self.labels is not None):
        assert len(self.labels) == len(self.points)
        lbl = self.labels[i_point_closest]
      else:
        lbl = "index %d" % i_point_closest
      txt = "pick point: %s" % lbl
      print txt
      self.parent.SetStatusText(txt)

class OpenGLSettingsToolbox (wx.MiniFrame) :
  def __init__ (self, parent) :
    wx.MiniFrame.__init__(self, parent, -1, title="OpenGL settings",
      pos=(100,100), style=wx.CAPTION|wx.CLOSE_BOX|wx.RAISED_BORDER)
    self.parent = parent
    self.widgets = {}
    panel = wx.Panel(self, -1)
    main_sizer = wx.BoxSizer(wx.VERTICAL)
    fog_box = wx.CheckBox(panel, -1, "Use fog")
    fog_box.SetValue(parent.flag_show_fog)
    main_sizer.Add(fog_box, 0, wx.ALL, 5)
    self.fog_box = fog_box
    szr = wx.FlexGridSizer(rows=0, cols=2, vgap=5, hgap=5)
    main_sizer.Add(szr, 0, 0, 0)
    slab_label = wx.StaticText(panel, -1, "Slab:")
    slab_slider = wx.Slider(panel, -1, int(parent.slab_scale * 100),
      minValue=1, maxValue=100)
    szr.Add(slab_label, 0, wx.ALL, 5)
    szr.Add(slab_slider, 0, wx.ALL, 5)
    fog_label = wx.StaticText(panel, -1, "Fog scale:")
    fog_slider = wx.Slider(panel, -1, int(parent.fog_scale_factor * 100),
      minValue=1, maxValue=100)
    szr.Add(fog_label, 0, wx.ALL, 5)
    szr.Add(fog_slider, 0, wx.ALL, 5)
    self.widgets['slab_scale'] = slab_slider
    self.widgets['fog_scale_factor'] = fog_slider
    self.SetSizer(main_sizer)
    main_sizer.Fit(panel)
    self.Fit()
    self.Bind(wx.EVT_SLIDER, self.OnUpdate)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdate)
    self.Bind(wx.EVT_CLOSE, self.OnClose)

  def OnUpdate (self, event=None) :
    for setting_name, widget in self.widgets.iteritems() :
      new_value = float(widget.GetValue()) / 100.0
      setattr(self.parent, setting_name, new_value)
    self.parent.flag_show_fog = self.fog_box.GetValue()
    self.parent.OnRedrawGL()

  def OnClose (self, event=None) :
    self.Destroy()
    self.parent._settings_widget = None

class App(wx.App):

  def __init__(self, title="gltbx.wx_viewer", default_size=(600,600)):
    self.title = title
    self.default_size = wx.Size(*default_size)
    wx.App.__init__(self, 0)

  def OnInit(self):
    self.frame = wx.Frame(
      None, -1,
      self.title,
      wx.DefaultPosition,
      self.default_size)
    self.frame.Bind(wx.EVT_CLOSE, self.OnFrameClose)

    self.frame.CreateStatusBar()

    tb = self.frame.CreateToolBar(
      style = wx.TB_HORIZONTAL | wx.NO_BORDER | wx.TB_FLAT | wx.TB_TEXT)
    isinstance(tb, wx.ToolBar)
    tb.SetToolBitmapSize((32,32))

    import gltbx.wx_viewers_images as images

    tb.Bind(wx.EVT_TOOL, self.OnToolClick)

    self.mcs_center_id = wx.NewId()
    tb.AddSimpleTool(self.mcs_center_id,
      images.mcs_img.as_wx_Bitmap(),
      "Rotate around MCS center",
      "Object are to be rotated around the Minimum Covering Sphere (MCS)"
      "centre."
      " Keyboard shortcut: m")

    self.center_on_screen_id = wx.NewId()
    tb.AddSimpleTool(self.center_on_screen_id,
      images.centre_img.as_wx_Bitmap(),
      "Center in window",
      "Shifts object so that centre of rotation is centered in window."
      " Keyboard shortcut: c")

    self.fit_on_screen_id = wx.NewId()
    tb.AddSimpleTool(self.fit_on_screen_id,
      images.fit_img.as_wx_Bitmap(),
      "Fit in window",
      "Resizes and shifts object to fit into window."
      " Keyboard shortcut: f")

    self.mark_snap_back_id = wx.NewId()
    tb.AddSimpleTool(self.mark_snap_back_id,
      images.mark_snap_back_img.as_wx_Bitmap(),
      "Mark orientation for snap-back",
      "Marks object orientation as the target of a subsequent snap-back."
      " Keyboard shortcut: k")

    self.snap_back_id = wx.NewId()
    tb.AddSimpleTool(self.snap_back_id,
      images.snap_back_img.as_wx_Bitmap(),
      "Snap back orientation",
      "Rotates object back to the last marked orientation."
      " Keyboard shortcut: a")

    self.toggle_spin_id = wx.NewId()
    tb.AddCheckTool(self.toggle_spin_id,
      images.spin_img.as_wx_Bitmap(),
      shortHelp="Spin on/off",
      longHelp="Turns auto-spin on/off. Keyboard shortcut: s")

    tb.Realize()

    menuBar = wx.MenuBar()
    file_menu = wx.Menu()
    item = file_menu.Append(-1, "E&xit\tAlt-X", "Exit demo")
    self.Bind(wx.EVT_MENU, self.OnExitApp, item)
    menuBar.Append(file_menu, "&File")

    self.frame.SetMenuBar(menuBar)
    self.init_view_objects()
    self.update_status_bar()
    self.view_objects.SetFocus()
    self.SetTopWindow(self.frame)
    self.frame.Show(True)
    return True

  def OnExitApp(self, event):
    self.frame.Close(True)

  def OnFrameClose(self, event):
    f = getattr(self.view_objects, "CleanupBeforeFrameClose", None)
    if (f is not None): f()
    self.frame.Destroy()

  def OnToolClick(self, event):
    id = event.GetId()
    if (id == self.mcs_center_id):
      self.view_objects.move_rotation_center_to_mcs_center()
    elif (id == self.center_on_screen_id):
      self.view_objects.move_to_center_of_viewport(
        self.view_objects.rotation_center)
    elif (id == self.fit_on_screen_id):
      self.view_objects.fit_into_viewport()
    elif id == self.mark_snap_back_id:
      self.view_objects.mark_rotation()
    elif (id == self.snap_back_id):
      self.view_objects.snap_back_rotation()
    elif (id == self.toggle_spin_id):
      self.view_objects.autospin_allowed \
        = not self.view_objects.autospin_allowed
      self.view_objects.autospin = False
      self.update_status_bar()
    else:
      raise RuntimeError("Unknown event Id: %d" % id)

  def update_status_bar(self):
    self.frame.SetStatusText("Auto Spin %s"
      % ["Off", "On"][int(self.view_objects.autospin_allowed)])
