
# TODO:
#  - cached scenes

from __future__ import absolute_import, division, print_function
from crys3d import hklview
#from crys3d.leapmotion import wxGLWindowLeapEnabled as wxGLWindow
from gltbx.wx_viewer import wxGLWindow
import gltbx.viewer_utils
import gltbx.gl_managed
import gltbx.quadrics
import gltbx.util
import gltbx.fonts
from gltbx.glu import *
from gltbx.gl import *
import wx

class hklview_3d (wxGLWindow) :
  def __init__ (self, *args, **kwds) :
    wxGLWindow.__init__(self, *args, **kwds)
    self.Bind(wx.EVT_LEFT_DCLICK, self.OnDoubleClick)
    self.Bind(wx.EVT_ERASE_BACKGROUND, lambda evt: None)
    # FIXME orthographic is definitely best for this application, but it isn't
    # working properly right now
    #self.orthographic = True
    parent = self.GetParent()
    if (parent is None) :
      parent = kwds.get("parent")
    assert (parent is not None)
    self.settings = parent.settings
    self.buffer_factor = 2.0
    self.min_slab = 4
    self.min_viewport_use_fraction = 0.1
    self.min_dist = 4.0
    self.flag_show_fog = True
    self.flag_use_lights = True
    self.flag_use_quadrics = False
    self.minimum_covering_sphere = None
    self.spheres_display_list = None
    self.points_display_list = None
    self.labels_display_list = None
    self.miller_array = None
    self.d_min = None
    self.scene = None
    self.animation_time = 0
    #self.fps = gltbx.viewer_utils.fps_monitor()
    # XXX prevent exception when no data are loaded
    from scitbx.math import minimum_covering_sphere
    from scitbx.array_family import flex
    points = flex.vec3_double([(0.0,0.0,0.0),(1.0,1.0,1.0)])
    mcs = minimum_covering_sphere(points=points, epsilon=0.1)
    self.minimum_covering_sphere = mcs
    #self.Bind(wx.EVT_SHOW, self.OnPaint, self)

  def set_miller_array (self, miller_array, zoom=False, merge=None) :
    self.miller_array = miller_array
    self.merge = merge
    self.d_min = miller_array.d_min()
    self.construct_reciprocal_space(merge=merge)
    # XXX for some reason InitGL does not get called until well after the
    # window is shown on (at least) Windows 7, so we need to check whether
    # it's safe to make OpenGL calls before zooming
    if (zoom) and (not self.GL_uninitialised) :
      self.fit_into_viewport()

  def construct_reciprocal_space (self, merge=None) :
    self.scene = hklview.scene(miller_array=self.miller_array,
      merge=merge,
      settings=self.settings)
    from scitbx.math import minimum_covering_sphere
    mcs = minimum_covering_sphere(points=self.scene.points,
                                  epsilon=0.1)
    self.minimum_covering_sphere = mcs
    self.spheres_display_list = None
    self.points_display_list = None
    self.labels_display_list = None
    self.rotation_center = (0,0,0)

  #--- OpenGL methods
  def InitGL(self):
    gltbx.util.handle_error()
    if (self.settings.black_background) :
      glClearColor(0.,0.,0.,0.)
    else :
      glClearColor(0.95,0.95,0.95,0.)
    self.minimum_covering_sphere_display_list = None
    glDepthFunc(GL_LESS)
    glEnable(GL_ALPHA_TEST)
    glEnable(GL_DEPTH_TEST)
    glEnable(GL_BLEND)
    glEnable(GL_LINE_SMOOTH)
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)
    self.initialize_modelview()
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
    self.proto_sphere = gltbx.quadrics.proto_ellipsoid(slices=8, stacks=8)
    gltbx.util.handle_error()

  def OnRedrawGL (self, event=None) :
    if (self.minimum_covering_sphere is None) :
      gltbx.util.handle_error()
      glClear(GL_COLOR_BUFFER_BIT)
      glClear(GL_DEPTH_BUFFER_BIT)
      glFlush()
      self.SwapBuffers()
      gltbx.util.handle_error()
    else :
      wxGLWindow.OnRedrawGL(self, event)

  def initialize_modelview (self) :
    if (self.minimum_covering_sphere is not None) :
      wxGLWindow.initialize_modelview(self)
    else :
      self.setup_lighting()

  def DrawGL(self):
    if (self.GL_uninitialised) or (self.miller_array is None) :
      return
    if (self.settings.show_axes) :
      self.draw_axes()
    if (self.settings.spheres) :
      self.draw_spheres()
    else :
      self.draw_points()
    if (self.settings.show_labels) :
      self.draw_labels()
    #self.fps.update()

  def draw_spheres (self) :
    glMatrixMode(GL_MODELVIEW)
    glShadeModel(GL_SMOOTH)
    glEnable(GL_DEPTH_TEST)
    glEnable(GL_LIGHTING)
    glEnable(GL_LIGHT0)
    glEnable(GL_NORMALIZE)
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
    if (self.spheres_display_list is None) :
      self.spheres_display_list = gltbx.gl_managed.display_list()
      self.spheres_display_list.compile()
      colors = self.scene.colors
      radii = self.scene.radii
      points = self.scene.points
      sphere = self.proto_sphere
      max_radius = self.scene.max_radius
      scale_factor = self.settings.sphere_detail / self.scene.max_radius
      assert (colors.size() == radii.size() == self.scene.points.size())
      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, [0.1, 0.1, 0.1, 1.0])
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, [0.1, 0.1, 0.1, 1.0])
      for i, hkl in enumerate(points) :
        col = list(colors[i]) + [1.0]
        #glColor3f(*colors[i])
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, col)
        if (self.flag_use_quadrics) :
          r = max(radii[i] / 10, 0.001)
          sphere.draw(hkl, (r, r, r, 0, 0, 0))
        else :
          detail = max(4, int(radii[i] * scale_factor))
          glPushMatrix()
          glTranslated(*hkl)
          gltbx.util.SolidSphere(radius=radii[i],
            slices=detail,
            stacks=detail)
          glPopMatrix()
      self.spheres_display_list.end()
    self.spheres_display_list.call()

  def draw_points (self) :
    glMatrixMode(GL_MODELVIEW)
    glDisable(GL_LIGHTING)
    if (self.points_display_list is None) :
      self.points_display_list = gltbx.gl_managed.display_list()
      self.points_display_list.compile()
      colors = self.scene.colors
      radii = self.scene.radii
      glLineWidth(2.0)
      gltbx.viewer_utils.draw_stars(
        points=self.scene.points,
        colors=colors,
        points_visible=self.scene.visible_points,
        radii=radii)
      self.points_display_list.end()
    self.points_display_list.call()

  def draw_axes (self) :
    h_axis = self.scene.axes[0]
    k_axis = self.scene.axes[1]
    l_axis = self.scene.axes[2]
    gltbx.fonts.ucs_bitmap_8x13.setup_call_lists()
    glDisable(GL_LIGHTING)
    if (self.settings.black_background) :
      glColor3f(1.0, 1.0, 1.0)
    else :
      glColor3f(0.,0.,0.)
    glLineWidth(1.0)
    glBegin(GL_LINES)
    glVertex3f(0.,0.,0.)
    glVertex3f(h_axis[0]*100, h_axis[1]*100, h_axis[2]*100)
    glEnd()
    glRasterPos3f(0.5+h_axis[0]*100, 0.2+h_axis[1]*100, 0.2+h_axis[2]*100)
    gltbx.fonts.ucs_bitmap_8x13.render_string("h")
    glBegin(GL_LINES)
    glVertex3f(0.,0.,0.)
    glVertex3f(k_axis[0]*100, k_axis[1]*100, k_axis[2]*100)
    glEnd()
    glRasterPos3f(0.2+k_axis[0]*100, 0.5+k_axis[1]*100, 0.2+k_axis[2]*100)
    gltbx.fonts.ucs_bitmap_8x13.render_string("k")
    glBegin(GL_LINES)
    glVertex3f(0.,0.,0.)
    glVertex3f(l_axis[0]*100, l_axis[1]*100, l_axis[2]*100)
    glEnd()
    glRasterPos3f(0.2+l_axis[0]*100, 0.5+l_axis[1]*100, 0.2+l_axis[2]*100)
    gltbx.fonts.ucs_bitmap_8x13.render_string("l")
    glEnable(GL_LINE_STIPPLE)
    glLineStipple(4, 0xAAAA)
    glBegin(GL_LINES)
    glVertex3f(0.,0.,0.)
    glVertex3f(-h_axis[0]*100, -h_axis[1]*100, -h_axis[2]*100)
    glEnd()
    glBegin(GL_LINES)
    glVertex3f(0.,0.,0.)
    glVertex3f(-k_axis[0]*100, -k_axis[1]*100, -k_axis[2]*100)
    glEnd()
    glBegin(GL_LINES)
    glVertex3f(0.,0.,0.)
    glVertex3f(-l_axis[0]*100, -l_axis[1]*100, -l_axis[2]*100)
    glEnd()
    glDisable(GL_LINE_STIPPLE)

  def draw_labels (self) :
    glMatrixMode(GL_MODELVIEW)
    glDisable(GL_LIGHTING)
    gltbx.fonts.ucs_bitmap_8x13.setup_call_lists()
    if (self.labels_display_list is None) :
      self.labels_display_list = gltbx.gl_managed.display_list()
      self.labels_display_list.compile()
      colors = self.scene.colors
      points = self.scene.points
      indices = self.scene.indices
      assert (colors.size() == indices.size() == points.size())
      for i_seq in self.scene.label_points :
        x, y, z = points[i_seq]
        glColor3f(*colors[i_seq])
        glRasterPos3f(x+0.5, y+0.5, z+0.5)
        gltbx.fonts.ucs_bitmap_8x13.render_string("%d,%d,%d" % indices[i_seq])
      self.labels_display_list.end()
    self.labels_display_list.call()

  #--- user input and settings
  def update_settings (self) :
    self.construct_reciprocal_space(merge=self.merge)
    if (self.settings.black_background) :
      glClearColor(0.,0.,0.,0.)
    else :
      glClearColor(0.95,0.95,0.95,0.)
    self.Refresh()

  def process_pick_points (self) :
    self.closest_point_i_seq = None
    if (self.pick_points is not None) and (self.scene is not None) :
      closest_point_i_seq = gltbx.viewer_utils.closest_visible_point(
        points=self.scene.points,
        atoms_visible=self.scene.visible_points,
        point0=self.pick_points[0],
        point1=self.pick_points[1])
      if (closest_point_i_seq is not None) :
        self.closest_point_i_seq = closest_point_i_seq
    if (self.closest_point_i_seq is not None) :
      self.scene.label_points.add(self.closest_point_i_seq)
      self.labels_display_list = None
      self.GetParent().update_clicked(index=self.closest_point_i_seq)
      #hkl, d_min, value = self.scene.get_reflection_info(
      #  self.closest_point_i_seq)
      #self.GetParent().update_clicked(hkl, d_min, value)
    else :
      self.GetParent().update_clicked(index=None)

  def clear_labels (self) :
    if (self.scene is not None) :
      self.scene.clear_labels()
      self.labels_display_list = None
      self.Refresh()

  def OnChar (self, event) :
    key = event.GetKeyCode()
    if (key == ord("m")) :
      self.OnScale(0.02)
    elif (key == ord("n")) :
      self.OnScale(-0.02)
    elif (key == ord("d")) :
      self.adjust_slab(-0.04)
    elif (key == ord("f")) :
      self.adjust_slab(0.04)
    elif (key == ord("Q")) :
      self.flag_use_quadrics = not self.flag_use_quadrics
      self.spheres_display_list = None
    elif (key == wx.WXK_UP) :
      self.rotate_view(0, 0, 0, 5, event.ShiftDown())
    elif (key == wx.WXK_DOWN) :
      self.rotate_view(0, 0, 0, -5, event.ShiftDown())
    elif (key == wx.WXK_LEFT) :
      self.rotate_view(0, 0, -5, 0, event.ShiftDown())
    elif (key == wx.WXK_RIGHT) :
      self.rotate_view(0, 0, 5, 0, event.ShiftDown())
    self.OnRedrawGL()

  def OnDoubleClick (self, event) :
    self.get_pick_points((event.GetX(), event.GetY()))
    self.process_pick_points()
    if self.closest_point_i_seq is not None :
      self.rotation_center = self.scene.points[self.closest_point_i_seq]
      self.move_to_center_of_viewport(self.rotation_center)

  def OnMouseWheel (self, event) :
    scale = event.GetWheelRotation()
    self.adjust_slab(0.01 * scale)
    self.OnRedrawGL()
