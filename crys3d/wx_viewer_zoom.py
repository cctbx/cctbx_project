from __future__ import division
from gltbx import wx_viewer
import gltbx.util
from gltbx.gl import *
from gltbx.glu import *

class viewer_with_automatic_zoom (wx_viewer.wxGLWindow) :
  def process_key_stroke (self, key) :
    if key == wx.WXK_UP :
      self.fog_start_offset += 1
    elif key == wx.WXK_DOWN :
      self.fog_start_offset -= 1
    elif key == wx.WXK_LEFT :
      self.clip_near -= 1
    elif key == wx.WXK_RIGHT :
      self.clip_near += 1
    self.update_scene = True

  def setup_distances (self) :
    s = self.minimum_covering_sphere
    r = self.buffer_factor*s.radius()
    #z = -gltbx.util.object_as_eye_coordinates(s.center())[2]
    z = -gltbx.util.object_as_eye_coordinates(self.rotation_center)[2]
    self.near = max(self.min_near, z-r)
    self.far = max(self.near*(1.e-6), z+r)

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
    self.setup_fog()

  def get_clipping_distances (self) :
    near = self.near + self.clip_near
    far = self.far + self.clip_far
    if near > far :
      near = far - 1
    if near < self.min_near :
      near = self.min_near
    return (near, far)

  def setup_fog (self) :
    if self.flag_show_fog :
      near, far = self.get_clipping_distances()
      # TODO: this needs work.
      fog_start = 0.25*(far - near) + near + self.fog_start_offset
      fog_end = far - self.fog_end_offset
      #print "%6.1f - %6.1f (%6.1f - %6.1f)" % (near, far, fog_start, fog_end)
      glMatrixMode(GL_MODELVIEW)
      glEnable(GL_FOG)
      glFogi(GL_FOG_MODE, GL_LINEAR)
      glFogf(GL_FOG_START, fog_start)
      glFogf(GL_FOG_END, fog_end)
      glFogfv(GL_FOG_COLOR, [self.r_back, self.g_back, self.b_back, 1.0])
    else :
      glDisable(GL_FOG)
