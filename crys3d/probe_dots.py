
import gltbx.gl_managed
from gltbx.gl import *
from libtbx import group_args

class probe_dots_mixin (object) :
  def __init__ (self) :
    self._dots = []
    self.probe_dots_dispaly_list = None
    self.flag_show_probe_dots = True
    self.flag_show_probe_hb = True
    self.flag_show_probe_wc = True
    self.flag_show_probe_cc = True
    self.flag_show_probe_so = True
    self.flag_show_probe_bo = True

  def read_probe_dots_unformatted (self, file_name) :
    lines = open(file_name).readlines()
    for line in lines :
      fields = line.split(":")
      gap = float(fields[6])
      if (fields[2] == "hb") : c = (0.2, 0.8, 0.4)
      elif (gap > 0.35) :  c = (0.2,0.2,1.0)
      elif (gap > 0.25) :  c = (0.2, 0.6, 1.0)
      elif (gap > 0.15) :  c = (0.2, 1.0, 1.0)
      elif (gap > 0.0) :   c = (0.2, 1.0, 0.2)
      elif (gap > -0.1) :  c = (0.8, 1.0, 0.2)
      elif (gap > -0.2) :  c = (1.0, 1.0, 0.0)
      elif (gap > -0.3) :  c = (1.0, 0.6, 0.0)
      elif (gap > -0.4) :  c = (1.0, 0.0, 0.0)
      else :               c = (1.0, 0.4, 0.4)
      dot = group_args(
        dot_type=fields[2],
        atom1=fields[3],
        atom2=fields[4],
        color=c,
        xyz1=(float(fields[7]), float(fields[8]), float(fields[9])),
        xyz2=(float(fields[14]), float(fields[15]), float(fields[16])))
      self._dots.append(dot)

  def draw_probe_dots (self) :
    glEnable(GL_LINE_SMOOTH)
    glEnable(GL_BLEND)
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)
    glDisable(GL_LIGHTING)
    if (self.probe_dots_display_list is None) :
      self.probe_dots_display_list = gltbx.gl_managed.display_list()
      self.probe_dots_display_list.compile()
      glLineWidth(2.0)
      for dot in self._dots :
        glColor3f(*(dot.color))
        glBegin(GL_LINES)
        glVertex3f(*(dot.xyz1))
        glVertex3f(*(dot.xyz2))
        glEnd()
      self.probe_dots_display_list.end()
    self.probe_dots_display_list.call()
