from __future__ import absolute_import, division, print_function

from gltbx import gl

class display_lists_owner:

  def __init__(self, range_values):
    self.range_values = range_values
    self.list = gl.glGenLists(range_values)

  def __del__(self):
    try: gl.glDeleteLists(self.list, self.range_values)
    except RuntimeError as e:
      if (str(e) != 'OpenGL: invalid operation'): raise
      # else: apparently the GL context was destroyed already

class display_list:

  def __init__(self, index=0, owner=None):
    if (owner is None):
      assert index == 0
      self.owner = display_lists_owner(range_values=1)
      self.gl_index = self.owner.list
    else:
      assert index < self.owner.range_values
      self.owner = owner
      self.gl_index = self.owner.list + index

  def compile(self, execute=False):
    if (execute):
      mode = gl.GL_COMPILE_AND_EXECUTE
    else:
      mode = gl.GL_COMPILE
    gl.glNewList(self.gl_index, mode)

  def end(self):
    gl.glEndList()

  def call(self):
    gl.glCallList(self.gl_index)

class display_lists:

  def __init__(self, range_values):
    self.owner = display_lists_owner(range_values=range_values)

  def __getitem__(self, index):
    return display_list(index=index, owner=self.owner)


class material_model(object):

  def __init__(self,
               ambient_front_colour,
               diffuse_front_colour,
               specular_front_colour=(1,1,1,1),
               ambient_back_colour=None,
               diffuse_back_colour=None,
               specular_back_colour=None,
               specular_focus=30):
    self.ambient_front_colour = ambient_front_colour
    self.diffuse_front_colour = diffuse_front_colour
    self.specular_front_colour = specular_front_colour
    if ambient_back_colour is None:
      ambient_back_colour = ambient_front_colour
    self.ambient_back_colour = ambient_back_colour
    if diffuse_back_colour is None:
      diffuse_back_colour = diffuse_front_colour
    self.diffuse_back_colour = diffuse_back_colour
    if specular_back_colour is None:
      specular_back_colour = specular_front_colour
    self.specular_back_colour = specular_back_colour
    self.specular_focus = specular_focus

  def execute(self, specular=True):
    from gltbx.gl import glMaterialfv, glMaterialf
    glMaterialfv(gl.GL_BACK, gl.GL_AMBIENT, self.ambient_back_colour)
    glMaterialfv(gl.GL_FRONT, gl.GL_AMBIENT, self.ambient_front_colour)
    glMaterialfv(gl.GL_BACK, gl.GL_DIFFUSE, self.diffuse_back_colour)
    glMaterialfv(gl.GL_FRONT, gl.GL_DIFFUSE, self.diffuse_front_colour)
    if specular:
      glMaterialfv(gl.GL_BACK, gl.GL_SPECULAR, self.specular_back_colour)
      glMaterialfv(gl.GL_FRONT, gl.GL_SPECULAR, self.specular_front_colour)
      glMaterialf(gl.GL_FRONT_AND_BACK, gl.GL_SHININESS, self.specular_focus)
