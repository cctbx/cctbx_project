from __future__ import division

from gltbx import gl

class display_lists_owner:

  def __init__(self, range):
    self.range = range
    self.list = gl.glGenLists(range=range)

  def __del__(self):
    try: gl.glDeleteLists(list=self.list, range=self.range)
    except RuntimeError, e:
      if (str(e) != 'OpenGL: invalid operation'): raise
      # else: apparently the GL context was destroyed already

class display_list:

  def __init__(self, index=0, owner=None):
    if (owner is None):
      assert index == 0
      self.owner = display_lists_owner(range=1)
      self.gl_index = self.owner.list
    else:
      assert index < self.owner.range
      self.owner = owner
      self.gl_index = self.owner.list + index

  def compile(self, execute=False):
    if (execute):
      mode = gl.GL_COMPILE_AND_EXECUTE
    else:
      mode = gl.GL_COMPILE
    gl.glNewList(list=self.gl_index, mode=mode)

  def end(self):
    gl.glEndList()

  def call(self):
    gl.glCallList(list=self.gl_index)

class display_lists:

  def __init__(self, range):
    self.owner = display_lists_owner(range=range)

  def __getitem__(self, index):
    return display_list(index=index, owner=self.owner)


class material_model(object):

  def __init__(self,
               front_colour=(102/255, 204/255, 1),
               back_colour=(1, 204/255, 102/255),
               ambient=0.5,
               diffuse=1.,
               specular=0.25,
               specular_focus=10):
    self.front_colour = front_colour
    self.back_colour = back_colour
    self.ambient = ambient
    self.diffuse = diffuse
    self.specular = specular
    self.specular_focus = specular_focus

  def ambient_colours(self):
    x = self.ambient
    return ([ c*x for c in self.front_colour ]+[1],
            [ c*x for c in self.back_colour ]+[1])
  ambient_colours = property(ambient_colours)

  def diffuse_colours(self):
    x = self.diffuse
    return ([ c*x for c in self.front_colour ]+[1],
            [ c*x for c in self.back_colour ]+[1])
  diffuse_colours = property(diffuse_colours)

  def specular_colour(self):
    x = self.specular
    return [x]*3 + [1]
  specular_colour = property(specular_colour)

  def execute(self, specular=True):
    from gl import glMaterialfv, glMaterialf
    fc, bc = self.ambient_colours
    glMaterialfv(gl.GL_BACK, gl.GL_AMBIENT, bc)
    glMaterialfv(gl.GL_FRONT, gl.GL_AMBIENT, fc)
    fc, bc = self.diffuse_colours
    glMaterialfv(gl.GL_BACK, gl.GL_DIFFUSE, bc)
    glMaterialfv(gl.GL_FRONT, gl.GL_DIFFUSE, fc)
    if specular:
      glMaterialfv(gl.GL_FRONT_AND_BACK, gl.GL_SPECULAR, self.specular_colour)
      glMaterialf(gl.GL_FRONT_AND_BACK, gl.GL_SHININESS, self.specular_focus)
