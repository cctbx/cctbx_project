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
