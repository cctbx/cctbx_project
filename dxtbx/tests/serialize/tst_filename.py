from __future__ import absolute_import, division
from __future__ import print_function
from builtins import object
from dxtbx.serialize.filename import temp_chdir, load_path

class Test(object):

  def __init__(self):
    pass

  def run(self):
    self.tst_temp_chdir()
    self.tst_load_path()

  def tst_temp_chdir(self):
    from os import getcwd
    from os.path import join, realpath
    import libtbx
    cwd = realpath(getcwd())
    new_path = realpath(join(libtbx.env.dist_path('dxtbx'), 'serialize'))
    with temp_chdir(new_path):
      assert(realpath(getcwd()) == new_path)
    assert(realpath(getcwd()) == cwd)
    print('OK')

  def tst_load_path(self):
    import os
    from os.path import join, abspath, expanduser
    os.environ['HELLO_WORLD'] = 'EXPANDED'
    new_path = join('~', '$HELLO_WORLD', 'path')
    path = load_path(new_path)
    assert(path == join(expanduser('~'), 'EXPANDED', 'path'))
    new_path = join('$HELLO_WORLD', 'path')
    path = load_path(new_path)
    assert(path == abspath(join('EXPANDED', 'path')))
    print('OK')

if __name__ == '__main__':
  test = Test()
  test.run()
