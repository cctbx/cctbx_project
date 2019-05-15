from __future__ import absolute_import, division, print_function
from gltbx import gl
from gltbx import glu
from libtbx.test_utils import Exception_expected, show_diff
import sys
from six.moves import range

def exercise_converter():
  textures = []
  gl.glGenTextures(3, textures)
  assert textures == [0,0,0]
  for i in range(10000):
    textures = []
    gl.glGenTextures(3, textures)
    assert textures == [0,0,0]
  for i in range(10000):
    textures = [9,3,5]
    gl.glGenTextures(3, textures)
    assert textures == [9,3,5]
  try:
    textures = [9,3,5]
    gl.glGenTextures(4, textures)
  except RuntimeError as e:
    assert not show_diff(str(e), """\
Argument "textures" has the wrong number of elements:
  expected size: 4
     given size: 3""")
  else: raise Exception_expected
  try:
    textures = [9,"foo",5]
    gl.glGenTextures(3, textures)
  except RuntimeError as e:
    assert not show_diff(str(e), """\
Argument "textures" has one or more elements of the wrong type.""")
  else: raise Exception_expected

def exercise_all():
  print("trying glGetError()...", end=' ')
  sys.stdout.flush()
  error = gl.glGetError()
  print("OK:", error)
  sys.stdout.flush()
  print("trying gluErrorString()...", end=' ')
  sys.stdout.flush()
  msg = glu.gluErrorString(error=error)
  print("OK:", msg)
  assert msg in ["no error", "invalid operation"]
  sys.stdout.flush()
  forever = "--forever" in sys.argv[1:]
  while True:
    exercise_converter()
    if (not forever): break

if (__name__ == "__main__"):
  exercise_all()
