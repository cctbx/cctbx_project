from __future__ import absolute_import, division, print_function

import libtbx.phil
from libtbx.phil import experimental

def exercise():
  master_phil = libtbx.phil.parse("""
first_prop = 0
  .type = int
second_prop = "default text"
  .type = str
first_scope
  .multiple = True
  .optional = True
{
  key = None
    .type = int
    .optional = False
  number = None
    .type = int
    .optional = False
  text = "Default first scope text"
    .type = str
    .optional = False
  second_scope
    .multiple = True
    .optional = True
  {
    key = None
      .type = int
      .optional = False
    flag = None
      .type = bool
    list = None
      .type = floats(size=2)
  }
}
""")

  default_str = """
first_scope {
  key = 1
  number = 1.0
  text = "First key 1 text"
  second_scope {
    key = 0
    flag = False
  }
  second_scope {
    key = 1
    list = [0, 0]
  }
  second_scope {
    key = 2
  }
}
first_scope {
  key = 2
  text = "First key 2 text"
  second_scope {
    key = 0
    flag = True
    list = [0, 1]
  }
}
"""

  overlay_str = """
first_prop = 1
first_scope {
  key = 1
  number = 2
  second_scope {
    key = 0
    list = [1, 0]
  }
  second_scope {
    key = 1
    flag = False
    list = [1, 1]
  }
  second_scope {
    key = 4
    list = [0, 2]
  }
}
"""

  # XXX Correct use of _phil and _params?
  default_phil = master_phil.fetch(sources=[libtbx.phil.parse(default_str)])
  default_params = default_phil.extract()

  overlay_phil = master_phil.fetch(sources=[libtbx.phil.parse(overlay_str)])
  overlay_params = overlay_phil.extract()

  experimental.merge_params_by_key(default_params, overlay_params, "key")

  assert default_params.first_prop == 1
  assert default_params.second_prop == "default text"
  assert len(default_params.first_scope) == 2

  for fs in default_params.first_scope:
    if (fs.key == 1):
      assert fs.number == 2 and fs.text == "Default first scope text"
      assert len(fs.second_scope) == 4
      for ss in fs.second_scope:
        if (ss.key == 0):
          assert ss.flag == False and ss.list == [1, 0]
        elif (ss.key == 1):
          assert ss.flag == False and ss.list == [1, 1]
        elif (ss.key == 2):
          assert ss.flag is None and ss.list is None
        elif (ss.key == 4):
          assert ss.flag is None and ss.list == [0, 2]

    elif (fs.key == 2):
      assert fs.number is None and fs.text == "First key 2 text"
      assert len(fs.second_scope) == 1
      assert fs.second_scope[0].key == 0 and \
          fs.second_scope[0].flag == True and \
          fs.second_scope[0].list == [0, 1]


if (__name__ == "__main__"):
  exercise()
  print("OK")
