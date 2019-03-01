
from __future__ import absolute_import, division, print_function

def run():
  import libtbx.load_env # import dependency
  for module in libtbx.env.module_list :
    print(module.name)

if (__name__ == "__main__"):
  run()
