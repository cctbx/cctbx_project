from __future__ import division
from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
def run(args):
  assert len(args) == 0
  import libtbx.load_env
  from io import StringIO
  import os
  op = os.path
  t_dir = libtbx.env.under_dist(
    module_name="fable", path="test/valid", test=op.isdir)
  excluded_file_names = set("""\
blockdata_unnamed.f
""".splitlines())
  from fable.command_line import show_calls
  for file_name in os.listdir(t_dir):
    if (not file_name.endswith(".f")): continue
    if (file_name in excluded_file_names): continue
    sys.stdout = StringIO()
    show_calls.run(args=[op.join(t_dir, file_name)])
    sys.stdout = sys.__stdout__
  print("OK")

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
