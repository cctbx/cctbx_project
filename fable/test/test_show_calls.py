from __future__ import absolute_import, division, print_function

import os

def test_that_show_calls_does_not_crash(testsdir):
  t_dir = os.path.join(testsdir, 'valid')
  files = [f for f in os.listdir(t_dir) if f.endswith('.f') and f != 'blockdata_unnamed.f']
  assert files
  import fable.command_line.show_calls
  for filename in files:
    fable.command_line.show_calls.run(args=[os.path.join(t_dir, filename)])
