from __future__ import absolute_import, division, print_function

def test_scan_to_string_does_not_crash_on_empty_scan():
  import dxtbx.model
  print(dxtbx.model.Scan())
