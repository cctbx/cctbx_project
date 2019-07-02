from __future__ import absolute_import, division, print_function
import os
from iotbx import mtz
import libtbx.load_env
import os.path

if (__name__ == "__main__"):
  fname = libtbx.env.find_in_repositories(
    relative_path="iotbx/regression/data/insulin_unmerged_cutted_from_ccp4.mtz",
    test=os.path.isfile)
  m = mtz.object(file_name=fname)

  miller_dict = m.as_miller_arrays_dict()

  keys = sorted(miller_dict.keys())

  assert keys == [('HKL_base', 'HKL_base', 'BATCH'),
                  ('HKL_base', 'HKL_base', 'BGPKRATIOS'),
                  ('HKL_base', 'HKL_base', 'FLAG'),
                  ('HKL_base', 'HKL_base', 'FRACTIONCALC'),
                  ('HKL_base', 'HKL_base', 'I'),
                  ('HKL_base', 'HKL_base', 'IPR'),
                  ('HKL_base', 'HKL_base', 'LP'),
                  ('HKL_base', 'HKL_base', 'MPART'),
                  ('HKL_base', 'HKL_base', 'M_ISYM'),
                  ('HKL_base', 'HKL_base', 'ROT'),
                  ('HKL_base', 'HKL_base', 'SIGI'),
                  ('HKL_base', 'HKL_base', 'SIGIPR'),
                  ('HKL_base', 'HKL_base', 'WIDTH'),
                  ('HKL_base', 'HKL_base', 'XDET'),
                  ('HKL_base', 'HKL_base', 'YDET')]
