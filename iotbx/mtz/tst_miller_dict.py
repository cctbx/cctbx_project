from __future__ import division
import os
from iotbx import mtz

if (__name__ == "__main__") :
  m = mtz.object(file_name=os.path.join(os.environ["CEXAM"],"data",
                           "insulin_unmerged.mtz"))
  
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
