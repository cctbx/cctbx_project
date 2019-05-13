from __future__ import division
from __future__ import print_function
import os
from iotbx import mtz
import libtbx.load_env
import os.path

if (__name__ == "__main__"):
  fname = libtbx.env.find_in_repositories(
    relative_path="iotbx/regression/data/insulin_unmerged_cutted_from_ccp4.mtz",
    test=os.path.isfile)
  m = mtz.object(file_name=fname)
  m.show_summary()

  h = m.extract_miller_indices()
  j = m.extract_original_index_miller_indices()
  misym = m.extract_integers("M_ISYM")

  for idx in xrange(len(h)):
    print("asu:%17s    orig:%17s    M/ISYM:%4d"%(h[idx],j[idx],misym.data[idx]))
