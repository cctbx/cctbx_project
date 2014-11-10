from __future__ import division
import os
from iotbx import mtz

if (__name__ == "__main__") :
  m = mtz.object(file_name=os.path.join(os.environ["CEXAM"],"data",
                           "insulin_unmerged.mtz"))
  m.show_summary()

  h = m.extract_miller_indices()
  j = m.extract_original_index_miller_indices()
  misym = m.extract_integers("M_ISYM")

  for idx in xrange(len(h)):
    print "asu:%17s    orig:%17s    M/ISYM:%4d"%(h[idx],j[idx],misym.data[idx])
