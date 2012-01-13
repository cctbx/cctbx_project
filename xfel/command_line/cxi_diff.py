# LIBTBX_SET_DISPATCHER_NAME cxi.diff
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

import sys

from libtbx import easy_pickle

from cxi_xdr_xes.cspad_ana import cspad_tbx

def run(args):
  assert len(args) == 3
  d1 = easy_pickle.load(args[0])
  d2 = easy_pickle.load(args[1])

  image_1 = d1["DATA"]
  image_2 = d2["DATA"]

  assert image_1.all() == image_2.all()
  diff_image = image_1 - image_2
  d = cspad_tbx.dpack(
    data=diff_image,
    timestamp=cspad_tbx.evt_timestamp(),
    distance=1,
  )
  easy_pickle.dump(args[2], d)


if __name__ == '__main__':
  run(sys.argv[1:])
