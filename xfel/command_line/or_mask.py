from __future__ import division
from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME cxi.or_mask
from libtbx import easy_pickle
from scitbx.array_family import flex
import sys,os
files = sys.argv[1:]
srcs = files[:-1]
dest = files[-1]
print(srcs,">",dest)

data1= easy_pickle.load(srcs[0])
ddata = data1["DATA"]
try:
  idata1 = ddata.iround()
except AttributeError:
  idata1 = ddata.as_double().iround()
discover_mask_pix_val = flex.sum(idata1.as_long())//(idata1!=0).count(True)
print("I think the mask pixel value is", discover_mask_pix_val)
bdata1 = (idata1!=0)
for item in srcs[1:]:
  dataN= easy_pickle.load(item)
  try:
    bdata1 |= (dataN["DATA"].iround()!=0)
  except AttributeError:
    bdata1 |= (dataN["DATA"].as_double().iround()!=0)
dirname = os.path.dirname(dest)
if dirname is not "" and not os.path.isdir(dirname):
  os.makedirs(dirname)
debug_fixer = flex.bool(list(bdata1)).as_int().as_double()*discover_mask_pix_val
debug_fixer.reshape(flex.grid(ddata.focus()))
data1["DATA"]=debug_fixer
easy_pickle.dump(dest,data1)
