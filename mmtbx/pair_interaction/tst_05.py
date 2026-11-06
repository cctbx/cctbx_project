from __future__ import absolute_import, division, print_function
"""
load_wfc: <el>_wfc_obj.pkl in the mmtbx/pair_interaction folder (or the
folder given) is loaded; a missing element is parsed from the yoink text
tables and its pickle written into that folder; nothing is written into the
current directory; objects are cached per process.
"""
import os, tempfile, stat
import libtbx.load_env
from mmtbx.pair_interaction import pair_interaction

def run():
  if libtbx.env.find_in_repositories("qrefine/plugin/yoink/dat") is None:
    return # the text tables are in qrefine
  tmp = tempfile.mkdtemp()
  cwd_before = sorted(os.listdir("."))
  pair_interaction._WFC_CACHE.clear()
  o1 = pair_interaction.load_wfc("xe", folder=tmp)        # parsed, written
  assert os.path.isfile(os.path.join(tmp, "xe_wfc_obj.pkl"))
  assert sorted(os.listdir(".")) == cwd_before              # nothing in cwd
  pair_interaction._WFC_CACHE.clear()
  o2 = pair_interaction.load_wfc("xe", folder=tmp)        # loaded
  assert o1.ngrid == o2.ngrid and o1.zz == o2.zz == 54, (o2.ngrid, o2.zz)
  assert pair_interaction.load_wfc("xe", folder=tmp) is o2 # cached
  # a folder that cannot be written: parsed, not written, no exception
  ro = os.path.join(tmp, "ro")
  os.mkdir(ro)
  os.chmod(ro, stat.S_IRUSR | stat.S_IXUSR)
  pair_interaction._WFC_CACHE.clear()
  assert pair_interaction.load_wfc("kr", folder=ro).zz == 36
  assert os.listdir(ro) == []
  os.chmod(ro, stat.S_IRWXU)
  # the repository pickles, one- and two-letter elements
  pair_interaction._WFC_CACHE.clear()
  assert pair_interaction.load_wfc("c_").zz == 6
  assert pair_interaction.load_wfc("c").zz == 6
  assert pair_interaction.load_wfc("cl").zz == 17
  assert pair_interaction.load_wfc("Cl").zz == 17
  assert sorted(os.listdir(".")) == cwd_before

if(__name__=="__main__"):
  run()
  print("OK")
