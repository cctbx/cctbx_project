
from iotbx import map_tools
from libtbx import easy_run
import libtbx.load_env # import dependency
import os

def exercise_maps () :
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/wizards/partial_refine_001.pdb",
    test=os.path.isfile)
  if (pdb_file is None) :
    print "phenix_regression not available, skipping"
    return
  refine_file = pdb_file[:-4] + "_map_coeffs.mtz"
  easy_run.call("cp %s ." % refine_file)
  easy_run.call("rm -f *.ccp4")
  server = map_tools.server(file_name="partial_refine_001_map_coeffs.mtz")
  files = server.convert_phenix_maps(file_base="refine")
  files = [ os.path.basename(file_name) for file_name in files ]
  assert (files == ['refine_mFo-DFc.ccp4', 'refine_mFo-DFc_2.ccp4',
                    'refine_2mFo-DFc.ccp4', 'refine_2mFo-DFc_no_fill.ccp4'])
  for fn in files :
    assert os.path.isfile(fn)
  resolve_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/wizards/resolve_1_offset.mtz",
    test=os.path.isfile)
  easy_run.call("cp %s ." % resolve_file)
  server = map_tools.server(file_name="resolve_1_offset.mtz")
  map_file = server.convert_resolve_map(pdb_file=None,
    force=False)
  assert (map_file == 'resolve_1_offset.ccp4')
  assert os.path.exists(map_file)

if (__name__ == "__main__") :
  exercise_maps()
  print "OK"
