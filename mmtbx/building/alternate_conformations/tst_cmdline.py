
from __future__ import division
from mmtbx.building.alternate_conformations import tst_build_simple
from mmtbx.command_line import build_alternate_conformations
from mmtbx.validation import rotalyze
from iotbx import file_reader
from libtbx.utils import null_out
import libtbx.load_env
import warnings
import os.path

def exercise () :
  tst_build_simple.prepare_inputs("tst_altconfs_cmdline")
  assert os.path.isfile("tst_altconfs_cmdline_start.pdb")
  assert os.path.isfile("tst_altconfs_cmdline.mtz")
  args = [
    "tst_altconfs_cmdline_start.pdb",
    "tst_altconfs_cmdline.mtz",
    "selection='chain A and resseq 3'",
    "output.file_name=tst_altconfs_cmdline_out.pdb",
    "expected_occupancy=0.4",
    "window_size=2",
    "nproc=1",
    "rsr_fofc_map_target=False",
    "create_dir=False",
    "prefix=tst_altconfs_cmdline_out",
    "--verbose",
  ]
  build_alternate_conformations.run(args=args, out=null_out())
  pdb_out = file_reader.any_file("tst_altconfs_cmdline_out.pdb")
  hierarchy =  pdb_out.file_object.construct_hierarchy()
  validate = rotalyze.rotalyze(pdb_hierarchy=hierarchy,
    data_version="8000",
    outliers_only=False)
  rota_out = [ (r.id_str(), r.rotamer_name) for r in validate.results ]
  #print rota_in
  assert (((' A   3 AASN', 't0') in rota_out) and
          ((' A   3 BASN', 'm-40') in rota_out))
  #assert (rota_out == rota_in)

if (__name__ == "__main__") :
  if (not libtbx.env.has_module("phenix")) :
    warnings.warn("phenix missing, skipping this test")
  else :
    exercise()
    print "OK"
