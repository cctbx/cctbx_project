from iotbx import pdb
import iotbx.pdb.interpretation
import iotbx.pdb.remark_3_interpretation
from cctbx.array_family import flex
from libtbx.utils import format_cpu_times
import sys, math, time, os
from mmtbx.tls.tls import *
from libtbx.test_utils import approx_equal
import libtbx.load_env
from cctbx import adptbx


def uaniso_from_tls_and_back():
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="regression/pdb/1OC2_tst.pdb", test=os.path.isfile)
  if (pdb_file is None):
    print "Skipping uaniso_from_tls_and_back(): input file not available"
    return
  xray_structure = pdb.as_xray_structure(file_name = pdb_file)
  stage_1 = pdb.interpretation.stage_1(file_name   = pdb_file)

  input_tls_data = iotbx.pdb.remark_3_interpretation.extract_tls_parameters(
                                                      stage_1.remark_3_records)

  input_tls_data.print_as_matrices()
  tls_obj = tls(tls_parameters = input_tls_data,
                xray_structure = xray_structure,
                stage_1        = stage_1)
  u_from_tls = tls_obj.u_from_tls()

  #i = 0
  for utls,atom in zip(u_from_tls, stage_1.atom_attributes_list):
    updb = atom.Ucart
    #i += 1
    #print "      ", i
    #print "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f "%(utls[0],utls[1],utls[2],utls[3],utls[4],utls[5])
    #print "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f "%(updb[0],updb[1],updb[2],updb[3],updb[4],updb[5])
    assert approx_equal(utls,updb, 1.e-4)

  tls_from_u = tls_obj.tls_from_u()

  assert approx_equal(input_tls_data.T, tls_from_u.T, 1.e-4)
  assert approx_equal(input_tls_data.L, tls_from_u.L, 1.e-4)
  assert approx_equal(input_tls_data.S, tls_from_u.S, 1.e-4)
  assert approx_equal(input_tls_data.origin, tls_from_u.origin, 1.e-6)

if (__name__ == "__main__"):
  uaniso_from_tls_and_back()
  print format_cpu_times()
