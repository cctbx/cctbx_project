from __future__ import division
import os

import libtbx.load_env
from libtbx.test_utils import approx_equal, show_diff
import iotbx.pdb


def exercise_mmcif_tls():
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/3orl.pdb",
    test=os.path.isfile)
  mmcif_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/3orl.cif",
    test=os.path.isfile)

  if pdb_file is None or mmcif_file is None:
    print "Skipping exercise_mmcif_tls(): missing phenix_regression directory."
    return

  pdb_input = iotbx.pdb.input(file_name=pdb_file)
  pdb_hierarchy = pdb_input.construct_hierarchy()
  cif_input = iotbx.pdb.input(file_name=mmcif_file)
  cif_hierarchy = cif_input.construct_hierarchy()

  pdb_tls_params = pdb_input.extract_tls_params(pdb_hierarchy).tls_params

  cif_block = cif_input.cif_block
  cif_tls_params = cif_input.extract_tls_params(cif_hierarchy).tls_params

  assert len(pdb_tls_params) == len(cif_tls_params) == 3
  for pdb_tls, cif_tls in zip(pdb_tls_params, cif_tls_params):
    assert approx_equal(pdb_tls.t, cif_tls.t)
    assert approx_equal(pdb_tls.l, cif_tls.l)
    assert approx_equal(pdb_tls.s, cif_tls.s)
    assert approx_equal(pdb_tls.origin, cif_tls.origin)
    assert not show_diff(pdb_tls.selection_string, cif_tls.selection_string)

  # this one has phenix selection strings
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/4g9h.pdb",
    test=os.path.isfile)
  mmcif_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/4g9h.cif",
    test=os.path.isfile)

  pdb_input = iotbx.pdb.input(file_name=pdb_file)
  pdb_hierarchy = pdb_input.construct_hierarchy()
  pdb_tls_params = pdb_input.extract_tls_params(pdb_hierarchy).tls_params

  cif_input = iotbx.pdb.input(file_name=mmcif_file)
  cif_hierarchy = cif_input.construct_hierarchy()
  cif_block = cif_input.cif_block
  cif_tls_params = cif_input.extract_tls_params(cif_hierarchy).tls_params

  assert len(pdb_tls_params) == len(cif_tls_params) == 18
  for pdb_tls, cif_tls in zip(pdb_tls_params, cif_tls_params):
    assert approx_equal(pdb_tls.t, cif_tls.t)
    assert approx_equal(pdb_tls.l, cif_tls.l)
    assert approx_equal(pdb_tls.s, cif_tls.s)
    assert approx_equal(pdb_tls.origin, cif_tls.origin)
    assert not show_diff(pdb_tls.selection_string, cif_tls.selection_string)

  # in this one the tls data items are not looped
  mmcif_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/2xw9.cif",
    test=os.path.isfile)
  cif_input = iotbx.pdb.input(file_name=mmcif_file)
  cif_hierarchy = cif_input.construct_hierarchy()

  cif_block = cif_input.cif_block
  cif_tls_params = cif_input.extract_tls_params(cif_hierarchy).tls_params

  assert len(cif_tls_params) == 1
  cif_tls = cif_tls_params[0]
  assert approx_equal(
    cif_tls.t, [0.0275, 0.0202, 0.0138, -0.0004, 0.0088, -0.0002])
  assert approx_equal(
    cif_tls.l, [0.0554, 0.0231, 0.0573, -0.0127, 0.0112, -0.017])
  assert approx_equal(
    cif_tls.s, [-0.0001, -0.0012, -0.0037, -0.0006, 0.001, 0.0007, -0.0023, -0.0001, -0.0009])
  assert approx_equal(cif_tls.origin, [-1.219, 1.557, 13.138])
  assert approx_equal(cif_tls.selection_string, "(chain A and resseq 1:228)")


def run():
  exercise_mmcif_tls()
  print "OK"

if __name__ == '__main__':
  run()
