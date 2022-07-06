from __future__ import absolute_import, division, print_function
from libtbx.test_utils import assert_lines_in_text
import libtbx.load_env
from six.moves import cStringIO as StringIO
import os, sys

def run_function_and_grab_stdout(func, params):
  saved_stdout = sys.stdout
  inp = StringIO()
  sys.stdout = inp
  func(params)
  stdout_text = inp.getvalue()
  sys.stdout = saved_stdout
  return stdout_text

def exercise_pdb_shake():
  from iotbx.examples.pdb_shake import run
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="iotbx/regression/data/2ERL.pdb",
    test=os.path.isfile)
  stdout_text = run_function_and_grab_stdout(run, [pdb_file])
  assert_lines_in_text(stdout_text, 'END')


def exercise_pdb_to_map_simple():
  from iotbx.examples.pdb_to_map_simple import run
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="iotbx/regression/data/pdb1zff.ent",
    test=os.path.isfile)
  stdout_text = run_function_and_grab_stdout(run, [pdb_file])
  # print(stdout_text)
  assert_lines_in_text(stdout_text, 'block last:  (26, 17, 26)')

def exercise_shelx_latt_sym_to_space_group_symbol():
  from iotbx.examples.shelx_latt_sym_to_space_group_symbol import run
  fname = libtbx.env.find_in_repositories(
    relative_path="iotbx/regression/data/shelx.ins",
    test=os.path.isfile)
  stdout_text = run_function_and_grab_stdout(run, [fname])
  # print(stdout_text)
  assert_lines_in_text(stdout_text, "P 1 21 1")

def exercise_direct_methods_light():
  from iotbx.examples.direct_methods_light import run
  data_fname = libtbx.env.find_in_repositories(
    relative_path="iotbx/regression/data/cu3182Isup2.hkl",
    test=os.path.isfile)
  model_fname = libtbx.env.find_in_repositories(
    relative_path="iotbx/regression/data/cu3182sup1.cif",
    test=os.path.isfile)
  run_function_and_grab_stdout(run, [data_fname, model_fname])

def exercise_pdb_symmetry_copies():
  from iotbx.examples.pdb_symmetry_copies import run
  fname = libtbx.env.find_in_repositories(
    relative_path="iotbx/regression/data/1yjp.pdb",
    test=os.path.isfile)
  stdout_text = run_function_and_grab_stdout(run, [fname])
  # print(stdout_text)
  assert_lines_in_text(stdout_text, "Writing file: symmetry_copies.pdb")
  assert os.path.isfile("symmetry_copies.pdb")

if (__name__ == "__main__"):
  exercise_pdb_shake()
  exercise_pdb_to_map_simple()
  exercise_shelx_latt_sym_to_space_group_symbol()
  exercise_direct_methods_light()
  exercise_pdb_symmetry_copies()
  print("OK")
