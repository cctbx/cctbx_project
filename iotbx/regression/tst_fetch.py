from __future__ import absolute_import, division, print_function

from libtbx.utils import format_cpu_times
import requests
from iotbx.cli_parser import run_program
from iotbx.pdb.fetch import fetch, get_link
from libtbx.test_utils import assert_lines_in_text
from mmtbx.command_line.fetch_pdb import custom_args_proc
from mmtbx.programs.fetch import Program
import sys, os

def exercise_1():
  string_1yjp = fetch(id = "1yjp").read().decode()
  # print("%s" % string_1yjp)
  assert_lines_in_text(string_1yjp, """\
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
""")

def exercise_2():
  string_1yjp_sf = fetch(id = "1yjp", entity='sf').read().decode()
  print("%s" % string_1yjp_sf)
  assert_lines_in_text(string_1yjp_sf, """\
loop_
_refln.crystal_id
_refln.wavelength_id
_refln.scale_group_code
_refln.index_h
_refln.index_k
_refln.index_l
_refln.status
_refln.F_meas_au
_refln.F_meas_sigma_au
1 1 1  -12    0    2 o     3.49   2.34
1 1 1  -12    0    3 o     3.31   2.43
1 1 1  -12    0    4 o     7.32   4.25 """)

def exercise_3():
  """
  pdb: 6yvd
  EMD-10944
  maps are only 2.8 Mb and has half-maps
  """
  result = run_program(program_class=Program,
                       custom_process_arguments=custom_args_proc,
                       args=["6yvd", "action=all"])
  for i in range(len(result[0])):
    result[0][i] = os.path.basename(result[0][i])
  for fn in [ "6yvd.pdb",
              "6yvd.cif",
              "emd_10944.map.gz",
              "emd_10944_half_map_1.map.gz",
              "emd_10944_half_map_2.map.gz",
              "6yvd.fa",
              ]:
    assert os.path.isfile(fn), "File %s not found" % fn
    assert fn in result[0]

def exercise_4():
    """
    pdb: 6yvd
    EMD-10944
    Testing with action=model
    """
    # Remove relevant files before running the command
    files_to_remove = [
        "6yvd.pdb",
        "6yvd.cif",
        "emd_10944.map.gz",
        "emd_10944_half_map_1.map.gz",
        "emd_10944_half_map_2.map.gz",
        "6yvd.fa",
    ]
    for file in files_to_remove:
        if os.path.exists(file):
            os.remove(file)
    result = run_program(program_class=Program,
                         custom_process_arguments=custom_args_proc,
                         args=["6yvd", "action=model"])
    for i in range(len(result[0])):
      result[0][i] = os.path.basename(result[0][i])
    # Files that should be present with action=model
    expected_files = [
        "6yvd.pdb",
        "6yvd.cif",
    ]
    # Files that should not be present with action=model
    unexpected_files = [
        "emd_10944.map.gz",
        "emd_10944_half_map_1.map.gz",
        "emd_10944_half_map_2.map.gz",
        "6yvd.fa",
    ]
    for fn in expected_files:
      assert os.path.isfile(fn), "File {0} not found, but it should be present with action=model".format(fn)
      assert fn in result[0], "File {0} not found, but it should be present with action=model".format(fn)
    for fn in unexpected_files:
      assert not os.path.isfile(fn), "File {0} found, but it should not be present with action=model".format(fn)
      assert fn not in result[0], "File {0} found, but it should not be present with action=model".format(fn)

def exercise_5():
    """
    Testing all possible actions in fetch_pdb
    """
    actions = ['model', 'data', 'half_maps', 'sequence', 'all']
    pdb_id = '6yvd'
    files_to_remove = [
        "%s.pdb" % pdb_id,
        "%s.cif" % pdb_id,
        "emd_10944.map.gz",
        "emd_10944_half_map_1.map.gz",
        "emd_10944_half_map_2.map.gz",
        "%s.fa" % pdb_id,
    ]
    for action in actions:
      # Remove relevant files before running the command
      for file in files_to_remove:
        if os.path.exists(file):
          os.remove(file)
      result = run_program(program_class=Program,
                           custom_process_arguments=custom_args_proc,
                           args=[pdb_id, "action=%s" % action])
      for i in range(len(result[0])):
        result[0][i] = os.path.basename(result[0][i])
      if action == 'model':
        assert os.path.isfile("%s.pdb" % pdb_id)
        assert "%s.pdb" % pdb_id in result[0]
        assert os.path.isfile("%s.cif" % pdb_id)
        assert "%s.cif" % pdb_id in result[0]
        assert not os.path.isfile("emd_10944.map.gz")
        assert "emd_10944.map.gz" not in result[0]
        assert not os.path.isfile("%s.fa" % pdb_id)
        assert "%s.fa" % pdb_id not in result[0]
        assert not os.path.isfile("emd_10944_half_map_1.map.gz")
        assert "emd_10944_half_map_1.map.gz" not in result[0]
        assert not os.path.isfile("emd_10944_half_map_2.map.gz")
        assert "emd_10944_half_map_2.map.gz" not in result[0]
      elif action == 'data':
        assert os.path.isfile("emd_10944.map.gz")
        assert "emd_10944.map.gz" in result[0]
        assert not os.path.isfile("%s.pdb" % pdb_id)
        assert "%s.pdb" % pdb_id not in result[0]
        assert not os.path.isfile("%s.cif" % pdb_id)
        assert "%s.cif" % pdb_id not in result[0]
        assert not os.path.isfile("emd_10944_half_map_1.map.gz")
        assert "emd_10944_half_map_1.map.gz" not in result[0]
        assert not os.path.isfile("emd_10944_half_map_2.map.gz")
        assert "emd_10944_half_map_2.map.gz" not in result[0]
      elif action == 'sequence':
        assert os.path.isfile("%s.fa" % pdb_id)
        assert "%s.fa" % pdb_id in result[0]
        assert not os.path.isfile("%s.pdb" % pdb_id)
        assert "%s.pdb" % pdb_id not in result[0]
        assert not os.path.isfile("%s.cif" % pdb_id)
        assert "%s.cif" % pdb_id not in result[0]
        assert not os.path.isfile("emd_10944.map.gz")
        assert "emd_10944.map.gz" not in result[0]
        assert not os.path.isfile("emd_10944_half_map_1.map.gz")
        assert "emd_10944_half_map_1.map.gz" not in result[0]
        assert not os.path.isfile("emd_10944_half_map_2.map.gz")
        assert "emd_10944_half_map_2.map.gz" not in result[0]
      elif action == 'all':
        assert os.path.isfile("%s.pdb" % pdb_id)
        assert "%s.pdb" % pdb_id in result[0]
        assert os.path.isfile("%s.cif" % pdb_id)
        assert "%s.cif" % pdb_id in result[0]
        assert os.path.isfile("emd_10944.map.gz")
        assert "emd_10944.map.gz" in result[0]
        assert os.path.isfile("emd_10944_half_map_1.map.gz")
        assert "emd_10944_half_map_1.map.gz" in result [0]
        assert os.path.isfile("emd_10944_half_map_2.map.gz")
        assert "emd_10944_half_map_2.map.gz" in result[0]
        assert os.path.isfile("%s.fa" % pdb_id)
        assert "%s.fa" % pdb_id in result[0]

def exercise_get_link():
  r = []
  for ft in ['model_pdb', 'model_cif', 'sequence', 'sf', 'em_map']:
    r.append(get_link('rcsb', ft, '1ab2', emdb_number="1111"))
  # print (r)
  assert r == ['https://files.rcsb.org/pub/pdb/data/structures/divided/pdb/ab/pdb1ab2.ent.gz',
               'https://files.rcsb.org/pub/pdb/data/structures/divided/mmCIF/ab/1ab2.cif.gz',
               'https://www.rcsb.org/fasta/entry/1ab2',
               'https://files.rcsb.org/download/1ab2-sf.cif.gz',
               'https://files.rcsb.org/pub/emdb/structures/EMD-1111/map/emd_1111.map.gz'], r

if (__name__ == "__main__"):
  exercise_get_link()
  if sys.version_info.major >= 3:
    exception_occured = False
    try:
      r = requests.get('https://search.rcsb.org/')
    except Exception:
      print("OK but exception.")
      exception_occured = True
    if not exception_occured and r.ok and len(r.text) > 100:
      exercise_1()
      exercise_2()
      exercise_3()
      exercise_4()
      exercise_5()
      print("OK")
    else:
      print("OK but skipped.")
  print(format_cpu_times())
