from __future__ import absolute_import, division, print_function

from scitbx.array_family import flex
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from libtbx.test_utils import show_diff
from mmtbx.model.utils import merge_models
import time

cryst_info = """\
CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1      2
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.045585  0.000000  0.014006        0.00000
SCALE2      0.000000  0.205508  0.000000        0.00000
SCALE3      0.000000  0.000000  0.044560        0.00000
"""

cryst_info_2 = """\
CRYST1   21.937    5.866   23.477  90.00 107.08  90.00 P 1 21 1      2
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.045585  0.000000  0.014006        0.00000
SCALE2      0.000000  0.205508  0.000000        0.00000
SCALE3      0.000000  0.000000  0.044560        0.00000
"""

pdb_str_1 = """\
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
"""

pdb_str_2 = """\
ATOM      1  N   GLY B   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY B   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY B   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY B   1      -7.523   2.521   5.381  1.00 16.78           O
"""

pdb_str_2_1 = """\
ATOM      1  N   GLY C   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY C   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY C   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY C   1      -7.523   2.521   5.381  1.00 16.78           O
"""

pdb_str_3 = """\
MODEL        1
ATOM      1  N   GLY B   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY B   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY B   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY B   1      -7.523   2.521   5.381  1.00 16.78           O
TER
ENDMDL
MODEL        2
ATOM      1  N   GLY B   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY B   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY B   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY B   1      -7.523   2.521   5.381  1.00 16.78           O
TER
ENDMDL
"""

pdb_str_3_1 = """\
MODEL        1
ATOM      1  N   GLY C   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY C   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY C   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY C   1      -7.523   2.521   5.381  1.00 16.78           O
TER
ENDMDL
MODEL        2
ATOM      1  N   GLY D   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY D   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY D   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY D   1      -7.523   2.521   5.381  1.00 16.78           O
TER
ENDMDL
"""

pdb_str_4 = """\
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  CA  GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
"""

def tst_1():
  """
  testing assertions
  """
  m1 = mmtbx.model.manager(
      model_input = iotbx.pdb.input(lines=cryst_info+pdb_str_1, source_info=None),
      log = null_out())
  m2 = mmtbx.model.manager(
      model_input = iotbx.pdb.input(lines=cryst_info+pdb_str_2, source_info=None),
      log = null_out())
  m2_1 = mmtbx.model.manager(
      model_input = iotbx.pdb.input(lines=cryst_info_2+pdb_str_2, source_info=None),
      log = null_out())
  m2.process(make_restraints=True)
  m3 = mmtbx.model.manager(
      model_input = iotbx.pdb.input(lines=cryst_info+pdb_str_3, source_info=None),
      log = null_out())
  m4 = mmtbx.model.manager(
      model_input = iotbx.pdb.input(lines=cryst_info+pdb_str_4, source_info=None),
      log = null_out())

  try:
    new_model = merge_models([m1, m3])
  except AssertionError as e:
    assert not show_diff(str(e), """\
Cannot merge models that have differing numbers of hierarchy 'models'""")
  try:
    new_model = merge_models([m1, m2])
  except AssertionError as e:
    assert not show_diff(str(e), """\
Model #2 has GRM initialized.""")
  try:
    new_model = merge_models([m1, m2_1])
  except AssertionError as e:
    assert not show_diff(str(e), """\
Cannot merge models with different crystal symmetries (#2)""")
  try:
    new_model = merge_models([m1, m1])
  except AssertionError as e:
    assert not show_diff(str(e), """\
Duplicated chain ids are present.""")
  try:
    new_model = merge_models([m1, m4])
  except AssertionError as e:
    assert not show_diff(str(e), """\
Model #2 has erros: ['### ERROR: duplicate atom labels ###']""")

def tst_2():
  """
  Normal operation
  """
  answer = """\
CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1
SCALE1      0.045585  0.000000  0.014006        0.00000
SCALE2      0.000000  0.205508  0.000000        0.00000
SCALE3      0.000000  0.000000  0.044560        0.00000
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
TER
ATOM      5  N   GLY B   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      6  CA  GLY B   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      7  C   GLY B   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      8  O   GLY B   1      -7.523   2.521   5.381  1.00 16.78           O
TER
END
"""
  m1 = mmtbx.model.manager(
      model_input = iotbx.pdb.input(lines=cryst_info+pdb_str_1, source_info=None),
      log = null_out())
  m2 = mmtbx.model.manager(
      model_input = iotbx.pdb.input(lines=cryst_info+pdb_str_2, source_info=None),
      log = null_out())
  new_model = merge_models([m1, m2])
  assert not show_diff(new_model.model_as_pdb(), answer)

  # make sure that new_model is detached from original models
  m1.get_hierarchy().atoms().set_b(flex.double([10]*4))
  m2.get_hierarchy().atoms().set_b(flex.double([20]*4))
  assert m1.get_hierarchy().atoms().extract_b() == flex.double([10]*4)
  assert m2.get_hierarchy().atoms().extract_b() == flex.double([20]*4)
  assert new_model.get_hierarchy().atoms().extract_b() == flex.double([16.77, 16.57, 16.16, 16.78]*2)
  # print(new_model.model_as_pdb())

def tst_3():
  """
  3 models
  """
  m1 = mmtbx.model.manager(
      model_input = iotbx.pdb.input(lines=cryst_info+pdb_str_1, source_info=None),
      log = null_out())
  m2 = mmtbx.model.manager(
      model_input = iotbx.pdb.input(lines=cryst_info+pdb_str_2, source_info=None),
      log = null_out())
  m3 = mmtbx.model.manager(
      model_input = iotbx.pdb.input(lines=cryst_info+pdb_str_2_1, source_info=None),
      log = null_out())
  new_model = merge_models([m1, m2, m3])
  assert len(new_model.get_hierarchy().only_model().chains()) == 3, len(new_model.get_hierarchy().only_model().chains())

def tst_4():
  """
  2 models with 2 hierarchy models
  """
  m1 = mmtbx.model.manager(
      model_input = iotbx.pdb.input(lines=cryst_info+pdb_str_3, source_info=None),
      log = null_out())
  m2 = mmtbx.model.manager(
      model_input = iotbx.pdb.input(lines=cryst_info+pdb_str_3_1, source_info=None),
      log = null_out())
  new_model = merge_models([m1, m2])
  assert len(new_model.get_hierarchy().models()) == 2
  assert len(new_model.get_hierarchy().models()[0].chains()) == 2
  assert len(new_model.get_hierarchy().models()[1].chains()) == 2
  assert [c.id for c in new_model.get_hierarchy().models()[0].chains()] == ['B', 'C']
  assert [c.id for c in new_model.get_hierarchy().models()[1].chains()] == ['B', 'D']
  print(new_model.model_as_pdb())

if (__name__ == "__main__"):
  t0 = time.time()
  tst_1()
  tst_2()
  tst_3()
  tst_4()
  print("Total time: %-8.4f"%(time.time()-t0))
  print("OK")
