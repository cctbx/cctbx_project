from __future__ import absolute_import, division, print_function
from mmtbx.validation import clashscore
from iotbx import pdb
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
from libtbx.easy_pickle import loads, dumps
import libtbx.load_env
import time
import os

#protein
pdb_str_1 = """
ATOM    556  N   LEU A  71      32.763  35.831  23.090  1.00 12.71           N
ATOM    557  CA  LEU A  71      34.145  35.472  23.481  1.00 16.06           C
ATOM    558  C   LEU A  71      34.239  35.353  24.979  1.00 18.09           C
ATOM    559  O   LEU A  71      33.707  36.197  25.728  1.00 19.26           O
ATOM    560  CB  LEU A  71      35.114  36.564  22.907  1.00 17.10           C
ATOM    561  CG  LEU A  71      35.926  35.979  21.737  1.00 19.37           C
ATOM    562  CD1 LEU A  71      35.003  35.084  20.920  1.00 17.51           C
ATOM    563  CD2 LEU A  71      36.533  37.087  20.917  1.00 19.57           C
ATOM    564  N   ARG A  72      34.930  34.384  25.451  1.00 21.47           N
ATOM    565  CA  ARG A  72      35.161  34.174  26.896  1.00 25.83           C
ATOM    566  C   ARG A  72      36.671  34.296  27.089  1.00 27.74           C
ATOM    567  O   ARG A  72      37.305  33.233  26.795  1.00 30.65           O
ATOM    568  CB  ARG A  72      34.717  32.760  27.286  1.00 28.49           C
ATOM    569  CG  ARG A  72      35.752  32.054  28.160  1.00 31.79           C
ATOM    570  CD  ARG A  72      35.612  30.577  28.044  1.00 34.05           C
ATOM    571  NE  ARG A  72      35.040  30.252  26.730  1.00 35.08           N
ATOM    572  CZ  ARG A  72      34.338  29.103  26.650  1.00 34.67           C
ATOM    573  NH1 ARG A  72      34.110  28.437  27.768  1.00 35.02           N
ATOM    574  NH2 ARG A  72      34.014  28.657  25.457  1.00 34.97           N
ATOM    575  N   LEU A  73      37.197  35.397  27.513  0.45 28.93           N
ATOM    576  CA  LEU A  73      38.668  35.502  27.680  0.45 30.76           C
ATOM    577  C   LEU A  73      39.076  34.931  29.031  0.45 32.18           C
ATOM    578  O   LEU A  73      38.297  34.946  29.996  0.45 32.31           O
ATOM    579  CB  LEU A  73      39.080  36.941  27.406  0.45 30.53           C
ATOM    580  CG  LEU A  73      39.502  37.340  26.002  0.45 30.16           C
ATOM    581  CD1 LEU A  73      38.684  36.647  24.923  0.45 29.57           C
ATOM    582  CD2 LEU A  73      39.337  38.854  25.862  0.45 29.11           C
ATOM    583  N   ARG A  74      40.294  34.412  29.045  0.45 33.82           N
ATOM    584  CA  ARG A  74      40.873  33.802  30.253  0.45 35.33           C
ATOM    585  C   ARG A  74      41.765  34.829  30.944  0.45 36.22           C
ATOM    586  O   ARG A  74      42.945  34.994  30.583  0.45 36.70           O
ATOM    587  CB  ARG A  74      41.651  32.529  29.923  0.45 36.91           C
ATOM    588  CG  ARG A  74      41.608  31.444  30.989  0.45 38.62           C
ATOM    589  CD  ARG A  74      41.896  30.080  30.456  0.45 39.75           C
ATOM    590  NE  ARG A  74      43.311  29.735  30.563  0.45 41.13           N
ATOM    591  CZ  ARG A  74      44.174  29.905  29.554  0.45 41.91           C
ATOM    592  NH1 ARG A  74      43.754  30.312  28.356  0.45 42.75           N
ATOM    593  NH2 ARG A  74      45.477  29.726  29.763  0.45 41.93           N
END
"""

def exercise_clashscore_old():
  if (not libtbx.env.has_module(name="probe")):
    print("Skipping exercise_clashscore(): probe not configured")
    return

  pdb_io = pdb.input(source_info=None, lines=pdb_str_1)
  pdb_hierarchy = pdb_io.construct_hierarchy()
  cs = clashscore.clashscore(pdb_hierarchy=pdb_hierarchy, out=null_out())
  for unpickle in [False, True] :
    if (unpickle):
      cs = loads(dumps(cs))
    c_score = cs.get_clashscore()
    assert approx_equal(c_score, 35.29, eps=0.01)
    bad_clashes_list = cs.results
    assert ([ c.format_old() for c in bad_clashes_list ] ==
      [' A  72  ARG  HG2  A  72  ARG  O   :1.038',
       ' A  72  ARG  CG   A  72  ARG  O   :0.465',
       ' A  71  LEU  HA   A  71  LEU HD12 :0.446']), [ c.format_old() for c in bad_clashes_list ]

  #test nuclear distances
  cs = clashscore.clashscore(pdb_hierarchy=pdb_hierarchy, nuclear=True)
  for unpickle in [False, True] :
    if (unpickle):
      cs = loads(dumps(cs))
    c_score = cs.get_clashscore()
    assert approx_equal(c_score, 58.82, eps=0.01)
    bad_clashes_list = cs.results
    assert ([ c.format_old() for c in bad_clashes_list ] ==
      [ ' A  72  ARG  HG2  A  72  ARG  O   :1.082',
        ' A  72  ARG  CG   A  72  ARG  O   :0.622',
        ' A  71  LEU  HA   A  71  LEU HD12 :0.535',
        ' A  72  ARG  HB3  A  72  ARG  HE  :0.475',
        ' A  72  ARG  HD3  A  72  ARG HH11 :0.451'])

  #test B factor cutoff
  cs = clashscore.clashscore(pdb_hierarchy=pdb_hierarchy, b_factor_cutoff=40)
  for unpickle in [False, True] :
    if (unpickle):
      cs = loads(dumps(cs))
    c_score = cs.get_clashscore()
    assert approx_equal(c_score, 35.29, eps=0.01)
    c_score_b_cutoff = cs.get_clashscore_b_cutoff()
    assert approx_equal(c_score_b_cutoff, 39.47, eps=0.01)
    bad_clashes_list = cs.results
    assert ([ c.format_old() for c in bad_clashes_list ] ==
      [' A  72  ARG  HG2  A  72  ARG  O   :1.038',
       ' A  72  ARG  CG   A  72  ARG  O   :0.465',
       ' A  71  LEU  HA   A  71  LEU HD12 :0.446'])

def exercise_clashscore():
  if (not libtbx.env.has_module(name="probe")):
    print("Skipping exercise_clashscore(): probe not configured")
    return

  pdb_io = pdb.input(source_info=None, lines=pdb_str_1)
  pdb_hierarchy = pdb_io.construct_hierarchy()
  cs = clashscore.clashscore(
      pdb_hierarchy=pdb_hierarchy,
      fast = False,
      condensed_probe=True,
      out=null_out())
  for unpickle in [False, True]:
    if unpickle:
      cs = loads(dumps(cs))
    c_score = cs.get_clashscore()
    assert approx_equal(c_score, 35.29, eps=0.01)
    bad_clashes_list = cs.results
    assert ([ c.format_old() for c in bad_clashes_list ] ==
      [' A  72  ARG  HG2  A  72  ARG  O   :1.048',
      ' A  71  LEU  HA   A  71  LEU HD12 :0.768',
      ' A  72  ARG  CG   A  72  ARG  O   :0.720']), [ c.format_old() for c in bad_clashes_list ]

  #test nuclear distances
  cs = clashscore.clashscore(
      pdb_hierarchy=pdb_hierarchy,
      fast = False,
      condensed_probe=True,
      nuclear=True)
  for unpickle in [False, True] :
    if (unpickle):
      cs = loads(dumps(cs))
    c_score = cs.get_clashscore()
    assert approx_equal(c_score, 58.82, eps=0.01)
    bad_clashes_list = cs.results
    assert ([ c.format_old() for c in bad_clashes_list ] ==
      [' A  72  ARG  HG2  A  72  ARG  O   :1.085',
      ' A  71  LEU  HA   A  71  LEU HD12 :0.793',
      ' A  72  ARG  CG   A  72  ARG  O   :0.720',
      ' A  72  ARG  HD3  A  72  ARG HH11 :0.669',
      ' A  72  ARG  HB3  A  72  ARG  HE  :0.647']), [ c.format_old() for c in bad_clashes_list ]

  #test B factor cutoff
  cs = clashscore.clashscore(
      pdb_hierarchy=pdb_hierarchy,
      fast = False,
      condensed_probe=True,
      b_factor_cutoff=40)
  for unpickle in [False, True] :
    if (unpickle):
      cs = loads(dumps(cs))
    c_score = cs.get_clashscore()
    assert approx_equal(c_score, 35.29, eps=0.01)
    c_score_b_cutoff = cs.get_clashscore_b_cutoff()
    assert approx_equal(c_score_b_cutoff, 39.47, eps=0.01)
    bad_clashes_list = cs.results
    assert ([ c.format_old() for c in bad_clashes_list ] ==
      [' A  72  ARG  HG2  A  72  ARG  O   :1.048',
      ' A  71  LEU  HA   A  71  LEU HD12 :0.768',
      ' A  72  ARG  CG   A  72  ARG  O   :0.720']), [ c.format_old() for c in bad_clashes_list ]

# TODO
def exercise_full_validation():
  from phenix.validation import analyze_all
  import iotbx.phil
  open("tmp_validation_neutron.pdb", "w").write(pdb_str_1)
  test_phil = iotbx.phil.parse("""
model_vs_data {
  pdb_file = tmp_validation_neutron.pdb
  scattering_table = *neutron
}
""")
  working_phil = analyze_all.model_vs_data_params.fetch(source=test_phil)
  params = working_phil.extract()
  validation = analyze_all.validation_result(
    params=params,
    tmp_dir=os.getcwd(),
    out=null_out(),
    quiet=True)
  c_score = validation.molprobity_result.get_clashscore()
  assert approx_equal(c_score, 58.82, eps=0.01)
  params.model_vs_data.scattering_table = "n_gaussian"
  validation = analyze_all.validation_result(
    params=params,
    tmp_dir=os.getcwd(),
    out=null_out(),
    quiet=True)
  c_score = validation.molprobity_result.get_clashscore()
  assert approx_equal(c_score, 35.29, eps=0.01)

if (__name__ == "__main__"):
  t0 = time.time()
  exercise_clashscore_old()
  exercise_clashscore()
  #exercise_full_validation()
  print("OK. Time: %8.3f"%(time.time()-t0))
