from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
import iotbx.pdb
from cctbx.maptbx.bcr import qmap
from cctbx import maptbx
import time

import boost_adaptbx.boost.python as bp
ext = bp.import_ext("cctbx_maptbx_bcr_bcr_ext")

pdb_str = """
remark CRYST1   20.000   25.000   30.000  90.00  90.00  90.00 P 1
CRYST1   20.000   25.000   30.000  96.06  99.11  91.14 P 1
ATOM      1  N   HIS A 109      10.000  14.000   9.000  1.10 30.10           N
ATOM      2  C   HIS A 109      13.000  10.000  10.000  0.90 20.00           C
ATOM      3  O   HIS A 109      12.000  12.000  11.000  0.50 10.00           O
TER
END
"""

def scale(x,y):
  scale = flex.sum(x*y)/flex.sum(y*y)
  return x, y*scale

def run(d_min, pdb_file, RadFact, RadAdd, table = "wk1995"):
  pdb_inp = iotbx.pdb.input(file_name=pdb_file)
  #pdb_inp = iotbx.pdb.input(source_info=None, lines = pdb_str)
  cs = pdb_inp.crystal_symmetry()
  #pdb_inp.write_pdb_file(file_name="m.pdb")
  xrs = pdb_inp.xray_structure_simple()
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell        = cs.unit_cell(),
    space_group_info = cs.space_group_info(),
    symmetry_flags   = maptbx.use_space_group_symmetry,
    step             = d_min/2.
    )
  n_real = crystal_gridding.n_real()
  xrs = pdb_inp.xray_structure_simple()
  xrs.scattering_type_registry(
    table = table,
    d_min = d_min,
    types_without_a_scattering_contribution=["?"])
  #
  # FFT
  #
  start = time.perf_counter()
  fc = xrs.structure_factors(d_min=d_min, algorithm="fft").f_calc()
  fft_map = fc.fft_map(
    crystal_gridding = crystal_gridding,
    f_000            = xrs.f_000())
  fft_map.apply_volume_scaling()
  m1 = fft_map.real_map_unpadded()
  ttimeFFT = time.perf_counter()-start
  #
  # VRM
  #
  o = qmap.compute(
    xray_structure = xrs,
    n_real         = n_real,
    resolution     = d_min,
    RadFact        = RadFact,
    RadAdd         = RadAdd,
    resolutions    = None,
    use_exp_table  = True, # XXX <<<<
    debug          = False,
    show_BCR       = True)
  OmegaMap = o.map_data()
  ttimeVRM = o.time_map
  #
  m1, m2 = m1.as_1d(),OmegaMap.as_1d()
  #
  cc_all = flex.linear_correlation(x=m1, y=m2).coefficient()
  #
  s = m2==0.0
  m1, m2 = m1.select(~s), m2.select(~s)
  cc_mol = flex.linear_correlation(x=m1, y=m2).coefficient()
  m1, m2 = scale(x=m1, y=m2)
  diff = m1-m2
  mi,ma,me = flex.min(diff), flex.max(diff), flex.mean(diff)
  miFFT,maFFT,meFFT = flex.min(m1), flex.max(m1), flex.mean(m1)
  miVRM,maVRM,meVRM = flex.min(m2), flex.max(m2), flex.mean(m2)
  #
  #
  its = ["CC all/mol: %5.3f %5.3f",
         "DIFF min/max/mean: %6.3f %6.3f %6.3f",
         "tFFT: %6.3f",
         "tVRM: %6.3f",
         "tFFT/tVRM: %6.3f",
         "mFFT min/max/mean: %6.3f %6.3f %6.3f",
         "mVRM min/max/mean: %6.3f %6.3f %6.3f",
         ]
  f = " ".join(its)
  return f%(cc_all, cc_mol, mi,ma,me, ttimeFFT, ttimeVRM, ttimeFFT/ttimeVRM,
            miFFT,maFFT,meFFT,
            miVRM,maVRM,meVRM)


if (__name__ == "__main__"):
  start = time.perf_counter()
  files = ["1crn_box_000H.pdb",]
  for it in [
      [[0.3, 0.5, 1.0], [1, 2, 3, 4, 5, 6, 7, 8]],
      [[2, 3, 4],       [1, 2, 3, 4, 5, 6, 7, 8]],
      [[5, 6],          [0.25, 0.5, 1, 2]]
      ]:
    d_mins, RadFact_s = it
    for pdb_file in files:
      print(pdb_file)
      for RadFact in RadFact_s:
        for RadAdd in [0, 0.5]:
          print("  RadFact, RadAdd:", RadFact, RadAdd)
          for d_min in d_mins:
            v = run(d_min    = d_min,
                    pdb_file = pdb_file,
                    RadFact  = RadFact,
                    RadAdd   = RadAdd)
            print("    d_min: %3s"%str(d_min), v)
  print("Time:", time.perf_counter()-start)
  print("OK")
