"""Test extract_rfactors_resolutions_sigma"""
from __future__ import absolute_import, division, print_function
from iotbx.pdb import extract_rfactors_resolutions_sigma
from libtbx.test_utils import approx_equal

example_cns = """
REMARK   2
REMARK   2 RESOLUTION. 2.20 ANGSTROMS.
REMARK   3
REMARK   3 REFINEMENT.
REMARK   3   PROGRAM     : CNS 0.9
REMARK   3   AUTHORS     : BRUNGER,ADAMS,CLORE,DELANO,GROS,GROSSE-
REMARK   3               : KUNSTLEVE,JIANG,KUSZEWSKI,NILGES, PANNU,
REMARK   3               : READ,RICE,SIMONSON,WARREN
REMARK   3
REMARK   3  REFINEMENT TARGET : ENGH & HUBER
REMARK   3
REMARK   3  DATA USED IN REFINEMENT.
REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : 2.20
REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : 99.00
REMARK   3   DATA CUTOFF            (SIGMA(F)) : 2.000
REMARK   3   DATA CUTOFF HIGH         (ABS(F)) : NULL
REMARK   3   DATA CUTOFF LOW          (ABS(F)) : NULL
REMARK   3   COMPLETENESS (WORKING+TEST)   (%) : 95.0
REMARK   3   NUMBER OF REFLECTIONS             : 22492
REMARK   3
REMARK   3  FIT TO DATA USED IN REFINEMENT.
REMARK   3   CROSS-VALIDATION METHOD          : NULL
REMARK   3   FREE R VALUE TEST SET SELECTION  : RANDOM
REMARK   3   R VALUE            (WORKING SET) : 0.171
REMARK   3   FREE R VALUE                     : 0.221
REMARK   3   FREE R VALUE TEST SET SIZE   (%) : NULL
REMARK   3   FREE R VALUE TEST SET COUNT      : 2116
REMARK   3   ESTIMATED ERROR OF FREE R VALUE  : NULL
REMARK   3
"""
example_refmac = """
REMARK   2
REMARK   2 RESOLUTION. 1.55 ANGSTROMS.
REMARK   3
REMARK   3 REFINEMENT.
REMARK   3   PROGRAM     : REFMAC
REMARK   3   AUTHORS     : MURSHUDOV,VAGIN,DODSON
REMARK   3
REMARK   3  DATA USED IN REFINEMENT.
REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : 1.55
REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : 18.00
REMARK   3   DATA CUTOFF            (SIGMA(F)) : 0.000
REMARK   3   COMPLETENESS FOR RANGE        (%) : 97.4
REMARK   3   NUMBER OF REFLECTIONS             : 10797
REMARK   3
REMARK   3  FIT TO DATA USED IN REFINEMENT.
REMARK   3   CROSS-VALIDATION METHOD          : AFTER RIGID BODY REFINEMENT
REMARK   3   FREE R VALUE TEST SET SELECTION  : RANDOM
REMARK   3   R VALUE     (WORKING + TEST SET) : 0.200
REMARK   3   R VALUE            (WORKING SET) : 0.198
REMARK   3   FREE R VALUE                     : 0.248
REMARK   3   FREE R VALUE TEST SET SIZE   (%) : 5.000
REMARK   3   FREE R VALUE TEST SET COUNT      : 520
REMARK   3
"""

example_shelx = """
REMARK   2
REMARK   2 RESOLUTION. 1.05 ANGSTROMS.
REMARK   3
REMARK   3 REFINEMENT.
REMARK   3   PROGRAM     : SHELXL-97
REMARK   3   AUTHORS     : G.M.SHELDRICK
REMARK   3
REMARK   3  DATA USED IN REFINEMENT.
REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : 1.05
REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : 20.00
REMARK   3   DATA CUTOFF            (SIGMA(F)) : 0.000
REMARK   3   COMPLETENESS FOR RANGE        (%) : 90.4
REMARK   3   CROSS-VALIDATION METHOD           : FREE R
REMARK   3   FREE R VALUE TEST SET SELECTION   : RANDOM
REMARK   3
REMARK   3  FIT TO DATA USED IN REFINEMENT (NO CUTOFF).
REMARK   3   R VALUE   (WORKING + TEST SET, NO CUTOFF) : 0.140
REMARK   3   R VALUE          (WORKING SET, NO CUTOFF) : 0.140
REMARK   3   FREE R VALUE                  (NO CUTOFF) : 0.166
REMARK   3   FREE R VALUE TEST SET SIZE (%, NO CUTOFF) : 5.400
REMARK   3   FREE R VALUE TEST SET COUNT   (NO CUTOFF) : 698
REMARK   3   TOTAL NUMBER OF REFLECTIONS   (NO CUTOFF) : 13042
REMARK   3
REMARK   3  FIT/AGREEMENT OF MODEL FOR DATA WITH F>4SIG(F).
REMARK   3   R VALUE   (WORKING + TEST SET, F>4SIG(F)) : 0.131
REMARK   3   R VALUE          (WORKING SET, F>4SIG(F)) : 0.131
REMARK   3   FREE R VALUE                  (F>4SIG(F)) : 0.158
REMARK   3   FREE R VALUE TEST SET SIZE (%, F>4SIG(F)) : 5.300
REMARK   3   FREE R VALUE TEST SET COUNT   (F>4SIG(F)) : 502
REMARK   3   TOTAL NUMBER OF REFLECTIONS   (F>4SIG(F)) : 9454
REMARK   3
"""

def run():
  result = extract_rfactors_resolutions_sigma.get_r_rfree_sigma(
    remark_2_and_3_records = example_cns.splitlines(), file_name = None)
  assert approx_equal(result.r_work      , 0.171)
  assert approx_equal(result.r_free      , 0.221)
  assert approx_equal(result.high        , 2.2  )
  assert approx_equal(result.low         , 99.0 )
  assert approx_equal(result.resolution  , 2.2  )
  assert approx_equal(result.sigma       , 2.0  )
  #
  result = extract_rfactors_resolutions_sigma.get_r_rfree_sigma(
    remark_2_and_3_records = example_refmac.splitlines(), file_name = None)
  assert approx_equal(result.r_work      , 0.198 )
  assert approx_equal(result.r_free      , 0.248 )
  assert approx_equal(result.high        , 1.55  )
  assert approx_equal(result.low         , 18.0  )
  assert approx_equal(result.resolution  , 1.55  )
  assert approx_equal(result.sigma       , 0.0   )
  #
  result = extract_rfactors_resolutions_sigma.get_r_rfree_sigma(
    remark_2_and_3_records = example_shelx.splitlines(), file_name = None)
  assert approx_equal(result.r_work      , 0.140 )
  assert approx_equal(result.r_free      , 0.166 )
  assert approx_equal(result.high        , 1.05  )
  assert approx_equal(result.low         , 20.0  )
  assert approx_equal(result.resolution  , 1.05  )
  assert approx_equal(result.sigma       , 0.0   )
  result = extract_rfactors_resolutions_sigma.extract(
    file_lines=example_shelx.splitlines(), file_name=None)
  assert approx_equal(result.r_work      , 0.140 )

if (__name__ == "__main__"):
  run()
  print("OK")
