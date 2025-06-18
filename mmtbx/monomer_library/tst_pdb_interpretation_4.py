from __future__ import absolute_import, division, print_function
from mmtbx.monomer_library import pdb_interpretation, server
from libtbx.utils import format_cpu_times
import sys


def exercise_long_vdw(mon_lib_srv, ener_lib):
  """
  Exercise breaking GRM when ion does not have charge, therefore
  assigned with 99 in pdb_interpretation:
  self.charges[i_seq] = 99 in def assign_from_monomer_mapping(self, conf_altloc, mm)

  At the same time in GRM's pair_proxies()
  if (self.nonbonded_distance_cutoff < max_vdw_dist):
  the above values are not consistently calculated. The
  self.nonbonded_distance_cutoff is essentially =
  self.nonbonded_params.find_max_vdw_distance(nonbonded_types=self.nonbonded_types)
  calculated without charges, while
  max_vdw_dist = self._pair_proxies.nonbonded_proxies.max_vdw_distance
  calculated in nonbonded_proxies with charges taken into account.
  Therefore the above self.charges[i_seq] = 99 makes difference and results in
  inconsistency.

  See the calls to the same function get_nonbonded_distance() in
  cctbx_project/cctbx/geometry_restraints/nonbonded.h from
    find_max_vdw_distance in the same file, no charges passed and
    make_nonbonded_asu_proxy in nonbonded_sorted.h, with charges.

  Example is from 7eax.
  """

  raw_records = """\
CRYST1   68.876   72.551   85.067  90.00  90.00  90.00 P 1 21 1
SCALE1      0.014519  0.000000  0.000000        0.00000
SCALE2      0.000000  0.013783  0.000000        0.00000
SCALE3      0.000000  0.000000  0.011755        0.00000
ATOM      1  CA  CYS A 124     -37.266  -1.426  -1.806  1.00 38.06           C
ATOM      2  C   CYS A 124     -38.616  -1.270  -2.515  1.00 38.16           C
ATOM      3  CB  CYS A 124     -36.251  -2.133  -2.694  1.00 40.42           C
ATOM      4  SG  CYS A 124     -36.398  -1.704  -4.443  1.00 43.45           S
ATOM      5  N   THR A 125     -39.009  -0.033  -2.825  1.00 37.76           N
ATOM      6  O   THR A 125     -39.229   2.353  -4.155  1.00 35.41           O
ATOM      7  C   MET A 133     -36.674   5.017  -4.524  1.00 39.62           C
ATOM      8  CB  MET A 133     -35.806   3.991  -6.649  1.00 44.00           C
ATOM      9  CG  MET A 133     -35.545   2.619  -6.073  1.00 46.53           C
ATOM     10  SD  MET A 133     -36.737   1.387  -6.638  1.00 52.32           S
ATOM     11  CE  MET A 133     -35.903   0.786  -8.101  1.00 51.94           C
ATOM     12  N   PHE A 134     -37.217   4.201  -3.618  1.00 39.78           N
ATOM     13  CA  PHE A 134     -37.019   4.279  -2.147  1.00 37.97           C
ATOM     14  C   PHE A 134     -36.414   2.955  -1.667  1.00 38.92           C
ATOM     15  O   PHE A 134     -37.089   1.920  -1.790  1.00 36.23           O
ATOM     16  N   CYS A 135     -35.177   2.986  -1.165  1.00 42.43           N
ATOM     17  CA  CYS A 135     -34.429   1.792  -0.686  1.00 44.07           C
ATOM     18  C   CYS A 135     -33.892   2.035   0.730  1.00 43.35           C
ATOM     19  CB  CYS A 135     -33.281   1.448  -1.628  1.00 44.78           C
ATOM     20  SG  CYS A 135     -32.962  -0.333  -1.746  1.00 48.37           S
ATOM     21  CA  CYS A 141     -32.516  -2.538  -5.168  1.00 40.78           C
ATOM     22  CB  CYS A 141     -31.888  -1.151  -5.100  1.00 42.43           C
ATOM     23  SG  CYS A 141     -32.764   0.085  -6.092  1.00 48.71           S
TER
HETATM   24 SB    SB A 402     -34.697   0.682  -3.911  1.00 56.40          SB
END
"""
  params = pdb_interpretation.master_params.extract()
  params.use_ncs_to_build_restraints = False

  processed_pdb_file = pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    params=params,
    raw_records=raw_records,
    force_symmetry=True)
  processed_pdb_file.geometry_restraints_manager()

def run(args):
  mon_lib_srv = server.server()
  ener_lib = server.ener_lib()
  exercise_long_vdw(mon_lib_srv, ener_lib)
  print(format_cpu_times())

if (__name__ == "__main__"):
  run(sys.argv[1:])
