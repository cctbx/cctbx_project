from __future__ import absolute_import, division, print_function
import time

import mmtbx.refinement.minimization_ncs_constraints
#import mmtbx.refinement.real_space.weight
from scitbx.array_family import flex
import mmtbx.ncs.ncs_utils as nu
import iotbx.ncs as ncs
import mmtbx.utils
import iotbx.pdb
import mmtbx.ncs
import math
from six.moves import range


pdb_str_answer = """\
CRYST1   26.628   30.419   28.493  90.00  90.00  90.00 P 1
ATOM      1  N   THR A   1      15.638  20.419  12.645  1.00 10.00           N
ATOM      2  CA  THR A   1      15.527  19.061  12.125  1.00 10.00           C
ATOM      3  C   THR A   1      15.602  18.033  13.249  1.00 10.00           C
ATOM      4  O   THR A   1      16.417  18.157  14.163  1.00 10.00           O
ATOM      5  CB  THR A   1      16.628  18.760  11.092  1.00 10.00           C
ATOM      6  OG1 THR A   1      16.528  19.682  10.000  1.00 10.00           O
ATOM      7  CG2 THR A   1      16.491  17.340  10.565  1.00 10.00           C
TER
ATOM      1  N   THR B   1      10.624  17.093  17.237  1.00 10.00           N
ATOM      2  CA  THR B   1      11.140  16.654  15.945  1.00 10.00           C
ATOM      3  C   THR B   1      12.494  15.968  16.097  1.00 10.00           C
ATOM      4  O   THR B   1      13.351  16.426  16.852  1.00 10.00           O
ATOM      5  CB  THR B   1      11.279  17.832  14.963  1.00 10.00           C
ATOM      6  OG1 THR B   1      10.000  18.446  14.765  1.00 10.00           O
ATOM      7  CG2 THR B   1      11.818  17.350  13.625  1.00 10.00           C
TER
ATOM      1  N   THR C   1      12.949  10.000  18.493  1.00 10.00           N
ATOM      2  CA  THR C   1      12.735  10.702  17.233  1.00 10.00           C
ATOM      3  C   THR C   1      13.937  11.570  16.875  1.00 10.00           C
ATOM      4  O   THR C   1      14.508  12.241  17.734  1.00 10.00           O
ATOM      5  CB  THR C   1      11.474  11.584  17.284  1.00 10.00           C
ATOM      6  OG1 THR C   1      10.328  10.770  17.559  1.00 10.00           O
ATOM      7  CG2 THR C   1      11.273  12.305  15.960  1.00 10.00           C
TER
END
"""

pdb_str_poor = """\
CRYST1   26.628   30.419   28.493  90.00  90.00  90.00 P 1
ATOM      1  N   THR A   1      15.886  19.796  13.070  1.00 10.00           N
ATOM      2  CA  THR A   1      15.489  18.833  12.050  1.00 10.00           C
ATOM      3  C   THR A   1      15.086  17.502  12.676  1.00 10.00           C
ATOM      4  O   THR A   1      15.739  17.017  13.600  1.00 10.00           O
ATOM      5  CB  THR A   1      16.619  18.590  11.033  1.00 10.00           C
ATOM      6  OG1 THR A   1      16.963  19.824  10.392  1.00 10.00           O
ATOM      7  CG2 THR A   1      16.182  17.583   9.980  1.00 10.00           C
TER       8      THR A   1
ATOM      1  N   THR B   1      10.028  17.193  16.617  1.00 10.00           N
ATOM      2  CA  THR B   1      11.046  16.727  15.681  1.00 10.00           C
ATOM      3  C   THR B   1      12.336  16.360  16.407  1.00 10.00           C
ATOM      4  O   THR B   1      12.772  17.068  17.313  1.00 10.00           O
ATOM      5  CB  THR B   1      11.356  17.789  14.609  1.00 10.00           C
ATOM      6  OG1 THR B   1      10.163  18.098  13.879  1.00 10.00           O
ATOM      7  CG2 THR B   1      12.418  17.281  13.646  1.00 10.00           C
TER      16      THR B   1
ATOM      1  N   THR C   1      12.121   9.329  18.086  1.00 10.00           N
ATOM      2  CA  THR C   1      12.245  10.284  16.991  1.00 10.00           C
ATOM      3  C   THR C   1      13.707  10.622  16.718  1.00 10.00           C
ATOM      4  O   THR C   1      14.493  10.814  17.645  1.00 10.00           O
ATOM      5  CB  THR C   1      11.474  11.584  17.284  1.00 10.00           C
ATOM      6  OG1 THR C   1      10.087  11.287  17.482  1.00 10.00           O
ATOM      7  CG2 THR C   1      11.619  12.563  16.129  1.00 10.00           C
TER      24      THR C   1
END
"""

pdb_str_poor2 = """\
CRYST1   26.628   30.419   28.493  90.00  90.00  90.00 P 1
ATOM      1  N   THR A   1      15.886  19.796  13.070  1.00 10.00           N
ATOM      2  CA  THR A   1      15.489  18.833  12.050  1.00 10.00           C
ATOM      3  C   THR A   1      15.086  17.502  12.676  1.00 10.00           C
ATOM      4  O   THR A   1      15.739  17.017  13.600  1.00 10.00           O
ATOM      5  CB  THR A   1      16.619  18.590  11.033  1.00 10.00           C
ATOM      6  OG1 THR A   1      16.963  19.824  10.392  1.00 10.00           O
ATOM      7  CG2 THR A   1      16.182  17.583   9.980  1.00 10.00           C
TER       8      THR A   1
"""

def run(prefix="tst", d_min=1.0):
  """
  NCS constraints: xyz, adp, and operators.
  """
  pdb_file_name_answer = "%s_answer.pdb"%prefix
  of=open(pdb_file_name_answer, "w")
  print(pdb_str_answer, file=of)
  of.close()
  #
  pdb_file_name_poor = "%s_poor.pdb"%prefix
  of=open(pdb_file_name_poor, "w")
  print(pdb_str_poor, file=of)
  of.close()
  #
  pdb_inp_answer = iotbx.pdb.input(file_name=pdb_file_name_answer)
  ph_answer = pdb_inp_answer.construct_hierarchy()
  ph_answer.atoms().reset_i_seq()
  xrs_answer = pdb_inp_answer.xray_structure_simple()
  #
  pdb_inp_poor = iotbx.pdb.input(file_name=pdb_file_name_poor)
  ph_poor = pdb_inp_poor.construct_hierarchy()
  ph_poor_obj = iotbx.pdb.input(source_info=None, lines=pdb_str_poor2)
  ph_poor.atoms().reset_i_seq()
  xrs_poor = pdb_inp_poor.xray_structure_simple()
  #
  ppf = mmtbx.utils.process_pdb_file_srv().process_pdb_files(
    raw_records=pdb_str_poor.splitlines())[0]
  mmtbx.utils.assert_xray_structures_equal(
    x1=ppf.xray_structure(show_summary = False),
    x2=xrs_poor)
  restraints_manager = mmtbx.restraints.manager(
    geometry      = ppf.geometry_restraints_manager(show_energies = False),
    normalization = True)
  #
  fc = xrs_answer.structure_factors(d_min=d_min).f_calc()
  fft_map = fc.fft_map(resolution_factor = 0.25)
  fft_map.apply_sigma_scaling()
  map_data = fft_map.real_map_unpadded()
  #
  transforms_obj = ncs.input(hierarchy=ph_answer)
  x = transforms_obj.get_ncs_restraints_group_list().concatenate_rot_tran()
  x = nu.shake_transformations(
    x = x,
    shake_angles_sigma=1*math.pi/180,
    shake_translation_sigma=0.1)
  nrgl = transforms_obj.get_ncs_restraints_group_list()
  nrgl.update_rot_tran(x=x)
  ncs_restraints_group_list = transforms_obj.get_ncs_restraints_group_list()
  refine_selection = flex.size_t(range(transforms_obj.truncated_hierarchy.atoms_size()))
  for i in range(5):
    data_weight = 1
    tfg_obj = mmtbx.refinement.minimization_ncs_constraints.\
      target_function_and_grads_real_space(
        map_data                   = map_data,
        xray_structure             = xrs_poor,
        ncs_restraints_group_list  = ncs_restraints_group_list,
        refine_selection           = refine_selection,
        real_space_gradients_delta = d_min/4,
        restraints_manager         = restraints_manager,
        data_weight                = data_weight,
        refine_sites               = True)
    minimized = mmtbx.refinement.minimization_ncs_constraints.lbfgs(
      target_and_grads_object      = tfg_obj,
      xray_structure               = xrs_poor,
      ncs_restraints_group_list    = ncs_restraints_group_list,
      refine_selection             = refine_selection,
      finite_grad_differences_test = False,
      max_iterations               = 60,
      refine_sites                 = True)
    xrs_poor = tfg_obj.xray_structure
    ph_poor.adopt_xray_structure(tfg_obj.xray_structure)
  ph_poor.write_pdb_file(file_name="refined.pdb")


if (__name__ == "__main__"):
  t0=time.time()
  run()
  print("Time: %6.4f"%(time.time()-t0))
  print("OK")
