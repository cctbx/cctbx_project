from __future__ import absolute_import, division, print_function
import mmtbx.monomer_library.pdb_interpretation
import cctbx.geometry_restraints.flags
from cctbx.array_family import flex
import scitbx.rigid_body
import scitbx.lbfgs
import time
from mmtbx import monomer_library
import mmtbx.refinement.rigid_body
import mmtbx.refinement.real_space.rigid_body
import scitbx.lbfgs
import mmtbx.geometry_restraints
import cctbx.geometry_restraints
import cctbx.geometry_restraints.flags
from cctbx import adp_restraints # import dependency
import mmtbx.utils
import scitbx.rigid_body
import iotbx


pdb_str = """\
CRYST1   14.914   14.930   14.600  90.00  90.00  90.00 P 1
ATOM     17  N   TYR A  20       5.877   9.930   9.600  1.00 23.81           N
ATOM     18  CA  TYR A  20       6.364   9.061   8.536  1.00 35.69           C
ATOM     19  C   TYR A  20       5.278   8.091   8.080  1.00 28.94           C
ATOM     20  O   TYR A  20       5.000   7.100   8.756  1.00 25.26           O
ATOM     21  CB  TYR A  20       7.603   8.290   8.999  1.00 33.37           C
ATOM     22  CG  TYR A  20       8.217   7.412   7.931  1.00 31.48           C
ATOM     23  CD1 TYR A  20       9.080   7.943   6.981  1.00 36.91           C
ATOM     24  CD2 TYR A  20       7.940   6.052   7.877  1.00 21.36           C
ATOM     25  CE1 TYR A  20       9.646   7.145   6.004  1.00 36.56           C
ATOM     26  CE2 TYR A  20       8.501   5.247   6.904  1.00 29.56           C
ATOM     27  CZ  TYR A  20       9.353   5.798   5.970  1.00 35.08           C
ATOM     28  OH  TYR A  20       9.914   5.000   5.000  1.00 38.73           O
TER
END
"""

def exercise(d_min = 2.0, resolution_factor = 0.1):
  params = iotbx.phil.parse(
    monomer_library.pdb_interpretation.grand_master_phil_str,
    process_includes=True).extract()
  params.pdb_interpretation.reference_coordinate_restraints.enabled = True
  params.pdb_interpretation.reference_coordinate_restraints.sigma = 1.0
  params.pdb_interpretation.reference_coordinate_restraints.selection = \
      "name C or name N"
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv              = monomer_library.server.server(),
    ener_lib                 = monomer_library.server.ener_lib(),
    params                   = params.pdb_interpretation,
    raw_records              = flex.std_string(pdb_str.splitlines()),
    strict_conflict_handling = True,
    force_symmetry           = True,
    log                      = None)
  pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy.deep_copy()
  reference_selection = flex.size_t()
  reference_sites = flex.vec3_double()
  for a in pdb_hierarchy.only_residue().atoms():
    if(a.name.strip().upper() in ["N", "C"]):
      reference_selection.append(a.i_seq)
      reference_sites.append(a.xyz)
  pdb_hierarchy.write_pdb_file(file_name = "answer.pdb")
  xrs_answer = processed_pdb_file.xray_structure()
  f_calc = xrs_answer.structure_factors(d_min = d_min).f_calc()
  fft_map = f_calc.fft_map(resolution_factor=resolution_factor)
  fft_map.apply_sigma_scaling()
  target_map = fft_map.real_map_unpadded()
  mtz_dataset = f_calc.as_mtz_dataset(column_root_label = "FCmap")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "answer.mtz")
  #
  rot_obj = scitbx.rigid_body.rb_mat_zyz(phi = 20, psi = 30, the = 40)
  xrs_poor = xrs_answer.deep_copy_scatterers()
  xrs_poor.apply_rigid_body_shift(
    rot = rot_obj.rot_mat().as_mat3(), trans = [1,2,3])
  pdb_hierarchy.adopt_xray_structure(xrs_poor)
  pdb_hierarchy.write_pdb_file(file_name = "poor.pdb")
  #
  restraints_manager = processed_pdb_file.geometry_restraints_manager()
  residue_poor = pdb_hierarchy.only_residue()
  residue_poor.atoms().set_xyz(xrs_poor.sites_cart())
  lbfgs_termination_params=scitbx.lbfgs.termination_parameters(
    max_iterations = 250)
  flags = cctbx.geometry_restraints.flags.flags(reference_coordinate=True, default=True)
  states_collector = mmtbx.utils.states(
    pdb_hierarchy  = pdb_hierarchy)
  states_collector.add(sites_cart = residue_poor.atoms().extract_xyz())
  for i in [1,2]:
    minimized = mmtbx.refinement.real_space.rigid_body.refine(
      residue                     = residue_poor,
      density_map                 = target_map,
      geometry_restraints_manager = restraints_manager,
      real_space_target_weight    = 1,
      real_space_gradients_delta  = d_min*resolution_factor,
      lbfgs_termination_params    = lbfgs_termination_params,
      unit_cell                   = xrs_poor.unit_cell(),
      cctbx_geometry_restraints_flags = flags,
      states_collector            = states_collector)
    xrs_poor = xrs_poor.replace_sites_cart(minimized.sites_cart_residue)
    pdb_hierarchy.adopt_xray_structure(xrs_poor)
    pdb_hierarchy.write_pdb_file(file_name = "refined.pdb")
  #
  dist = xrs_answer.max_distance(other = xrs_poor)
  assert dist < 0.15, dist
  dist = xrs_answer.mean_distance(other = xrs_poor)
  assert dist < 0.08, dist
  states_collector.write(file_name = "all.pdb")

if(__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print("Time: %6.2f" % (time.time()-t0))
