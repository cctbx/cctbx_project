from libtbx import test_utils
import libtbx.load_env

def run():
  tst_list = (
  "$D/mmtbx/refinement/tst_rigid_body.py",
  "$D/mmtbx/tst_model.py",
  "$D/mmtbx/tst_fmodel.py",
  "$D/mmtbx/tst_fmodel_fd.py",
  "$D/ncs/tst_restraints.py",
  "$D/mmtbx/ncs/tst_restraints.py",
  "$D/mmtbx/regression/tst_adp_restraints.py",
  "$D/mmtbx/scaling/tst_scaling.py",
  "$D/mmtbx/monomer_library/tst_cif_types.py",
  "$D/mmtbx/monomer_library/tst_motif.py",
  "$D/mmtbx/monomer_library/tst_selection.py",
  "$D/mmtbx/monomer_library/tst_tyr_from_gly_and_bnz.py",
  "$D/mmtbx/hydrogens/build_hydrogens.py",
  "$D/mmtbx/max_lik/tst_maxlik.py",
  "$D/masks/tst_masks.py",
  "$D/max_lik/tst_max_lik.py",
  "$D/mmtbx/dynamics/tst_cartesian_dynamics.py",
  "$D/mmtbx/tls/tst_tls.py",
  "$D/mmtbx/tls/tst_tls_refinement_fft.py",
  "$D/mmtbx/examples/f_model_manager.py",
  "$D/mmtbx/bulk_solvent/tst_bulk_solvent_and_scaling_ls.py",
  "$D/mmtbx/bulk_solvent/tst_bulk_solvent_and_scaling_ml.py",
  )

  build_dir = libtbx.env.under_build("mmtbx")
  dist_dir = libtbx.env.dist_path("mmtbx")

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
