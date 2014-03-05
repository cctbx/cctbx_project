from __future__ import division
from iotbx.pdb.multimer_reconstruction import multimer
import iotbx.pdb
import mmtbx.f_model
import mmtbx.refinement.minimization_ncs_constraints
from libtbx.test_utils import approx_equal
import os
from libtbx import adopt_init_args
from scitbx.array_family import flex
import mmtbx.utils
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.monomer_library.server
import mmtbx.geometry_restraints

ncs_1_copy="""\
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  0.496590 -0.643597  0.582393        0.00000
MTRIX2   2  0.867925  0.376088 -0.324443        0.00000
MTRIX3   2 -0.010221  0.666588  0.745356        0.00000
MTRIX1   3 -0.317946 -0.173437  0.932111        0.00000
MTRIX2   3  0.760735 -0.633422  0.141629        0.00000
MTRIX3   3  0.565855  0.754120  0.333333        0.00000
ATOM      1  N   THR A   1       8.111  11.079  10.645  1.00 20.00           N
ATOM      2  CA  THR A   1       8.000   9.721  10.125  1.00 20.00           C
ATOM      3  C   THR A   1       8.075   8.693  11.249  1.00 20.00           C
ATOM      4  O   THR A   1       8.890   8.817  12.163  1.00 20.00           O
ATOM      5  CB  THR A   1       9.101   9.420   9.092  1.00 20.00           C
ATOM      6  OG1 THR A   1       9.001  10.342   8.000  1.00 20.00           O
ATOM      7  CG2 THR A   1       8.964   8.000   8.565  1.00 20.00           C
TER
"""

class ncs_minimization_test(object):

  def __init__(self,
               n_macro_cycle,
               sites,
               u_iso,
               finite_grad_differences_test,
               use_restraints=False):
    """ create temp test files and data for tests """
    adopt_init_args(self, locals())
    # 1 NCS copy: starting template to generate whole asu; place into P1 box
    pdb_inp = iotbx.pdb.input(source_info=None, lines=ncs_1_copy)
    mtrix_object = pdb_inp.process_mtrix_records()
    ph = pdb_inp.construct_hierarchy()
    xrs = pdb_inp.xray_structure_simple()
    xrs_one_ncs = xrs.orthorhombic_unit_cell_around_centered_scatterers(
      buffer_size=8)
    ph.adopt_xray_structure(xrs_one_ncs)
    of = open("one_ncs_in_asu.pdb", "w")
    print >> of, mtrix_object.format_MTRIX_pdb_string()
    print >> of, ph.as_pdb_string(crystal_symmetry=xrs_one_ncs.crystal_symmetry())
    of.close()
    # 1 NCS copy -> full asu (expand NCS). This is the answer-structure
    m = multimer(file_name="one_ncs_in_asu.pdb",
                 round_coordinates=False,
                 reconstruction_type='cau',error_handle=True,eps=1e-2)
    assert m.number_of_transforms == 2, m.number_of_transforms
    xrs_asu = m.assembled_multimer.extract_xray_structure(
      crystal_symmetry = xrs_one_ncs.crystal_symmetry())
    m.write("full_asu.pdb")
    # force ASU none-rounded coordinates into xray structure
    xrs_asu.set_sites_cart(m.sites_cart())
    assert xrs_asu.crystal_symmetry().is_similar_symmetry(
      xrs_one_ncs.crystal_symmetry())
    # Generate Fobs from answer structure
    f_obs = abs(xrs_asu.structure_factors(d_min=2, algorithm="direct").f_calc())
    r_free_flags = f_obs.generate_r_free_flags()
    mtz_dataset = f_obs.as_mtz_dataset(column_root_label="F-obs")
    mtz_dataset.add_miller_array(
      miller_array=r_free_flags,
      column_root_label="R-free-flags")
    mtz_object = mtz_dataset.mtz_object()
    mtz_object.write(file_name = "data.mtz")
    # Shake structure - subject to refinement input
    xrs_shaken = xrs_one_ncs.deep_copy_scatterers()
    if sites: xrs_shaken.shake_sites_in_place(mean_distance=0.5)
    if(self.u_iso):
      u_random = flex.random_double(xrs_shaken.scatterers().size())
      xrs_shaken = xrs_shaken.set_u_iso(values=u_random)
    ph.adopt_xray_structure(xrs_shaken)
    of = open("one_ncs_in_asu_shaken.pdb", "w")
    print >> of, mtrix_object.format_MTRIX_pdb_string()
    print >> of, ph.as_pdb_string(crystal_symmetry=xrs.crystal_symmetry())
    of.close()
    self.f_obs = f_obs
    self.r_free_flags = r_free_flags
    self.xrs_one_ncs = xrs_one_ncs

  def get_restraints_manager(self,pdb_string=None):
    if not self.use_restraints:
      restraints_manager = None
    else:
      assert pdb_string
      processed_pdb_files_srv = mmtbx.utils.process_pdb_file_srv(
        crystal_symmetry=self.xrs_one_ncs.crystal_symmetry())
      processed_pdb_file, pdb_inp = processed_pdb_files_srv.\
        process_pdb_files(raw_records=pdb_string.splitlines())
      geometry = processed_pdb_file.geometry_restraints_manager()
      restraints_manager = mmtbx.restraints.manager(
        geometry      = geometry,
        normalization = True)
      # restraints_manager.crystal_symmetry = self.xrs_one_ncs.crystal_symmetry()
    return restraints_manager

  def run_test(self):
    ### Refinement
    params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
    params.algorithm = "direct"
    # Get the xray_structure of the shaken ASU
    m_shaken = multimer(
      file_name="one_ncs_in_asu_shaken.pdb",
      round_coordinates=False,
      reconstruction_type='cau',error_handle=True,eps=1e-2)
    xrs_shaken_asu = m_shaken.assembled_multimer.extract_xray_structure(
      crystal_symmetry=self.xrs_one_ncs.crystal_symmetry())
    # force non-rounded coordinates into xray structure
    xrs_shaken_asu.set_sites_cart(m_shaken.sites_cart())
    # Save the shaken ASU for inspection
    m_shaken.write(
      pdb_output_file_name='asu_shaken.pdb',
      crystal_symmetry=self.xrs_one_ncs.crystal_symmetry())
    # Create a boolean selection string for selecting chains in NCS
    selection_str = 'chain A'
    ncs_selection = m_shaken.assembled_multimer.\
      atom_selection_cache().selection(selection_str)
    assert ncs_selection.count(True) > 0
    fmodel = mmtbx.f_model.manager(
      f_obs                        = self.f_obs,
      r_free_flags                 = self.r_free_flags,
      xray_structure               = xrs_shaken_asu,
      sf_and_grads_accuracy_params = params,
      target_name                  = "ls_wunit_k1")
    r_start = fmodel.r_work()
    assert r_start > 0.15
    print "start r_factor: %6.4f" % r_start
    pdb_str = m_shaken.assembled_multimer.as_pdb_string(
      crystal_symmetry=self.xrs_one_ncs.crystal_symmetry())
    restraints_manager = self.get_restraints_manager(pdb_string=pdb_str)
    for macro_cycle in xrange(self.n_macro_cycle):
      minimized = mmtbx.refinement.minimization_ncs_constraints.lbfgs(
        fmodel                       = fmodel,
        restraints_manager           = restraints_manager,
        rotation_matrices            = m_shaken.rotation_matrices,
        translation_vectors          = m_shaken.translation_vectors,
        ncs_atom_selection           = ncs_selection,
        finite_grad_differences_test = self.finite_grad_differences_test,
        refine_sites                 = self.sites,
        refine_u_iso                 = self.u_iso)
      refine_type = 'adp'*self.u_iso + 'sites'*self.sites
      print "  macro_cycle %3d (%s)   r_factor: %6.4f"%(macro_cycle,
        refine_type, fmodel.r_work())
      assert approx_equal(fmodel.r_work(), minimized.fmodel.r_work())
    # check results
    if(self.u_iso):
      assert approx_equal(fmodel.r_work(), 0, 1.e-5)
    elif(self.sites):
      assert approx_equal(fmodel.r_work(), 0, 1.e-5)
    else: assert 0
    # output refined model
    xrs_refined = fmodel.xray_structure
    m_shaken.assembled_multimer.adopt_xray_structure(fmodel.xray_structure)
    output_file_name = "refined_u_iso%s_sites%s.pdb"%(str(self.u_iso),
      str(self.sites))
    m_shaken.write(output_file_name)
    # check final model
    pdb_inp_answer = iotbx.pdb.input(source_info=None, lines=ncs_1_copy)
    pdb_inp_refined = iotbx.pdb.input(file_name=output_file_name)
    xrs1 = pdb_inp_answer.xray_structure_simple()
    xrs2 = pdb_inp_refined.xray_structure_simple().select(ncs_selection)
    mmtbx.utils.assert_xray_structures_equal(
      x1 = xrs1,
      x2 = xrs2,
      sites = False)
    delta = flex.vec3_double([xrs1.center_of_mass()]*xrs2.scatterers().size())-\
            flex.vec3_double([xrs2.center_of_mass()]*xrs2.scatterers().size())
    xrs2.set_sites_cart(sites_cart = xrs2.sites_cart()+delta)
    mmtbx.utils.assert_xray_structures_equal(
      x1 = xrs1,
      x2 = xrs2)

  def clean_up_temp_test_files(self):
    """delete temporary test files """
    files_to_delete = ['one_ncs_in_asu.pdb','full_asu.pdb',
                       'one_ncs_in_asu_shaken.pdb',
                       'asu_shaken.pdb','data.mtz']
    for fn in files_to_delete:
      if os.path.isfile(fn): os.remove(fn)

if __name__ == "__main__":
  for sites, u_iso, n_macro_cycle in [(True, False, 100), (False, True, 50)]:
    t = ncs_minimization_test(
      n_macro_cycle=n_macro_cycle,
      sites=sites,
      u_iso=u_iso,
      finite_grad_differences_test = False,
      use_restraints = False)
    t.run_test()
    t.clean_up_temp_test_files()
