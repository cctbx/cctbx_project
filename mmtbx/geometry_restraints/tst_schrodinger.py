from __future__ import print_function, division

import os
import random
import unittest

import numpy as np

import iotbx.cif
import iotbx.pdb
import iotbx.phil
import mmtbx.model
from cctbx.array_family import flex
from cctbx.geometry_restraints.lbfgs import lbfgs
from cctbx.geometry_restraints.flags import flags as Flags
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal
from mmtbx.geometry_restraints import external
from mmtbx.monomer_library import pdb_interpretation
from scitbx.lbfgs import termination_parameters


if True:
  random.seed(0)
  flex.set_random_seed(0)


MASTER_PARAMS_STR = pdb_interpretation.master_params_str \
                   + external.external_energy_params_str


class TestSchrodingerLigand(unittest.TestCase):

  MAE_FN = 'rq3.mae'
  PDB_FN = 'rq3.pdb'
  CIF_FN = 'rq3.cif'


  def setUp(self):
    self.params = iotbx.phil.parse(
        input_string=MASTER_PARAMS_STR, process_includes=True).extract()
    self.params.schrodinger.use_schrodinger = True
    self.params.schrodinger.maestro_file = self.MAE_FN
    self.params.schrodinger.debug = False

    # Build model manager
    pdb_inp = iotbx.pdb.input(file_name=self.PDB_FN, source_info=None)
    cif_reader = iotbx.cif.reader(file_path=self.CIF_FN, strict=True)
    cif_object = (self.CIF_FN, cif_reader.model())
    restraint_objects = [cif_object]
    self.model = mmtbx.model.manager(
        model_input=pdb_inp, pdb_interpretation_params=self.params,
        restraint_objects=restraint_objects, log=null_out())
    self.grm = self.model.get_restraints_manager().geometry
    self.sites_cart = self.model.get_sites_cart()

  def tearDown(self):
    try:
      self.grm.stop_server()
    except AttributeError:
      pass

  def update_grm(self):
    self.model.set_pdb_interpretation_params(self.params)
    self.grm = self.model.get_restraints_manager().geometry

  def test_simple(self):
    """Perform a simple energy and gradient calculation."""
    print("GETTING ENERGIES")
    energies = self.grm.energies_sites(
        sites_cart=self.sites_cart, compute_gradients=True)
    assert energies.gradients is not None

    print("GETTING ENERGIES2")
    # Check for translation invariance
    self.sites_cart += [10, 10, 10]
    energies_trans = self.grm.energies_sites(
        sites_cart=self.sites_cart, compute_gradients=True)
    assert approx_equal(energies.target, energies_trans.target)
    assert approx_equal(energies.gradients, energies_trans.gradients)

    print("GETTING ENERGIES3")
    # Next move an atom away and check if energy and gradient for that atom is
    # increasing.
    self.sites_cart[0] += np.asarray([1, 0 ,0], dtype=np.float64)
    energies_move = self.grm.energies_sites(
        sites_cart=self.sites_cart, compute_gradients=True)
    assert energies.target < energies_move.target
    gnorm = np.linalg.norm(energies.gradients[0])
    gnorm_move = np.linalg.norm(energies_move.gradients[0])
    assert gnorm < gnorm_move

  def test_select(self):
    """Test selection method of the geometry restraints manager."""
    # Do energies tests
    energies = self.grm.energies_sites(
        sites_cart=self.sites_cart, compute_gradients=True)

    # Select atom 1
    iselection = flex.size_t(1, 0)
    sites_cart_select = self.sites_cart.select(iselection)
    grm_select = self.grm.select(iselection=iselection)
    energies_select = grm_select.energies_sites(
        sites_cart=sites_cart_select, compute_gradients=True)
    assert energies.gradients[0] == energies_select.gradients[0]
    assert energies_select.gradients.size() == 1
    # The grm of PHENIX forgets about its environment after selection, so atoms
    # at the border are a combination of PHENIX and Schrodinger force field.
    # TODO

  def test_minimize(self):
    """Perform L-BFGS minimization."""
    lbfgs_termination_params = termination_parameters(max_iterations=3)
    # Get initial energies
    energies_start = self.grm.energies_sites(
        sites_cart=self.sites_cart, compute_gradients=False)
    # Perform first minimization round
    lbfgs(self.sites_cart, self.grm, 
        lbfgs_termination_params=lbfgs_termination_params)
    energies_middle = self.grm.energies_sites(
        sites_cart=self.sites_cart, compute_gradients=False)
    # Perform second minimization round
    lbfgs(self.sites_cart, self.grm, 
        lbfgs_termination_params=lbfgs_termination_params)
    energies_end = self.grm.energies_sites(
        sites_cart=self.sites_cart, compute_gradients=False)
    assert energies_end.target < energies_middle.target < energies_start.target

  def test_flags(self):
    energies = self.grm.energies_sites(self.sites_cart)
    # Remove bond energies
    self.params.schrodinger.flags.bond = False
    self.update_grm()
    energies_no_bond = self.grm.energies_sites(self.sites_cart)
    # Remove angle energies
    self.params.schrodinger.flags.angle = False
    self.update_grm()
    energies_no_angle = self.grm.energies_sites(self.sites_cart)

    assert energies_no_bond.target != energies.target
    assert energies_no_bond.target != energies_no_angle.target
    assert energies.target != energies_no_angle.target

  def _test_mix_term(self, term):

    energies_schrodinger = self.grm.energies_sites(
        self.sites_cart, compute_gradients=True)
    setattr(self.params.schrodinger.mix, term, True)
    self.update_grm()
    energies_mix_term = self.grm.energies_sites(
        self.sites_cart, compute_gradients=True)

    setattr(self.params.schrodinger.mix, term, False)
    setattr(self.params.schrodinger.flags, term, False)
    self.update_grm()
    energies_no_term = self.grm.energies_sites(
        self.sites_cart, compute_gradients=True)

    self.params.schrodinger.use_schrodinger = False
    self.update_grm()
    flags = Flags()
    setattr(flags, term, True)
    energies_phenix = self.grm.energies_sites(
        self.sites_cart, flags=flags, compute_gradients=True)

    # Test energy difference
    phenix_residual = getattr(energies_phenix, term + '_residual_sum')
    delta_mix = energies_mix_term.target - energies_no_term.target
    assert approx_equal(phenix_residual, delta_mix)
    # Test gradient difference
    delta_gradients = energies_mix_term.gradients - energies_no_term.gradients
    assert np.allclose(energies_phenix.gradients, delta_gradients)

  def test_mix_bond(self):
    self._test_mix_term('bond')

  def test_mix_angle(self):
    self._test_mix_term('angle')

  def test_mix_dihedral(self):
    self._test_mix_term('dihedral')

  def test_mix_nonbonded(self):
    self._test_mix_term('nonbonded')

  #def test_selection(self):
  #  """Test force field selection."""
  #    pass

  #def test_scaling_factor(self):
  #    pass

  #def test_forcefield(self):
  #    pass


class TestSchrodingerProtein(TestSchrodingerLigand):

  MAE_FN = '1g9v_final-prepped.maegz'
  PDB_FN = '1g9v_final-prepped.pdb'



if (__name__ == "__main__") :
  if external.schrodinger_installed:
      unittest.main()
      print("OK")
