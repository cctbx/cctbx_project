from __future__ import division
from libtbx import adopt_init_args
from cctbx import geometry_restraints
from cctbx.array_family import flex
import scitbx.restraints
from stdlib import math
import sys

class energies(scitbx.restraints.energies):

  def __init__(self, sites_cart,
               unit_cell=None,
               bond_proxies=None,
               nonbonded_proxies=None,
               nonbonded_function=None,
               angle_proxies=None,
               dihedral_proxies=None,
               reference_coordinate_proxies=None,
               reference_dihedral_proxies=None,
               ncs_dihedral_manager=None,
               chirality_proxies=None,
               planarity_proxies=None,
               parallelity_proxies=None,
               bond_similarity_proxies=None,
               generic_restraints_manager=None,
               ramachandran_manager=None,
               external_energy_function=None,
               compute_gradients=True,
               gradients=None,
               disable_asu_cache=False,
               normalization=False,
               extension_objects=[]):
    adopt_init_args(self, locals())
    scitbx.restraints.energies.__init__(self,
                                        compute_gradients=compute_gradients,
                                        gradients=gradients,
                                        gradients_size=sites_cart.size(),
                                        gradients_factory=flex.vec3_double,
                                        normalization=normalization)
    self.n_dihedral_restraints = None
    self.dihedral_restraints_residual_sum = 0
    if (nonbonded_proxies is not None): assert nonbonded_function is not None
    if (compute_gradients):
      if (self.gradients is None):
        self.gradients = flex.vec3_double(sites_cart.size(), [0,0,0])
      else:
        assert self.gradients.size() == sites_cart.size()
    if (bond_proxies is None):
      self.n_bond_proxies = None
      self.bond_residual_sum = 0
    else:
      self.n_bond_proxies = bond_proxies.n_total()
      self.bond_residual_sum = geometry_restraints.bond_residual_sum(
        sites_cart=sites_cart,
        sorted_asu_proxies=bond_proxies,
        gradient_array=self.gradients,
        disable_cache=disable_asu_cache)
      self.number_of_restraints += self.n_bond_proxies
      self.residual_sum += self.bond_residual_sum
    if (nonbonded_proxies is None):
      self.n_nonbonded_proxies = None
      self.nonbonded_residual_sum = 0
    else:
      self.n_nonbonded_proxies = nonbonded_proxies.n_total()
      self.nonbonded_residual_sum = geometry_restraints.nonbonded_residual_sum(
        sites_cart=sites_cart,
        sorted_asu_proxies=nonbonded_proxies,
        gradient_array=self.gradients,
        function=nonbonded_function,
        disable_cache=False)
      self.number_of_restraints += self.n_nonbonded_proxies
      self.residual_sum += self.nonbonded_residual_sum
    if (angle_proxies is None):
      self.n_angle_proxies = None
      self.angle_residual_sum = 0
    else:
      self.n_angle_proxies = len(angle_proxies)
      if unit_cell is None: # ignore proxy.i_seqs
        self.angle_residual_sum = geometry_restraints.angle_residual_sum(
          sites_cart=sites_cart,
          proxies=angle_proxies,
          gradient_array=self.gradients)
      else:
        self.angle_residual_sum = geometry_restraints.angle_residual_sum(
          unit_cell=unit_cell,
          sites_cart=sites_cart,
          proxies=angle_proxies,
          gradient_array=self.gradients)
      self.number_of_restraints += self.n_angle_proxies
      self.residual_sum += self.angle_residual_sum
    if (dihedral_proxies is None):
      self.n_dihedral_proxies = None
      self.dihedral_residual_sum = 0
    else:
      self.n_dihedral_proxies = len(dihedral_proxies)
      if unit_cell is None: # ignore proxy.i_seqs
        self.dihedral_residual_sum = geometry_restraints.dihedral_residual_sum(
          sites_cart=sites_cart,
          proxies=dihedral_proxies,
          gradient_array=self.gradients)
      else:
        self.dihedral_residual_sum = geometry_restraints.dihedral_residual_sum(
          unit_cell=unit_cell,
          sites_cart=sites_cart,
          proxies=dihedral_proxies,
          gradient_array=self.gradients)
      self.number_of_restraints += self.n_dihedral_proxies
      self.residual_sum += self.dihedral_residual_sum

    if reference_coordinate_proxies is None:
      self.n_reference_coordinate_proxies = None
      self.reference_coordinate_residual_sum = 0
    else:
      import boost.python
      ext = boost.python.import_ext("mmtbx_reference_coordinate_ext")
      self.n_reference_coordinate_proxies = reference_coordinate_proxies.size()
      self.reference_coordinate_residual_sum = \
          ext.reference_coordinate_residual_sum(
              sites_cart=sites_cart,
              proxies=reference_coordinate_proxies,
              gradient_array=self.gradients)
      self.number_of_restraints += self.n_reference_coordinate_proxies
      self.residual_sum += self.reference_coordinate_residual_sum

    if (reference_dihedral_proxies is None):
      self.n_reference_dihedral_proxies = None
      self.reference_dihedral_residual_sum = 0
    else:
      self.n_reference_dihedral_proxies = len(reference_dihedral_proxies)
      if unit_cell is None: # ignore proxy.i_seqs
        self.reference_dihedral_residual_sum = geometry_restraints.dihedral_residual_sum(
          sites_cart=sites_cart,
          proxies=reference_dihedral_proxies,
          gradient_array=self.gradients)
      else:
        self.reference_dihedral_residual_sum = geometry_restraints.dihedral_residual_sum(
          unit_cell=unit_cell,
          sites_cart=sites_cart,
          proxies=reference_dihedral_proxies,
          gradient_array=self.gradients)
      self.number_of_restraints += self.n_reference_dihedral_proxies
      self.residual_sum += self.reference_dihedral_residual_sum

    if (ncs_dihedral_manager is None):
      self.n_ncs_dihedral_proxies = None
      self.ncs_dihedral_residual_sum = 0
    else:
      self.n_ncs_dihedral_proxies = ncs_dihedral_manager.get_n_proxies()
      if unit_cell is None: # ignore proxy.i_seqs
        self.ncs_dihedral_residual_sum = geometry_restraints.dihedral_residual_sum(
          sites_cart=sites_cart,
          proxies=ncs_dihedral_manager.ncs_dihedral_proxies,
          gradient_array=self.gradients)
      else:
        self.ncs_dihedral_residual_sum = geometry_restraints.dihedral_residual_sum(
          unit_cell=unit_cell,
          sites_cart=sites_cart,
          proxies=ncs_dihedral_manager.ncs_dihedral_proxies,
          gradient_array=self.gradients)
      self.number_of_restraints += self.n_ncs_dihedral_proxies
      self.residual_sum += self.ncs_dihedral_residual_sum

    if (chirality_proxies is None):
      self.n_chirality_proxies = None
      self.chirality_residual_sum = 0
    else:
      self.n_chirality_proxies = len(chirality_proxies)
      self.chirality_residual_sum = geometry_restraints.chirality_residual_sum(
        sites_cart=sites_cart,
        proxies=chirality_proxies,
        gradient_array=self.gradients)
      self.number_of_restraints += self.n_chirality_proxies
      self.residual_sum += self.chirality_residual_sum
    if (planarity_proxies is None):
      self.n_planarity_proxies = None
      self.planarity_residual_sum = 0
    else:
      self.n_planarity_proxies = len(planarity_proxies)
      if unit_cell is None: # ignore proxy.i_seqs
        self.planarity_residual_sum = geometry_restraints \
          .planarity_residual_sum(
            sites_cart=sites_cart,
            proxies=planarity_proxies,
            gradient_array=self.gradients)
      else:
        self.planarity_residual_sum = geometry_restraints \
          .planarity_residual_sum(
            unit_cell=unit_cell,
            sites_cart=sites_cart,
            proxies=planarity_proxies,
            gradient_array=self.gradients)
      self.number_of_restraints += self.n_planarity_proxies
      self.residual_sum += self.planarity_residual_sum
    if parallelity_proxies is None:
      self.n_parallelity_proxies = None
      self.parallelity_residual_sum = 0
    else:
      self.n_parallelity_proxies = len(parallelity_proxies)
      if unit_cell is None: # ignore proxy.i_seqs
        self.parallelity_residual_sum = geometry_restraints \
          .parallelity_residual_sum(
            sites_cart=sites_cart,
            proxies=parallelity_proxies,
            gradient_array=self.gradients)
      else:
        self.parallelity_residual_sum = geometry_restraints \
          .parallelity_residual_sum(
            unit_cell=unit_cell,
            sites_cart=sites_cart,
            proxies=parallelity_proxies,
            gradient_array=self.gradients)
      self.number_of_restraints += self.n_parallelity_proxies
      self.residual_sum += self.parallelity_residual_sum

    if (bond_similarity_proxies is None):
      self.n_bond_similarity_proxies = None
      self.bond_similarity_residual_sum = 0
    else:
      self.n_bond_similarity_proxies = len(bond_similarity_proxies)
      if unit_cell is None: # ignore proxy.i_seqs
        self.bond_similarity_residual_sum = \
          geometry_restraints.bond_similarity_residual_sum(
            sites_cart=sites_cart,
            proxies=bond_similarity_proxies,
            gradient_array=self.gradients)
      else:
        self.bond_similarity_residual_sum = \
          geometry_restraints.bond_similarity_residual_sum(
            unit_cell=unit_cell,
            sites_cart=sites_cart,
            proxies=bond_similarity_proxies,
            gradient_array=self.gradients)
      self.number_of_restraints += self.n_bond_similarity_proxies
      self.residual_sum += self.bond_similarity_residual_sum

    if (generic_restraints_manager is None) :
      self.n_generic_proxies = 0
      self.generic_restraint_residual_sum = 0
    else :
      self.n_generic_proxies = generic_restraints_manager.get_n_proxies()
      self.generic_restraint_residual_sum = \
        generic_restraints_manager.restraints_residual_sum(
          sites_cart=sites_cart,
          gradient_array=self.gradients)
      self.number_of_restraints += self.n_generic_proxies
      self.residual_sum += self.generic_restraint_residual_sum

    if ramachandran_manager is None:
      self.n_ramachandran_proxies = 0
      self.ramachandran_restraints_residual_sum = 0
    else:
      self.n_ramachandran_proxies = self.ramachandran_manager.get_n_proxies()
      self.ramachandran_restraints_residual_sum = self.ramachandran_manager.restraints_residual_sum(
          sites_cart=sites_cart,
          gradient_array=self.gradients)
      self.number_of_restraints += self.n_ramachandran_proxies
      self.residual_sum += self.ramachandran_restraints_residual_sum
    if (external_energy_function is not None) :
      self.external_energy = external_energy_function(
        sites_cart=sites_cart,
        gradient_array=self.gradients)
      self.residual_sum += self.external_energy
    else :
      self.external_energy = 0
    for extension_obj in self.extension_objects:
      extension_obj.energies_add(energies_obj=self)
    self.finalize_target_and_gradients()

  # Not used anymore? -- used in model_statistics.py
  def get_filtered_n_bond_proxies(self):
    correct_id_proxies = [i.origin_id for i in self.bond_proxies.simple].count(0)
    return correct_id_proxies

  def get_filtered_n_angle_proxies(self):
    correct_id_proxies = [i.origin_id for i in self.angle_proxies].count(0)
    return correct_id_proxies

  def bond_deviations_z(self):
    '''
    Calculate rmsz of bond deviations

    Compute rmsz, the Root-Mean-Square of the z-scors for a set of data
    using z_i = {x_i - mu / sigma}  and rmsz = sqrt(mean(z*z))
    x_i: atcual bond length
    mu: geometry restraints mean
    sigma:  geometry restraints standard deviation
    z_i: z-score for bond i
    z: array of z_i

    The sigma and the (x_i - mu) are model constrains, geometry restraints. They function extracts
    from self, not calculated from data.

    :returns:
    b_rmsz: rmsz, root mean square of the z-scors of all bonds
    b_z_min/max: min/max abolute values of z-scors
    '''
    if(self.n_bond_proxies is not None):
      bond_deltas = self.bond_proxies.deltas(
          sites_cart=self.sites_cart, origin_id=0)
      if len(bond_deltas) >0:
        sigmas = [geometry_restraints.weight_as_sigma(x.weight) for x in self.bond_proxies.simple]
        z_scores = flex.double([(bond_delta/sigma) for bond_delta,sigma in zip(bond_deltas,sigmas)])
        b_rmsz = math.sqrt(flex.mean_default(z_scores*z_scores,0))
        b_z_max = flex.max_default(flex.abs(z_scores), 0)
        b_z_min = flex.min_default(flex.abs(z_scores), 0)
        return b_z_min, b_z_max, b_rmsz
      else:
        return 0,0,0

  def bond_deviations(self):
    if(self.n_bond_proxies is not None):
      bond_deltas = self.bond_proxies.deltas(
          sites_cart=self.sites_cart, origin_id=0)
      if len(bond_deltas) >0:
        b_sq  = bond_deltas * bond_deltas
        b_ave = math.sqrt(flex.mean_default(b_sq, 0))
        b_max = math.sqrt(flex.max_default(b_sq, 0))
        b_min = math.sqrt(flex.min_default(b_sq, 0))
        return b_min, b_max, b_ave
      else:
        return 0,0,0

  def angle_deviations_z(self):
    '''
    Calculate rmsz of angles deviations

    Compute rmsz, the Root-Mean-Square of the z-scors for a set of data
    using z_i = {x_i - mu / sigma}  and rmsz = sqrt(mean(z*z))

    Compute rmsz, the Root-Mean-Square of the z-scors for a set of data
    using z_i = {x_i - mu / sigma}  and rmsz = sqrt(mean(z*z))
    x_i: atcual bond angle
    mu: geometry restraints mean
    sigma:  geometry restraints standard deviation
    z_i: z-score for bond i
    z: array of z_i

    The sigma and the (x_i - mu) are model constrains, geometry restraints. They function extracts
    from self, not calculated from data.

    :returns:
    a_rmsz: rmsz, root mean square of the z-scors of all angles
    a_z_min/max: min/max values of z-scors
    '''
    if(self.n_angle_proxies is not None):
      angle_deltas = self.angle_proxies.proxy_select(origin_id=0).deltas(
          sites_cart=self.sites_cart)
      if len(angle_deltas) > 0:
        sigmas = [geometry_restraints.weight_as_sigma(x.weight) for x in self.angle_proxies]
        z_scores = flex.double([(angle_delta/sigma) for angle_delta,sigma in zip(angle_deltas,sigmas)])
        a_rmsz = math.sqrt(flex.mean_default(z_scores*z_scores,0))
        a_z_max = flex.max_default(flex.abs(z_scores), 0)
        a_z_min = flex.min_default(flex.abs(z_scores), 0)
        return a_z_min, a_z_max, a_rmsz
      else:
        return 0,0,0

  def angle_deviations(self):
    if(self.n_angle_proxies is not None):
      angle_deltas = self.angle_proxies.proxy_select(origin_id=0).deltas(
          sites_cart=self.sites_cart)
      if len(angle_deltas) > 0:
        a_sq  = angle_deltas * angle_deltas
        a_ave = math.sqrt(flex.mean_default(a_sq, 0))
        a_max = math.sqrt(flex.max_default(a_sq, 0))
        a_min = math.sqrt(flex.min_default(a_sq, 0))
        return a_min, a_max, a_ave
      else:
        return 0,0,0

  def nonbonded_distances(self):
    return geometry_restraints.nonbonded_deltas(
      sites_cart         = self.sites_cart,
      sorted_asu_proxies = self.nonbonded_proxies)

  def nonbonded_deviations(self):
    if(self.n_nonbonded_proxies is not None):
      nonbonded_deltas = self.nonbonded_distances()
      r_sq  = nonbonded_deltas * nonbonded_deltas
      r_ave = math.sqrt(flex.mean_default(r_sq, 0))
      r_max = math.sqrt(flex.max_default(r_sq, 0))
      r_min = math.sqrt(flex.min_default(r_sq, 0))
      return r_min, r_max, r_ave

  def dihedral_deviations(self):
    if(self.n_dihedral_proxies is not None):
      dihedral_deltas = geometry_restraints.dihedral_deltas(
        sites_cart = self.sites_cart,
        proxies    = self.dihedral_proxies)
      d_sq  = dihedral_deltas * dihedral_deltas
      d_ave = math.sqrt(flex.mean_default(d_sq, 0))
      d_max = math.sqrt(flex.max_default(d_sq, 0))
      d_min = math.sqrt(flex.min_default(d_sq, 0))
      return d_min, d_max, d_ave

  def reference_dihedral_deviations(self):
    if(self.n_reference_dihedral_proxies is not None):
      reference_dihedral_deltas = geometry_restraints.reference_dihedral_deltas(
        sites_cart = self.sites_cart,
        proxies    = self.reference_dihedral_proxies)
      d_sq  = reference_dihedral_deltas * reference_dihedral_deltas
      d_ave = math.sqrt(flex.mean_default(d_sq, 0))
      d_max = math.sqrt(flex.max_default(d_sq, 0))
      d_min = math.sqrt(flex.min_default(d_sq, 0))
      return d_min, d_max, d_ave

  def ncs_dihedral_deviations(self):
    if(self.n_ncs_dihedral_proxies is not None):
      ncs_dihedral_deltas = geometry_restraints.ncs_dihedral_deltas(
        sites_cart = self.sites_cart,
        proxies    = self.ncs_dihedral_manager.ncs_dihedral_proxies)
      d_sq  = ncs_dihedral_deltas * ncs_dihedral_deltas
      d_ave = math.sqrt(flex.mean_default(d_sq, 0))
      d_max = math.sqrt(flex.max_default(d_sq, 0))
      d_min = math.sqrt(flex.min_default(d_sq, 0))
      return d_min, d_max, d_ave

  def chirality_deviations(self):
    if(self.n_chirality_proxies is not None):
      chirality_deltas = geometry_restraints.chirality_deltas(
        sites_cart = self.sites_cart,
        proxies    = self.chirality_proxies)
      c_sq  = chirality_deltas * chirality_deltas
      c_ave = math.sqrt(flex.mean_default(c_sq, 0))
      c_max = math.sqrt(flex.max_default(c_sq, 0))
      c_min = math.sqrt(flex.min_default(c_sq, 0))
      return c_min, c_max, c_ave

  def planarity_deviations(self):
    if(self.n_planarity_proxies is not None):
      planarity_deltas = geometry_restraints.planarity_deltas_rms(
        sites_cart = self.sites_cart,
        proxies    = self.planarity_proxies)
      p_sq  = planarity_deltas * planarity_deltas
      p_ave = math.sqrt(flex.mean_default(p_sq, 0))
      p_max = math.sqrt(flex.max_default(p_sq, 0))
      p_min = math.sqrt(flex.min_default(p_sq, 0))
      return p_min, p_max, p_ave

  def parallelity_deviations(self):
    if self.n_parallelity_proxies is not None:
      parallelity_deltas = geometry_restraints.parallelity_deltas(
        sites_cart = self.sites_cart,
        proxies    = self.parallelity_proxies)
      p_sq  = parallelity_deltas * parallelity_deltas
      p_ave = math.sqrt(flex.mean_default(p_sq, 0))
      p_max = math.sqrt(flex.max_default(p_sq, 0))
      p_min = math.sqrt(flex.min_default(p_sq, 0))
      return p_min, p_max, p_ave

  def show(self, f=None, prefix=""):
    if (f is None): f = sys.stdout
    print >> f, prefix+"target: %.6g" % self.target
    if (self.n_bond_proxies is not None):
      print >> f, prefix+"  bond_residual_sum (n=%d): %.6g" % (
        self.n_bond_proxies, self.bond_residual_sum)
    if (self.n_nonbonded_proxies is not None):
      print >> f, prefix+"  nonbonded_residual_sum (n=%d): %.6g" % (
        self.n_nonbonded_proxies, self.nonbonded_residual_sum)
    if (self.n_angle_proxies is not None):
      print >> f, prefix+"  angle_residual_sum (n=%d): %.6g" % (
        self.n_angle_proxies, self.angle_residual_sum)
    if (self.n_dihedral_proxies is not None):
      print >> f, prefix+"  dihedral_residual_sum (n=%d): %.6g" % (
        self.n_dihedral_proxies, self.dihedral_residual_sum)
    if (self.n_reference_dihedral_proxies is not None):
      print >> f, prefix+"  reference_dihedral_residual_sum (n=%d): %.6g" % (
        self.n_reference_dihedral_proxies, self.reference_dihedral_residual_sum)

    if (self.n_reference_coordinate_proxies is not None):
      print >> f, prefix+"  reference_coordinate_residual_sum (n=%d): %.6g" % (
        self.n_reference_coordinate_proxies, self.reference_coordinate_residual_sum)

    if (self.n_ncs_dihedral_proxies is not None):
      print >> f, prefix+"  ncs_dihedral_residual_sum (n=%d): %.6g" % (
        self.n_ncs_dihedral_proxies, self.ncs_dihedral_residual_sum)
    if (self.n_chirality_proxies is not None):
      print >> f, prefix+"  chirality_residual_sum (n=%d): %.6g" % (
        self.n_chirality_proxies, self.chirality_residual_sum)
    if (self.n_planarity_proxies is not None):
      print >> f, prefix+"  planarity_residual_sum (n=%d): %.6g" % (
        self.n_planarity_proxies, self.planarity_residual_sum)
    if (self.n_parallelity_proxies is not None):
      print >> f, prefix+"  parallelity_residual_sum (n=%d): %.6g" % (
        self.n_parallelity_proxies, self.parallelity_residual_sum)
    if (self.n_ramachandran_proxies is not None):
      print >> f, prefix+"  ramachandran_residual_sum (n=%d): %.6g" % (
        self.n_ramachandran_proxies, self.ramachandran_restraints_residual_sum)
    if (self.n_bond_similarity_proxies is not None):
      print >> f, prefix+"  bond_similarity_residual_sum (n=%d): %.6g" % (
        self.n_bond_similarity_proxies, self.bond_similarity_residual_sum)
    for extension_obj in self.extension_objects:
      extension_obj.energies_show(energies_obj=self, f=f, prefix=prefix)
    if (self.gradients is not None):
      print >> f, prefix+"  norm of gradients: %.6g" % self.gradients.norm()
