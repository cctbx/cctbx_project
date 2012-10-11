from __future__ import division
from cctbx.array_family import flex
import scitbx.rigid_body
import scitbx.lbfgs
from mmtbx.utils import rotatable_bonds
from libtbx import adopt_init_args
import iotbx.pdb
import mmtbx.refinement.real_space
import cctbx.geometry_restraints.flags

def get_rotamer_iterator(mon_lib_srv, residue):
  rotamer_iterator = mon_lib_srv.rotamer_iterator(
    fine_sampling = True,
    comp_id=residue.resname,
    atom_names=residue.atoms().extract_name(),
    sites_cart=residue.atoms().extract_xyz())
  if (rotamer_iterator is None) :
    return None
  if (rotamer_iterator.problem_message is not None):
    return None
  if (rotamer_iterator.rotamer_info is None):
    return None
  return rotamer_iterator

class cluster(object):
  def __init__(self,
               axis,
               atoms_to_rotate,
               vector=None,
               selection=None):
    adopt_init_args(self, locals())

class manager(object):
  def __init__(self,
               target_map,
               residue,
               special_position_settings,
               mon_lib_srv,
               rotamer_manager,
               sites_cart_all=None,
               use_clash_filter = False,
               use_slope=True,
               debug=False,
               use_torsion_search=True,
               use_rotamer_iterator=True):
    adopt_init_args(self, locals())
    get_class = iotbx.pdb.common_residue_names_get_class
    assert get_class(residue.resname) == "common_amino_acid", residue.resname
    self.vector_selections = None
    self.clusters = None
    self.axes_and_atoms_aa_specific = None
    if(use_torsion_search): self.clusters = self.define_clusters()
    if(use_slope): self.get_vector()
    self.fit(debug=debug)

  def define_clusters(self):
    # XXX could be a plase to exclude H
    result = []
    backrub_axis  = []
    backrub_atoms_to_rotate = []
    backrub_atoms_to_evaluate = []
    counter = 0 # XXX DOES THIS RELY ON ORDER
    for atom in self.residue.atoms():
      if(atom.name.strip().upper() in ["N", "C"]):
        backrub_axis.append(counter)
      else:
        backrub_atoms_to_rotate.append(counter)
      if(atom.name.strip().upper() in ["CA", "O", "CB"]):
        backrub_atoms_to_evaluate.append(counter)
      counter += 1
    result.append(cluster(
      axis=backrub_axis,
      atoms_to_rotate=backrub_atoms_to_rotate,
      selection=backrub_atoms_to_evaluate))
    self.axes_and_atoms_aa_specific = \
      rotatable_bonds.axes_and_atoms_aa_specific(residue = self.residue,
        mon_lib_srv = self.mon_lib_srv)
    if(self.axes_and_atoms_aa_specific is not None):
      for i_aa, aa in enumerate(self.axes_and_atoms_aa_specific):
        if(i_aa == len(self.axes_and_atoms_aa_specific)-1):
          result.append(cluster(
            axis=aa[0],
            atoms_to_rotate=aa[1],
            selection=flex.size_t(aa[1]),
            vector=None)) # XXX
        else:
          result.append(cluster(
            axis=aa[0],
            atoms_to_rotate=aa[1],
            selection=flex.size_t([aa[1][0]]),
            vector=None)) # XXX
    return result

  def get_vector(self):
    self.vector_selections = []
    if(self.axes_and_atoms_aa_specific is not None and self.clusters is not None):
      for i_aa, aa in enumerate(self.axes_and_atoms_aa_specific):
        for aa_ in aa[0]:
          if(not aa_ in self.vector_selections):
            self.vector_selections.append(aa_)
      self.vector_selections.append(
        self.clusters[len(self.clusters)-1].atoms_to_rotate)
      #if(self.vector_selections.size()>1):
      #  self.vector_selections = self.vector_selections[1:]
      for cl in self.clusters:
        cl.vector = self.vector_selections

  def get_scorer(self, vector=None, use_binary=False, use_clash_filter=False):
    return mmtbx.refinement.real_space.score(
      target_map                = self.target_map,
      sites_cart_all            = self.sites_cart_all,
      special_position_settings = self.special_position_settings,
      residue                   = self.residue,
      vector                    = vector,
      use_binary                = use_binary,
      use_clash_filter          = use_clash_filter,
      rotamer_eval              = self.rotamer_manager)

  def fit(self, debug=False):
    if(self.use_rotamer_iterator):
      rotamer_iterator = get_rotamer_iterator(
        mon_lib_srv = self.mon_lib_srv,
        residue     = self.residue)
      if(rotamer_iterator is not None):
        score_rotamers = self.get_scorer(
          vector     = self.vector_selections,
          use_clash_filter = self.use_clash_filter,
          use_binary = True)
        cntr=0
        for rotamer, rotamer_sites_cart in rotamer_iterator:
          if(debug): print rotamer.id,"-"*50
          if(cntr==0):
            score_rotamers.reset_with(sites_cart = rotamer_sites_cart.deep_copy())
            if(debug): print score_rotamers.target, 0, rotamer.id
          else:
            score_rotamers.update(sites_cart = rotamer_sites_cart.deep_copy())
            if(debug): print score_rotamers.target, 0, rotamer.id
          cntr+=1
          if(self.use_torsion_search):
            score_rotamer = self.get_scorer()
            score_rotamers.update(sites_cart = rotamer_sites_cart.deep_copy())
            mmtbx.refinement.real_space.torsion_search(
              clusters   = self.clusters,
              sites_cart = rotamer_sites_cart.deep_copy(),
              scorer     = score_rotamer)
            if(debug): print score_rotamers.target, 1
            score_rotamers.update(sites_cart = score_rotamer.sites_cart.deep_copy(), tmp=rotamer.id)
            if(debug): print score_rotamers.target, 2
          else:
            if(debug): print score_rotamers.target, 3
            score_rotamers.update(sites_cart = rotamer_sites_cart.deep_copy())
            if(debug): print score_rotamers.target, 4
          if(debug): print rotamer.id, score_rotamers.target, score_rotamers.tmp, score_rotamer.target
        if(debug): print "final:", score_rotamers.target, score_rotamers.tmp
        if(score_rotamers.sites_cart is not None):
          self.residue.atoms().set_xyz(new_xyz = score_rotamers.sites_cart)
          if(self.sites_cart_all is not None):
            self.sites_cart_all = self.sites_cart_all.set_selected(
              self.residue.atoms().extract_i_seq(), score_rotamers.sites_cart)
          #
          if 0:
            score_rotamer = self.get_scorer()
            score_rotamer.update(sites_cart = score_rotamers.sites_cart.deep_copy())
            mmtbx.refinement.real_space.torsion_search_nested(
                clusters   = self.clusters,
                sites_cart = score_rotamers.sites_cart,
                scorer     = score_rotamer)
            self.residue.atoms().set_xyz(new_xyz = score_rotamer.sites_cart)
            print "final-final:", score_rotamer.target
          #
    elif(self.use_torsion_search):
      score_residue = self.get_scorer()
      mmtbx.refinement.real_space.torsion_search(
        clusters   = self.clusters,
        sites_cart = self.residue.atoms().extract_xyz(),
        scorer = score_residue,
        start = 0,
        stop  = 360,
        step  = 1)
      self.residue.atoms().set_xyz(new_xyz=score_residue.sites_cart)

  def rigid_body_refine(self, max_iterations = 250):
    import mmtbx.geometry_restraints
    import cctbx.geometry_restraints.manager
    import mmtbx.refinement.real_space.rigid_body
    reference_sites = flex.vec3_double()
    reference_selection = flex.size_t()
    cntr = 0
    for atom in self.residue.atoms():
      if(atom.name.strip().upper() in ["N", "C"]):
        reference_sites.append(atom.xyz)
        reference_selection.append(cntr)
        cntr += 1
    generic_restraints_manager = mmtbx.geometry_restraints.manager()
    restraints_manager = cctbx.geometry_restraints.manager.manager(
      generic_restraints_manager = generic_restraints_manager)
    restraints_manager.generic_restraints_manager.reference_manager.\
      add_coordinate_restraints(
        sites_cart = reference_sites,
        selection  = reference_selection,
        sigma      = 0.5)
    flags = cctbx.geometry_restraints.flags.flags(generic_restraints=True)
    lbfgs_termination_params=scitbx.lbfgs.termination_parameters(
      max_iterations = max_iterations)
    minimized = mmtbx.refinement.real_space.rigid_body.refine(
      residue                     = self.residue,
      density_map                 = self.target_map,
      geometry_restraints_manager = restraints_manager,
      real_space_target_weight    = 1,
      real_space_gradients_delta  = 2.0*0.25,
      lbfgs_termination_params    = lbfgs_termination_params,
      unit_cell                   = self.unit_cell,
      cctbx_geometry_restraints_flags = flags)
    return minimized.sites_cart_residue
