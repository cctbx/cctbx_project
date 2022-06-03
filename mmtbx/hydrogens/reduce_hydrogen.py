from __future__ import absolute_import, division, print_function
import six
import sys, time
import mmtbx.model
import iotbx.pdb
import boost_adaptbx.boost.python as bp
from libtbx.utils import null_out
from libtbx import group_args
from cctbx.array_family import flex
from mmtbx.ligands.ready_set_utils import add_n_terminal_hydrogens_to_residue_group
from cctbx.geometry_restraints.linking_class import linking_class
#
from cctbx.maptbx.box import shift_and_box_model

ext = bp.import_ext("cctbx_geometry_restraints_ext")
get_class = iotbx.pdb.common_residue_names_get_class
# For development
print_time = False

# ==============================================================================

def mon_lib_query(residue, mon_lib_srv):
    md, ani = mon_lib_srv.get_comp_comp_id_and_atom_name_interpretation(
      residue_name=residue.resname,
      atom_names=residue.atoms().extract_name())
    return md

# ==============================================================================

class place_hydrogens():
  '''
  Add H atoms to a model

  Parameters
  ----------
  use_neutron_distances : bool
    use neutron distances instead of X-ray

  adp_scale : float
    scale factor for isotropic B of H atoms.
    B(H-atom) = adp_scale * B(parent non-H atom)

  keep_existing_H : bool
    keep existing H atoms in model, only place missing H
  '''

# ------------------------------------------------------------------------------

  def __init__(self,
               model,
               use_neutron_distances = False,
               n_terminal_charge     = 'residue_one',
               adp_scale             = 1,
               exclude_water         = True,
               stop_for_unknowns     = False,
               keep_existing_H       = False,
               validate_e            = False):

    self.model                 = model
    self.use_neutron_distances = use_neutron_distances
    self.n_terminal_charge     = n_terminal_charge
    self.adp_scale             = adp_scale
    self.exclude_water         = exclude_water
    self.stop_for_unknowns     = stop_for_unknowns
    self.keep_existing_H       = keep_existing_H
    self.validate_e            = validate_e
    #
    self.no_H_placed_mlq        = list()
    self.site_labels_disulfides = list()
    self.site_labels_no_para    = list()
    self.charged_atoms          = list()
    self.sl_removed             = list()
    self.n_H_initial            = 0
    self.n_H_final              = 0

# ------------------------------------------------------------------------------

  def run(self):
    '''
    Function that places H atoms
    '''
    model_has_bogus_cs = False

    # TODO temporary fix until the code is moved to model class
    # check if box cussion of 5 A is enough to prevent symm contacts
    cs = self.model.crystal_symmetry()
    if (cs is None) or (cs.unit_cell() is None):
      self.model = shift_and_box_model(model = self.model)
      model_has_bogus_cs = True

    # Remove existing H if requested
    self.n_H_initial = self.model.get_hd_selection().count(True)
    if not self.keep_existing_H:
      self.model = self.model.select(~self.model.get_hd_selection())

    t0 = time.time()
    # Add H atoms and place them at center of coordinates
    pdb_hierarchy = self.add_missing_H_atoms_at_bogus_position()
    if print_time:
      print("add_missing_H_atoms_at_bogus_position:", round(time.time()-t0, 2))

    #print(pdb_hierarchy.composition().n_hd)

#    f = open("intermediate1.pdb","w")
#    f.write(self.model.model_as_pdb())

    # place N-terminal propeller hydrogens
#    if self.n_terminal_charge != 'no_charge':
#      for m in pdb_hierarchy.models():
#        for chain in m.chains():
#          rgs = chain.residue_groups()[0]
#          # by default, place NH3 only at residue with resseq 1
#          if (self.n_terminal_charge == 'residue_one' and rgs.resseq_as_int() != 1):
#            continue
#          elif (self.n_terminal_charge == 'first_in_chain'):
#            pass
#          for ag in rgs.atom_groups():
#            if (get_class(name=ag.resname) in
#                ['common_amino_acid', 'modified_amino_acid', 'd_amino_acid']):
#              if ag.get_atom('H'):
#                ag.remove_atom(ag.get_atom('H'))
#            rc = add_n_terminal_hydrogens_to_residue_group(rgs)
#            # rc is always empty list?

    pdb_hierarchy.sort_atoms_in_place()
    pdb_hierarchy.atoms().reset_serial()
#    f = open("intermediate2.pdb","w")
#    f.write(self.model.model_as_pdb())

    p = mmtbx.model.manager.get_default_pdb_interpretation_params()
    p.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
    p.pdb_interpretation.use_neutron_distances = self.use_neutron_distances
    p.pdb_interpretation.proceed_with_excessive_length_bonds=True
    #p.pdb_interpretation.automatic_linking.link_metals = True

    t0 = time.time()
    #p.pdb_interpretation.restraints_library.cdl=False # XXX this triggers a bug !=360
    ro = self.model.get_restraint_objects()
    self.model = mmtbx.model.manager(
      model_input       = None,
      pdb_hierarchy     = pdb_hierarchy,
      stop_for_unknowns = self.stop_for_unknowns,
      crystal_symmetry  = self.model.crystal_symmetry(),
      restraint_objects = ro,
      log               = null_out())
    self.model.process(pdb_interpretation_params=p,
      make_restraints=True)
    if print_time:
      print("get new model obj and grm:", round(time.time()-t0, 2))

    #f = open("intermediate3.pdb","w")
    #f.write(self.model.model_as_pdb())

    # Only keep H that have been parameterized in riding H procedure
    sel_h = self.model.get_hd_selection()
    if sel_h.count(True) == 0:
      return

    # get rid of isolated H atoms.
    #For example when heavy atom is missing, H needs not to be placed
    sel_isolated = self.model.isolated_atoms_selection()
    self.sel_lone_H = sel_h & sel_isolated
    self.model = self.model.select(~self.sel_lone_H)

    t0 = time.time()
    # get riding H manager --> parameterize all H atoms
    sel_h = self.model.get_hd_selection()
    self.model.setup_riding_h_manager(use_ideal_dihedral = True)
    riding_h_manager = self.model.riding_h_manager
    if riding_h_manager is None:
      return
    sel_h_in_para = flex.bool(
      [bool(x) for x in riding_h_manager.h_parameterization])
    sel_h_not_in_para = sel_h_in_para.exclusive_or(sel_h)
    self.site_labels_no_para = [atom.id_str().replace('pdb=','').replace('"','')
      for atom in self.model.get_hierarchy().atoms().select(sel_h_not_in_para)]
    #
    self.model = self.model.select(~sel_h_not_in_para)
    if print_time:
      print("set up riding H manager and some cleanup:", round(time.time()-t0, 2))

  #  f = open("intermediate4.pdb","w")
  #  f.write(model.model_as_pdb())

    if self.validate_e:
      t0 = time.time()
      self.validate_electrons()
      if print_time:
        print("validate electrons:", round(time.time()-t0, 2))

    t0 = time.time()
    # Reset occupancies, ADPs and idealize H atom positions
    self.model.reset_adp_for_hydrogens(scale = self.adp_scale)
    self.model.reset_occupancy_for_hydrogens_simple()
    self.model.idealize_h_riding()
    if print_time:
      print("reset adp, occ; idealize:", round(time.time()-t0, 2))

    t0 = time.time()
    self.exclude_H_on_links()
    if print_time:
      print("all links:", round(time.time()-t0, 2))

    # place N-terminal propeller hydrogens
    if self.n_terminal_charge != 'no_charge':
      hierarchy = self.model.get_hierarchy()
      for m in hierarchy.models():
        for chain in m.chains():
          rgs = chain.residue_groups()[0]
          # by default, place NH3 only at residue with resseq 1
          if (self.n_terminal_charge == 'residue_one' and rgs.resseq_as_int() != 1):
            continue
          elif (self.n_terminal_charge == 'first_in_chain'):
            pass
          add_charge = True
          for ag in rgs.atom_groups():
            if (get_class(name=ag.resname) in
                ['common_amino_acid', 'modified_amino_acid', 'd_amino_acid']):
              if ag.get_atom('N'):
                N = ag.get_atom('N')
                if N.i_seq in self.exclusion_iseqs:
                  add_charge = False
              if ag.get_atom('H'):
                H = ag.get_atom('H')
                ag.remove_atom(H)
                H_label = H.id_str().replace('pdb=','').replace('"','')
                if H_label in self.site_labels_no_para:
                  self.site_labels_no_para.remove(H_label)
            if add_charge:
              rc = add_n_terminal_hydrogens_to_residue_group(rgs)
      hierarchy.sort_atoms_in_place()
      hierarchy.atoms().reset_serial()
      self.model = mmtbx.model.manager(
        model_input       = None,
        pdb_hierarchy     = hierarchy,
        stop_for_unknowns = self.stop_for_unknowns,
        crystal_symmetry  = self.model.crystal_symmetry(),
        restraint_objects = ro,
        log               = null_out())

    self.n_H_final = self.model.get_hd_selection().count(True)

# ------------------------------------------------------------------------------

  def add_missing_H_atoms_at_bogus_position(self):
    '''Add missing H atoms at bogus positions to the pdb_hierarchy

    This procedure changes the hierarchy in place.
    All H atoms are placed at center of coordinates + (0.5, 0.5, 0.5)
    The translation is necessary because sometimes the center of coordinates
    coincides with the position of a heavy atom.

    In one residue/entity, all newly placed H are superposed, they will be
    moved to their expected position later.
    '''
    # TODO temporary fix until v3 names are in mon lib
    alternative_names = [
      ('HA1', 'HA2', 'HA3'),
      ('HB1', 'HB2', 'HB3'),
      ('HG1', 'HG2', 'HG3'),
      ('HD1', 'HD2', 'HD3'),
      ('HE1', 'HE2', 'HE3'),
      ('HG11', 'HG12', 'HG13')
      ]
    # end TODO
    pdb_hierarchy = self.model.get_hierarchy()
    mon_lib_srv = self.model.get_mon_lib_srv()
    #XXX This breaks for 1jxt, residue 2, TYR
    #get_class = iotbx.pdb.common_residue_names_get_class
    #no_H_placed_resnames = list()
    for m in pdb_hierarchy.models():
      for chain in m.chains():
        for rg in chain.residue_groups():
          n_atom_groups = len(rg.atom_groups())
          for ag in rg.atom_groups():
            if n_atom_groups == 3 and ag.altloc == '':
              continue
            #print list(ag.atoms().extract_name())
            if(get_class(name=ag.resname) == "common_water"): continue
            actual = [a.name.strip().upper() for a in ag.atoms()]
            #
            mlq = mon_lib_query(residue=ag, mon_lib_srv=mon_lib_srv)
            #if (get_class(name=ag.resname) in ['modified_rna_dna', 'other']):
            if mlq is None:
              self.no_H_placed_mlq.append(ag.resname)
              continue

            expected_h = list()
            atom_dict = mlq.atom_dict()
            for k, v in six.iteritems(atom_dict):
              if(v.type_symbol=="H"):
                expected_h.append(k)
            # TODO start: temporary fix until v3 names are in mon lib
            for altname in alternative_names:
              if (altname[0] in expected_h and altname[1] in expected_h):
                if (atom_dict[altname[0]].type_energy == 'HCH2' and
                    atom_dict[altname[1]].type_energy == 'HCH2'):
                  expected_h.append(altname[2])
                  expected_h.remove(altname[0])
            # TODO end
            missing_h = list(set(expected_h).difference(set(actual)))
            if 0: print(ag.resname, missing_h)
            #new_xyz = ag.atoms().extract_xyz().mean()
            new_xyz = flex.double(ag.atoms().extract_xyz().mean()) + \
              flex.double([0.5,0.5,0.5])
            new_xyz = tuple(new_xyz)

            hetero = ag.atoms()[0].hetero
            segid = ag.atoms()[0].segid

            for mh in missing_h:
              # TODO: this should be probably in a central place
              if len(mh) < 4: mh = (' ' + mh).ljust(4)
              a = (iotbx.pdb.hierarchy.atom()
                .set_name(new_name=mh)
                .set_element(new_element="H")
                .set_xyz(new_xyz=new_xyz)
                .set_hetero(new_hetero=hetero)
                .set_segid(new_segid=segid))

              ag.append_atom(a)

    return pdb_hierarchy

# ------------------------------------------------------------------------------

  def validate_electrons(self):
    from elbow.quantum import electrons
    atom_valences = electrons.electron_distribution(
      # XXX How do we get this working on models with alternate locations?
      self.model.get_hierarchy(), # needs to be altloc free
      self.model.get_restraints_manager().geometry,
      verbose=False,
    )
    atom_valences.validate(ignore_water=True,
                           raise_if_error=False)
    self.charged_atoms = atom_valences.get_charged_atoms()

# ------------------------------------------------------------------------------

  def exclude_H_on_links(self):
    """Remove H atoms bound to heavy atoms that form a link

    An exception are HD1 and HE2 of HIS. The mover functionality in reduce will
    take care of those.
    """
    origin_ids = linking_class()
    grm = self.model.get_restraints_manager()
    bond_proxies_simple, asu = grm.geometry.get_all_bond_proxies(
      sites_cart = self.model.get_sites_cart())
    elements = self.model.get_hierarchy().atoms().extract_element()
    exclusion_iseqs = list()
    exclusion_dict = dict()
    all_proxies = [p for p in bond_proxies_simple]
    for proxy in asu:
      all_proxies.append(proxy)
    # Loop through bond proxies to find links (i.e. proxies with origin_id != 0)
    for proxy in all_proxies:
      if(  isinstance(proxy, ext.bond_simple_proxy)): i,j=proxy.i_seqs
      elif(isinstance(proxy, ext.bond_asu_proxy)):    i,j=proxy.i_seq,proxy.j_seq
      else: assert 0 # never goes here
      if proxy.origin_id != 0:
        exclusion_iseqs.extend([i,j])
        exclusion_dict[i] = proxy.origin_id
        exclusion_dict[j] = proxy.origin_id
    sel_remove = flex.size_t()

    atoms = self.model.get_atoms()
    # Find H atoms bound to linked atoms
    removed_dict = dict()
    for proxy in all_proxies:
      if(  isinstance(proxy, ext.bond_simple_proxy)): i,j=proxy.i_seqs
      elif(isinstance(proxy, ext.bond_asu_proxy)):    i,j=proxy.i_seq,proxy.j_seq
      else: assert 0 # never goes here
      # Exception for HIS HD1 and HE2
      if (atoms[i].parent().resname == 'HIS' and
        atoms[i].name.strip() in ['HD1','DD1', 'HE2', 'DE2']): continue
      if (atoms[j].parent().resname == 'HIS' and
        atoms[j].name.strip() in ['HD1','DD1', 'HE2', 'DE2']): continue
      if(elements[i] in ["H","D"] and j in exclusion_iseqs):
        sel_remove.append(i)
        removed_dict[i] = exclusion_dict[j]
      if(elements[j] in ["H","D"] and i in exclusion_iseqs):
        sel_remove.append(j)
        removed_dict[j] = exclusion_dict[i]
    #
    sl_removed = [(atom.id_str().replace('pdb=','').replace('"',''),
                   origin_ids.get_origin_key(removed_dict[atom.i_seq]))
        for atom in self.model.get_hierarchy().atoms().select(sel_remove)]
    #
    self.model = self.model.select(~flex.bool(self.model.size(), sel_remove))
    self.sl_removed = sl_removed
    self.exclusion_iseqs = exclusion_iseqs


# ------------------------------------------------------------------------------

  def show(self, log):
    '''
    Informative output
    '''
    if log is None: log = sys.stdout
    #
    if (not self.keep_existing_H and self.n_H_initial):
      msg = 'Number of hydrogen atoms trimmed from input model: %s \n'
      print(msg % self.n_H_initial, file=log)
    #
    msg = 'Number of hydrogen atoms added to the input model: %s \n'
    print(msg % self.n_H_final, file=log)
    #
    if self.no_H_placed_mlq:
      msg = '''
No H atoms were placed on the following residues because no restraints
were found:'''
      print(msg, file=log)
      for resname in self.no_H_placed_mlq:
        print(resname, file=log)
    #
    if self.site_labels_disulfides:
      msg = '''
The following cysteine HG atoms were not placed because the sulfur atom is
involved in a disulfide bond'''
      print(msg, file=log)
      for label in self.site_labels_disulfides:
        print(label, file=log)
    #
    if self.site_labels_no_para:
      msg = '''
The following H atoms were not placed because they could not be parameterized
(not enough restraints information)'''
      print(msg, file=log)
      for label in self.site_labels_no_para:
        print(label)
    if self.charged_atoms:
      msg = '''
The following heavy atom have an unusual electron count. This could be because
heavy atoms or H atoms are missing.'''
      print(msg, file=log)
      for item in self.charged_atoms:
        idstr = item[0].id_str().replace('pdb=','').replace('"','')
        if 'HOH' in idstr: continue
        print(idstr, item[1])

    if self.sl_removed:
      print()
      msg = '''Atom %s was not placed because it is involved in %s'''
      for item in self.sl_removed:
        print(msg % (item[0], item[1]), file=log)

# ------------------------------------------------------------------------------

  def get_model(self):
    return self.model

  def get_counts(self):
    return group_args(
      number_h_final  = self.n_H_final,
      no_H_placed_mlq = self.no_H_placed_mlq,
      site_labels_disulfides = self.site_labels_disulfides,
      site_labels_no_para = self.site_labels_no_para)

# ==============================================================================

# stub for reduce parameters
# TODO can be parameters or phil, depending on how many options are really needed
reduce_master_params_str = """
flip_NQH = True
  .type = bool
  .help = add H and rotate and flip NQH groups
search_time_limit = 600
  .type = int
  .help = max seconds to spend in exhaustive search (default=600)
"""

def optimize(model):
  """
  Carry out reduce optimization

  Parameters
  ----------
  model
      mmtbx model object that contains H atoms
      H atoms should be at approprite distances

  Returns
  -------
  model
      mmtbx model object with optimized H atoms
  """
  # hierarchy object --> has hierarchy of structure
  pdb_hierarchy = model.get_hierarchy()
  # geometry restraints manager --> info about ideal bonds, angles; what atoms are bonded, etc.
  grm = model.get_restraints_manager()

  print("Reduce optimization happens here")

  return model
