from __future__ import absolute_import, division, print_function
import time
from cctbx import sgtbx
from cctbx import geometry_restraints

from mmtbx.monomer_library import linking_setup
from mmtbx.monomer_library import linking_utils
from mmtbx.monomer_library import glyco_utils
from mmtbx.monomer_library import bondlength_defaults
from libtbx.utils import Sorry
from functools import cmp_to_key

origin_ids = geometry_restraints.linking_class.linking_class()

hydrogens = ["H", "D", "T"]

class ResidueLinkClass(dict):
  def remove_link(self, residues, atom_names):
    if residues in self:
      if atom_names in self[residues]:
        self[residues].remove(atom_names)

def _apply_link_using_proxies(link,
                              atom_group1,
                              atom_group2,
                              bond_params_table,
                              bond_asu_table,
                              geometry_proxy_registries,
                        #      distance,
                              rt_mx_ji,
                              origin_id=None,
                              ):
  assert origin_id
  ######################################
  def _get_restraint_i_seqs(atom_group1,
                            atom_group2,
                            restraint,
                            ):
    i_seqs = []
    keys = restraint.cif_keywords()
    if "value_dist" in keys:
      attrs = [
        "atom_1_comp_id",
        "atom_id_1",
        "atom_2_comp_id",
        "atom_id_2",
        ]
    elif "period" in keys:
      attrs = [
        "atom_1_comp_id",
        "atom_id_1",
        "atom_2_comp_id",
        "atom_id_2",
        "atom_3_comp_id",
        "atom_id_3",
        "atom_4_comp_id",
        "atom_id_4",
        ]
    elif "value_angle" in keys:
      attrs = [
        "atom_1_comp_id",
        "atom_id_1",
        "atom_2_comp_id",
        "atom_id_2",
        "atom_3_comp_id",
        "atom_id_3",
        ]
    elif "volume_sign" in keys:
      attrs = [
        "atom_centre_comp_id",
        "atom_id_centre",
        "atom_1_comp_id",
        "atom_id_1",
        "atom_2_comp_id",
        "atom_id_2",
        "atom_3_comp_id",
        "atom_id_3",
        ]
    elif "plane_id" in keys:
      attrs = [
        "atom_comp_id",
        "atom_id",
        ]
    else:
      assert 0
    for i, attr in enumerate(attrs):
      if i%2:
        # name
        name = getattr(restraint, attr)
        for atom in atoms:
          # uses names to confirm link
          if atom.name.strip()==name.strip():
            i_seqs.append(atom.i_seq)
            break
        else:
          # name not found, could be hydrogen or ...
          return None
      else:
        # atoms
        if getattr(restraint, attr)==1:
          atoms = atom_group1.atoms()
        else:
          atoms = atom_group2.atoms()
    return i_seqs
  ###############
  def _check_i_seqs(atom_group1, atom_group2, i_seqs):
    atoms = []
    for i_seq in i_seqs:
      for atom in list(atom_group1.atoms())+list(atom_group2.atoms()):
        if atom.i_seq==i_seq:
          atoms.append(atom)
          break
    d2 = linking_utils.get_distance2(*atoms) # XXXX needs to be sym aware
    if d2>9: return False
    return True
  #############
  assert link
  count = 0
  #
  bond_i_seqs = []
  for bond in link.bond_list:
    i_seqs = _get_restraint_i_seqs(atom_group1,
                                   atom_group2,
                                   bond,
      )
    if i_seqs is None: continue
    if not _check_i_seqs(atom_group1, atom_group2, i_seqs): # check distances
      tmp = atom_group2
      atom_group2 = atom_group1
      atom_group1 = tmp
      i_seqs = _get_restraint_i_seqs(atom_group1,
                                     atom_group2,
                                     bond,
        )
      if i_seqs is None: continue
    value = "value_dist"
    assert origin_id
    proxy = geometry_restraints.bond_simple_proxy(
      i_seqs=i_seqs,
      distance_ideal=getattr(bond, value),
      weight=1/bond.value_dist_esd**2,
      origin_id=origin_id,
      )
    bond_params_table.update(i_seq=i_seqs[0],
                             j_seq=i_seqs[1],
                             params=proxy)
    #if rt_mx_ji is None: continue
    bond_asu_table.add_pair(
      i_seq=i_seqs[0],
      j_seq=i_seqs[1],
      rt_mx_ji=rt_mx_ji,
      )
    count+=1
    bond_i_seqs.append(i_seqs)
  #
  for angle in link.angle_list:
    i_seqs = _get_restraint_i_seqs(atom_group1,
                                   atom_group2,
                                   angle,
        )
    if i_seqs is None: continue
    proxy = geometry_restraints.angle_proxy(
      i_seqs=i_seqs,
      angle_ideal=angle.value_angle,
      weight=1/angle.value_angle_esd**2,
      origin_id=origin_id,
      )
    geometry_proxy_registries.angle.add_if_not_duplicated(proxy=proxy)
  #
  for tor in link.tor_list:
    i_seqs = _get_restraint_i_seqs(atom_group1,
                                   atom_group2,
                                   tor,
        )
    if i_seqs is None: continue
    proxy = geometry_restraints.dihedral_proxy(
      i_seqs=i_seqs,
      angle_ideal=tor.value_angle,
      weight=1/tor.value_angle_esd**2,
      periodicity=tor.period,
      origin_id=origin_id,
      )
    geometry_proxy_registries.dihedral.add_if_not_duplicated(proxy=proxy)
  #
  for chir in link.chir_list:
    i_seqs = _get_restraint_i_seqs(atom_group1,
                                   atom_group2,
                                   chir,
        )
    if i_seqs is None: continue
    volume_ideal = 2.4
    if chir.volume_sign[:4].lower()=="nega":
      volume_ideal = -2.4
    elif chir.volume_sign[:4].lower()=="zero":
      volume_ideal = 0.
    both_signs=False
    if chir.volume_sign=='both': both_signs=True
    proxy = geometry_restraints.chirality_proxy(
      i_seqs=i_seqs,
      volume_ideal=volume_ideal,
      both_signs=both_signs,
      weight=25.,
      origin_id=origin_id,
      )
    geometry_proxy_registries.chirality.add_if_not_duplicated(proxy=proxy)
  #
  planes = {}
  weights = {}
  for plane in link.plane_list:
    i_seqs = _get_restraint_i_seqs(atom_group1,
                                   atom_group2,
                                   plane,
        )
    if i_seqs is None: continue
    planes.setdefault(plane.plane_id, [])
    planes[plane.plane_id]+=i_seqs
    weights.setdefault(plane.plane_id, [])
    weights[plane.plane_id].append(1/plane.dist_esd**2)
  if planes:
    for plane_id in planes:
      if len(planes[plane_id])<4: continue
      proxy = geometry_restraints.planarity_proxy(
        i_seqs=planes[plane_id],
        weights=weights[plane_id],
        origin_id=origin_id,
        )
      geometry_proxy_registries.planarity.add_if_not_duplicated(proxy=proxy)
  return count, bond_i_seqs

def possible_cyclic_peptide(atom1,
                            atom2,
                            atoms_in_first_last_rgs,
                            verbose=False,
                            ):
  if verbose:
    print(atom1.quote(),atom2.quote())
  names=[atom1.name, atom2.name]
  names.sort()
  if verbose: print('names',names)
  if names!=[' C  ', ' N  ']: return False
  len_fl = 0
  len_fl += atoms_in_first_last_rgs.get(atom1.i_seq, -1)
  len_fl += atoms_in_first_last_rgs.get(atom2.i_seq, -1)
  if len_fl == 1:
    chain1 = atom1.parent().parent().parent()
    chain2 = atom2.parent().parent().parent()
    if chain1.id == chain2.id:
      return True
    elif verbose: print('chain id differs', chain1.id, chain2.id)
  return False

def check_for_peptide_links(atom1,
                            atom2,
                            classes1,
                            classes2,
                            ):
  atom_group1 = atom1.parent()
  atom_group2 = atom2.parent()
  if classes1.common_amino_acid or classes2.common_amino_acid:
    if classes1.common_amino_acid:
      other = atom_group2
    else:
      other = atom_group1
    #
    # only for AA to other
    #
    if linking_utils.get_class(other.resname) not in ["other"]:
      return None
    # sulfur bridge
    key = "SS"
    if atom1.name.strip()=="SG" and atom2.name.strip()=="SG": # too strict?
      return key, False
    # test for backbond atoms
    count = 0
    for atom in other.atoms():
      if atom.name.strip() in ["C", "N", "O"]: count+=1
    if count!=3: return None
    key = "TRANS"
    if atom2.name.strip()=="C" and atom1.name.strip()=="N":
      return key, False
    elif atom1.name.strip()=="C" and atom2.name.strip()=="N":
      return key, True
  return False

def check_all_classes(pdb_hierarchy, class_type):
  """ This one is not used anywhere.
  If start using it, beware of getting classes for all atoms
  """
  assert 0
  found = False
  for residue_group in pdb_hierarchy.residue_groups():
    for atom in residue_group.atoms():
      classes = linking_utils.get_classes(atom)
      if getattr(classes, class_type):
        found = True
      break
    if found: break
  return found

class linking_mixins(object):
  def process_nonbonded_for_links(self,
                                  bond_params_table,
                                  bond_asu_table,
                                  geometry_proxy_registries,
                                  link_metals                 = True,
                                  link_residues               = True,
                                  link_carbohydrates          = True,
                                  link_amino_acid_rna_dna     = False,
                                  link_ligands                = False,
                                  link_small_molecules        = False,
                                  max_bonded_cutoff           = None,
                                  metal_coordination_cutoff   = 3.,
                                  amino_acid_bond_cutoff      = 2.,
                                  inter_residue_bond_cutoff   = 2.,
                                  second_row_buffer           = 0.5,
                                  carbohydrate_bond_cutoff    = 2.,
                                  ligand_bond_cutoff          = 2.,
                                  small_molecule_bond_cutoff  = 2.,
                                  include_selections          = None,
                                  exclude_selections          = None,
                                  exclude_hydrogens_from_bonding_decisions = False,
                                  log                         = None,
                                  verbose                     = False,
                                  ):
    assert hasattr(self, "_cif")
    if max_bonded_cutoff is None:
      max_bonded_cutoff = max(metal_coordination_cutoff,
                              amino_acid_bond_cutoff,
                              carbohydrate_bond_cutoff,
                              ligand_bond_cutoff,
                              small_molecule_bond_cutoff,
                              inter_residue_bond_cutoff+second_row_buffer,
                              )
    max_bonded_cutoff_standard = max_bonded_cutoff
    if include_selections:
      for selection_1, selection_2, cutoff in include_selections:
        max_bonded_cutoff = max(max_bonded_cutoff, cutoff)

    if max_bonded_cutoff > 15:
      raise Sorry("One of the following parameters: \nmetal_coordination_"+
          "cutoff, amino_acid_bond_cutoff,"+
          "inter_residue_bond_cutoff, \ncarbohydrate_bond_cutoff,"
          "bonds.bond_distance_cutoff \nis greater than 15A. Please check and"+
          " correct these parameters.")
    if verbose and log is not None:
      print("""
      metal_coordination_cutoff %s
      amino_acid_bond_cutoff    %s
      carbohydrate_bond_cutoff  %s
      inter_residue_bond_cutoff %s
      second_row_buffer         %s
      """ % ( metal_coordination_cutoff,
              amino_acid_bond_cutoff,
              carbohydrate_bond_cutoff,
              inter_residue_bond_cutoff,
              second_row_buffer,
              ), file=log)
    from cctbx import crystal
    from cctbx.array_family import flex
    #
    def _nonbonded_pair_objects(max_bonded_cutoff=3., i_seqs=None):
      if i_seqs is None:
        atoms = self.pdb_hierarchy.atoms()
        i_seqs = flex.size_t()
        for atom in atoms:
          i_seqs.append(atom.i_seq)
      if (self.model_indices is not None):
        model_indices = self.model_indices.select(i_seqs)
      conformer_indices = self.conformer_indices.select(i_seqs)
      sym_excl_indices = self.sym_excl_indices.select(i_seqs)
      donor_acceptor_excl_groups = self.donor_acceptor_excl_groups.select(i_seqs)
      asu_mappings = self.special_position_settings.asu_mappings(
        buffer_thickness=max_bonded_cutoff)
      sites_cart = self.sites_cart.select(i_seqs)
      asu_mappings.process_sites_cart(
        original_sites=sites_cart,
        site_symmetry_table=self.site_symmetry_table().select(i_seqs))
      pair_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
      nonbonded_proxies = geometry_restraints.nonbonded_sorted_asu_proxies(
        model_indices=model_indices,
        conformer_indices=conformer_indices,
        sym_excl_indices=sym_excl_indices,
        donor_acceptor_excl_groups=donor_acceptor_excl_groups,
        nonbonded_params=geometry_restraints.nonbonded_params(
          default_distance=1),
        nonbonded_types=flex.std_string(conformer_indices.size()),
        nonbonded_charges=flex.int(conformer_indices.size(), 0),
        nonbonded_distance_cutoff_plus_buffer=max_bonded_cutoff,
        min_cubicle_edge=5,
        shell_asu_tables=[pair_asu_table])
      return nonbonded_proxies, sites_cart, pair_asu_table, asu_mappings, i_seqs
    #
    if(log is not None):
      print("""
  Automatic linking
    Parameters for automatic linking
      Linking & cutoffs
        Metal                : %-5s - %0.2f
        Amino acid           : %-5s - %0.2f
        Carbohydrate         : %-5s - %0.2f
        Ligands              : %-5s - %0.2f
        Small molecules      : %-5s - %0.2f
        Amino acid - RNA/DNA : %-5s
      """ % (link_metals,
             metal_coordination_cutoff,
             link_residues,
             amino_acid_bond_cutoff,
             link_carbohydrates,
             carbohydrate_bond_cutoff,
             link_ligands,
             ligand_bond_cutoff,
             link_small_molecules,
             small_molecule_bond_cutoff,
             link_amino_acid_rna_dna,
             ), file=log)
    t0=time.time()
    atoms = self.pdb_hierarchy.atoms()
    bond_data = []
    bond_data_i_seqs = {}
    simple_bonds = 0
    sym_bonds = 0
    link_data = []
    simple_links = 0
    sym_links = 0
    done = ResidueLinkClass()
    links = {}
    custom_links = {}
    exclude_out_lines = {}

    atom_classes = [linking_utils.get_classes(a) for a in atoms]
    atoms_in_first_last_rgs = {}
    for c in self.pdb_hierarchy.chains():
      for i in [0,-1]:
        if not c.residue_groups(): continue
        rg = c.residue_groups()[i]
        abs_i = abs(i)
        for a in rg.atoms():
          atoms_in_first_last_rgs[a.i_seq] = abs_i
    skip_if_longer = linking_setup.update_skip_if_longer(amino_acid_bond_cutoff,
                                                        rna_dna_bond_cutoff=3.4,
                                                        intra_residue_bond_cutoff=inter_residue_bond_cutoff,
                                                        saccharide_bond_cutoff=carbohydrate_bond_cutoff,
                                                        metal_coordination_cutoff=metal_coordination_cutoff,
                                                        sulfur_bond_cutoff=2.5,
                                                        other_bond_cutoff=2,
                                                        )

    # main loop
    nonbonded_proxies, sites_cart, pair_asu_table, asu_mappings, nonbonded_i_seqs = \
        _nonbonded_pair_objects(max_bonded_cutoff=max_bonded_cutoff,
          )
    initial_pair_asu_table_table = bond_asu_table.table().deep_copy()
    for ii, item in enumerate(nonbonded_proxies.sorted_value_proxies_generator(
        by_value="delta",
        sites_cart=sites_cart,
        cutoff=max_bonded_cutoff)):
      i_seq, j_seq, distance, sym_op, rt_mx_ji, proxy = item
      #
      # include & exclude selection
      #
      origin_id = None
      updated_skip_if_longer = None
      if ( include_selections and
           distance>=max_bonded_cutoff_standard
           ):
        for selection_1, selection_2, bond_cutoff in include_selections:
          if (( i_seq in selection_1 and j_seq in selection_2 ) or
              ( i_seq in selection_2 and j_seq in selection_1 )
              ):
            metal_coordination_cutoff=bond_cutoff
            amino_acid_bond_cutoff=bond_cutoff
            carbohydrate_bond_cutoff=bond_cutoff
            ligand_bond_cutoff=bond_cutoff
            small_molecule_bond_cutoff=bond_cutoff
            inter_residue_bond_cutoff=bond_cutoff
            saccharide_bond_cutoff=bond_cutoff
            link_residues=True
            updated_skip_if_longer = linking_setup.update_skip_if_longer(
                amino_acid_bond_cutoff,
                rna_dna_bond_cutoff=3.4,
                intra_residue_bond_cutoff=inter_residue_bond_cutoff,
                saccharide_bond_cutoff=saccharide_bond_cutoff,
                metal_coordination_cutoff=metal_coordination_cutoff,
                sulfur_bond_cutoff=2.5,
                other_bond_cutoff=2,
                )
            break
        else:
          continue # exclude this nonbond from consideration
      exclude_this_nonbonded = False
      if exclude_selections:
        for selection_1, selection_2 in exclude_selections:
          if selection_2: # check both
            if (( i_seq in selection_1 and j_seq in selection_2 ) or
                ( i_seq in selection_2 and j_seq in selection_1 )
                ):
              exclude_this_nonbonded = True
              break
          else:
            if i_seq in selection_1 or j_seq in selection_1:
              exclude_this_nonbonded = True
      # XXX this is a poor job!!!
      if bond_asu_table.contains(i_seq, j_seq, 0): continue
      if bond_asu_table.contains(i_seq, j_seq, 1): continue
      atom1 = atoms[i_seq]
      atom2 = atoms[j_seq]
      if atom1.element_is_hydrogen() or atom2.element_is_hydrogen():
        continue
      if exclude_this_nonbonded:
        key = (selection_1,selection_2)
        if key not in exclude_out_lines:
          exclude_out_lines[key] = \
              '    bond %s\n%s%s\n%smodel %.2f' % (
                atom1.id_str(),
                ' '*9,
                atom2.id_str(),
                ' '*6,
                distance)
        continue
      #
      # moving sections of this to outside the loop
      #  - SF4
      #  - ZN-CYS-HIS
      #
      moved = ['SF4', 'F3S', 'FES']
      if ( atom1.parent().resname in moved or
           atom2.parent().resname in moved
           ): continue
      moved = ['ZN', 'CYS']
      if ( atom1.parent().resname.strip() in moved and
           atom2.parent().resname.strip() in moved
           ): continue
      moved = ['ZN', 'HIS']
      if ( atom1.parent().resname.strip() in moved and
           atom2.parent().resname.strip() in moved
           ): continue
      #
      #
      if verbose:
        print('='*80)
        print('nonbonded', i_seq, j_seq, atom1.quote(), end=' ')
        print(atom2.quote(), end=' ')
        print("Distance: %0.2f" % distance, rt_mx_ji, sym_op)

      # don't link atoms not in the same conformer (works for models also)...
      # They will never be in different conformers, we are looping over
      # nonbonded interactions here! This check takes 14% of this whole function.
      # if not atom1.is_in_same_conformer_as(atom2):
      #   assert 0
      #   continue
      # don't link atoms in same residue group
      if atom1.parent().parent()==atom2.parent().parent(): continue
      atom_group1 = atom1.parent()
      atom_group2 = atom2.parent()
      # dont't like atom groups in different altloc expect " "
      if atom_group1.altloc.strip()==atom_group2.altloc.strip(): pass
      elif atom_group1.altloc.strip()=="": pass
      elif atom_group2.altloc.strip()=="": pass
      else: continue
      # don't link some classes
      classes1 = atom_classes[i_seq]
      classes2 = atom_classes[j_seq]
      use_only_bond_cutoff = False
      if verbose:
        print("""
Residue classes

%s
%s """ % (classes1, classes2))
      # why was this commented out???
      if not link_ligands and (classes1.other or classes2.other): continue
      if (not link_small_molecules and
          (classes1.common_small_molecule or classes2.common_small_molecule)): continue
      # is_proxy_set between any of the atoms ????????
      if linking_utils.allow_cis_trans(classes1, classes2):
        if not link_residues:
          continue
        # special amino acid linking
        #  - cyclic
        #  - beta, delta ???
        if possible_cyclic_peptide(atom1, atom2, atoms_in_first_last_rgs): # first & last peptide
          if verbose: print('possible_cyclic_peptide')
          use_only_bond_cutoff = True
      if sym_op:
        if classes1.common_amino_acid and classes2.common_saccharide: continue
        if classes2.common_amino_acid and classes1.common_saccharide: continue
      #
      # bonded atoms can't link to same atom, eg MG-PG and MG-O1P
      #
      if bond_data:
        bonded = False
        if i_seq in bond_data_i_seqs:
          for t_i in bond_data_i_seqs[i_seq]:
            if bond_asu_table.contains(j_seq, t_i, 0): bonded = True
        if j_seq in bond_data_i_seqs:
          for t_i in bond_data_i_seqs[j_seq]:
            if bond_asu_table.contains(i_seq, t_i, 0): bonded = True
        if bonded: continue
      #
      key = [
        atom1.id_str()[9:-1],
        atom2.id_str()[9:-1],
        ]
      key.sort()
      if sym_op:
        key.append(str(sym_op))
      key = tuple(key)
      # hydrogens
      if not exclude_hydrogens_from_bonding_decisions:
        if atom1.element.strip() in hydrogens:
          done[atom2.id_str()] = atom1.id_str()
        if atom2.element.strip() in hydrogens:
          done[atom1.id_str()] = atom2.id_str()
      # bond length cutoff & some logic
      if not linking_utils.is_atom_pair_linked(
          atom1,
          atom2,
          classes1.important_only,
          classes2.important_only,
          distance=distance,
          skip_if_longer=updated_skip_if_longer if updated_skip_if_longer is not None else skip_if_longer,
          second_row_buffer=second_row_buffer,
          saccharide_bond_cutoff=carbohydrate_bond_cutoff,
          metal_coordination_cutoff=metal_coordination_cutoff,
          use_only_bond_cutoff=use_only_bond_cutoff,
          link_metals=link_metals,
          verbose=verbose,
          ):
        if verbose:
          print("is not linked", atom1.quote(),atom2.quote(),key)
          print('link_metals',link_metals)
        if ( atom1.element.strip().upper() in hydrogens or
             atom2.element.strip().upper() in hydrogens):
          if verbose: print('hydrogens')
          pass
        else:
          done.setdefault(key, [])
          done[key].append([atom1.name, atom2.name])
        continue
      # check some valences...
      if not (classes1.common_element or classes2.common_element):
        if not linking_utils.check_valence(self.pdb_hierarchy, atom1):
          print("  Atom %s rejected from bonding due to valence issues." % atom1.quote(), file=log)
          continue
        if not linking_utils.check_valence(self.pdb_hierarchy, atom2):
          print("  Atom %s rejected from bonding due to valence issues." % atom2.quote(), file=log)
          continue
      # got a link....

      class1 = classes1.important_only
      class2 = classes2.important_only
      class_key = [class1, class2]
      class_key.sort()
      class_key = tuple(class_key)
      if verbose: print('class_key',class_key)
      #
      if link_metals!=True and "metal" in class_key: continue
      #atoms_must_be = {}
      if not link_residues:
        if class_key in [
            ("common_amino_acid", "common_amino_acid"),
            ("common_amino_acid", "d_amino_acid"),
            ("common_amino_acid", "uncommon_amino_acid"),
           ]:
          continue
      else:
        if len(key)>2: continue
      #else:
      #  atoms_must_be.setdefault(("common_amino_acid",
      #                            "common_amino_acid"),["C", "N"])
      #  atoms_must_be.setdefault(("common_amino_acid", "other"),["C", "N"])
      if not link_carbohydrates and "common_saccharide" in class_key: continue
      if not link_amino_acid_rna_dna:
        if "common_amino_acid" in class_key and "common_rna_dna" in class_key:
          continue
      #
      names = [atom1.name, atom2.name]
      if verbose: print('names',names)
      names.sort()
      atom1_key = None
      atom2_key = None
      if class1 in linking_setup.maximum_per_atom_links: # is one
        atom1_key = atom1.id_str()
      if class2 in linking_setup.maximum_per_atom_links: # is one
        atom2_key = atom2.id_str()
      if verbose:
        print('-'*80)
        print('class_key',class_key)
        print('done')
        for k, item in done.items():
          print("> %s : %s" % (k, item))
        print('key',key)
        print('atom keys',atom1_key, atom2_key)
      # exclude duplicate symmetry op.
      if key in done:
        if names in done[key]: continue
      if atom1.parent().altloc==atom2.parent().altloc:
        if atom1_key:
          if atom1_key in done: continue
          done[atom1_key] = key
        if atom2_key:
          if atom2_key in done: continue
          done[atom2_key] = key
      if verbose: print(done)
      #
      current_number_of_links = len(done.setdefault(key, []))
      if(current_number_of_links >=
         linking_setup.maximum_inter_residue_links.get(class_key, 1)
         ):
        if verbose:
          print("too many links:",current_number_of_links,linking_setup.maximum_inter_residue_links.get(class_key, 1), class_key)
        continue
      #
      done[key].append(names)
      done_key = key
      # get all possible links
      i_seqs = []
      for atom in atom_group1.atoms():
        i_seqs.append(atom.i_seq)
      j_seqs = []
      for atom in atom_group2.atoms():
        j_seqs.append(atom.i_seq)
      ij_seqs = []
      for i in i_seqs:
        for j in j_seqs:
          tmp = [i,j]
          tmp.sort()
          ij_seqs.append(tuple(tmp))
      # check that a link not already made
      link_found = False
      if verbose:
        print('len simple bond proxies',len(geometry_proxy_registries.bond_simple.proxies))

      # Consistency check - debugging only
      # for bsp in geometry_proxy_registries.bond_simple.proxies:
      #   if bsp.i_seqs[1] not in initial_pair_asu_table_table[bsp.i_seqs[0]].keys():
      #     print "ERROR!!!", bsp.i_seqs
      # STOP()

      # Proposed fast loop: Time building additional restraints for
      # ribosome went from 5272 to 204 seconds.
      for p in ij_seqs:
        if p[1] in initial_pair_asu_table_table[p[0]].keys():
          link_found = True
          break
      # VERY SLOW !!! - original loop
      # for bond_simple_proxy in geometry_proxy_registries.bond_simple.proxies:
      #   if bond_simple_proxy.i_seqs in ij_seqs:
      #     link_found = True
      #     break
      if link_found: continue
      # check for any link between atom groups based on residue name, eg ASN-NAG
      # get predefined link
      link, swap, key = linking_utils.is_atom_group_pair_linked(
        atom_group1,
        atom_group2,
        self.mon_lib_srv,
        )
      if link is None:
        link, swap, key = linking_utils.is_atom_pair_linked_tuple(atom1,
                                                                  atom2,
                                                                  classes1.important_only,
                                                                  classes2.important_only,
                                                                  self.mon_lib_srv,
                                                                  )
        if link=='TRANS':
          key=link
          link = self.mon_lib_srv.link_link_id_dict[link]
      if verbose:
        print('link',link)
        print('swap',swap)
        print('key',key)
      if swap:
        tmp = atom_group2
        atom_group2 = atom_group1
        atom_group1 = tmp
      space_group = self.special_position_settings.space_group()
      #
      if len(done_key)==2:
        link_rt_mx_ji = sgtbx.rt_mx(symbol="x,y,z", t_den=space_group.t_den())
      else:
        link_rt_mx_ji = sgtbx.rt_mx(symbol=done_key[2], t_den=space_group.t_den())
      #
      if link:
        # apply a standard link
        origin_id = origin_ids.get_origin_id('link_%s' % key,
                                             return_none_if_absent=True,
                                             )
        if verbose:
          print('apply standard link', key, origin_id)
          assert origin_id, 'origin_id for "link_%s" not found' % key
        if origin_id is None:
          # user defined links should not be applied here
          continue
        count, bond_i_seqs = _apply_link_using_proxies(
          link,
          atom_group1,
          atom_group2,
          bond_params_table,
          bond_asu_table,
          geometry_proxy_registries,
          rt_mx_ji=link_rt_mx_ji,
          origin_id=origin_id,
        )
        origin_id = None
        if len(bond_i_seqs)==0:
          if verbose:
            print('failed to link using %s' % key)
          continue
        links.setdefault(key, [])
        links[key].append([atom_group1, atom_group2])
        links[key][-1]+=bond_i_seqs[0] # odd?
        if verbose: print("predefined residue named link",key)
        continue
      #
      #if atoms_must_be:
      #  # this could be fancier...
      #  # link_residues is peptide and SG links
      #  atoms_must_be_key = [atom1.element.strip(), atom2.element.strip()]
      #  #atoms_must_be_key = [atom1.name.strip(), atom2.name.strip()]
      #  atoms_must_be_key.sort()
      #  if class_key in atoms_must_be and "S" not in atoms_must_be_key:
      #    if atoms_must_be[class_key]!=atoms_must_be_key:
      #      continue
      rc = linking_utils.process_atom_groups_for_linking_single_link(
        self.pdb_hierarchy,
        atom1,
        atom2,
        verbose=verbose,
        )
      if not rc:
        done.remove_link(done_key, names)
        continue
      pdbres, link_key, link_atoms = rc
      assert len(link_key)==1
      key = link_key[0]
      link = self.mon_lib_srv.link_link_id_dict.get(key, None)
      if verbose:
        print('pdbres',pdbres)
        print('link',link)
        print('link_key',link_key)
        print('link_atoms',link_atoms)
        for tatoms in link_atoms:
          for atom in tatoms: print(atom.quote())
      is_glyco_link = (key.find("ALPHA1")>-1 or
                       key.find("BETA1")>-1 or
                       key.find("ALPHA2")>-1 or
                       key.find("BETA2")>-1
                       )
      if is_glyco_link: # is handled in elif
        key, cif, bond_i_seqs = \
          glyco_utils.apply_glyco_link_using_proxies_and_atoms(
            atom_group2,
            atom_group1,
            bond_params_table,
            bond_asu_table,
            geometry_proxy_registries,
            rt_mx_ji=link_rt_mx_ji,
            link_carbon_dist=carbohydrate_bond_cutoff,
            origin_id=None,
          )
        links.setdefault(key, [])
        links[key].append([atom_group1, atom_group2])
        links[key][-1]+=bond_i_seqs
        continue
      elif link:
        origin_id = origin_ids['link_%s' % key]
        count, bond_i_seqs = _apply_link_using_proxies(
          link,
          atom_group1,
          atom_group2,
          bond_params_table,
          bond_asu_table,
          geometry_proxy_registries,
          rt_mx_ji=link_rt_mx_ji,
          origin_id=origin_id,
        )
        origin_id=None
        links.setdefault(key, [])
        links[key].append([atom_group1, atom_group2])
        if not bond_i_seqs:
          raise Sorry('Failed to find atoms and/or bond for "%s" in "%s" and "%s"' % (
            key,
            atom_group1.id_str(),
            atom_group2.id_str()))
        links[key][-1]+=bond_i_seqs[0]
        continue
      else:
        # possible peptide or rna/dna link
        rc = check_for_peptide_links(atom1, atom2, classes1, classes2)
        # no peptide links across symmetry
        if len(done_key)==3:
          rc = None
        if rc:
          key, swap = rc
          link = self.mon_lib_srv.link_link_id_dict.get(key)
          if swap:
            tmp = atom_group2
            atom_group2 = atom_group1
            atom_group1 = tmp
          origin_id = origin_ids['link_%s' % key]
          rc = _apply_link_using_proxies(link,
                                         atom_group1,
                                         atom_group2,
                                         bond_params_table,
                                         bond_asu_table,
                                         geometry_proxy_registries,
                                         rt_mx_ji=link_rt_mx_ji,
                                         origin_id=origin_id,
            )
          if not rc:
            tmp = atom_group2
            atom_group2 = atom_group1
            atom_group1 = tmp
            rc = _apply_link_using_proxies(link,
                                           atom_group1,
                                           atom_group2,
                                           bond_params_table,
                                           bond_asu_table,
                                           geometry_proxy_registries,
                                           rt_mx_ji=link_rt_mx_ji,
                                           origin_id=origin_id,
              )
          origin_id=None
          # not added to links so not LINK record
          if sym_op:
            sym_links += 1
            link_data.append( (atoms[i_seq].id_str(),
                               atoms[j_seq].id_str(),
                               rt_mx_ji,
                               key,
                               )
              )
          else:
            simple_links += 1
            link_data.append( (atoms[i_seq].id_str(),
                               atoms[j_seq].id_str(),
                               None, #rt_mx_ji,
                               key,
                               )
              )
          continue
      #
      custom_links.setdefault(ii, [])
      custom_links[ii].append([atom_group1, atom_group2, atom1, atom2])
      # simple
      origin_id=origin_ids['Misc. bond']
      if ((classes1.common_rna_dna or
        classes1.ccp4_mon_lib_rna_dna) and
       (classes2.common_rna_dna or classes2.ccp4_mon_lib_rna_dna)):
        bond_name = "h-dna"
        assert 0
      elif (class1=="metal" or
            class2=="metal"):
        origin_id = origin_ids['metal coordination']
      if sym_op:
        sym_bonds += 1
        bond_data.append( (atoms[i_seq].id_str(),
                           atoms[j_seq].id_str(),
                           rt_mx_ji,
                           origin_id,
            )
          )
        bond_data_i_seqs.setdefault(i_seq, [])
        bond_data_i_seqs.setdefault(j_seq, [])
        bond_data_i_seqs[i_seq].append(j_seq)
        bond_data_i_seqs[j_seq].append(i_seq)
        pair_asu_table.add_pair(proxy)
      else:
        simple_bonds += 1
        bond_data.append( (atoms[i_seq].id_str(),
                           atoms[j_seq].id_str(),
                           None, #rt_mx,
                           origin_id,
            )
          )
        bond_data_i_seqs.setdefault(i_seq, [])
        bond_data_i_seqs.setdefault(j_seq, [])
        bond_data_i_seqs[i_seq].append(j_seq)
        bond_data_i_seqs[j_seq].append(i_seq)
        pair_asu_table.add_pair(proxy.i_seqs)

    # END MAIN LOOP for ii, item in enumerate(nonbonded?)
    #
    #
    if verbose:
      for key in sorted(custom_links):
        print('-'*80)
        print(key)
        for pair in custom_links[key]:
          for atom in pair:
            try: print(atom.quote())
            except Exception: print(atom)

    pair_sym_table = pair_asu_table.extract_pair_sym_table()
    n_simple, n_symmetry = 0, 0
    self.pdb_link_records.setdefault("LINK", [])
    retain = []
    for ijk, sym_pair in enumerate(pair_sym_table.iterator()):
      i_seq, j_seq = sym_pair.i_seqs()
      origin_id = bond_data[ijk][-1]
      assert i_seq == nonbonded_i_seqs[i_seq]
      assert j_seq == nonbonded_i_seqs[j_seq]
      atom1 = atoms[i_seq]
      atom2 = atoms[j_seq]
      # check for NA linkage
      classes1 = atom_classes[i_seq]
      classes2 = atom_classes[j_seq]
      ans = bondlength_defaults.run(atom1, atom2)
      equil = 2.3
      weight = 0.02
      slack = 0.
      if len(ans) >0:
        equil = ans[0]
        if len(ans) > 1:
          weight = ans[1]
          if len(ans) > 2:
            slack = ans[2]
      if equil is None:
        equil = 2.3
      added_to_asu_table = False
      try:
        #bond_asu_table.add_pair([i_seq, j_seq])
        bond_asu_table.add_pair(
          i_seq=i_seq,
          j_seq=j_seq,
          rt_mx_ji=sym_pair.rt_mx_ji)
        added_to_asu_table = True
      except RuntimeError as e:
        error = """
    Difficulties linking atoms
      %s
      %s
    Suggestions include providing restraints for any unknown residues.
        """ % (atom1.quote(), atom2.quote())
        print(error, file=log)
      if added_to_asu_table:
        retain.append(ijk)
        if (sym_pair.rt_mx_ji.is_unit_mx()): n_simple += 1
        else:                                n_symmetry += 1
        assert origin_id
        bond_params_table.update(
          i_seq=i_seq,
          j_seq=j_seq,
          params=geometry_restraints.bond_params(
            distance_ideal=equil,
            weight=1.0/weight**2,
            slack=slack,
            origin_id=origin_id,
          ))
        # adding link to PDB
        self.pdb_link_records["LINK"].append([self.pdb_atoms[i_seq],
                                              self.pdb_atoms[j_seq],
                                              sym_pair.rt_mx_ji
                                              ])

    # output
    if link_data:
      print("  Number of additional links: simple=%d, symmetry=%d" % (
        simple_links,
        sym_bonds,
        ), file=log)
      for label1, label2, sym_op, link_name in sorted(link_data):
        if sym_op is None:
          print("    Simple link:   %s - %s" % (label1, label2), file=log)
      for label1, label2, sym_op, bond_type in sorted(bond_data):
        if sym_op:
          print("    Symmetry link: %s - %s sym. op: %s" % (label1,
                                                                    label2,
                                                                    sym_op,
            ), file=log)
    if(log is not None):
      print("  Number of custom bonds: simple=%d, symmetry=%d" % (
        n_simple, n_symmetry), file=log)
    if (n_symmetry == 0):
      blanks = ""
    else:
      blanks = "  "
    #
    def _sort_on_id(ag1, ag2):
      ag1=ag1[0]
      ag2=ag2[0]
      if ag1.id_str()[4:]==ag2.id_str()[4:]:
        if ag1.altloc<ag2.altloc:
          return -1
        return 1
      elif ag1.id_str()[4:]<ag2.id_str()[4:]:
        return -1
      return 1
    #
    if exclude_out_lines:
      print("  Excluded links - shortest distance candidate listed for each exclusion", file=log)
      for key, item in exclude_out_lines.items():
        print(item, file=log)
    if links:
      explained = []
      print("  Links applied", file=log)
      for key in sorted(links):
        print("    %s" % key, file=log)
        links[key].sort(key=cmp_to_key(_sort_on_id))
        for ag1, ag2, i_seq, j_seq in links[key]:
          self.pdb_link_records["LINK"].append([self.pdb_atoms[i_seq],
                                                self.pdb_atoms[j_seq],
                                                "x,y,z", #link_rt_mx_ji,
                                                ])
          if ag1.altloc or ag2.altloc:
            print('      "%s" - "%s" : altloc "%s" - "%s"' % (
              ag1.id_str(),
              ag2.id_str(),
              ag1.altloc,
              ag2.altloc,
              ), file=log)
          else:
            print('      "%s" - "%s"' % (ag1.id_str(),
                                                 ag2.id_str(),
              ), file=log)
          explain=""
          if key.find("ALPHA")==0 or key.find("BETA")==0:
            true_alpha_beta = glyco_utils.get_alpha_beta(ag2.resname,
                                                         fake=False,
            )
            if true_alpha_beta and key.find(true_alpha_beta.upper())==-1:
              one = "a beta"
              two = "an alpha"
              if true_alpha_beta=="alpha":
                one = "an alpha"
                two = "a beta"
              explain = "%s~> Even though %s is %s isomer," % (' '*7,
                                                               ag2.resname,
                                                               one)
              explain += " %s linkage is required..." % (two)
              if explain not in explained:
                print(explain, file=log)
              explained.append(explain)
    if bond_data:
      print("  Number of additional bonds: simple=%d, symmetry=%d" % (
        simple_bonds,
        sym_bonds,
        ), file=log)
      for caption, bond_type in [
          ("Coordination",'metal coordination'),
          ("Other bonds",'bond'),
        ]:
        print("  %s:" % caption, file=log)
        for ijk, (label1, label2, sym_op, bt) in enumerate(sorted(bond_data)):
          if sym_op is None and bt == bond_type and ijk in retain:
            print("    Simple bond:   %s - %s" % (label1, label2), file=log)
        for ijk, (label1, label2, sym_op, bt) in enumerate(sorted(bond_data)):
          if sym_op and bt == bond_type and ijk in retain:
            print("    Symmetry bond: %s - %s sym. op: %s" % (label1,
                                                                      label2,
                                                                      sym_op,
            ), file=log)
    if(log is not None):
      print('  Time building additional restraints: %0.2f' % (
        time.time()-t0), file=log)
