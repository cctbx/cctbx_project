from __future__ import division
import time
from cctbx import sgtbx
from cctbx import geometry_restraints

from mmtbx.monomer_library import linking_setup
from mmtbx.monomer_library import linking_utils
from mmtbx.monomer_library import glyco_utils
import bondlength_defaults

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
                              ):
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

    # what what what
    # print i_seqs

    value = "value_dist"
    proxy = geometry_restraints.bond_simple_proxy(
      i_seqs=i_seqs,
      distance_ideal=getattr(bond, value),
      weight=1/bond.value_dist_esd**2)
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
      weight=1/angle.value_angle_esd**2)
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
    proxy = geometry_restraints.chirality_proxy(
      i_seqs=i_seqs,
      volume_ideal=volume_ideal,
      both_signs=False,
      weight=25.,
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
        )
      geometry_proxy_registries.planarity.add_if_not_duplicated(proxy=proxy)
  return count

def possible_cyclic_peptide(atom1,
                            atom2,
                            ):
  if 0:
    print atom1.quote(),atom2.quote()
  chain1 = atom1.parent().parent().parent()
  chain2 = atom1.parent().parent().parent()
  if not chain1.id == chain2.id: return False
  fl = {}
  rgs = chain1.residue_groups()
  for i in range(0,-2,-1):
    for atom in rgs[i].atoms():
      if atom.quote()==atom1.quote():
        fl[i]=atom1
        break
      elif atom.quote()==atom2.quote():
        fl[i]=atom2
        break
  return len(fl)==2

def check_for_peptide_links(atom1,
                            atom2,
                            ):
  classes1 = linking_utils.get_classes(atom1)
  classes2 = linking_utils.get_classes(atom2)
  atom_group1 = atom1.parent()
  atom_group2 = atom2.parent()
  if classes1.common_amino_acid or classes2.common_amino_acid:
    if classes1.common_amino_acid:
      other = atom_group2
    else:
      other = atom_group1
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

class linking_mixins(object):
  def process_intra_chain_links(self,
                                model,
                                mon_lib_srv,
                                log,
                                residue_group_cutoff2=400.,
                                #bond_cutoff=2.75,
                                amino_acid_bond_cutoff=1.9,
                                rna_dna_bond_cutoff=3.5,
                                intra_residue_bond_cutoff=1.99,
                                verbose=False,
                                ):
    assert 0
    t0 = time.time()
    ########################################
    # must be after process_apply_cif_link #
    ########################################
    import linking_utils
    #
    def generate_first_atom_of_residue_groups(model, chain_id, verbose=False):
      assert 0
      for i_chain, chain in enumerate(model.chains()):
        if chain.id!=chain_id: continue
        for i_residue_group, residue_group in enumerate(chain.residue_groups()):
          if verbose: print '  residue_group: resseq="%s" icode="%s"' % (
            residue_group.resseq, residue_group.icode)
          yield residue_group.atoms()[0]
    #
    chain_ids = []
    outls = {}
    for chain in model.chains():
      if chain.id not in chain_ids: chain_ids.append(chain.id)
    for chain_id in chain_ids:
      for i_atom, atom1 in enumerate(
        generate_first_atom_of_residue_groups(model, chain_id, verbose=verbose)
        ):
        classes1 = linking_utils.get_classes(atom1)
        for j_atom, atom2 in enumerate(
          generate_first_atom_of_residue_groups(model, chain_id, verbose=False)
          ):
          if i_atom>=j_atom: continue
          classes2 = linking_utils.get_classes(atom2)
          # what about gamma linking?
          if classes1.common_water: continue
          if classes2.common_water: continue
          if classes1.common_amino_acid and classes2.common_amino_acid: continue
          if classes1.common_rna_dna and classes2.common_rna_dna: continue
          d2 = linking_utils.get_distance2(atom1, atom2)
          if d2>residue_group_cutoff2: continue
          if verbose:
            print atom1.quote(), atom2.quote(), d2,residue_group_cutoff2
          #
          rc = linking_utils.process_atom_groups_for_linking(
            self.pdb_hierarchy,
            atom1,
            atom2,
            classes1,
            classes2,
            #bond_cutoff=bond_cutoff,
            amino_acid_bond_cutoff=amino_acid_bond_cutoff,
            rna_dna_bond_cutoff=rna_dna_bond_cutoff,
            intra_residue_bond_cutoff=intra_residue_bond_cutoff,
            verbose=verbose,
            )
          if rc is None: continue
          pdbres_pairs, data_links, atomss = rc
          for pdbres_pair, data_link, atoms in zip(pdbres_pairs,
                                                   data_links,
                                                   atomss,
                                                   ):
            if verbose:
              print '-'*80
              print pdbres_pairs
              print "data_link",data_link
              for atom in atoms:
                print atom.quote()
            outls.setdefault(data_link, "")
            tmp = self.process_custom_links(mon_lib_srv,
                                            pdbres_pair,
                                            data_link,
                                            atoms,
                                            verbose=verbose,
                                            )
            if verbose:
              print "rc :%s:" % tmp
            outls[data_link] = tmp
    # log output about detection
    remove=[]
    if self.apply_cif_links:
      outl = ""
      for i, apply in enumerate(self.apply_cif_links):
        if(not getattr(apply, "automatic", False)): continue
        if apply.data_link in mon_lib_srv.link_link_id_dict:
          outl += "%sLinking %s to %s using %s\n" % (
            " "*8,
            apply.pdbres_pair[0][7:],
            apply.pdbres_pair[1][7:],
            apply.data_link,
            )
        else:
          if getattr(apply, "possible_peptide_link", False):
            if 1:
              pass
            else:
              print >> log, "%sPossible peptide link %s to %s" % (
                " "*8,
                apply.pdbres_pair[0][7:],
                apply.pdbres_pair[1][7:],
                )
          elif getattr(apply, "possible_rna_dna_link", False):
            pass
          else:
            outl += "%sLinking %s to %s\n" % (
              " "*8,
              apply.pdbres_pair[0][7:],
              apply.pdbres_pair[1][7:],
              )
            outl += '%sCreating link for "%s"\n' % (" "*10, apply.data_link)
            mon_lib_srv.link_link_id_dict[apply.data_link] = None
        if outls.get(apply.data_link, ""):
          outl += outls[apply.data_link]
      if outl:
        print >> log, "%sAdding automatically detected intra-chain links" % (
          " "*6,
          )
        print >> log, outl[:-1]
    if remove:
      remove.sort()
      remove.reverse()
      for r in remove: del self.apply_cif_links[r]
    print >> log, "%sTime to detect intra-chain links : %0.1fs" % (
      " "*6,
      time.time()-t0)

  # def process_nonbonded_for_linking(pdb_inp,
  #                                   pdb_hierarchy,
  #                                   geometry_restaints_manager,
  #                                   verbose=False,
  #                                   ):
  #   assert 0
  #   sorted_nonbonded_proxies = get_nonbonded(pdb_inp,
  #                                            pdb_hierarchy,
  #                                            geometry_restaints_manager,
  #     )
  #   atoms = pdb_hierarchy.atoms()
  #   result = []
  #   for item in sorted_nonbonded_proxies:
  #     labels, i_seq, j_seq, distance, vdw_distance, sym_op, rt_mx_ji = item
  #     item = empty()
  #     item.labels = labels
  #     item.i_seq = i_seq
  #     item.j_seq = j_seq
  #     item.distance = distance
  #     if item.distance>2.75: break
  #     item.sym_op = sym_op
  #     item.rt_mx_ji = rt_mx_ji
  #     atom1 = atoms[i_seq]
  #     atom2 = atoms[j_seq]
  #     if verbose:
  #       print " Nonbonded: %s %s %0.3f %s %s" % (atoms[item.i_seq].id_str(),
  #                                                atoms[item.j_seq].id_str(),
  #                                                item.distance,
  #                                                item.sym_op,
  #                                                item.rt_mx_ji,
  #                                                ),
  #
  #     if is_atom_pair_linked(atom1, atom2):
  #       print " Linking?"
  #       result.append(item)
  #     else:
  #       print
  #   return result

  def process_nonbonded_for_links(self,
                                  bond_params_table,
                                  bond_asu_table,
                                  geometry_proxy_registries,
                                  na_params,
                                  link_metals                 = True,
                                  link_residues               = True,
                                  link_carbohydrates          = True,
                                  max_bonded_cutoff           = None,
                                  metal_coordination_cutoff   = 3.,
                                  amino_acid_bond_cutoff      = 2.,
                                  inter_residue_bond_cutoff   = 2.,
                                  second_row_buffer           = 0.5,
                                  carbohydrate_bond_cutoff    = 2.,
                                  log                         = None,
                                  verbose                     = False,
                                  ):
    assert hasattr(self, "cif")
    if max_bonded_cutoff is None:
      max_bonded_cutoff = max(metal_coordination_cutoff,
                              amino_acid_bond_cutoff,
                              na_params.bonds.bond_distance_cutoff,
                              carbohydrate_bond_cutoff,
                              inter_residue_bond_cutoff+second_row_buffer,
                              )
    hbonds_in_bond_list = []
    if verbose:
      print """
      metal_coordination_cutoff %s
      amino_acid_bond_cutoff    %s
      rna_dna_bond_cutoff       %s
      carbohydrate_bond_cutoff  %s
      inter_residue_bond_cutoff %s
      second_row_buffer         %s
      """ % ( metal_coordination_cutoff,
              amino_acid_bond_cutoff,
              na_params.bonds.bond_distance_cutoff,
              carbohydrate_bond_cutoff,
              inter_residue_bond_cutoff,
              second_row_buffer,
              )
    from cctbx import crystal
    from cctbx.array_family import flex
    #
    # def _nonbonded_pair_generator_from_pair_asu_table(max_bonded_cutoff=3.):
    #   assert 0
    #   xray_structure_simple=self.pdb_inp.xray_structure_simple()
    #   pair_asu_table = xray_structure_simple.pair_asu_table(
    #     distance_cutoff=10, #max_bonded_cutoff,
    #     ) #self.clash_threshold)
    #   bonded_i_seqs = []
    #   for bp in self.geometry_proxy_registries.bond_simple.proxies:
    #     bonded_i_seqs.append(bp.i_seqs)
    #   pair_sym_table = pair_asu_table.extract_pair_sym_table()
    #   atom_pairs_i_seqs, sym_atom_pairs_i_seqs = pair_sym_table.both_edge_list()
    #   nonbonded_pairs = list(set(atom_pairs_i_seqs).difference(set(bonded_i_seqs)))
    #   for ii, (i_seq, j_seq) in enumerate(nonbonded_pairs):
    #     yield (i_seq, j_seq, None)
    #     if ii>5: break
    #   nonbonded_pairs = list(set(sym_atom_pairs_i_seqs).difference(set(bonded_i_seqs)))
    #   for ii, (i_seq, j_seq, sym) in enumerate(nonbonded_pairs):
    #     yield (i_seq, j_seq, sym)
    #     if ii>5: break
    # #
    def _nonbonded_pair_objects(max_bonded_cutoff=3.,
                                i_seqs=None,
                                ):
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
    def _nonbonded_pair_generator_geometry_restraints_sort(
        nonbonded_proxies,
        max_bonded_cutoff=3.):
      rc = nonbonded_proxies.get_sorted(by_value="delta",
                                        sites_cart=sites_cart,
                                        include_proxy=True,
        )
      if rc is None: return
      rc, junk = rc
      for item in rc:
        yield item
    #
    # def _nonbonded_proxies_generator_geometry_restraints_sort(
    #     nonbonded_proxies,
    #     max_bonded_cutoff=3.,
    #     ):
    #   assert 0
    #   rc = nonbonded_proxies.get_sorted_proxies(by_value="delta",
    #                                             sites_cart=sites_cart,
    #     )
    #   for item in rc:
    #     yield item
    #
    print >> log, """
  Automatic linking
    Parameters for automatic linking
      Linking & cutoffs
        Metal        : %-5s - %0.2f
        Amimo acid   : %-5s - %0.2f
        RNA/DNA      : %-5s - %0.2f
        Carbohydrate : %-5s - %0.2f
      """ % (link_metals,
             metal_coordination_cutoff,
             link_residues,
             amino_acid_bond_cutoff,
             (na_params.enabled and na_params.bonds.enabled),
             na_params.bonds.bond_distance_cutoff,
             link_carbohydrates,
             carbohydrate_bond_cutoff,
             )
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

    # main loop
    nonbonded_proxies, sites_cart, pair_asu_table, asu_mappings, nonbonded_i_seqs = \
        _nonbonded_pair_objects(max_bonded_cutoff=max_bonded_cutoff,
          )
    for ii, item in enumerate(
        _nonbonded_pair_generator_geometry_restraints_sort(
          nonbonded_proxies=nonbonded_proxies,
          max_bonded_cutoff=max_bonded_cutoff,
          )
        ):
      labels, i_seq, j_seq, distance, vdw_distance, sym_op, rt_mx_ji, proxy = item
      # XXX this is a poor job!!!
      if bond_asu_table.contains(i_seq, j_seq, 0): continue
      if bond_asu_table.contains(i_seq, j_seq, 1): continue
      atom1 = atoms[i_seq]
      atom2 = atoms[j_seq]
      if verbose:
        print i_seq, j_seq, atom1.quote(),
        print atom2.quote(),
        print "Distance: %0.2f" % distance, rt_mx_ji, sym_op

      # don't link atoms not in the same conformer (works for models also)...
      if not atom1.is_in_same_conformer_as(atom2):
        assert 0
        continue
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
      classes1 = linking_utils.get_classes(atom1)
      classes2 = linking_utils.get_classes(atom2)
      use_only_bond_cutoff = False
      if verbose:
        print """
Residue classes

%s
%s """ % (classes1, classes2)
      # is_proxy_set between any of the atoms ????????
      if classes1.common_amino_acid and classes2.common_amino_acid:
        if not link_residues:
          continue
        # special amino acid linking
        #  - cyclic
        #  - beta, delta ???
        if possible_cyclic_peptide(atom1, atom2): # first & last peptide
          use_only_bond_cutoff = True
      if classes1.common_rna_dna and classes2.common_rna_dna:
        if not na_params.enabled:
          continue
      #else:
      #  continue # Hard hookup to disable all but DNA/RNA basepair linking
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
        key.append(str(rt_mx_ji))
      key = tuple(key)
      # bond length cutoff & some logic
      if not linking_utils.is_atom_pair_linked(
          atom1,
          atom2,
          distance=distance,
          max_bonded_cutoff=max_bonded_cutoff,
          amino_acid_bond_cutoff=amino_acid_bond_cutoff,
          rna_dna_bond_cutoff=na_params.bonds.bond_distance_cutoff,
          rna_dna_angle_cutoff=na_params.bonds.\
              angle_between_bond_and_nucleobase_cutoff,
          inter_residue_bond_cutoff=inter_residue_bond_cutoff,
          second_row_buffer=second_row_buffer,
          saccharide_bond_cutoff=carbohydrate_bond_cutoff,
          metal_coordination_cutoff=metal_coordination_cutoff,
          use_only_bond_cutoff=use_only_bond_cutoff,
          verbose=verbose,
          ):
        if verbose:
          print "is not linked", atom1.quote(),atom2.quote(),key
        continue
      # got a link....

      if classes1.common_rna_dna and classes2.common_rna_dna:
        hbonds_in_bond_list.append(tuple(sorted([atom1.i_seq, atom2.i_seq])))

      class1 = linking_utils.get_classes(atom1, #_group1.resname,
                                         important_only=True,
        )
      class2 = linking_utils.get_classes(atom2, #_group2.resname,
                                         important_only=True,
        )
      class_key = [class1, class2]
      class_key.sort()
      class_key = tuple(class_key)
      #
      if not link_metals and "metal" in class_key: continue
      if (not na_params.enabled and na_params.bonds.enabled
          and "common_rna_dna" in class_key): continue
      if not link_residues and "common_amino_acid" in class_key: continue
      if not link_carbohydrates and "common_saccharide" in class_key: continue
      #
      names = [atom1.name, atom2.name]
      names.sort()
      # exclude duplicate symmetry op.
      if key in done:
        if names in done[key]:
          continue
      #
      current_number_of_links = len(done.setdefault(key, []))
      if(current_number_of_links >=
         linking_setup.maximum_inter_residue_links.get(class_key, 1)
         ):
        if verbose:
          print "too many links:",current_number_of_links,linking_setup.maximum_inter_residue_links.get(class_key, 1), class_key
        continue
      #
      done[key].append(names)
      done_key = key
      # check for any link between atom groups based on residue name, eg ASN-NAG
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

      # !!! For every possible link we are looping over _all_ bond proxies?
      # Impossible for DNA/RNA links due to enormous runtime (hundreds links).
      link_found = False
      if not (classes1.common_rna_dna and classes2.common_rna_dna):
        for bond_simple_proxy in geometry_proxy_registries.bond_simple.proxies:
          if bond_simple_proxy.i_seqs in ij_seqs:
            link_found = True
            break
      if link_found: continue
      # get predefined link
      link, swap, key = linking_utils.is_atom_group_pair_linked(
        atom_group1,
        atom_group2,
        self.mon_lib_srv,
        )
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
        links.setdefault(key, [])
        links[key].append([atom_group1, atom_group2])
        _apply_link_using_proxies(link,
                                  atom_group1,
                                  atom_group2,
                                  bond_params_table,
                                  bond_asu_table,
                                  geometry_proxy_registries,
#                                  distance=distance,
                                  rt_mx_ji=link_rt_mx_ji,
                                  )
        if verbose: print "predefined residue named link",key
        self.cif["link_%s" % key] = link.cif_object
        continue
      #
      rc = linking_utils.process_atom_groups_for_linking_single_link(
        self.pdb_hierarchy,
        atom1,
        atom2,
        )
      if not rc:
        done.remove_link(done_key, names)
        continue
      pdbres, link_key, link_atoms = rc
      assert len(link_key)==1
      key = link_key[0]
      link = self.mon_lib_srv.link_link_id_dict.get(key, None)
      if key.find("ALPHA")>-1 or key.find("BETA")>-1:
        key, cif = glyco_utils.apply_glyco_link_using_proxies_and_atoms(
          atom_group2,
          atom_group1,
          bond_params_table,
          bond_asu_table,
          geometry_proxy_registries,
          rt_mx_ji=link_rt_mx_ji,
          )
        links.setdefault(key, [])
        links[key].append([atom_group1, atom_group2])
        self.cif["link_%s" % key] = cif
        continue
      elif link:
        links.setdefault(key, [])
        links[key].append([atom_group1, atom_group2])
        _apply_link_using_proxies(link,
                                  atom_group1,
                                  atom_group2,
                                  bond_params_table,
                                  bond_asu_table,
                                  geometry_proxy_registries,
                                  rt_mx_ji=link_rt_mx_ji,
                                  )
        continue
      else:
        # possible peptide or rna/dna link
        rc = check_for_peptide_links(atom1, atom2)
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
          rc = _apply_link_using_proxies(link,
                                         atom_group1,
                                         atom_group2,
                                         bond_params_table,
                                         bond_asu_table,
                                         geometry_proxy_registries,
                                         rt_mx_ji=link_rt_mx_ji,
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
              )
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
      if sym_op:
        sym_bonds += 1
        bond_data.append( (atoms[i_seq].id_str(),
                           atoms[j_seq].id_str(),
                           rt_mx_ji,
                           "h-dna" if (classes1.common_rna_dna and
                             classes2.common_rna_dna) else "bond",
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
                           "h-dna" if (classes1.common_rna_dna and
                             classes2.common_rna_dna) else "bond",
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
        print '-'*80
        print key
        for pair in custom_links[key]:
          for atom in pair:
            try: print atom.quote()
            except Exception: print atom

    pair_sym_table = pair_asu_table.extract_pair_sym_table()
    #for sym_pair in pair_sym_table.iterator():
    #  print sym_pair.i_seq, sym_pair.j_seq,
    #  print atoms[sym_pair.i_seq].quote(), atoms[sym_pair.j_seq].quote()
    n_simple, n_symmetry = 0, 0

    for sym_pair in pair_sym_table.iterator():
      i_seq, j_seq = sym_pair.i_seqs()
      assert i_seq == nonbonded_i_seqs[i_seq]
      assert j_seq == nonbonded_i_seqs[j_seq]
      atom1 = atoms[i_seq]
      atom2 = atoms[j_seq]
      if (sym_pair.rt_mx_ji.is_unit_mx()): n_simple += 1
      else:                                n_symmetry += 1
      # check for NA linkage
      classes1 = linking_utils.get_classes(atom1)
      classes2 = linking_utils.get_classes(atom2)
      if (classes1.common_rna_dna and classes2.common_rna_dna
          and na_params.enabled and na_params.bonds.enabled):
        ans = [na_params.bonds.target_value,
               na_params.bonds.sigma, na_params.bonds.slack]
      else:
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
      bond_params_table.update(
        i_seq=i_seq,
        j_seq=j_seq,
        params=geometry_restraints.bond_params(
          distance_ideal=equil,
          weight=1.0/weight**2,
          slack=slack))
      bond_asu_table.add_pair(
        i_seq=i_seq,
        j_seq=j_seq,
        rt_mx_ji=sym_pair.rt_mx_ji)
    # output
    if link_data:
      print >> log, "  Number of additional links: simple=%d, symmetry=%d" % (
        simple_links,
        sym_bonds,
        )
      for label1, label2, sym_op, link_name in sorted(link_data):
        if sym_op is None:
          print >> log, "    Simple link:   %s - %s" % (label1, label2)
      for label1, label2, sym_op, bond_type in sorted(bond_data):
        if sym_op:
          print >> log, "    Symmetry link: %s - %s sym. op: %s" % (label1,
                                                                    label2,
                                                                    sym_op,
            )
    print >> log,"  Number of custom bonds: simple=%d, symmetry=%d" % (
      n_simple, n_symmetry)
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
    if links:
      print >> log,  "  Links applied"
      for key in sorted(links):
        print >> log,  "    %s" % key
        links[key].sort(_sort_on_id)
        for ag1, ag2 in links[key]:
          if ag1.altloc or ag2.altloc:
            print >> log, '      "%s" - "%s" : altloc "%s" - "%s"' % (
              ag1.id_str(),
              ag2.id_str(),
              ag1.altloc,
              ag2.altloc,
              )
          else:
            print >> log, '      "%s" - "%s"' % (ag1.id_str(),
                                                 ag2.id_str(),
              )
    if bond_data:
      print >> log, "  Number of additional bonds: simple=%d, symmetry=%d" % (
        simple_bonds,
        sym_bonds,
        )
      for caption, bond_type in [("Nucleic acid basepair bonds",'h-dna'),
                                 ("Other bonds",'bond')]:
        print >> log, "  %s:" % caption
        for label1, label2, sym_op, bt in sorted(bond_data):
          if sym_op is None and bt == bond_type:
            print >> log, "    Simple bond:   %s - %s" % (label1, label2)
        for label1, label2, sym_op, bt in sorted(bond_data):
          if sym_op and bt == bond_type:
            print >> log, "    Symmetry bond: %s - %s sym. op: %s" % (label1,
                                                                      label2,
                                                                      sym_op,
            )
    print >> log, '  Time building additional restraints: %0.2f' % (time.time()-t0)
    return hbonds_in_bond_list
