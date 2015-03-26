from iotbx.pdb import common_residue_names_get_class
import sys
from scitbx.array_family import flex
import math
from iotbx.pdb import get_one_letter_rna_dna_name
from libtbx.utils import Sorry
from cctbx import geometry_restraints
import iotbx.phil

dna_rna_params_str = """
hbond_distance_cutoff = 3.4
  .type = float
  .short_caption = Distance cutoff for hydrogen bonds
angle_between_bond_and_nucleobase_cutoff = 35.0
  .type = float
  .short_caption = Angle between bond and nucleobase cutoff for \
    hydrogen bonds
  .help = If angle between supposed hydrogen bond and \
    basepair plane (defined by C4, C5, C6 atoms) is less than this \
    value (in degrees), the bond will not be established.

base_pair
  .multiple = True
  .style = noauto
{
  enabled = True
    .type = bool
  base1 = None
    .type = atom_selection
  base2 = None
    .type = atom_selection
  saenger_class = 0
    .type = int
    .optional = True
    .caption = Saenger number of basepairing type, 0 if unknown
    .help = reference
  restrain_planarity = False
    .type = bool
  planarity_sigma = 0.176
    .type = float
  restrain_hbonds = True
    .type = bool
  restrain_hb_angles = True
    .type = bool
  restrain_parallelity = True
    .type = bool
  parallelity_target = 0
    .type = float
  parallelity_sigma = 0.0335
    .type = float
}
stacking_pair
  .multiple = True
  .style = noauto
{
  base1 = None
    .type = atom_selection
  base2 = None
    .type = atom_selection
  enabled = True
    .type = bool
  angle = 0
    .type = float
  sigma = 0.027
    .type = float
}
"""

def create_user_defined_NA_basepair_restraints(self, log):
  # temporarily transfered from pdb_interpretation, just in case
  # user-defined
  na_params = self.params.nucleic_acid_restraints
  sel_cache = self.pdb_hierarchy.atom_selection_cache()
  new_hbond_proxies = []
  new_hbonds = []
  max_distance = na_params.bonds.bond_distance_cutoff
  planarity_proxies = []
  parallelity_proxies = []
  if na_params.basepair_planarity.sigma < 1e-5:
    raise Sorry("Sigma for basepair planarity should be > 1e-5")
  if na_params.basepair_parallelity.sigma < 1e-5:
    raise Sorry("Sigma for basepair parallelity should be > 1e-5")
  weight = 1./na_params.basepair_planarity.sigma**2
  for bp in na_params.base_pair:
    if bp.base1 is not None and bp.base2 is not None:
      a1 = self.pdb_atoms[sel_cache.iselection(bp.base1)[0]]
      a2 = self.pdb_atoms[sel_cache.iselection(bp.base2)[0]]
      r1 = a1.parent()
      r2 = a2.parent()
      if na_params.basepair_planarity.enabled:
        seqs = self.get_i_seqs_from_na_planes(r1, r2)
        for i_seqs, j_seqs in seqs:
          planarity_proxies.append(geometry_restraints.planarity_proxy(
            i_seqs=i_seqs+j_seqs,
            weights=[weight]*len(i_seqs+j_seqs)))
      if na_params.basepair_parallelity.enabled:
        seqs = self.get_i_seqs_from_na_planes(r1, r2)
        for i_seqs, j_seqs in seqs:
          parallelity_proxies.append(geometry_restraints.parallelity_proxy(
            i_seqs=i_seqs,
            j_seqs=j_seqs,
            weight=1./na_params.basepair_parallelity.sigma**2))
      if na_params.bonds.enabled:
        bonds = self.get_h_bonds_for_basepair(a1,a2)
        for b in bonds:
          a1 = self.pdb_atoms[b[0]]
          a2 = self.pdb_atoms[b[1]]
          if a1 is not None and a2 is not None:
            dist = a1.distance(a2)
            if dist <= na_params.bonds.bond_distance_cutoff:
              new_hbonds.append(b)
              max_distance = max(max_distance, dist)
            else:
              print >> log, "    Bond between", a1.id_str(), "and", a2.id_str(),
              print >> log, "length=%4.2f" % dist, "is rejected because of",
              print >> log, "nucleic_acid_restraints.bonds.bond_distance_cutoff=",
              print >> log, na_params.bonds.bond_distance_cutoff
  return new_hbonds, planarity_proxies, parallelity_proxies, max_distance




def output_hbonds(hbond_proxies, pdb_atoms):
  # actually just to save the piece of code.
  if hbond_proxies is not None:
    dashes = open('dashes.pml', 'w')
    for p in hbond_proxies:
      s1 = pdb_atoms[p.i_seqs[0]].id_str()
      s2 = pdb_atoms[p.i_seqs[1]].id_str()
      ps = "dist chain \"%s\" and resi %s and name %s, chain \"%s\" and resi %s and name %s\n" % (
        s1[14:15], s1[15:19], s1[5:8], s2[14:15], s2[15:19], s2[5:8])
      dashes.write(ps)
    dashes.close()


def additional_check_Gendron_2001(r_i, r_j):
  # distance between rings < 5.5A
  ring_i = []
  ring_j = []
  center_i = flex.vec3_double([[0,0,0]])
  center_j = flex.vec3_double([[0,0,0]])
  for ring, res, center in [(ring_i, r_i, center_i),
                          (ring_j, r_j, center_j)]:
    for atom_name in [" N1 ", " C2 ", " N3 ", " C4 ", " C5 ", " C6 ",
                      " N9 ", " C8 ", " N7 "]:
      atom = res.find_atom_by(name=atom_name)
      if atom is not None:
        center += flex.vec3_double([atom.xyz])
        ring.append(flex.vec3_double([atom.xyz]))
    if len(ring) > 2:
      center *= 1./float(len(ring))
    else:
      return False
  d = center_j-center_i
  dn = d.norm()
  # angle between 2 normals is < 30 degrees
  r_i1 = ring_i[0]-center_i
  r_i2 = ring_i[1]-center_i
  n_i = r_i1.cross(r_i2)
  r_j1 = ring_j[0]-center_j
  r_j2 = ring_j[1]-center_j
  n_j = r_j1.cross(r_j2)
  cos_ni_nj = n_i.dot(n_j)/n_i.norm()/n_j.norm()
  angle_ni_nj_degrees = math.degrees(math.acos(abs(cos_ni_nj[0])))
  # angle between center line and normal
  cos_d_ni =  d.dot(n_i)/dn/n_i.norm()
  angle_d_ni_degrees = math.degrees(math.acos(abs(cos_d_ni[0])))
  result = (dn < 5.5 and angle_ni_nj_degrees < 30 and
    angle_d_ni_degrees < 40)
  return result

def make_phil_stacking_pair_record(residue1, residue2, params=None,
    add_segid=None, nesting_depth=1):
  segid_add = ""
  if add_segid is not None:
    segid_add = "and %s" % add_segid
  res = "%sstacking_pair {\n" % ("  "*nesting_depth)
  res += "%sbase1 = chain %s and resid %s %s\n" % ("  "*(nesting_depth+1),
      residue1.parent().parent().id, residue1.resid(), segid_add)
  res += "%sbase2 = chain %s and resid %s %s\n" % ("  "*(nesting_depth+1),
      residue2.parent().parent().id, residue2.resid(), segid_add)
  if params is not None and len(params.stacking_pair) > 0:
    master_phil = iotbx.phil.parse(dna_rna_params_str)
    actual_params = master_phil.format(params)
    w_phil = master_phil.fetch_diff(actual_params)
    w_phil_ex = w_phil.extract()
    if hasattr(w_phil_ex, 'stacking_pair'):
      for k, v in w_phil_ex.stacking_pair[0].__dict__.iteritems():
        if not k.startswith('_'):
          res += "%s%s = %s\n" % ("  "*(nesting_depth+1), k, str(v))
  res += "%s}\n" % ("  "*nesting_depth)
  return res

def make_phil_base_pair_record(residue1, residue2, params=None,
    saenger_class=None, add_segid=None, nesting_depth=1):
  segid_add = ""
  if add_segid is not None:
    segid_add = "and %s" % add_segid
  res = "%sbase_pair {\n" % ("  "*nesting_depth)
  res += "%sbase1 = chain %s and resid %s %s\n" % ("  "*(nesting_depth+1),
      residue1.parent().parent().id, residue1.parent().resid(), segid_add)
  res += "%sbase2 = chain %s and resid %s %s\n" % ("  "*(nesting_depth+1),
      residue2.parent().parent().id, residue2.parent().resid(), segid_add)
  if saenger_class is not None:
    res += "%ssaenger_class = %d\n" % ("  "*(nesting_depth+1), saenger_class)
  if params is not None and len(params.base_pair) > 0:
    master_phil = iotbx.phil.parse(dna_rna_params_str)
    actual_params = master_phil.format(params)
    w_phil = master_phil.fetch_diff(actual_params).extract()
    if hasattr(w_phil, 'base_pair'):
      for k, v in w_phil.base_pair[0].__dict__.iteritems():
        if not k.startswith('_'):
          res += "%s%s = %s\n" % ("  "*(nesting_depth+1), k, str(v))
  res += "%s}\n" % ("  "*nesting_depth)
  return res

def get_phil_stacking_pairs(pdb_hierarchy, skip_gendron_check=False,
    prefix=None, params=None, log=sys.stdout, add_segid=None):
  pairs = []
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for conformer in chain.conformers():
        if conformer.is_na():
          prev_res = None
          cur_res = None
          for i_residue, residue in enumerate(conformer.residues()):
            prev_res = cur_res
            cur_res = residue
            if prev_res is not None and cur_res is not None:
              if (not skip_gendron_check and
                  additional_check_Gendron_2001(prev_res, cur_res)):
                p = make_phil_stacking_pair_record(prev_res, cur_res, params, add_segid)
                pairs.append(p)
  phil_str = ""
  for pair_phil in pairs:
    phil_str += pair_phil
  result = ""
  if prefix is not None:
    result = "%s {\n%s}" % (prefix, phil_str)
  else:
    result = phil_str
  return result

def consecutive_residues(atom1, atom2):
  atom1_idstr = atom1.id_str()
  atom2_idstr = atom2.id_str()
  try:
    if ((atom1_idstr[14:15] == atom2_idstr[14:15]) and
        abs(int(atom1_idstr[16:19]) - int(atom2_idstr[16:19])) < 2 ):
      return True
  except ValueError:
    pass
  return False

def get_h_bonds_for_basepair(a1, a2, distance_cutoff=100):
  # a1, a2.parent().atom_group_size have to  == 1
  new_hbonds = []
  r1 = a1.parent()
  r2 = a2.parent()
  r1n = get_one_letter_rna_dna_name(r1.resname)
  r2n = get_one_letter_rna_dna_name(r2.resname)

  message_template = "Residue with name '%s' cannot be processed "
  message_template += "automatically \nfor nucleic acid basepair restraints "
  message_template += "because of non-standard residue name."
  if r1n is None:
    raise Sorry(message_template % r1.resname)
  if r2n is None:
    raise Sorry(message_template % r2.resname)
  # Translate DNA resname to RNA for unification
  # RNA
  if r1n > r2n:
    t = r1
    r1 = r2
    r2 = t
    t = r1n
    r1n = r2n
    r2n = t
  r1n = 'U' if r1n == "T" else r1n
  r2n = 'U' if r2n == "T" else r2n
  from mmtbx.monomer_library import bondlength_defaults
  best_possible_link_list = []
  best_score = 100.
  best_class_number = None
  for class_number, data in bondlength_defaults.basepairs_lengths.iteritems():
    d1 = 0
    if (r1n, r2n) == data[0]:
      for l in data[1:]:
        a1 = r1.get_atom(l[0])
        a2 = r2.get_atom(l[1])
        d1 += abs(a1.distance(a2)-2.89)
      n_links_in_data = len(data[1:])
      if n_links_in_data < 1:
        raise Sorry("Corrupted dictionary in bondlength_defaults.py")
      d1 /=n_links_in_data
      if best_score > d1:
        best_possible_link_list = [(x[:2]) for x in data[1:]]
        best_score = d1
        best_class_number = class_number
  for n1, n2 in best_possible_link_list:
    a1 = r1.get_atom(n1)
    a2 = r2.get_atom(n2)
    if a1 is not None and a2 is not None and a1.distance(a2)<distance_cutoff:
      new_hbonds.append(tuple(sorted([a1.i_seq, a2.i_seq])))
  return new_hbonds, best_class_number

def final_link_direction_check(atom1, atom2, rna_dna_angle_cutoff=35):
  import math
  a1p = atom1.parent().get_atom('C4')
  a2p = atom1.parent().get_atom('C5')
  a3p = atom1.parent().get_atom('C6')
  v1p = flex.vec3_double([(a1p.xyz)]) - flex.vec3_double([(a2p.xyz)])
  v2p = flex.vec3_double([(a2p.xyz)]) - flex.vec3_double([(a3p.xyz)])
  vn = v1p.cross(v2p)
  vl = flex.vec3_double([(atom2.xyz)]) - flex.vec3_double([(atom1.xyz)])
  cos_phi = vn.dot(vl)/vn.norm()/vl.norm()
  #print "cos_phi:", cos_phi[0], "phi:", math.acos(cos_phi[0])*360/math.pi, abs(cos_phi[0]) < 0.55
  # we have cosine between normal to plane group and link, and want this angle
  # to be around 90 degrees
  return 90 - math.degrees(math.acos(abs(cos_phi[0]))) < rna_dna_angle_cutoff

def get_phil_base_pairs(pdb_hierarchy, nonbonded_proxies,
    hbond_distance_cutoff=3.4, prefix=None, params=None,
    log=sys.stdout, add_segid=None):
  hbonds = []
  result = ""
  atoms = pdb_hierarchy.atoms()
  sites_cart = atoms.extract_xyz()
  get_sorted_result = nonbonded_proxies.get_sorted(
      by_value="delta",
      sites_cart=sites_cart)
  if get_sorted_result is None:
    return result
  sorted_nonb, n_not_shown = get_sorted_result

  # print "n_not_shown", n_not_shown
  # Get potential hbonds
  n_nonb = len(sorted_nonb)
  i = 0
  while i < n_nonb and sorted_nonb[i][3] < hbond_distance_cutoff:
    (labels, i_seq, j_seq, dist, vdw_distance, sym_op_j, rt_mx) = sorted_nonb[i]
    a1 = atoms[i_seq]
    ag1 = a1.parent()
    a2 = atoms[j_seq]
    ag2 = a2.parent()
    if (common_residue_names_get_class(ag1.resname, consider_ccp4_mon_lib_rna_dna=True) in \
          ["common_rna_dna", "ccp4_mon_lib_rna_dna"] and
        common_residue_names_get_class(ag2.resname, consider_ccp4_mon_lib_rna_dna=True) in \
          ["common_rna_dna", "ccp4_mon_lib_rna_dna"] and
        (a1.element in ["N", "O"] and a2.element in ["N", "O"]) and
        not consecutive_residues(a1, a2) and
        (ag1.altloc.strip() == ag2.altloc.strip()) and
        final_link_direction_check(a1, a2)):
      hbonds.append((i_seq, j_seq))
    i += 1
  # check and define basepairs
  pairs = []
  for hb in hbonds:
    new_hbonds, class_number = get_h_bonds_for_basepair(atoms[hb[0]],
        atoms[hb[1]], distance_cutoff=hbond_distance_cutoff)
    if len(new_hbonds) > 1:
      p = make_phil_base_pair_record(atoms[hb[0]].parent(), atoms[hb[1]].parent(),
          params, saenger_class=class_number, add_segid=add_segid)
      pairs.append(p)

  phil_str = ""
  # print "N basepairs:", len(pairs)
  for pair_phil in pairs:
    phil_str += pair_phil
  if prefix is not None:
    result = "%s {\n%s}" % (prefix, phil_str)
  else:
    result = phil_str
  return result

def get_plane_i_seqs_from_residues(r1, r2, grm):
  # Return [([i_seqA],[j_seqA]),([i_seqB],[j_seqB])] of
  # atoms in planar groups for given atom_group r. Handles up to 2 alternative
  # conformations (the lenght of resulting array will be appropriate: A and B
  # in i_seqA i_seqB - altlocs) If there is no alt confs - the list will
  # only one tuple.
  # Uses geometry_restraints_manager that should already have basic planarity
  # restraints to find the largest plane group in residue, which should be the
  # nucleobase in case of nucleic acids.
  i_seqs = []
  result = []
  r1_i_seqs = {}
  r2_i_seqs = {}
  for r, r_i_seqs in [(r1, r1_i_seqs), (r2, r2_i_seqs)]:
    for conf in r.parent().conformers():
      r_i_seqs[conf.altloc] = []
      conf_iseqs = set(conf.atoms().extract_i_seq())
      for p in grm.planarity_proxies:
        if (conf_iseqs.issuperset(p.i_seqs)
            and len(p.i_seqs) > len(r_i_seqs[conf.altloc])):
          r_i_seqs[conf.altloc] = list(p.i_seqs)
  if len(r1_i_seqs) > len(r2_i_seqs):
    t = r1_i_seqs
    r1_i_seqs = r2_i_seqs
    r2_i_seqs = t
  assert len(r1_i_seqs) > 0 and len(r1_i_seqs) > 0
  if len(r1_i_seqs) == 1:
    if len(r2_i_seqs) == 1:
      if (('' in r1_i_seqs or '' in r2_i_seqs)
          or (r1_i_seqs.keys()[0] == r2_i_seqs.keys()[0])):
        result.append((r1_i_seqs[r1_i_seqs.keys()[0]],
                       r2_i_seqs[r2_i_seqs.keys()[0]]))
    else:
      if ('' in r1_i_seqs):
        for k,v in r2_i_seqs.iteritems():
          result.append((r1_i_seqs[''], v))
  else:
    for k, v in r1_i_seqs.iteritems():
      if k in r2_i_seqs:
        result.append((v, r2_i_seqs[k]))
  # check whether sets of iseqs are different
  if len(result) > 1 and result[0] == result[1]:
    # STOP()
    return [result[0]]
  return result

def get_stacking_proxies(pdb_hierarchy, stacking_phil_params, grm):
  result = []
  assert pdb_hierarchy is not None
  if len(stacking_phil_params) < 1:
    return result
  if grm is None:
    return result
  selection_cache = pdb_hierarchy.atom_selection_cache()
  pdb_atoms = pdb_hierarchy.atoms()
  for stacking_pair in stacking_phil_params:
    if stacking_pair.base1 is not None and stacking_pair.base2 is not None:
      selected_atoms_1 = selection_cache.iselection(stacking_pair.base1)
      selected_atoms_2 = selection_cache.iselection(stacking_pair.base2)
      if len(selected_atoms_1) == 0:
        raise Sorry("Selection %s in stacking_pair retusulted in 0 atoms." % (
            stacking_pair.base1))
      if len(selected_atoms_2) == 0:
        raise Sorry("Selection %s in stacking_pair retusulted in 0 atoms." % (
            stacking_pair.base2))
      a1 = pdb_atoms[selected_atoms_1[0]]
      a2 = pdb_atoms[selected_atoms_2[0]]
      r1 = a1.parent()
      r2 = a2.parent()
      seqs = get_plane_i_seqs_from_residues(r1, r2, grm)
      for i_seqs, j_seqs in seqs:
        if len(i_seqs) > 2 and len(j_seqs) > 2:
          if stacking_pair.sigma < 1e-5:
            raise Sorry("Sigma for stacking restraints should be > 1e-5")
          proxy=geometry_restraints.parallelity_proxy(
            i_seqs=flex.size_t(i_seqs),
            j_seqs=flex.size_t(j_seqs),
            weight=1/(stacking_pair.sigma**2),
            target_angle_deg=0,
            slack=0,
            top_out=False,
            limit=1)
          result.append(proxy)
  return result

def get_basepair_plane_proxies(
    pdb_hierarchy,
    bp_phil_params,
    grm):
  assert pdb_hierarchy is not None
  result_planarities = []
  result_parallelities = []
  if len(bp_phil_params) < 1:
    return result_planarities, result_parallelities
  if grm is None:
    return result_planarities, result_parallelities
  selection_cache = pdb_hierarchy.atom_selection_cache()
  pdb_atoms = pdb_hierarchy.atoms()
  for base_pair in bp_phil_params:
    if base_pair.base1 is not None and base_pair.base2 is not None:
      selected_atoms_1 = selection_cache.iselection(base_pair.base1)
      selected_atoms_2 = selection_cache.iselection(base_pair.base2)
      if len(selected_atoms_1) == 0:
        raise Sorry("Selection %s in base_pair retusulted in 0 atoms." % (
            base_pair.base1))
      if len(selected_atoms_2) == 0:
        raise Sorry("Selection %s in base_pair retusulted in 0 atoms." % (
            base_pair.base2))
      a1 = pdb_atoms[selected_atoms_1[0]]
      a2 = pdb_atoms[selected_atoms_2[0]]
      r1 = a1.parent()
      r2 = a2.parent()
      seqs = get_plane_i_seqs_from_residues(r1, r2, grm)
      for i_seqs, j_seqs in seqs:
        if len(i_seqs) > 2 and len(j_seqs) > 2:
          if base_pair.restrain_parallelity:
            if base_pair.parallelity_sigma < 1e-5:
              raise Sorry("Sigma for parallelity basepair restraints should be > 1e-5")
            proxy=geometry_restraints.parallelity_proxy(
              i_seqs=flex.size_t(i_seqs),
              j_seqs=flex.size_t(j_seqs),
              weight=1/(base_pair.parallelity_sigma**2),
              target_angle_deg=0,
              slack=0,
              top_out=False,
              limit=1)
            result_parallelities.append(proxy)
          if base_pair.restrain_planarity:
            if base_pair.planarity_sigma < 1e-5:
              raise Sorry("Sigma for planarity basepair restraints should be > 1e-5")
            w = 1./(base_pair.planarity_sigma**2)
            proxy=geometry_restraints.planarity_proxy(
              i_seqs=flex.size_t(i_seqs+j_seqs),
              weights=[w]*len(i_seqs+j_seqs))
            result_planarities.append(proxy)
  return result_planarities, result_parallelities

def get_hb_lenght_targets(atoms):
  restraint_values = { 'N6 O4' : (3.00, 0.11),
                       'O6 N4' : (2.93, 0.10),
                       'N2 O2' : (2.78, 0.10)
  }
  anames = [atoms[0].name.strip(),
            atoms[1].name.strip()]
  rnames = [atoms[0].parent().resname,
            atoms[1].parent().resname]
  if sorted(anames) == ['N1', 'N3']:
    if get_one_letter_rna_dna_name(rnames[0]) in ['G', 'C']:
      return (2.88, 0.07) # G-C
    else:
      return (2.82, 0.08) # A-T
  else:
    key1 = "%s %s" % (anames[0], anames[1])
    key2 = "%s %s" % (anames[1], anames[0])
    vals = restraint_values.get(key1, None)
    if vals is not None:
      return vals
    vals = restraint_values.get(key2, None)
    if vals is not None:
      return vals
  return (2.91, 0.15) # values for all other bonds

def get_angle_proxies_for_bond(atoms):
  angle_values = {'O6 N4': [(122.8, 3.00), (117.3, 2.86)],
                  'N4 O6': [(117.3, 2.86), (122.8, 3.00)],
                  'N2 O2': [(122.2, 2.88), (120.7, 2.20)],
                  'O2 N2': [(120.7, 2.20), (122.2, 2.88)],
                  'N6 O4': [(115.6, 8.34), (121.2, 4.22)],
                  'O4 N6': [(121.2, 4.22), (115.6, 8.34)]}
  proxies = []
  anames = [atoms[0].name.strip(),
            atoms[1].name.strip()]
  rnames = [atoms[0].id_str().split()[1],
            atoms[1].id_str().split()[1]]
  if sorted(anames) == ['N1', 'N3']:
    if get_one_letter_rna_dna_name(rnames[0]) in ['G', 'C']:
      if anames[0] == 'N1':
        vals = [(119.1, 2.59), (116.3, 2.66)]
      else:
        vals = [(116.3, 2.66), (119.1, 2.59)]
    else:
      if anames[0] == 'N1':
        vals = [(116.2, 3.46), (115.8, 2.88)]
      else:
        vals = [(115.8, 2.88), (116.2, 3.46)]
  else:
    key = "%s %s" % (anames[0], anames[1])
    vals = angle_values.get(key, None)
  if vals is not None:
    for i in range(2):
      atoms_for_angle = [None, None, None]
      aname = anames[i]
      if (aname == 'N1' or aname == 'N2' or aname == 'N3' or aname == 'O2'):
        atoms_for_angle[0] = atoms[i].parent().get_atom('C2')
      elif (aname == 'N4' or aname == 'O4'):
        atoms_for_angle[0] = atoms[i].parent().get_atom('C4')
      elif (aname == 'N6' or aname == 'O6'):
        atoms_for_angle[0] = atoms[i].parent().get_atom('C6')
      if atoms_for_angle[0] is not None:
        atoms_for_angle[1] = atoms[i]
        atoms_for_angle[2] = atoms[1-i]
      if atoms_for_angle.count(None) == 0:
        i_seqs_for_angle = [x.i_seq for x in atoms_for_angle]
        p = geometry_restraints.angle_proxy(
          i_seqs=i_seqs_for_angle,
          angle_ideal=vals[i][0],
          weight=1./vals[i][1]**2)
      proxies.append(p)
  return proxies

def get_h_bonds_for_particular_basepair(atoms, saenger_class=0):
  result = []
  if saenger_class == 0:
    return result
  assert len(atoms) == 2
  new_hbonds = []
  r1 = atoms[0].parent()
  r2 = atoms[1].parent()
  r1n = get_one_letter_rna_dna_name(r1.resname)
  r2n = get_one_letter_rna_dna_name(r2.resname)

  message_template = "Residue with name '%s' cannot be processed "
  message_template += "automatically \nfor nucleic acid basepair restraints "
  message_template += "because of non-standard residue name."
  if r1n is None:
    raise Sorry(message_template % r1.resname)
  if r2n is None:
    raise Sorry(message_template % r2.resname)
  # Translate DNA resname to RNA for unification
  # RNA
  if r1n > r2n:
    t = r1
    r1 = r2
    r2 = t
    t = r1n
    r1n = r2n
    r2n = t
  r1n = 'U' if r1n == "T" else r1n
  r2n = 'U' if r2n == "T" else r2n
  from mmtbx.monomer_library import bondlength_defaults
  if bondlength_defaults.basepairs_lengths[saenger_class][0] != (r1n, r2n):
    raise Sorry("Saenger class does not match residue names")
  hbonds = []
  for b in bondlength_defaults.basepairs_lengths[saenger_class][1:]:
    hba1 = r1.get_atom(b[0])
    hba2 = r2.get_atom(b[1])
    hbonds.append((hba1, hba2))
  return hbonds

def get_basepair_hbond_proxies(
    pdb_hierarchy,
    bp_phil_params):
  assert pdb_hierarchy is not None
  bond_proxies_result = []
  angle_proxies_result = []
  if len(bp_phil_params) < 1:
    return bond_proxies_result, angle_proxies_result
  selection_cache = pdb_hierarchy.atom_selection_cache()
  pdb_atoms = pdb_hierarchy.atoms()
  # dashes = open('dashes.pml', 'w')
  pdb_atoms = pdb_hierarchy.atoms()
  for base_pair in bp_phil_params:
    if base_pair.base1 is not None and base_pair.base2 is not None:
      selected_atoms_1 = selection_cache.iselection(base_pair.base1)
      selected_atoms_2 = selection_cache.iselection(base_pair.base2)
      if len(selected_atoms_1) == 0:
        raise Sorry("Selection %s in base_pair retusulted in 0 atoms." % (
            base_pair.base1))
      if len(selected_atoms_2) == 0:
        raise Sorry("Selection %s in base_pair retusulted in 0 atoms." % (
            base_pair.base2))
      a1 = pdb_atoms[selected_atoms_1[0]]
      a2 = pdb_atoms[selected_atoms_2[0]]
      hbonds = get_h_bonds_for_particular_basepair((a1, a2), base_pair.saenger_class)
      for hb in hbonds:
        if base_pair.restrain_hbonds:
          hb_target, hb_sigma = get_hb_lenght_targets(hb)
          p = geometry_restraints.bond_simple_proxy(
            i_seqs=[hb[0].i_seq, hb[1].i_seq],
            distance_ideal=hb_target,
            weight=1.0/hb_sigma**2,
            slack=0,
            top_out=False,
            limit=1)
          bond_proxies_result.append(p)
          # print "bond:", hb[0].id_str(), hb[1].id_str(), "(%4.2f, %4.2f)" % (hb_target, hb_sigma)
        # s1 = pdb_atoms[hb[0].i_seq].id_str()
        # s2 = pdb_atoms[hb[1].i_seq].id_str()
        # ps = "dist chain \"%s\" and resi %s and name %s, chain \"%s\" and resi %s and name %s\n" % (
        #   s1[14:15], s1[15:19], s1[5:8], s2[14:15], s2[15:19], s2[5:8])
        # dashes.write(ps)
        if base_pair.restrain_hb_angles:
          angle_proxies_result += get_angle_proxies_for_bond(hb)
  # dashes.close()
  return bond_proxies_result, angle_proxies_result
