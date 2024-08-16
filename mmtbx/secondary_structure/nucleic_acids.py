from __future__ import absolute_import, division, print_function
from iotbx.pdb import common_residue_names_get_class
import sys
from scitbx.array_family import flex
import math
from iotbx.pdb import get_one_letter_rna_dna_name
from libtbx.utils import Sorry
from cctbx import geometry_restraints
from mmtbx.monomer_library import bondlength_defaults
import iotbx.phil

import six
from six.moves import range

origin_ids = geometry_restraints.linking_class.linking_class()

dna_rna_params_str = """
hbond_distance_cutoff = 3.4
  .type = float
  .short_caption = Distance cutoff for hydrogen bonds
  .help = Hydrogen bonds with length exceeding this limit will not be \
    established
scale_bonds_sigma = 1.
  .type = float
  .short_caption = Scale h-bond sigma
  .help = All sigmas for h-bond length will be multiplied \
    by this number. The smaller number is tighter restraints.
  .expert_level = 3

base_pair
  .multiple = True
  .style = noauto
{
  enabled = True
    .type = bool
    .help = Restraint this particular base-pair
  base1 = None
    .type = atom_selection
    .help = Selection string selecting at least one atom in the desired residue
  base2 = None
    .type = atom_selection
    .help = Selection string selecting at least one atom in the desired residue
  saenger_class = 0
    .type = int
    .optional = True
    .caption = Saenger number of basepairing type, 0 if unknown
    .help = Type of base-pairing
  restrain_planarity = False
    .type = bool
    .help = Apply planarity restraint to this base-pair
  planarity_sigma = 0.176
    .type = float
  restrain_hbonds = True
    .type = bool
    .help = Restrain hydrogen bonds length for this base-pair
  restrain_hb_angles = True
    .type = bool
    .help = Restrain angles around hydrogen bonds for this base-pair
  restrain_parallelity = True
    .type = bool
    .help = Apply parallelity restraint to this base-pair
  parallelity_target = 0
    .type = float
  parallelity_sigma = 0.0335
    .type = float
}
stacking_pair
  .multiple = True
  .style = noauto
{
  enabled = True
    .type = bool
    .help = Restraint this particular base-pair
  base1 = None
    .type = atom_selection
    .help = Selection string selecting at least one atom in the desired residue
  base2 = None
    .type = atom_selection
    .help = Selection string selecting at least one atom in the desired residue
  angle = 0
    .type = float
  sigma = 0.027
    .type = float
}
"""

def output_hbonds(hbond_proxies, pdb_atoms):
  # actually just to save the piece of code.
  if hbond_proxies is not None:
    dashes = open('dashes.pml', 'w')
    for p in hbond_proxies:
      awl1 = pdb_atoms[p.i_seqs[0]].fetch_labels()
      awl2 = pdb_atoms[p.i_seqs[1]].fetch_labels()
      ps = "dist chain \"%s\" and resi %s and name %s, chain \"%s\" and resi %s and name %s\n" % (
        awl1.chain_id, awl1.resseq, awl1.name, awl2.chain_id, awl2.resseq, awl2.name)
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

def format_base_string(base_str, residue, segid=None):
  segid_add = "and segid '%s'" % segid if segid is not None else ""
  resid = "%s" % (
      residue.resid() if hasattr(residue, "resid") else residue.parent().resid())
  chain_add = "chain '%s' and " % residue.parent().parent().id
  base = "%s = %s resid %s %s\n" % (base_str, chain_add,
      resid, segid_add)
  return base

def make_phil_stacking_pair_record(residue1, residue2, params=None,
    add_segid=None, nesting_depth=1):
  res = "%sstacking_pair {\n" % ("  "*nesting_depth)
  res += "%s%s" % ("  "*(nesting_depth+1), format_base_string(
      "base1", residue1, add_segid))
  res += "%s%s" % ("  "*(nesting_depth+1), format_base_string(
      "base2", residue2, add_segid))
  # add non-defaults!
  if params is not None and len(params.stacking_pair) > 0:
    master_phil = iotbx.phil.parse(dna_rna_params_str)
    actual_params = master_phil.format(params)
    w_phil = master_phil.fetch_diff(actual_params).extract()
    if hasattr(w_phil, 'stacking_pair'):
      for k, v in six.iteritems(w_phil.stacking_pair[0].__dict__):
        if not k.startswith('_'):
          res += "%s%s = %s\n" % ("  "*(nesting_depth+1), k, str(v))
  res += "%s}\n" % ("  "*nesting_depth)
  return res

def make_phil_base_pair_record(residue1, residue2, params=None,
    saenger_class=None, add_segid=None, nesting_depth=1):
  res = "%sbase_pair {\n" % ("  "*nesting_depth)
  res += "%s%s" % ("  "*(nesting_depth+1), format_base_string(
      "base1", residue1, add_segid))
  res += "%s%s" % ("  "*(nesting_depth+1), format_base_string(
      "base2", residue2, add_segid))
  if saenger_class is not None:
    res += "%ssaenger_class = %d\n" % ("  "*(nesting_depth+1), saenger_class)
  # add non-defaults!
  if params is not None and len(params.base_pair) > 0:
    master_phil = iotbx.phil.parse(dna_rna_params_str)
    actual_params = master_phil.format(params)
    w_phil = master_phil.fetch_diff(actual_params).extract()
    if hasattr(w_phil, 'base_pair'):
      for k, v in six.iteritems(w_phil.base_pair[0].__dict__):
        if not k.startswith('_'):
          res += "%s%s = %s\n" % ("  "*(nesting_depth+1), k, str(v))
  res += "%s}\n" % ("  "*nesting_depth)
  return res

def get_phil_stacking_pairs(pdb_hierarchy, skip_gendron_check=False,
    prefix=None, params=None, log=sys.stdout, add_segid=None):
  pairs = []
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      if chain.is_na():
        for conformer in chain.conformers():
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
  awl1 = atom1.fetch_labels()
  awl2 = atom2.fetch_labels()
  if ((awl1.chain_id == awl2.chain_id) and
      abs(awl1.resseq_as_int() - awl2.resseq_as_int()) < 2 ):
    return True
  return False

def unify_residue_names_and_order(r1, r2):
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
  return r1, r2, r1n, r2n

def get_h_bonds_for_basepair(a1, a2, distance_cutoff=100, log=sys.stdout, verbose=-1):
  # a1, a2.parent().atom_group_size have to  == 1
  new_hbonds = []
  r1 = a1.parent()
  r2 = a2.parent()
  r1, r2, r1n, r2n = unify_residue_names_and_order(a1.parent(), a2.parent())
  best_possible_link_list = []
  best_score = 100.
  best_class_number = None
  for class_number, data in six.iteritems(bondlength_defaults.basepairs_lengths):
    d1 = 0
    if (r1n, r2n) == data[0]:
      for l in data[1:]:
        a1 = r1.get_atom(l[0])
        a2 = r2.get_atom(l[1])
        # Bug notification: a2 could be None: bug report form
        # ms@mrc-lmb.cam.ac.uk on Nov 30, 2015, file was not provided.
        # UPD. a1 also could be None, investigated and improved potential
        # hbond filtering procedure. Hopefully, this will resolve the issue.
        if a1 is None or a2 is None:
          missing_atom = l[0] if a1 is None else l[1]
          missing_res = r1 if a1 is None else r2
          print("Warning! %s atom is missing from residue %s " % (
              missing_atom, missing_res.id_str()), file=log)
          if verbose > 1:
            print("Atoms present in the residue:", file=log)
            for a in missing_res.atoms():
              print(a.id_str(), file=log)
          print("  Was trying to link: %s%s with %s%s, Saenger class: %d" % (
              r1.id_str(), l[0], r2.id_str(), l[1], class_number), file=log)
          # a1_id = a1.id_str() if a1 is not None else "None"
          # a2_id = a2.id_str() if a2 is not None else "None"
          # msg = "Something is wrong in .pdb file around '%s' or '%s'.\n" % (
          #     a1_id, a2_id)
          # msg += "If it is not clear to you, please contact developers and "
          # msg += "supply above message and .pdb file used."
          continue
          # raise Sorry(msg)
        d1 += abs(a1.distance(a2)-2.89)
      n_links_in_data = len(data[1:])
      if n_links_in_data < 1:
        raise Sorry("Corrupted dictionary in bondlength_defaults.py")
      d1 /=n_links_in_data
      if verbose > 2:
        print("  Class %d penalty=%.3f" % (class_number, d1), file=log)
      if best_score > d1:
        best_possible_link_list = [(x[:2]) for x in data[1:]]
        best_score = d1
        best_class_number = class_number
  for n1, n2 in best_possible_link_list:
    a1 = r1.get_atom(n1)
    a2 = r2.get_atom(n2)
    if verbose > 2:
      print("    %s --> %s distance = %.3f" % (
          a1.id_str(), a2.id_str(), a1.distance(a2)), file=log)
    if a1 is not None and a2 is not None and a1.distance(a2)<distance_cutoff:
      new_hbonds.append(tuple([a1, a2] if a1.i_seq<a2.i_seq else [a2, a1]))
      # new_hbonds.append(tuple(sorted([a1.i_seq, a2.i_seq])))
  return new_hbonds, best_class_number

def final_link_direction_check(atom1, atom2, rna_dna_angle_cutoff=35):
  import math
  a1p = atom1.parent().get_atom('C4')
  a2p = atom1.parent().get_atom('C5')
  a3p = atom1.parent().get_atom('C6')
  if a1p is None or a2p is None or a3p is None:
    # this is something strange, but a user managed to get here and sent
    # bug report
    return False
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
    prefix=None, params=None,
    log=sys.stdout, add_segid=None, verbose=-1):
  hbond_distance_cutoff = 3.4
  if params is not None:
    hbond_distance_cutoff = params.hbond_distance_cutoff
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
        a1.name.find("P") < 0 and a2.name.find("P") < 0 and
        a1.name.find("'") < 0 and a2.name.find("'") < 0 and
        not consecutive_residues(a1, a2) and
        (ag1.altloc.strip() == ag2.altloc.strip()) and
        final_link_direction_check(a1, a2)):
      hbonds.append((i_seq, j_seq))
    i += 1
  # check and define basepairs
  pairs = []
  for hb in hbonds:
    if verbose > 1:
      print("Making pair with", atoms[hb[0]].id_str(), atoms[hb[1]].id_str(), file=log)
    new_hbonds, class_number = get_h_bonds_for_basepair(
        atoms[hb[0]],
        atoms[hb[1]],
        distance_cutoff=hbond_distance_cutoff,
        log=log,
        verbose=verbose)
    if verbose > 1:
      print("  Picked class: %d, number of h-bonds under cutoff:%d" % (class_number, len(new_hbonds)), end=' ', file=log)
    if len(new_hbonds) > 1:
      p = make_phil_base_pair_record(atoms[hb[0]].parent(), atoms[hb[1]].parent(),
          params, saenger_class=class_number, add_segid=add_segid)
      if verbose > 1:
        print("  OK", file=log)
      pairs.append(p)
    else:
      if verbose > 0:
        s = " ".join(["Basepairing for residues '%s' and '%s'" % (
            atoms[hb[0]].id_str()[10:-1], atoms[hb[1]].id_str()[10:-1]),
          "was rejected because only 1 h-bond was found"])
        if verbose > 1:
          print("Rejected", file=log)

  phil_str = ""
  # print "N basepairs:", len(pairs)
  for pair_phil in pairs:
    phil_str += pair_phil
  if prefix is not None:
    result = "%s {\n%s}" % (prefix, phil_str)
  else:
    result = phil_str
  return result

def get_plane_i_seqs_from_residues(r1, r2, grm,mon_lib_srv, plane_cache):
  # Return [([i_seqA],[j_seqA]),([i_seqB],[j_seqB])] of
  # atoms in planar groups for given atom_group r. Handles up to 2 alternative
  # conformations (the lenght of resulting array will be appropriate: A and B
  # in i_seqA i_seqB - altlocs) If there is no alt confs - the list will
  # only one tuple.
  # Uses geometry_restraints_manager that should already have basic planarity
  # restraints to find the largest plane group in residue, which should be the
  # nucleobase in case of nucleic acids.
  def convert_to_good_name(rn):
    result = rn
    if len(rn)>1 and rn[0] == 'D':
      result = rn[1]+rn[0]
    return result
  def print_warning_msg(rn):
    # print resname
    # print r.resname
    # print new_res.resname.strip()
    print("Warning, Cannot make NA restraints for %s residue (no planarity definition)" % rn)
  i_seqs = []
  result = []
  r1_i_seqs = {}
  r2_i_seqs = {}
  new_r1 = {}
  new_r2 = {}
  for r, r_i_seqs, new_r in [(r1, r1_i_seqs, new_r1), (r2, r2_i_seqs, new_r2)]:
    for conf in r.parent().conformers():
      for c_residue in conf.residues():
        if c_residue.resseq == r.parent().resseq:
          new_res = c_residue
          break
      r_i_seqs[conf.altloc] = []
      new_r[conf.altloc] = []

      # new getting i_seqs in plane
      resname = convert_to_good_name(new_res.resname.strip())
      if resname not in plane_cache:
        libdef = mon_lib_srv.get_comp_comp_id_direct(resname)
        if libdef is None:
          print_warning_msg(resname)
          continue
          # raise Sorry('Cannot make NA restraints for %s residue' % resname)
        planes = libdef.get_planes()
        if planes is not None and len(planes) > 0:
          best_index = 0
          best_len = len(planes[0].plane_atoms)
          for i in range(1, len(planes)):
            if len(planes[i].plane_atoms) > best_len:
              best_len = len(planes[i].plane_atoms)
              best_index = i
          plane_cache[resname] = planes[best_index].plane_atoms
        else:
          print_warning_msg(resname)
          continue
      # now this residue is in cache even if it wasn't before
      for plane_atom in plane_cache.get(resname,[]):
        good_atom_id = plane_atom.atom_id.replace("*","'")
        if good_atom_id == "C5M":
          good_atom_id = "C7"
        good_atom_id = " %s" % good_atom_id
        if len(good_atom_id) == 3:
          good_atom_id += " "
        a = new_res.find_atom_by(name=good_atom_id)
        if a is not None:
          new_r[conf.altloc].append(a.i_seq)
      new_r[conf.altloc] = sorted(new_r[conf.altloc])
      r_i_seqs[conf.altloc] = new_r[conf.altloc]
  if len(r1_i_seqs) > len(r2_i_seqs):
    t = r1_i_seqs
    r1_i_seqs = r2_i_seqs
    r2_i_seqs = t
  assert len(r1_i_seqs) > 0 and len(r1_i_seqs) > 0
  if len(r1_i_seqs) == 1:
    if len(r2_i_seqs) == 1:
      if (('' in r1_i_seqs or '' in r2_i_seqs)
          or (list(r1_i_seqs.keys())[0] == list(r2_i_seqs.keys())[0])):  # FIXME: indexing keys breaks compat py2/3 if more than 1 key
        result.append((r1_i_seqs[list(r1_i_seqs.keys())[0]],
                       r2_i_seqs[list(r2_i_seqs.keys())[0]]))
    else:
      if ('' in r1_i_seqs):
        for k,v in six.iteritems(r2_i_seqs):
          result.append((r1_i_seqs[''], v))
  else:
    for k, v in six.iteritems(r1_i_seqs):
      if k in r2_i_seqs:
        result.append((v, r2_i_seqs[k]))
  # check whether sets of iseqs are different
  if len(result) > 1 and result[0] == result[1]:
    # STOP()
    return [result[0]]
  return result

def get_stacking_proxies(pdb_hierarchy, stacking_phil_params, grm,
    mon_lib_srv, plane_cache):
  result = []
  assert pdb_hierarchy is not None
  if len(stacking_phil_params) < 1:
    return result
  if grm is None:
    return result
  selection_cache = pdb_hierarchy.atom_selection_cache()
  pdb_atoms = pdb_hierarchy.atoms()
  for stacking_pair in stacking_phil_params:
    if (stacking_pair.base1 is not None and stacking_pair.base2 is not None
        and stacking_pair.enabled):
      selected_atoms_1 = selection_cache.iselection(stacking_pair.base1)
      selected_atoms_2 = selection_cache.iselection(stacking_pair.base2)
      if len(selected_atoms_1) == 0:
        raise Sorry("Selection %s in stacking_pair resulted in 0 atoms." % (
            stacking_pair.base1))
      if len(selected_atoms_2) == 0:
        raise Sorry("Selection %s in stacking_pair resulted in 0 atoms." % (
            stacking_pair.base2))
      a1 = pdb_atoms[selected_atoms_1[0]]
      a2 = pdb_atoms[selected_atoms_2[0]]
      r1 = a1.parent()
      r2 = a2.parent()
      seqs = get_plane_i_seqs_from_residues(r1, r2, grm, mon_lib_srv, plane_cache)
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
            limit=1,
            origin_id=origin_ids.get_origin_id('basepair stacking'))
          result.append(proxy)
  return result

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
  rnames = [atoms[0].parent().resname,
            atoms[1].parent().resname]
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
          weight=1./vals[i][1]**2,
          origin_id=origin_ids.get_origin_id('hydrogen bonds'))
        proxies.append(p)
  return proxies

def get_h_bonds_for_particular_basepair(atoms, saenger_class=0):
  result = []
  if saenger_class == 0:
    return result
  if saenger_class == 16:
    # We don't have values for this one.
    return result
  assert len(atoms) == 2
  new_hbonds = []
  r1, r2, r1n, r2n = unify_residue_names_and_order(
      atoms[0].parent(), atoms[1].parent())
  if bondlength_defaults.basepairs_lengths.get(saenger_class, None) == None:
    raise Sorry("""
Bad or unknown Saenger class. Presently we don't have enough
data to support Saenger class #16. """)
  if bondlength_defaults.basepairs_lengths[saenger_class][0] != (r1n, r2n):
    print(bondlength_defaults.basepairs_lengths[saenger_class][0], r1n, r2n,saenger_class)
    print(r1.id_str(), r2.id_str())
    raise Sorry("Saenger class does not match residue names")
  hbonds = []
  for b in bondlength_defaults.basepairs_lengths[saenger_class][1:]:
    hba1 = r1.get_atom(b[0])
    hba2 = r2.get_atom(b[1])
    hbonds.append((hba1, hba2))
  return hbonds

def get_basepair_proxies(
    pdb_hierarchy,
    bp_phil_params,
    grm,
    mon_lib_srv,
    plane_cache,
    hbond_distance_cutoff=3.4,
    scale_bonds_sigma=1.):
  assert pdb_hierarchy is not None
  bond_proxies_result = []
  angle_proxies_result = []
  result_planarities = []
  result_parallelities = []
  if len(bp_phil_params) < 1:
    return bond_proxies_result, angle_proxies_result, result_planarities, result_parallelities
  if grm is None:
    return bond_proxies_result, angle_proxies_result, result_planarities, result_parallelities
  selection_cache = pdb_hierarchy.atom_selection_cache()
  pdb_atoms = pdb_hierarchy.atoms()
  for base_pair in bp_phil_params:
    if (base_pair.base1 is not None and base_pair.base2 is not None
        and base_pair.enabled):
      selected_atoms_1 = selection_cache.iselection(base_pair.base1)
      selected_atoms_2 = selection_cache.iselection(base_pair.base2)
      if len(selected_atoms_1) == 0:
        raise Sorry("Selection %s in base_pair resulted in 0 atoms." % (
            base_pair.base1))
      if len(selected_atoms_2) == 0:
        raise Sorry("Selection %s in base_pair resulted in 0 atoms." % (
            base_pair.base2))
      a1 = pdb_atoms[selected_atoms_1[0]]
      a2 = pdb_atoms[selected_atoms_2[0]]
      # get hbonds
      bp_proxies, ap_proxies = get_bp_hbond_proxies(
          a1, a2, base_pair, hbond_distance_cutoff, scale_bonds_sigma)
      bond_proxies_result += bp_proxies
      angle_proxies_result += ap_proxies
      # get planarity/parallelity
      plan_p, parr_p = get_bp_plan_proxies(
          a1, a2, base_pair, grm, mon_lib_srv, plane_cache)
      result_planarities += plan_p
      result_parallelities += parr_p
  return bond_proxies_result, angle_proxies_result, result_planarities, result_parallelities

def get_bp_hbond_proxies(a1, a2, base_pair, hbond_distance_cutoff,
    scale_bonds_sigma):
  bp_result = []
  ap_result = []
  if base_pair.saenger_class == 0:
    hbonds, saenger_class = get_h_bonds_for_basepair(
      a1, a2, distance_cutoff=hbond_distance_cutoff,
      log=sys.stdout, verbose=-1)
    base_pair.saenger_class = saenger_class
  hbonds = get_h_bonds_for_particular_basepair((a1, a2), base_pair.saenger_class)
  for hb in hbonds:
    if hb[0] is None or hb[1] is None:
      print("NA hbond rejected because one of the atoms is absent")
      continue
    dist = hb[0].distance(hb[1])
    if dist < hbond_distance_cutoff:
      if base_pair.restrain_hbonds:
        hb_target, hb_sigma = get_hb_lenght_targets(hb)
        p = geometry_restraints.bond_simple_proxy(
          i_seqs=[hb[0].i_seq, hb[1].i_seq],
          distance_ideal=hb_target,
          weight=1.0/(hb_sigma*scale_bonds_sigma)**2,
          slack=0,
          top_out=False,
          limit=1,
          origin_id=origin_ids.get_origin_id('hydrogen bonds'))
        bp_result.append(p)
        # print "bond:", hb[0].id_str(), hb[1].id_str(), "(%4.2f, %4.2f)" % (hb_target, hb_sigma)
        # s1 = pdb_atoms[hb[0].i_seq].id_str()
        # s2 = pdb_atoms[hb[1].i_seq].id_str()
        # ps = "dist chain \"%s\" and resi %s and name %s, chain \"%s\" and resi %s and name %s\n" % (
        #   s1[14:15], s1[15:19], s1[5:8], s2[14:15], s2[15:19], s2[5:8])
        # dashes.write(ps)
      if base_pair.restrain_hb_angles:
        ap_result += get_angle_proxies_for_bond(hb)
    else:
      # Use log channel!
      pass
      #print("NA hbond rejected:",hb[0].id_str(), hb[1].id_str(), "distance=%.2f" % dist)
  return bp_result, ap_result

def get_bp_plan_proxies(a1, a2, base_pair, grm, mon_lib_srv, plane_cache):
  result_plan_p = []
  result_parr_p = []
  seqs = get_plane_i_seqs_from_residues(
      a1.parent(), a2.parent(), grm,mon_lib_srv, plane_cache)
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
          limit=1,
          origin_id=origin_ids.get_origin_id('basepair parallelity'))
        result_parr_p.append(proxy)
      if base_pair.restrain_planarity:
        if base_pair.planarity_sigma < 1e-5:
          raise Sorry("Sigma for planarity basepair restraints should be > 1e-5")
        w = 1./(base_pair.planarity_sigma**2)
        proxy=geometry_restraints.planarity_proxy(
          i_seqs=flex.size_t(i_seqs+j_seqs),
          weights=[w]*len(i_seqs+j_seqs),
          origin_id=origin_ids.get_origin_id('basepair planarity'))
        result_plan_p.append(proxy)
  return result_plan_p, result_parr_p
