
from __future__ import division
from libtbx import easy_run
from libtbx.utils import Sorry
from libtbx import Auto
from mmtbx import utils
import os, sys, re

dna_rna_params_str = """
base_pair
  .multiple = True
  .optional = True
  .style = noauto
{
  base1 = None
    .type = str
  base2 = None
    .type = str
  saenger_class = None
    .type = str
    .optional = True
    .caption = Saenger number if applicable
    .help = reference
  leontis_westhof_class = *Auto wwt
    .type = choice
    .caption = Watson-Crick
    .help = reference
# TODO
#  planar = False
#    .type = bool
}
"""

saenger_list = [ "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
                 "XI", "XII", "XIII", "XIV", "XV", "XVI", "XVII", "XVIII",
                 "XIX", "XX", "XXI", "XXII", "XXIII", "XXIV", "XXV", "XXVI",
                 "XXVII", "XXVIII" ]

class pair_database (object) :
  def __init__ (self) :
    self._h_bond_pairs = {}
    self._pseudo_bond_pairs = {}
    base, ext = os.path.splitext(__file__)
    dblines = open("%s.data" % base).readlines()
    for i, line in enumerate(dblines) :
      if line.startswith("#") :
        continue
      fields = line.strip().split()
      assert len(fields) >= 5
      lw_class = fields[0]
      if lw_class == '_':
        lw_class = 'other'
      saenger_class = fields[1]
      if saenger_class == '_':
        saenger_class = None
      paired_bases = fields[2]
      hydrogen_flag = fields[3]
      atom_pairs = [ (p.split(",")[0],
                      p.split(",")[1]) for p in fields[4:] ]
      distances = [ [p.split(",")[2],
                     p.split(",")[3],
                     p.split(",")[4]] for p in fields[4:] ]
      for i in range(len(distances)) :
        for j in range(3) :
          if (distances[i][j] == '_') :
            distances[i][j] = None
        distances[i] = tuple(distances[i])
      if hydrogen_flag == '+' :
        db = self._h_bond_pairs
      else :
        db = self._pseudo_bond_pairs
      for key in [saenger_class, lw_class]:
        if key != None and key != 'other':
          if key in db :
            if paired_bases in db[key] :
              raise RuntimeError("Duplicate entry in base pair dictionary, line %d"
                % i)
          else :
            db[key] = {}
          db[key][paired_bases] = [atom_pairs,distances]

  def get_pair_saenger_classes (self, base_pair, use_hydrogens=False) :
    if use_hydrogens :
      db = self._h_bond_pairs
    else :
      db = self._pseudo_bond_pairs
    allowed_types = []
    for pair_type, pair_rules in db.iteritems() :
      if (not pair_type in saenger_list) :
        continue
      if base_pair in pair_rules :
        atom_pairs, distances = pair_rules[base_pair]
        allowed_types.append((pair_type, atom_pairs, distances))
    return allowed_types

  def invert_bases (self, base_pair) :
    bases = base_pair.split('-')
    assert len(bases) == 2
    return bases[1] + '-' + bases[0]

  def get_atoms (self, base_pair, pair_type, use_hydrogens=False) :
    if use_hydrogens :
      db = self._h_bond_pairs
    else :
      db = self._pseudo_bond_pairs
    pair_rules = db[pair_type]
    if base_pair in pair_rules :
      return pair_rules[base_pair][0]
    elif self.invert_bases(base_pair) in pair_rules :
      return invert_pairs(pair_rules[self.invert_bases(base_pair)][0])
    else :
      raise RuntimeError("No entry for base pair %s with type %s (H=%s)." %
        (base_pair, pair_type, str(use_hydrogens)))

  def get_distances (self, base_pair, pair_type, use_hydrogens=False) :
    if use_hydrogens :
      db = self._h_bond_pairs
    else :
      db = self._pseudo_bond_pairs
    pair_rules = db[pair_type]
    if base_pair in pair_rules :
      return pair_rules[base_pair][1]
    elif self.invert_bases(base_pair) in pair_rules :
      return pair_rules[self.invert_bases(base_pair)][1]
    else :
      raise RuntimeError("No entry for base pair %s with type %s (H=%s)." %
        (base_pair, pair_type, str(use_hydrogens)))

  def get_pair_type (self, base_pair, atom_pairs, use_hydrogens=False) :
    return_pair_type = None
    inverted_atom_pairs = invert_pairs(atom_pairs)
    if use_hydrogens :
      db = self._h_bond_pairs
    else:
      db = self._pseudo_bond_pairs
    for pair_type in db:
      if base_pair in db[pair_type]:
        db[pair_type][base_pair][0].sort(key=sort_tuple)
        atom_pairs.sort(key=sort_tuple)
        if db[pair_type][base_pair][0] == atom_pairs:
          if return_pair_type is None:
            return_pair_type = pair_type
          elif (not is_saenger(return_pair_type)) and (is_saenger(pair_type)):
            return_pair_type = pair_type
          elif (is_saenger(return_pair_type)) and (not is_saenger(pair_type)):
            pass
          else:
            print return_pair_type, pair_type
            raise RuntimeError("Redundant entries found for base pair %s." % base_pair)
      elif self.invert_bases(base_pair) in db[pair_type]:
        db[pair_type][self.invert_bases(base_pair)][0].sort(key=sort_tuple)
        inverted_atom_pairs.sort(key=sort_tuple)
        if db[pair_type][self.invert_bases(base_pair)][0] == \
           inverted_atom_pairs:
          if return_pair_type is None:
            return_pair_type = pair_type
          else:
            raise RuntimeError("Redundant entries found for base pair %s." % base_pair)
    return return_pair_type

def is_saenger(key):
  if key in saenger_list:
    return True
  else:
    return False

def invert_pairs(h_bond_atoms):
  inverted_h_bond_atoms = []
  for pair in h_bond_atoms:
    inverted_h_bond_atoms.append((pair[1],pair[0]))
  return inverted_h_bond_atoms

def sort_tuple(tuple):
  return tuple[0]

db = pair_database()

def get_h_bond_atoms(residues,
                     saenger_class,
                     leontis_westhof_class,
                     use_hydrogens=False):
  base_pair = clean_single_base_name(residues[0]).strip() + '-' + \
              clean_single_base_name(residues[1]).strip()
  if (saenger_class is not None) :
    pair_type = saenger_class.upper()
  elif (leontis_westhof_class is not None) :
    assert (leontis_westhof_class != "Auto")
    pair_type = leontis_westhof_class.upper()
  else :
    raise Sorry("Base pair type not specified.")
  return db.get_atoms(base_pair, pair_type, use_hydrogens)

def get_distances(residues,
                  saenger_class,
                  leontis_westhof_class,
                  use_hydrogens=False):
  base_pair = clean_single_base_name(residues[0]).strip() + '-' + \
              clean_single_base_name(residues[1]).strip()
  if (saenger_class is not None) :
    pair_type = saenger_class.upper()
  elif (leontis_westhof_class is not None) :
    pair_type = leontis_westhof_class.upper()
  else :
    raise Sorry("Base pair type not specified.")
  return db.get_distances(base_pair, pair_type, use_hydrogens)

def run_probe(pdb_hierarchy, flags=None, add_hydrogens=True):
  bare_chains = utils.find_bare_chains_with_segids(
                  pdb_hierarchy=pdb_hierarchy)
  if bare_chains:
    tmp_pdb_hierarchy = pdb_hierarchy.deep_copy()
    tmp_pdb_hierarchy.atoms().reset_i_seq()
    seg_dict = \
      utils.seg_id_to_chain_id(pdb_hierarchy=tmp_pdb_hierarchy)
    rename_txt = utils.assign_chain_ids(pdb_hierarchy=tmp_pdb_hierarchy,
                                        seg_dict=seg_dict)
  else:
    tmp_pdb_hierarchy = pdb_hierarchy
  reduce_output = run_reduce(tmp_pdb_hierarchy,
                             remove_hydrogens=False).stdout_lines
  if bare_chains:
    reverted_reduce_output = ""
    for line in reduce_output:
      if line.startswith("ATOM"):
        segid = line[72:76]
        tmp_line = ""
        if seg_dict.get(segid):
          tmp_line = line[0:20] + "  " + line[22:]
          reverted_reduce_output += tmp_line+'\n'
        else:
          reverted_reduce_output += line+'\n'
      else:
        reverted_reduce_output += line+'\n'
  else:
    reverted_reduce_output = reduce_output
  cmd = "phenix.probe " + flags + " -"
  probe_output = easy_run.fully_buffered(cmd,
           stdin_lines=reverted_reduce_output)
  return probe_output

def clean_base_names(pdb_line):
  clean_line = pdb_line
  #clean RNA lines
  clean_line = re.sub(r'^((ATOM  |HETATM).{11}) ([ACGU])([Rr])',r'\1  \3',clean_line)
  clean_line = re.sub(r'^((ATOM  |HETATM).{11})([ACGU])  ',r'\1  \3',clean_line)
  #clean DNA lines
  clean_line = re.sub(r'^((ATOM  |HETATM).{11}) ([ACGT])([Dd])',r'\1 D\3',clean_line)
  return clean_line

def clean_single_base_name(base):
  clean_base = base
  #clean RNA base names
  clean_base = re.sub(r' ([ACGU])([Rr])',r'  \1',clean_base)
  clean_base = re.sub(r'([ACGU])  ',r'  \1',clean_base)
  #clean DNA base names
  clean_base = re.sub(r' ([ACGT])([Dd])',r' D\1',clean_base)
  return clean_base


def run_reduce(hierarchy, remove_hydrogens=True):
  #log = sys.stderr
  trim = "phenix.reduce -quiet -trim -"
  build = "phenix.reduce -quiet -build -allalt -"
  input_str = ""
  pdb_string = hierarchy.as_pdb_string()
  clean_lines = []
  for line in pdb_string.splitlines() :
    # XXX for some reason the ANISOU records are confusing Reduce - maybe
    # because residue names aren't being capitalized there?  fortunately we
    # don't need ANISOUs here anyway.
    if not line.startswith("ANISOU") :
      # *'s in atom names don't impact base, so leaving alone for now
      #line = re.sub(r'^((ATOM  |HETATM).{9})\*',r'\1'+"\'",line)
      line = clean_base_names(line)
      clean_lines.append(line)
  clean_pdb_string = "\n".join(clean_lines)
  if(remove_hydrogens):
    clean = easy_run.fully_buffered(trim,
                                    stdin_lines=clean_pdb_string)
    input_str = clean.stdout_lines
  else:
    input_str = clean_pdb_string
  output = easy_run.fully_buffered(build,
                                   stdin_lines=input_str)
  return output

def get_base_pairs(pdb_hierarchy, probe_flags=None):
  #db = mmtbx.base_pairing.pair_database()
  if probe_flags is None:
    probe_flags = "-Both -Unformated -NOCLASHOUT -NOVDWOUT \"BASE\" \"BASE\""
  probe_out = run_probe(pdb_hierarchy=pdb_hierarchy,
                                  flags=probe_flags)
  hbond_hash={}
  pair_hash={}
  reduced_pair_hash={}
  probe_iter = iter(probe_out.stdout_lines)
  for line in probe_iter:
    temp = line.split(":")
    bases = [temp[3], temp[4]]
    #sort so C before G, A before U
    bases.sort(key=get_resname)
    base1 = bases[0][0:10]
    base2 = bases[1][0:10]
    atom1 = bases[0][11:15]
    atom2 = bases[1][11:15]
    base_key = base1+base2
    bond_key = atom1.strip()+','+atom2.strip()
    #filter out carbon weak h-bonds for the time being
    if atom1.strip()[0] == 'C' or atom2.strip()[0] == 'C':
      continue
    if (not base_key in hbond_hash):
      hbond_hash[base_key] = {}
    if (not bond_key in hbond_hash[base_key]) :
      hbond_hash[base_key][bond_key] = 0
    hbond_hash[base_key][bond_key]+=1

  for key in hbond_hash:
    counter = 0
    pair_hash[key]=[]
    for bond in hbond_hash[key] :
      if (hbond_hash[key][bond] > 0) :
        pair_hash[key].append(tuple(bond.split(',')))
    if (len(pair_hash[key]) >= 2) :
      pair_hash[key].sort(key=sort_tuple)
      reduced_pair_hash[key] = pair_hash[key]
  base_pair_list = []
  for pair in reduced_pair_hash :
    bases = (pair[:10], pair[10:])
    base_pair = pair[7:10].strip() + '-' + pair[17:20].strip()
    pair_type = db.get_pair_type(base_pair, reduced_pair_hash[pair], use_hydrogens=True)
    if (pair_type is not None) :
      base_pair_list.append([bases, pair_type])
  return base_pair_list

def get_resname(probe_field):
  resname = probe_field[6:10]
  return resname

def get_phil_base_pairs (pdb_hierarchy, probe_flags=None, prefix=None,
    log=sys.stderr, add_segid=None) :
  print >> log, "  Running PROBE to identify base pairs"
  base_pair_list = get_base_pairs(pdb_hierarchy, probe_flags)
  if len(base_pair_list) == 0 :
    return None
  phil_strings = []
  segid_extra = ""
  if add_segid is not None :
    segid_extra = """and segid "%s" """ % add_segid
  for (bases, pair_type) in base_pair_list :
    chains = []
    for base in bases :
      chain = base[0:2].strip()
      if (chain == "") :
        chain = " "
      chains.append(chain)
    if (is_saenger(pair_type)):
      type_key = "saenger_class"
    else:
      type_key = "leontis_westhof_class"
    distances = db.get_distances(bases[0][-3:].strip()+'-'+bases[1][-3:].strip(),
                                 pair_type,
                                 use_hydrogens=False)
    phil_strings.append("""base_pair {
  base1 = \"\"\"chain "%s" %sand resseq %s\"\"\"
  base2 = \"\"\"chain "%s" %sand resseq %s\"\"\"
  %s = %s
}""" % (chains[0], segid_extra, bases[0][2:6], chains[1], segid_extra,
        bases[1][2:6], type_key, pair_type))
  phil_str = """nucleic_acids {\n%s\n}""" % ("\n".join(phil_strings))
  if prefix is not None :
    return """%s {\n%s\n}""" % (prefix, phil_str)
  return phil_str

def create_hbond_proxies (
    build_proxies,
    base_pairs,
    pdb_hierarchy,
    restraint_type,
    hbond_counts,
    distance_ideal,
    distance_cut,
    remove_outliers,
    use_hydrogens,
    sigma,
    slack,
    use_db_values=True,
    raise_error_if_missing_selection=False,
    log=sys.stderr) :
  assert (restraint_type in ["simple", "lennard_jones"])
  assert (slack >= 0) and (sigma >= 0)
  from mmtbx.secondary_structure import utils as ss_utils
  selection_cache = pdb_hierarchy.atom_selection_cache()
  atoms = pdb_hierarchy.atoms()
  n_proxies = 0
  for base_pair in base_pairs :
    if (base_pair.base1 is None) or (base_pair.base2 is None) :
      msg = ("  One or more selections missing for base pair:\n" +
          "    base 1: %s\n" +
          "    base 2: %s") % (str(base_pair.base1),str(base_pair.base2))
      if (raise_error_if_missing_selection) :
        raise Sorry("%s\n  Please define both bases." % msg)
      else :
        print >> log, msg
        continue
    try :
      resname1 = ss_utils.get_residue_name_from_selection(
        resi_sele=base_pair.base1,
        selection_cache=selection_cache,
        atoms=atoms)
      resname2 = ss_utils.get_residue_name_from_selection(
        resi_sele=base_pair.base2,
        selection_cache=selection_cache,
        atoms=atoms)
      atom_pairs = get_h_bond_atoms(
        residues=(resname1,resname2),
        saenger_class=base_pair.saenger_class,
        leontis_westhof_class=base_pair.leontis_westhof_class,
        use_hydrogens=use_hydrogens)
      distance_values = get_distances(
        residues=(resname1,resname2),
        saenger_class=base_pair.saenger_class,
        leontis_westhof_class=base_pair.leontis_westhof_class,
        use_hydrogens=use_hydrogens)
      assert (len(atom_pairs) == len(distance_values))
      for i, (name1, name2) in enumerate(atom_pairs) :
        sele1 = """name %s and %s""" % (name1, base_pair.base1)
        sele2 = """name %s and %s""" % (name2, base_pair.base2)
        # XXX these aren't necessarily in donor/acceptor order, but it
        # doesn't really matter here since we're only using the simple
        # potential type.
        (i_seq,j_seq) = ss_utils.hydrogen_bond_from_selection_pair(sele1,
          sele2, selection_cache)
        if (hbond_counts[i_seq] > 2) \
           or (hbond_counts[j_seq] > 2):
          print >> log, "One or more atoms already bonded:"
          print >> log, "  %s" % atoms[i_seq].fetch_labels().id_str()
          print >> log, "  %s" % atoms[j_seq].fetch_labels().id_str()
          continue
        hbond_counts[i_seq] += 1
        hbond_counts[j_seq] += 1
        if (use_db_values) :
          bp_distance_cut = -1
          if (distance_values[i][0] is not None) :
            bp_distance = float(distance_values[i][0])
          else:
            bp_distance = distance_ideal
          if (distance_values[i][1] is not None) :
            bp_sigma = float(distance_values[i][1])
          else:
            bp_sigma = sigma
          if (distance_values[i][2] is not None) :
            bp_slack = float(distance_values[i][2])
          else:
            bp_slack = slack
        else :
          bp_distance = distance_ideal
          bp_distance_cut = distance_cut
          bp_sigma = sigma
          bp_slack = slack
        if (remove_outliers) and (distance_cut > 0) :
          dist = atoms[i_seq].distance(atoms[j_seq])
          if (dist > distance_cut) :
            continue
        build_proxies.add_proxy(
          i_seqs=[i_seq,j_seq],
          distance_ideal=bp_distance,
          distance_cut=bp_distance_cut,
          weight=1/(bp_sigma ** 2),
          slack=bp_slack)
        build_proxies.add_nonbonded_exclusion(i_seq, j_seq)
        n_proxies += 1
    except RuntimeError, e :
      print >> log, "  %s" % str(e)
  return n_proxies

def identify_base_pairs (base_pairs,
                         pdb_hierarchy,
                         use_hydrogens,
                         distance_ideal,
                         use_db_values=True) :
  import boost.python
  ext = boost.python.import_ext("mmtbx_secondary_structure_ext")
  atoms = pdb_hierarchy.atoms()
  atom_names = atoms.extract_name().strip()
  sites_cart = atoms.extract_xyz()
  sel_cache = pdb_hierarchy.atom_selection_cache()
  for base_pair in base_pairs :
    if ((base_pair.saenger_class is not None) or
        (base_pair.leontis_westhof_class != "Auto")) :
      continue
    base_1 = sel_cache.selection(base_pair.base1).iselection()
    base_2 = sel_cache.selection(base_pair.base2).iselection()
    basecode_1 = atoms[base_1[0]].fetch_labels().resname.strip()[0]
    basecode_2 = atoms[base_2[0]].fetch_labels().resname.strip()[0]
    bp = "%s-%s" % (basecode_1, basecode_2)
    classes = db.get_pair_saenger_classes(bp)
    class_dist_sq = []
    for (class_name, atom_pairs, distances) in classes :
      dist_sq = _get_distance_score_for_class(
        atom_pairs=atom_pairs,
        distances=distances,
        base_1=base_1,
        base_2=base_2,
        atom_names=atom_names,
        sites_cart=sites_cart,
        distance_ideal=distance_ideal,
        use_db_values=use_db_values)
      class_dist_sq.append((class_name, dist_sq))
    classes2 = db.get_pair_saenger_classes(bp[::-1])
    for (class_name, atom_pairs, distances) in classes2 :
      dist_sq = _get_distance_score_for_class(
        atom_pairs=atom_pairs,
        distances=distances,
        base_1=base_2, # note switch in order
        base_2=base_1,
        atom_names=atom_names,
        sites_cart=sites_cart,
        distance_ideal=distance_ideal,
        use_db_values=use_db_values)
      class_dist_sq.append((class_name, dist_sq))
    if (len(class_dist_sq) == 0) :
      base_pair.leontis_westhof_class = "wwt"
      continue
    class_dist_sq.sort(lambda x,y: cmp(x[1], y[1]))
    base_pair.saenger_class = class_dist_sq[0][0]
    #print base_pair.base1, base_pair.base2
    #for class_name, d_sq in class_dist_sq :
    #  print "%s %.3f" % (class_name, d_sq)

def _get_distance_score_for_class (atom_pairs,
                                   distances,
                                   base_1,
                                   base_2,
                                   atom_names,
                                   sites_cart,
                                   distance_ideal,
                                   use_db_values) :
  import boost.python
  ext = boost.python.import_ext("mmtbx_secondary_structure_ext")
  sum_dist_sq = 0
  for k, (atm1, atm2) in enumerate(atom_pairs) :
    if (use_db_values) :
      pair_dist_ = distances[k][0]
      if (pair_dist_ is not None) :
        pair_dist = float(pair_dist_)
      else :
        pair_dist = distance_ideal
    sum_dist_sq += ext.delta_distance_squared(
      base_1=base_1,
      base_2=base_2,
      name_1=atm1,
      name_2=atm2,
      atom_names=atom_names,
      sites_cart=sites_cart,
      distance_ideal=pair_dist)
  return sum_dist_sq / len(atom_pairs)

########################################################################
def exercise () :
  import libtbx.load_env
  gu_classes = db.get_pair_saenger_classes("G-U")
  assert (gu_classes[0][0] == "XXVII")
  if libtbx.env.has_module("probe") and libtbx.env.has_module("reduce"):
    assert (db.get_atoms("A-U", "WWT", True) == [('H61', 'O4'), ('N1', 'H3')])
    assert (db.get_atoms("DA-DT", "WWT", False) == [('N6', 'O4'), ('N1', 'N3')])
    assert (db.get_atoms("DT-DA", "WWT", False) == [('O4', 'N6'), ('N3', 'N1')])
    assert (db.get_atoms("G-C", "WWT", False) == [('O6', 'N4'), ('N1', 'N3'),
      ('N2', 'O2')])
    assert (db.get_atoms("C-G", "WWT", False) == [('N4', 'O6'), ('N3', 'N1'),
      ('O2', 'N2')])
    assert db.get_pair_type("A-U", [('H61', 'O4'), ('N1', 'H3')], True) == "XX"
    assert db.get_pair_type("A-U", [('N1', 'H3'), ('H61', 'O4')], True) == "XX"
    assert db.get_pair_type("C-G", [('N4', 'O6'), ('N3', 'N1'), ('O2', 'N2')], False) == "XIX"
  else:
    print "Skipping: probe and/or reduce not available"
  if libtbx.env.has_module("phenix_regression") :
    pdb_file = libtbx.env.find_in_repositories(
      relative_path="phenix_regression/pdb/1u8d.pdb",
      test=os.path.isfile)
    from iotbx.file_reader import any_file
    import libtbx.phil
    pdb_in = any_file(pdb_file).file_object
    hierarchy = pdb_in.construct_hierarchy()
    hierarchy.atoms().reset_i_seq()
    bp_phil = libtbx.phil.parse(dna_rna_params_str)
    params = bp_phil.fetch(source=libtbx.phil.parse("""
      base_pair {
        base1 = chain A and resseq 18
        base2 = chain A and resseq 78
        leontis_westhof_class = *Auto
      }
      base_pair {
        base1 = chain A and resseq 21
        base2 = chain A and resseq 75
        leontis_westhof_class = *Auto
      }
      base_pair {
        base1 = chain A and resseq 33
        base2 = chain A and resseq 66
        leontis_westhof_class = *Auto
      }""")).extract()
    identify_base_pairs(
      pdb_hierarchy=hierarchy,
      base_pairs=params.base_pair,
      use_hydrogens=False,
      distance_ideal=3.0)
    classes = [ bp.saenger_class for bp in params.base_pair ]
    assert (classes == ["XIX", "XX", "V"])
    # and now in reverse order...
    params = bp_phil.fetch(source=libtbx.phil.parse("""
      base_pair {
        base1 = chain A and resseq 78
        base2 = chain A and resseq 18
        leontis_westhof_class = *Auto
      }
      base_pair {
        base1 = chain A and resseq 75
        base2 = chain A and resseq 21
        leontis_westhof_class = *Auto
      }
      base_pair {
        base1 = chain A and resseq 66
        base2 = chain A and resseq 33
        leontis_westhof_class = *Auto
      }""")).extract()
    identify_base_pairs(
      pdb_hierarchy=hierarchy,
      base_pairs=params.base_pair,
      use_hydrogens=False,
      distance_ideal=3.0)
    classes = [ bp.saenger_class for bp in params.base_pair ]
    assert (classes == ["XIX", "XX", "V"])
    # TODO: test h-bond extraction in either order???
  print "OK"

if __name__ == "__main__" :
  exercise()
