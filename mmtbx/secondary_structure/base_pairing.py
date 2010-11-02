
from libtbx import easy_run
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
  pair_type = *wwt
    .type = choice
    .caption = Watson-Crick
# TODO
#  planar = False
#    .type = bool
}
"""

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
      assert len(fields) >= 4
      pair_type = fields[0]
      paired_bases = fields[1]
      hydrogen_flag = fields[2]
      atom_pairs = [ (p.split(",")[0], p.split(",")[1]) for p in fields[3:] ]
      if hydrogen_flag == '+' :
        db = self._h_bond_pairs
      else :
        db = self._pseudo_bond_pairs
      if pair_type in db :
        if paired_bases in db[pair_type] :
          raise RuntimeError("Duplicate entry in base pair dictionary, line %d"
            % i)
      else :
        db[pair_type] = {}
      db[pair_type][paired_bases] = atom_pairs

  def get_atoms (self, base_pair, pair_type, use_hydrogens=False) :
    if use_hydrogens :
      db = self._h_bond_pairs
    else :
      db = self._pseudo_bond_pairs
    pair_rules = db[pair_type]
    if base_pair in pair_rules :
      return pair_rules[base_pair]
    elif base_pair[::-1] in pair_rules :
      return invert_pairs(pair_rules[base_pair[::-1]])
    else :
      raise RuntimeError("No entry for base pair %s with type %s (H=%s)." %
        (base_pair, pair_type, str(use_hydrogens)))

  def get_pair_type (self, base_pair, atom_pairs, use_hydrogens=False) :
    return_pair_type = None
    if use_hydrogens :
      db = self._h_bond_pairs
    else:
      db = self._pseudo_bond_pairs
    for pair_type in db:
      if base_pair in db[pair_type]:
        if db[pair_type][base_pair].sort(key=sort_tuple) == atom_pairs.sort(key=sort_tuple):
          if return_pair_type is None:
            return_pair_type = pair_type
          else:
            raise RuntimeError("Redundant entries found for base pair %s." % base_pair)
      elif base_pair[::-1] in db[pair_type]:
        if db[pair_type][base_pair[::-1]].sort(key=sort_tuple) == \
           invert_pairs(atom_pairs).sort(key=sort_tuple):
          if return_pair_type is None:
            return_pair_type = pair_type
          else:
            raise RuntimeError("Redundant entries found for base pair %s." % base_pair)
    return return_pair_type

def invert_pairs(h_bond_atoms):
  inverted_h_bond_atoms = []
  for pair in h_bond_atoms:
    inverted_h_bond_atoms.append((pair[1],pair[0]))
  return inverted_h_bond_atoms

def sort_tuple(tuple):
  return tuple[0]

db = pair_database()

def get_h_bond_atoms(residues, pair_type, use_hydrogens=False):
  base_pair = residues[0].strip()[0] + residues[1].strip()[0]
  return db.get_atoms(base_pair, pair_type, use_hydrogens)

def run_probe(pdb_hierarchy, flags=None, add_hydrogens=True):
  reduce_output = run_reduce(pdb_hierarchy, remove_hydrogens=False).stdout_lines
  cmd = "phenix.probe " + flags + " -"
  probe_output = easy_run.fully_buffered(cmd,
           stdin_lines=reduce_output)
  return probe_output

def clean_base_names(pdb_line):
  clean_line = pdb_line
  #clean RNA lines
  clean_line = re.sub(r'^((ATOM  |HETATM).{11}) ([ACGU])([Rr])',r'\1  \3',clean_line)
  #clean DNA lines
  clean_line = re.sub(r'^((ATOM  |HETATM).{11}) ([ACGT])([Dd])',r'\1 D\3',clean_line)
  return clean_line

def run_reduce(hierarchy, remove_hydrogens=True):
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
    base_pair = pair[7:10].strip()+pair[17:20].strip()
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
    phil_strings.append("""base_pair {
  base1 = \"\"\"chain "%s" %sand resseq %s\"\"\"
  base2 = \"\"\"chain "%s" %sand resseq %s\"\"\"
  pair_type = %s
}""" % (chains[0], segid_extra, bases[0][2:6], chains[1], segid_extra,
        bases[1][2:6], pair_type))
  phil_str = """nucleic_acids {\n%s\n}""" % ("\n".join(phil_strings))
  if prefix is not None :
    return """%s {\n%s\n}""" % (prefix, phil_str)
  return phil_str

def exercise () :
  import libtbx.load_env
  if libtbx.env.has_module("probe") and libtbx.env.has_module("reduce"):
    assert (db.get_atoms("AU", "WWT", True) == [('H61', 'O4'), ('N1', 'H3')])
    assert (db.get_atoms("GC", "WWT", False) == [('O6', 'N4'), ('N1', 'N3'),
      ('N2', 'O2')])
    assert db.get_pair_type("AU", [('H61', 'O4'), ('N1', 'H3')], True) == "WWT"
    assert db.get_pair_type("AU", [('N1', 'H3'), ('H61', 'O4')], True) == "WWT"
    assert db.get_pair_type("CG", [('N4', 'O6'), ('N3', 'N1'), ('O2', 'N2')], False) == "WWT"
  else:
    print "Skipping: probe and/or reduce not available"
  print "OK"

if __name__ == "__main__" :
  exercise()
