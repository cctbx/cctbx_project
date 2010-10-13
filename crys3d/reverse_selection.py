# XXX: internal data for mouse selections in wx_selection_editor.py

from libtbx.utils import Sorry
from libtbx import adopt_init_args
import re

single_quote = re.compile(r"'")

class chain_selection_info (object) :
  def __init__ (self, chain_id) :
    adopt_init_args(self, locals())
    self.resseqs = set()
    self.resseqs_by_altloc = {}
    self.remove_resseqs = set()
    self.remove_resseqs_by_altloc = {}

  def add_range (self, start, end, altloc=None) :
    resseqs = set([ x for x in range(start, end+1)])
    if altloc is None :
      if len(self.remove_resseqs) != 0 :
        self.remove_resseqs -= resseqs
      else :
        self.resseqs |= resseqs
    else :
      remove_resseqs = self.remove_resseqs_by_altloc.get(altloc, set())
      if len(remove_resseqs) != 0 :
        remove_resseqs -= resseqs
        if len(remove_resseqs) == 0 :
          self.remove_resseqs_by_altloc.pop(altloc)
        else :
          self.remove_resseqs_by_altloc[altloc] = remove_resseqs
      else :
        old_resseqs = self.resseqs_by_altloc.get(altloc, set())
        old_resseqs |= resseqs
        self.resseqs_by_altloc[altloc] = old_resseqs

  def remove_range (self, start, end, altloc=None) :
    resseqs = set([ x for x in range(start, end+1)])
    if altloc is None :
      if len(self.resseqs) != 0 :
        self.resseqs -= resseqs
      else :
        self.remove_resseqs |= resseqs
    else :
      selected_resseqs = self.resseqs_by_altloc.get(altloc, set())
      if len(selected_resseqs) > 0 :
        selected_resseqs -= resseqs
        if len(selected_resseqs) == 0 :
          self.resseqs_by_altloc.pop(altloc)
        else :
          self.resseqs_by_altloc[altloc] = selected_resseqs
      else :
        old_resseqs = self.remove_resseqs_by_altloc.get(altloc, set())
        old_resseqs |= resseqs
        self.remove_resseqs_by_altloc[altloc] = old_resseqs
    if (len(self.resseqs) == 0 and len(self.remove_resseqs) == 0 and
        len(self.resseqs_by_altloc) == 0 and
        len(self.remove_resseqs_by_altloc) == 0) :
      return True
    return False

  def __str__ (self) :
    sele_str = "chain '%s'" % self.chain_id
    clauses = []
    if len(self.resseqs) > 0 :
      ranges = get_set_ranges(self.resseqs)
      clauses.extend(assemble_resseq_ranges(ranges))
    if len(self.resseqs_by_altloc) > 0 :
      for altloc in self.resseqs_by_altloc :
        ranges = get_set_ranges(self.resseqs_by_altloc[altloc])
        clauses.extend(assemble_resseq_ranges_with_altloc(ranges, altloc))
    if len(clauses) > 0 :
      sele_str += " and ((" + ") or (".join(clauses) + "))"
    clauses = []
    if len(self.remove_resseqs) > 0 :
      ranges = get_set_ranges(self.remove_resseqs)
      clauses.extend(assemble_resseq_ranges(ranges))
    if len(self.remove_resseqs_by_altloc) > 0 :
      for altloc in self.remove_resseqs_by_altloc :
        ranges = get_set_ranges(self.remove_resseqs_by_altloc[altloc])
        clauses.extend(assemble_resseq_ranges_with_altloc(ranges, altloc))
    if len(clauses) > 0 :
      sele_str += " and not ((" + ") or (".join(clauses) + "))"
    return sele_str

class residue_selection_info (object) :
  def __init__ (self, chain_id, resid, altloc=None) :
    adopt_init_args(self, locals())

  def __str__ (self) :
    if self.altloc is None :
      return "chain '%s' and resid '%s'" % (self.chain_id, self.resid)
    else :
      return "chain '%s' and resid '%s' and altloc '%s'" % (self.chain_id,
        self.resid, self.altloc)

class atom_selection_info (object) :
  def __init__ (self, i_seq, atom) : #chain_id, resid, atom_name, altloc) :
    adopt_init_args(self, locals())

  def __str__ (self) :
    atom = self.atom
    altloc = atom.altloc
    if altloc == "" :
      altloc = " "
    return ("chain '%s' and resid '%s' and name '%s' and altloc '%s'" %
      (atom.chain_id, atom.resid(), single_quote.sub("\\\'",atom.name),
       altloc))

#-----------------------------------------------------------------------
class mouse_selection_manager (object) :
  def __init__ (self) :
    self.saved_selection = "none"
    self.selection_string = "none"
    self.flag_overwrite_mode = False
    self._selection_callback = None
    self.start_i_seq = None
    self.end_i_seq = None

  def set_overwrite_mode (self, overwrite=True) :
    self.flag_overwrite_mode = overwrite

  def selection_callback (self, selection_string, atom_selection) :
    if self._selection_callback is not None :
      self._selection_callback(selection_string, atom_selection)

  def set_selection_callback (self, callback) :
    assert hasattr(callback, "__call__") or callback is None
    self._selection_callback = callback

  def set_mmtbx_selection_function (self, mmtbx_selection_function) :
    self.mmtbx_selection_function = mmtbx_selection_function

  def update_selection_handlers (self, pdb_hierarchy,
      mmtbx_selection_function) :
    from scitbx.array_family import flex
    self.mmtbx_selection_function = mmtbx_selection_function
    self.selection_cache = pdb_hierarchy.atom_selection_cache()
    self.atom_selection = flex.bool()
    #--- XXX: for testing only
    if not hasattr(self, "atom_index") :
      self.atom_index = []
      for atom in pdb_hierarchy.atoms_with_labels() :
        self.atom_index.append(atom)
    if not hasattr(self, "pdb_hierarchy") :
      self.pdb_hierarchy = pdb_hierarchy
    #---
    self.clear_selection()

  def clear_selection (self, apply_empty=True) :
    self.selected_chains = {}
    self.deselected_chains = {}
    self.selected_atoms = []
    self.deselected_atoms = []
    self.selected_residues = []
    self.deselected_residues = []
    self.selected_pair = []
    if apply_empty :
      self.apply_selection("none") #self.saved_selection)

  def selection_size (self) :
    return self.selection_i_seqs.size()

  def apply_selection (self, selection_string) :
    if selection_string is None or selection_string == "" :
      selection_string = "none"
    atom_selection = self.get_atom_selection(selection_string)
    error = False
    if atom_selection is None :
      error = True
      self.selection_string = "none"
      atom_selection = self.selection_cache.selection("none")
    else :
      self.selection_string = selection_string
    self.atom_selection = atom_selection
    self.selection_i_seqs = atom_selection.iselection()
    if self.selection_callback is not None :
      self.selection_callback(self.selection_string, self.atom_selection)
    if error : # wait until the end to do this
      raise Sorry("Invalid selection '%s'."%selection_string)

  def get_atom_selection (self, selection_string) :
    try :
      if self.mmtbx_selection_function is not None :
        atom_selection = self.mmtbx_selection_function(
          string=selection_string,
          cache=self.selection_cache)
      else :
        atom_selection = self.selection_cache.selection(selection_string)
    except KeyboardInterrupt :
      raise
    except Exception, e :
      atom_selection =None
    return atom_selection

  def revert_selection (self) :
    self.apply_selection(self.saved_selection)

  def set_initial_selection (self, selection_string) :
    self.saved_selection = selection_string
    self.apply_selection(selection_string)

  def reset_range_selection (self) :
    self.start_i_seq = None

  def start_range_selection (self, i_seq) :
    if self.flag_overwrite_mode :
      self.clear_selection(apply_empty=False)
    self.start_i_seq = i_seq

  def end_range_selection (self, i_seq, deselect, ignore_altloc) :
    success = False
    if self.start_i_seq is not None :
      start_atom = self.atom_index[self.start_i_seq]
      end_atom = self.atom_index[i_seq]
      if (start_atom.chain_id == end_atom.chain_id and
          (ignore_altloc or start_atom.altloc == end_atom.altloc)) :
        chain_id = start_atom.chain_id
        if self.flag_overwrite_mode and not deselect :
          if chain_id in self.selected_chains :
            del self.selected_chains[chain_id]
        altloc = None
        if not ignore_altloc :
          altloc = start_atom.altloc
        start_range = start_atom.resseq_as_int()
        end_range = end_atom.resseq_as_int()
        if start_range > end_range :
          start_range = end_atom.resseq_as_int()
          end_range = start_atom.resseq_as_int()
        chain_info = self.selected_chains.get(chain_id)
        if deselect :
          deselect_str = "chain '%s' and resseq %d:%d" % (chain_id,
            start_range, end_range)
          if altloc is not None :
            deselect_str += " and altloc '%s'" % altloc
          if chain_info is not None :
            remove_chain = chain_info.remove_range(start_range, end_range,
              altloc)
            if remove_chain :
              self.selected_chains.pop(chain_id)
            success = True
          deselection = self.get_atom_selection(deselect_str)
          self.remove_redundant_residues(deselection, self.selected_residues)
          self.remove_redundant_atoms(deselection, self.selected_atoms)
          self.remove_redundant_residues(deselection.__invert__(),
            self.deselected_residues)
          self.remove_redundant_atoms(deselection.__invert__(),
            self.deselected_atoms)
        else :
          if chain_info is None :
            chain_info = chain_selection_info(chain_id)
          chain_info.add_range(start_range, end_range, altloc)
          self.selected_chains[chain_id] = chain_info
          chain_sel = self.selection_cache.selection(str(chain_info))
          self.remove_redundant_residues(chain_sel, self.selected_residues)
          self.remove_redundant_atoms(chain_sel, self.selected_atoms)
          success = True
        self.construct_selection()
    self.reset_range_selection()
    return success

  def toggle_chain_selection (self, i_seq) :
    atom = self.atom_index[i_seq]
    chain_id = atom.chain_id
    if chain_id in self.selected_chains :
      chain_info = self.selected_chains.pop(chain_id)
      chain_sel = self.selection_cache.selection(str(chain_info))
      self.remove_redundant_residues(chain_sel, self.deselected_residues)
      self.remove_redundant_atoms(chain_sel, self.deselected_atoms)
    else :
      chain_info = chain_selection_info(chain_id)
      self.selected_chains[chain_id] = chain_info
      chain_sel = self.selection_cache.selection(str(chain_info))
      self.remove_redundant_residues(chain_sel, self.selected_residues)
      self.remove_redundant_atoms(chain_sel, self.selected_atoms)
    self.construct_selection()

  def remove_redundant_residues (self, main_selection, residue_list) :
    i = 0
    while i < len(residue_list) :
      residue_info = residue_list[i]
      resi_sel = self.selection_cache.selection(str(residue_info))
      if main_selection.is_super_set(resi_sel) :
        residue_list.pop(i)
      else :
        i += 1

  def remove_redundant_atoms (self, main_selection, atom_list) :
    i = 0
    while i < len(atom_list) :
      i_seq = atom_list[i]
      if main_selection[i_seq] :
        atom_list.pop(i)
      else :
        i += 1

  def toggle_residue_selection (self, i_seq, ignore_altloc=True) :
    atom = self.atom_index[i_seq]
    if ignore_altloc :
      resi_info = residue_selection_info(atom.chain_id, atom.resid())
    else :
      resi_info = residue_selection_info(atom.chain_id, atom.resid(),
        atom.altloc)
    resi_sel = self.selection_cache.selection(str(resi_info))
    if self.flag_overwrite_mode :
      self.clear_selection(apply_empty=False)
      if not self.atom_selection.all_eq(resi_sel) :
        self.selected_residues.append(resi_info)
    else :
      if self.atom_selection.is_super_set(resi_sel) :
        self.deselected_residues.append(resi_info)
        self.remove_redundant_residues(resi_sel, self.selected_residues)
      else :
        self.selected_residues.append(resi_info)
        self.remove_redundant_residues(resi_sel, self.deselected_residues)
    self.remove_redundant_atoms(resi_sel, self.selected_atoms)
    self.remove_redundant_atoms(resi_sel, self.deselected_atoms)
    self.construct_selection()

  def select_pair (self, i_seq, selection_type="residue",
      allow_duplicate=False) :
    if len(self.selected_pair) == 2 :
      self.selected_pair = []
    atom = self.atom_index[i_seq]
    if selection_type == "residue" :
      selected_object = residue_selection_info(atom.chain_id, atom.resid())
    else :
      selected_object = atom_selection_info(i_seq, atom)
    n_selected = len(self.selected_pair)
    if n_selected == 1 and not allow_duplicate :
      # force items to be unique (no self-pairing)
      if str(selected_object) == str(self.selected_pair[0]) :
        print "skipping"
        return
    self.selected_pair.append(selected_object)
    self.construct_selection()

  def get_pair_selections (self) :
    pass

  def toggle_atom_selection (self, i_seq) :
    atom = self.atom_index[i_seq]
    if self.atom_selection[i_seq] : #and not i_seq in self.deselected_atoms :
      self.deselected_atoms.append(i_seq)
      if i_seq in self.selected_atoms :
        self.selected_atoms.remove(i_seq)
    elif not i_seq in self.selected_atoms :
      self.selected_atoms.append(i_seq)
      if i_seq in self.deselected_atoms :
        self.deselected_atoms.remove(i_seq)
    self.construct_selection()

  def select_single_residue (self, i_seq) :
    atom = self.atom_index[i_seq]
    resi_info = residue_selection_info(atom.chain_id, atom.resid())
    resi_sel = self.selection_cache.selection(str(resi_info))
    is_current_selection = ((resi_sel==self.atom_selection).count(False) == 0)
    self.clear_selection()
    if is_current_selection :
      return False
    self.toggle_residue_selection(i_seq)

  def construct_selection (self) :
    final_selection = ""
    # Part 1: stuff we want
    clauses = []
    for chain_id, chain_info in self.selected_chains.iteritems() :
      clauses.append(str(chain_info))
    chains_selection_str = assemble_selection_clauses(clauses)
    selection1 = self.get_atom_selection(chains_selection_str)
    deselection1 = selection1.__invert__()
    self.remove_redundant_residues(selection1, self.selected_residues)
    self.remove_redundant_residues(deselection1, self.deselected_residues)
    for residue_info in self.selected_residues :
      residue_selection_str = str(residue_info)
      clauses.append(residue_selection_str)
    chains_and_resi_sel_str = assemble_selection_clauses(clauses)
    selection2 = self.get_atom_selection(chains_and_resi_sel_str)
    deselection2 = selection2.__invert__()
    self.remove_redundant_atoms(selection2, self.selected_atoms)
    self.remove_redundant_atoms(deselection2, self.deselected_atoms)
    for other_info in self.selected_pair :
      clauses.append(str(other_info))
    for i_seq in self.selected_atoms :
      atom_info = atom_selection_info(i_seq, self.atom_index[i_seq])
      clauses.append(str(atom_info))
    positive_selection = assemble_selection_clauses(clauses)
    # Part 2: stuff we don't want
    clauses = []
    for residue_info in self.deselected_residues :
      residue_selection_str = str(residue_info)
      clauses.append(residue_selection_str)
    for i_seq in self.deselected_atoms :
      atom_info = atom_selection_info(i_seq, self.atom_index[i_seq])
      clauses.append(str(atom_info))
    negative_selection = assemble_selection_clauses(clauses)
    # assemble final selection
    if positive_selection != "" :
      if negative_selection != "" :
        final_selection = "(%s) and not (%s)" % (positive_selection,
                                                 negative_selection)
      else :
        final_selection = positive_selection
    else :
      final_selection = "none"
    self.apply_selection(final_selection)

def assemble_selection_clauses (clauses) :
  if len(clauses) == 0 :
    return ""
  elif len(clauses) == 1 : # not necessary, but looks nicer
    return clauses[0]
  else :
    return "(" + ") or (".join(clauses) + ")"

def get_set_ranges (residue_set) :
  resseqs = sorted(residue_set)
  previous_resseq = resseqs[0]
  ranges = []
  current_start = resseqs[0]
  for resseq in resseqs[1:-1] :
    if resseq > (previous_resseq + 1) :
      ranges.append((current_start, previous_resseq))
      current_start = resseq
    previous_resseq = resseq
  ranges.append((current_start, resseqs[-1]))
  return ranges

def assemble_resseq_ranges (ranges) :
  clauses = []
  for (start, end) in ranges :
    if start == end :
      clauses.append("resseq %d" % start)
    else :
      clauses.append("resseq %d:%d" % (start, end))
  return clauses

def assemble_resseq_ranges_with_altloc (ranges, altloc) :
  clauses = []
  for (start, end) in ranges :
    if start == end :
      clauses.append("resseq %d and altloc '%s'" % (start, altloc))
    else :
      clauses.append("resseq %d:%d and altloc '%s'" % (start,end,altloc))
  return clauses

########################################################################
def exercise () :
  from mmtbx.monomer_library import pdb_interpretation
  import cStringIO
  open("tmp.pdb", "w").write("""\
CRYST1   50.800   50.800  155.300  90.00  90.00  90.00 P 43 21 2     8
ATOM      4  N   SER A   1       8.753  29.755  61.685  1.00 49.13
ATOM      5  CA  SER A   1       9.242  30.200  62.974  1.00 46.62
ANISOU    5  CA  SER A   1    343    490   2719    -45   -169    617
ATOM      6  C   SER A   1      10.453  29.500  63.579  1.00 41.99
ATOM      7  O   SER A   1      10.593  29.607  64.814  1.00 43.24
ANISOU    7  O   SER A   1    343    490   2719    -45   -169    617
ATOM      8  CB  SER A   1       8.052  30.189  63.974  1.00 53.00
ATOM      9  OG  SER A   1       7.294  31.409  63.930  1.00 57.79
ATOM     10  N   ARG A   2      11.360  28.819  62.827  1.00 36.48
ATOM     11  CA  ARG A   2      12.548  28.316  63.532  1.00 30.20
ATOM     12  C   ARG A   2      13.502  29.501  63.500  1.00 25.54
ATOM     13  O   ARG A   2      13.730  30.037  62.407  1.00 23.86
ATOM     14  CB  ARG A   2      13.241  27.119  62.861  1.00 27.44
ATOM     15  CG  ARG A   2      12.412  25.849  62.964  1.00 23.66
ATOM     16  CD  ARG A   2      13.267  24.651  63.266  1.00 23.98
ATOM     17  NE  ARG A   2      13.948  24.115  62.135  1.00 22.71
ATOM     18  CZ  ARG A   2      15.114  23.487  62.201  1.00 21.38
ATOM     19  NH1 ARG A   2      15.845  23.331  63.301  1.00 19.34
ATOM     20  NH2 ARG A   2      15.575  23.030  61.051  1.00 26.66
ATOM     21  N   PRO A   3J     13.947  29.997  64.680  1.00 22.94
ATOM     22  CA  PRO A   3J     14.902  31.100  64.827  1.00 20.19
ATOM     23  C   PRO A   3J     16.195  30.718  64.086  1.00 18.44
ATOM     24  O   PRO A   3J     16.545  29.521  64.086  1.00 19.76
ATOM     25  CB  PRO A   3J     15.133  31.218  66.313  1.00 19.17
ATOM     26  CG  PRO A   3J     14.065  30.364  66.951  1.00 15.12
ATOM     27  CD  PRO A   3J     13.816  29.289  65.966  1.00 19.56
ATOM     28  N  AILE A   4      16.953  31.648  63.512  1.00 15.29
ATOM     29  CA AILE A   4      18.243  31.372  62.859  1.00 14.32
ATOM     30  C  AILE A   4      19.233  32.112  63.743  1.00 13.54
ATOM     31  O  AILE A   4      19.105  33.315  64.009  1.00 11.84
ATOM     32  CB AILE A   4      18.298  31.951  61.406  1.00 13.62
ATOM     33  CG1AILE A   4      17.157  31.300  60.620  1.00 18.39
ATOM     34  CG2AILE A   4      19.661  31.747  60.743  1.00 13.64
ATOM     35  CD1AILE A   4      16.879  32.102  59.355  1.00 16.69
ATOM     28  N  BILE A   4      16.953  31.648  63.512  1.00 15.29
ATOM     29  CA BILE A   4      18.243  31.372  62.859  1.00 14.32
ATOM     30  C  BILE A   4      19.233  32.112  63.743  1.00 13.54
ATOM     31  O  BILE A   4      19.105  33.315  64.009  1.00 11.84
ATOM     32  CB BILE A   4      18.298  31.951  61.406  1.00 13.62
ATOM     33  CG1BILE A   4      17.157  31.300  60.620  1.00 18.39
ATOM     34  CG2BILE A   4      19.661  31.747  60.743  1.00 13.64
ATOM1200035  CD1BILE A   4      16.879  32.102  59.355  1.00 16.69
HETATM 1475  S   SO4 S 188      31.424  42.923  60.396  1.00 55.69           S4+
HETATM 1476  O1  SO4 S 188      31.631  41.513  60.336  1.00 59.84           O1-
HETATM 1477  O2  SO4 S 188      32.533  43.699  59.932  1.00 49.98           O1-
HETATM 1478  O3  SO4 S 188      31.128  43.217  61.738  1.00 59.44           O1-
HETATM 1479  O4  SO4 S 188      30.353  43.201  59.539  1.00 60.54           O1-
HETATM 1480  O   HOH W 200      29.478  23.354  61.364  1.00  8.67      WATE
END""")
  out = cStringIO.StringIO()
  processed_pdb_file = pdb_interpretation.run(args=["tmp.pdb"], log=out)
  m = mouse_selection_manager()
  pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  pdb_hierarchy.atoms().reset_i_seq()
  m.update_selection_handlers(
    pdb_hierarchy=pdb_hierarchy,
    mmtbx_selection_function=processed_pdb_file.all_chain_proxies.selection)
  assert m.selection_size() == 0
  m.apply_selection("chain A")
  assert m.selection_size() == 40
  m.clear_selection()
  m.toggle_chain_selection(5)
  assert m.selection_size() == 40
  m.toggle_residue_selection(10)
  assert m.selection_size() == 29
  m.toggle_atom_selection(10) # XXX: doesn't work!
  assert m.selection_size() == 29
  m.toggle_atom_selection(20)
  assert m.selection_size() == 28

  from iotbx import pdb
  from scitbx.array_family import flex
  pdb_hierarchy = pdb.input(source_info=None, lines=flex.split_lines("""\
HETATM 4049  O   HOH W   1       2.954  13.042  11.632  1.00 37.53           O
HETATM 4050  O   HOH W   2       5.539  14.595  10.951  1.00 31.25           O
HETATM 4051  O   HOH W   3      -2.971  14.661  14.669  1.00 38.68           O
HETATM 4052  O   HOH W   4       6.281  34.000   7.684  1.00 39.58           O
HETATM 4053  O   HOH W   5      16.004   9.039  10.335  1.00 37.31           O
HETATM 4054  O   HOH W   6       2.144   5.718  20.447  1.00 49.77           O
HETATM 4055  O   HOH W   7      -1.180  10.517  14.630  1.00 32.95           O
HETATM 4056  O  AHOH W   8       9.227   8.636  12.535  1.00 32.52           O
HETATM 4056  O  BHOH W   8       9.227   8.636  12.535  1.00 32.52           O
HETATM 4057  O  AHOH W   9      11.070  -0.570  15.047  1.00 30.24           O
HETATM 4057  O  BHOH W   9      11.070  -0.570  15.047  1.00 30.24           O
HETATM 4058  O  AHOH W  10      15.630  -6.169  12.853  1.00 31.08           O
HETATM 4058  O  BHOH W  10      15.630  -6.169  12.853  1.00 31.08           O
HETATM 4059  O   HOH W  11      14.854  -8.299  16.887  1.00 32.65           O
HETATM 4060  O   HOH W  12      27.586   0.391  24.184  1.00 31.29           O
HETATM 4061  O   HOH W  13       3.240   7.801  38.401  1.00 32.09           O
END""")).construct_hierarchy()
  m = mouse_selection_manager()
  pdb_hierarchy.atoms().reset_i_seq()
  m.update_selection_handlers(pdb_hierarchy=pdb_hierarchy,
    mmtbx_selection_function=None)
  assert m.selection_size() == 0
  m.apply_selection("chain W")
  assert m.selection_size() == 16
  m.start_range_selection(5)
  m.end_range_selection(15, deselect=False, ignore_altloc=True)
  assert m.selection_size() == 11
  m.start_range_selection(6)
  m.end_range_selection(10, deselect=False, ignore_altloc=True)
  assert m.selection_size() == 11 # no change because of deselect=False
  m.start_range_selection(6)
  m.end_range_selection(9, deselect=True, ignore_altloc=True)
  assert m.selection_size() == 6
  m.start_range_selection(10)
  m.end_range_selection(12, deselect=True, ignore_altloc=False)
  assert m.selection_size() == 5
  print "OK"

if __name__ == "__main__" :
  exercise()

#---end
