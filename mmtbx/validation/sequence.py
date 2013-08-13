from __future__ import division
from libtbx import easy_mp
from libtbx import easy_pickle
from libtbx import str_utils
import libtbx.phil
from libtbx.utils import Sorry, null_out
from libtbx import adopt_init_args, Auto
import sys
import os

master_phil = libtbx.phil.parse("""
  similarity_matrix =  blosum50  dayhoff *identity
    .type = choice
  min_allowable_identity = 0.5
    .type = float
""")

PROTEIN = 0
NUCLEIC_ACID = 1
UNK_AA = "U"
UNK_NA = "K"

def get_mean_coordinate (sites) :
  if (len(sites) == 0) :
    return None
  elif (len(sites) == 1) :
    return sites[0]
  else :
    from scitbx.array_family import flex
    v = flex.vec3_double(sites)
    return v.mean()

class chain (object) :
  """
  Stores information on a protein or nucleic acid chain, its alignment to
  the target sequence, and the coordinates of each residue.  For command-line
  use, much of this is irrelevant.  In the PHENIX GUI, much of this data is
  fed to the sequence/alignment viewer (wxtbx.sequence_view), which controls
  the graphics window(s).
  """
  def __init__ (self, chain_id, sequence, resids, chain_type, sec_str=None,
                resnames=None) :
    adopt_init_args(self, locals())
    assert (chain_type in [PROTEIN, NUCLEIC_ACID])
    assert (len(sequence) == len(resids))
    self.alignment = None
    self.sequence_name = None
    self.sequence_id = None
    self.identity = 0.
    self.n_missing = 0
    self.n_missing_start = 0
    self.n_missing_end = 0
    self.n_gaps = 0
    self.extra = []
    self.unknown = []
    self.mismatch = []
    self._xyz = []
    self._table = None
    self._flag_indices = []

  def set_alignment (self, alignment, sequence_name, sequence_id) :
    assert (len(alignment.a) == len(alignment.b)) and (len(alignment.a) > 0)
    self.alignment = alignment
    self.sequence_name = sequence_name
    self.sequence_id = sequence_id
    if (self.sec_str is not None) :
      raw_sec_str = self.sec_str
    else :
      raw_sec_str = "-" * len(self.sequence)
    self.sec_str = ""
    if (self.sequence_name in [None, ""]) :
      self.sequence_name = "(unnamed)"
    if (alignment is not None) :
      self.identity = alignment.calculate_sequence_identity(skip_chars=['X'])
    i = j = k = 0
    chain_started = False
    prev_char = None
    while (i < len(alignment.a)) :
      a = alignment.a[i]
      b = alignment.b[i]
      resid = None
      if (a != '-') :
        resid = self.resids[j]
        self.sec_str += raw_sec_str[j]
        j += 1
      else :
        self.sec_str += "-"
      i_resid = i - self.n_missing_start
      if (a == 'X') :
        if (not b in ["X","-"]) :
          self.n_missing += 1
          if (not chain_started and
              (self.resnames is None or self.resnames[i_resid] is None)):
            self.n_missing_start += 1
          elif (prev_char != 'X') :
            self.n_gaps += 1
      elif (a == '-') :
        if chain_started:
          if (i_resid) < len(self.resids):
            self.resids.insert(i_resid, None)
            if self.resnames is not None:
              self.resnames.insert(i_resid, None)
        elif (not b in ['-']) :
          self.n_missing += 1
          self.n_missing_start += 1
      if (not a in ["X","-"]) :
        chain_started = True
        if (b == "-") :
          self.extra.append(resid)
          self._flag_indices.append(i)
        elif (not b in ["X","-"]) and (b != a) :
          if (((a == UNK_AA) and (self.chain_type == PROTEIN)) or
              ((a == UNK_NA) and (self.chain_type == NUCLEIC_ACID))) :
            self.unknown.append(resid)
          else :
            self.mismatch.append(resid)
          self._flag_indices.append(i)
      prev_char = a
      i += 1
    i = len(alignment.a) - 1
    matches = alignment.matches()
    while (alignment.a[i] in ["X", "-"]) :
      i_resid = i - self.n_missing_start
      if (matches[i] == " " and
          (self.resnames is None or
           (i_resid) >= len(self.resnames) or
           self.resnames[(i_resid)] is None)):
        self.n_missing_end += 1
      i -= 1
    assert (len(self.sec_str) == len(alignment.a))

  def extract_coordinates (self, pdb_chain) :
    """
    Collect the coordinate of the central atom (CA or P) in each residue,
    padding the array with None so it matches the sequence and resid arrays.
    """
    assert (self.chain_id == pdb_chain.id)
    self._table = []
    k = 0
    for residue_group in pdb_chain.residue_groups() :
      resid = residue_group.resid()
      while (self.resids[k] != resid) :
        self._xyz.append(None)
        k += 1
      res_class = None
      if (resid in self.extra) :
        res_class = "not in sequence"
      elif (resid in self.mismatch) :
        res_class = "mismatch to sequence"
      elif (resid in self.unknown) :
        res_class = "special residue"
      xyz = None
      for atom in residue_group.atoms() :
        if (atom.name == " CA ") or (atom.name == " P  ") :
          xyz = atom.xyz
          break
      else :
        print "WARNING: can't find center of residue %s" % resid
        xyz = residue_group.atoms()[0].xyz
      self._xyz.append(xyz)
      if (res_class is not None) :
        self._table.append([self.chain_id, resid, res_class,
          "chain '%s' and resid %s" % (self.chain_id, resid), xyz])
      k += 1
    while (k < len(self.resids)) :
      self._xyz.append(None)
      k += 1
    assert (len(self.resids) == len(self._xyz))

  def iter_residue_groups(self, pdb_chain):
    assert (self.chain_id == pdb_chain.id)
    k = 0
    for residue_group in pdb_chain.residue_groups():
      resid = residue_group.resid()
      while (self.resids[k] != resid):
        k += 1
        yield None
      yield residue_group
      k += 1
    while k < len(self.resids):
      k += 1
      yield None

  def extract_residue_groups(self, pdb_chain):
    self.residue_groups = list(self.iter_residue_groups(pdb_chain))

  def get_outliers_table (self) :
    """Used in PHENIX validation GUI"""
    return self._table

  def get_highlighted_residues (self) :
    """Used for wxtbx.sequence_view to highlight mismatches, etc."""
    return self._flag_indices

  def get_coordinates_for_alignment_range (self, i1, i2) :
    assert (len(self._xyz) > 0)
    k = 0
    sites = []
    for j, a in enumerate(self.alignment.a) :
      if (j > i2) :
        break
      elif (j >= i1) :
        sites.append(self._xyz[k])
      if (a != '-') :
        k += 1
    return sites

  def get_mean_coordinate_for_alignment_range (self, *args, **kwds) :
    sites = self.get_coordinates_for_alignment_range(*args, **kwds)
    return get_mean_coordinate(sites)

  def get_coordinates_for_alignment_ranges (self, ranges) :
    sites = []
    for i1, i2 in ranges :
      sites.extend(self.get_coordinates_for_alignment_range(i1,i2))
    return sites

  def get_mean_coordinate_for_alignment_ranges (self, *args, **kwds) :
    sites = self.get_coordinates_for_alignment_ranges(*args, **kwds)
    return get_mean_coordinate(sites)

  def get_alignment (self, include_sec_str=False) :
    if (include_sec_str) :
      return [self.alignment.a, self.alignment.b, self.sec_str]
    else :
      return [self.alignment.a, self.alignment.b]

  def show_summary (self, out, verbose=True) :
    def print_resids (resids) :
      assert (len(resids) > 0)
      w = 0
      line = " ".join([ resid.strip() for resid in resids ])
      lines = str_utils.wordwrap(line, 60).split("\n")
      print >> out,   "    residue IDs: %s" % lines[0]
      for line in lines[1:] :
        print >> out, "                 %s" % line
    print >> out, "Chain '%s':" % self.chain_id
    if (self.alignment is None) :
      print >> out, "  No appropriate sequence match found!"
    else :
      print >> out, "  best matching sequence: %s" % self.sequence_name
      print >> out, "  sequence identity: %.2f%%" % (self.identity*100)
      if (self.n_missing > 0) :
        print >> out, "  %d residue(s) missing from PDB chain (%d at start, %d at end)" % (self.n_missing, self.n_missing_start, self.n_missing_end)
      if (self.n_gaps > 0) :
        print >> out, "  %d gap(s) in chain" % self.n_gaps
      if (len(self.mismatch) > 0) :
        print >> out, "  %d mismatches to sequence" % len(self.mismatch)
        print_resids(self.mismatch)
      if (len(self.extra) > 0) :
        print >> out, "  %d residues not found in sequence" % len(self.extra)
        print_resids(self.extra)
      if (len(self.unknown) > 0) :
        print >> out, "  %d residues of unknown type" % len(self.unknown)
        print_resids(self.unknown)
      if (verbose) :
        self.alignment.pretty_print(out=out,
          block_size=60,
          top_name="PDB file",
          bottom_name="sequence",
          show_ruler=False)

class validation (object) :
  def __init__ (self, pdb_hierarchy, sequences, params=None, log=None,
      nproc=Auto, include_secondary_structure=False,
      extract_coordinates=False, extract_residue_groups=False,
      minimum_identity=0) :
    assert (len(sequences) > 0)
    for seq_object in sequences :
      assert (seq_object.sequence != "")
    if (log is None) :
      log = sys.stdout
    if (params is None) :
      params = master_phil.extract()
    self.n_protein = 0
    self.n_rna_dna = 0
    self.n_other = 0
    self.chains = []
    self.minimum_identity = minimum_identity
    self.sequences = sequences
    self.sequence_mappings = [ None ] * len(sequences)
    for i_seq in range(1, len(sequences)) :
      seq_obj1 = sequences[i_seq]
      for j_seq in range(0, len(sequences)) :
        if (j_seq == i_seq) :
          break
        else :
          seq_obj2 = sequences[j_seq]
          if (seq_obj1.sequence == seq_obj2.sequence) :
            self.sequence_mappings[i_seq ] = j_seq
            break
    if (len(pdb_hierarchy.models()) > 1) :
      raise Sorry("Multi-model PDB files not supported.")
    helix_selection = sheet_selection = None
    if (include_secondary_structure) :
      import mmtbx.secondary_structure
      ssm = mmtbx.secondary_structure.manager(
        pdb_hierarchy=pdb_hierarchy,
        sec_str_from_pdb_file=None,
        assume_hydrogens_all_missing=None)
      ssm.find_automatically(log=null_out())
      helix_selection = ssm.alpha_selection()
      sheet_selection = ssm.beta_selection()
    pdb_chains = []
    for pdb_chain in pdb_hierarchy.models()[0].chains() :
      unk = UNK_AA
      chain_id = pdb_chain.id
      main_conf = pdb_chain.conformers()[0]
      if (main_conf.is_na()) :
        self.n_rna_dna += 1
        unk = UNK_NA
        chain_type = NUCLEIC_ACID
      elif (main_conf.is_protein()) :
        self.n_protein += 1
        chain_type = PROTEIN
      else :
        self.n_other += 1
        print >> log, "Skipping non-polymer chain '%s'" % chain_id
        continue
      pad = True
      pad_at_start = False
      seq = pdb_chain.as_padded_sequence(
        substitute_unknown=unk, pad=pad, pad_at_start=pad_at_start)
      resids = pdb_chain.get_residue_ids(pad=pad, pad_at_start=pad_at_start)
      resnames = pdb_chain.get_residue_names_padded(
        pad=pad, pad_at_start=pad_at_start)
      assert (len(seq) == len(resids) == len(resnames))
      sec_str = None
      if (helix_selection is not None) and (main_conf.is_protein()) :
        sec_str = main_conf.as_sec_str_sequence(helix_selection,
          sheet_selection, pad=pad, pad_at_start=pad_at_start)
        assert (len(sec_str) == len(seq))
      c = chain(chain_id=chain_id,
        sequence=seq,
        resids=resids,
        chain_type=chain_type,
        sec_str=sec_str,
        resnames=resnames)
      self.chains.append(c)
      pdb_chains.append(pdb_chain)
    if len(self.chains) == 0:
      raise Sorry("Could not find any polymer chains to align.")
    debug = False
    if debug:
      alignments_and_names = []
      for i in range(len(self.chains)):
        alignments_and_names.append(self.align_chain(i))
    else:
      alignments_and_names = easy_mp.pool_map(
        fixed_func=self.align_chain,
        args=range(len(self.chains)),
        processes=nproc)
    assert (len(alignments_and_names) == len(self.chains) == len(pdb_chains))
    for i, c in enumerate(self.chains) :
      alignment, seq_name, seq_id = alignments_and_names[i]
      if (alignment is None) : continue
      pdb_chain = pdb_chains[i]
      try :
        c.set_alignment(alignment, seq_name, seq_id)
      except Exception, e :
        print "Error processing chain %s" % c.chain_id
        print e
      else :
        if (extract_coordinates) :
          c.extract_coordinates(pdb_chain)
        if extract_residue_groups:
          c.extract_residue_groups(pdb_chain)
    self.sequences = None

  def align_chain (self, i) :
    import mmtbx.alignment
    chain = self.chains[i]
    best_alignment = None
    best_sequence = None
    best_seq_id = None
    best_identity = self.minimum_identity
    best_width = sys.maxint
    for i_seq, seq_object in enumerate(self.sequences) :
      alignment = mmtbx.alignment.align(
        seq_a=chain.sequence,
        seq_b=seq_object.sequence).extract_alignment()
      identity = alignment.calculate_sequence_identity(skip_chars=['X'])
      # if the identities of two alignments are equal, then we prefer the
      # alignment that has the narrowest range for the match
      width = alignment.match_codes.rfind('m') - alignment.match_codes.find('m')
      if ((identity > best_identity) or
          (identity == best_identity and width < best_width)):
        best_identity = identity
        best_alignment = alignment
        best_sequence = seq_object.name
        best_seq_id = i_seq
        best_width = width
    return best_alignment, best_sequence, best_seq_id

  def get_table_data (self) :
    table = []
    for c in self.chains :
      outliers = c.get_outliers_table()
      if (outliers is not None) :
        table.extend(outliers)
    return table

  def get_missing_chains (self) :
    missing = []
    for c in self.chains :
      if (c.alignment is None) :
        missing.append((c.chain_id, c.sequence))
    return missing

  def show (self, out=None) :
    if (out is None) :
      out = sys.stdout
    for chain in self.chains :
      chain.show_summary(out)

  def get_relative_sequence_copy_number (self) :
    """
    Count the number of copies of each sequence within the model, used for
    adjusting the input settings for Phaser-MR in Phenix.  This should
    automatically account for redundancy: only the first matching sequence is
    considered, and sequences which are non-unique will have a copy number of
    -1.  Thus if we have 4 copies of a sequence, and the PDB hierarchy has
    2 matching chains, the copy numbers will be [0.5, -1, -1, -1].
    """
    n_seq  = len(self.sequence_mappings)
    counts = [ 0 ] * n_seq
    for c in self.chains :
      if (c.sequence_id is not None) :
        counts[c.sequence_id] += 1
    redundancies = [ 1 ] * n_seq
    for i_seq in range(n_seq) :
      j_seq = self.sequence_mappings[i_seq]
      if (j_seq is not None) :
        redundancies[j_seq] += 1
        redundancies[i_seq] = 0
    counts_relative = [ 0 ] * n_seq
    for i_seq in range(n_seq) :
      if (redundancies[i_seq] == 0) :
        counts_relative[i_seq] = -1
      else :
        counts_relative[i_seq] = counts[i_seq] / redundancies[i_seq]
    return counts_relative

  def as_cif_block(self, cif_block=None):
    import iotbx.cif.model
    if cif_block is None:
      cif_block = iotbx.cif.model.block()

    struct_ref_loop = iotbx.cif.model.loop(header=(
      "_struct_ref.id",
      "_struct_ref.db_name",
      "_struct_ref.db_code",
      "_struct_ref.pdbx_db_accession",
      "_struct_ref.entity_id",
      "_struct_ref.pdbx_seq_one_letter_code",
      "_struct_ref.pdbx_align_begin",
      "_struct_ref.biol_id",
      #"_struct_ref.seq_align",
      #"_struct_ref.seq_dif",
      #"_struct_ref.details"
    ))

    # maybe we can find this information out somehow and pass it along?
    db_name = ""
    db_code = ""
    db_accession = ""
    for i_chain, chain in enumerate(self.chains):
      struct_ref_loop.add_row((
        i_chain+1, db_name, db_code, db_accession, "?",
        chain.alignment.b, chain.n_missing_start+1, "" ))

    struct_ref_seq_loop = iotbx.cif.model.loop(header=(
      "_struct_ref_seq.align_id",
      "_struct_ref_seq.ref_id",
      "_struct_ref_seq.pdbx_PDB_id_code",
      "_struct_ref_seq.pdbx_strand_id",
      "_struct_ref_seq.seq_align_beg",
      "_struct_ref_seq.pdbx_seq_align_beg_ins_code",
      "_struct_ref_seq.seq_align_end",
      "_struct_ref_seq.pdbx_seq_align_end_ins_code",
      "_struct_ref_seq.pdbx_db_accession",
      "_struct_ref_seq.db_align_beg",
      "_struct_ref_seq.db_align_end",
      "_struct_ref_seq.pdbx_auth_seq_align_beg",
      "_struct_ref_seq.pdbx_auth_seq_align_end"
    ))

    import re
    prog = re.compile("\d+")

    def decode_resid(resid):
      resid = resid.strip()
      s = prog.search(resid)
      assert s is not None
      resseq = resid[s.start():s.end()]
      ins_code = resid[s.end():]
      if len(ins_code) == 0: ins_code = "?"
      return resseq, ins_code

    align_id = 1
    for i_chain, chain in enumerate(self.chains):
      matches = chain.alignment.matches()
      #i_range_begin = 0
      i_range_end = -1
      while True:
        i_range_begin = matches.find("|", i_range_end+1)
        if i_range_begin == -1: break
        i_range_end = matches.find(" ", i_range_begin)
        if i_range_end == -1:
          i_range_end = len(matches)
        i_range_end -= 1
        i_a_range_begin = chain.alignment.i_seqs_a[i_range_begin]
        i_a_range_end = chain.alignment.i_seqs_a[i_range_end]
        i_b_range_begin = chain.alignment.i_seqs_b[i_range_begin]
        i_b_range_end = chain.alignment.i_seqs_b[i_range_end]
        resseq_begin, ins_code_begin = decode_resid(chain.resids[i_a_range_begin])
        resseq_end, ins_code_end = decode_resid(chain.resids[i_a_range_end])
        struct_ref_seq_loop.add_row((
          align_id, i_chain+1, "?", chain.chain_id, i_a_range_begin+1, ins_code_begin,
          i_a_range_end+1, ins_code_end, "?", i_b_range_begin+1, i_b_range_end+1,
          resseq_begin, resseq_end
        ))
        align_id +=1

    cif_block.add_loop(struct_ref_loop)
    cif_block.add_loop(struct_ref_seq_loop)

    return cif_block

# XXX I am not particularly proud of this code
def get_sequence_n_copies (
    sequences,
    pdb_hierarchy,
    force_accept_composition=False,
    copies_from_xtriage=None,
    copies_from_user=Auto,
    minimum_identity=0.3, # okay for MR, may not be suitable in all cases
    assume_xtriage_copies_from_sequence_file=None,
    out=sys.stdout,
    nproc=1) :
  """
  Utility function for reconciling the contents of the sequence file, the
  chains in the model, and the ASU.  Returns the number of copies of the
  sequence file to tell Phaser are present, or raises an error if this is
  either ambiguous or in conflict with the search model multiplicity.
  This is intended to allow the user to specify any combination of inputs -
  for instance, given a tetrameric search model, 2 copies in the ASU, and a
  monomer sequence file, the ASU contains 8 copies of the sequence(s).
  """
  print >> out, "Guessing the number of copies of the sequence file in the ASU"
  v = validation(
    pdb_hierarchy=pdb_hierarchy,
    sequences=sequences,
    log=out,
    nproc=1,
    include_secondary_structure=True,
    extract_coordinates=False,
    minimum_identity=0.3)
  missing = v.get_missing_chains()
  def raise_sorry (msg) :
    raise Sorry(msg + " (Add the parameter force_accept_composition=True to "+
      "disable this error message and continue assuming 1 copy of sequence "+
      "file.)")
  if (len(missing) > 0) :
    error = [
      "%d chain(s) not found in sequence file.  If the sequence does " % \
        len(missing),
      "not accurately represent the contents of the crystal after adjusting ",
      "for copy number, molecular replacement and building may fail.", ]
    if (force_accept_composition) :
      print >> out, "WARNING: %s" % error
    else :
      raise_sorry(error)
  counts = v.get_relative_sequence_copy_number()
  unique_counts = set(counts)
  if (-1 in unique_counts) : unique_counts.remove(-1)
  if (0 in unique_counts) :
    unique_counts.remove(0)
    print >> out, "WARNING: %d sequence(s) not found in model.  This is not" %\
      counts.count(0)
    print >> out, "usually a problem, but it may indicate an incomplete model."
  error = freq = None
  model_copies = copies_from_user
  if (model_copies is Auto) :
    model_copies = copies_from_xtriage
  if (model_copies is None) :
    model_copies = 1
    print >> out, "WARNING: assuming 1 copy of model"
  if (model_copies < 1) :
    raise Sorry("Must have at least one copy of model.")
  if (len(unique_counts) == 0) :
    error = "The sequence file provided does not appear to match the " +\
            "search model."
  elif (len(unique_counts) > 1) :
    error = "The sequence file provided does not map evenly to the "+\
      "search model; this usually means an error in the sequence file (or "+\
      "model) such as accidental duplication or deletion of a chain. "+\
      "[counts: %s]" % list(unique_counts)
  else :
    # XXX float is_integer introduced in Python 2.6
    def is_integer (x) : return (int(x) == x)
    seq_freq = list(unique_counts)[0]
    assert seq_freq > 0
    # 1-1 mapping between PDB hierarchy and sequence
    if (seq_freq == 1.0) :
      print >> out, "Assuming %d copies of sequence file (as well as model)" %\
        model_copies
      return model_copies
    elif (seq_freq > 1) :
      # hierarchy contains N copies of sequence
      if is_integer(seq_freq) :
        freq = int(seq_freq) * model_copies
        print >> out, "Assuming %d copies of sequence file present" % freq
        return freq
      # hierarchy contains X.Y copies of sequence
      else :
        error = "The search model does not appear to contain a round "+\
          "number of copies of the sequence; this usually means "+\
          "that your sequence file has too many or too few sequences.  " +\
          "[freq=%g]" % seq_freq
    # more copies of sequence than in model
    else :
      inverse_freq = 1 / seq_freq
      if is_integer(inverse_freq) :
        # too much sequence
        if (inverse_freq > model_copies) :
          error = "The number of copies of the model expected is less than "+\
            "the number indicated by the sequence file (%d versus %d).  " % \
            (model_copies, int(inverse_freq)) + \
            "This may mean that either the sequence file or the predicted "+\
            "number of copies is wrong."
          if ((copies_from_user is Auto) and
              (copies_from_xtriage is not None) and
              (assume_xtriage_copies_from_sequence_file)) :
            print >> out, error
            print >> out, ""
            print >> out, "Since the number of copies was guessed by Xtriage"
            print >> out, "based on the sequence file, it will be scaled"
            print >> out, "by %g to be appropriate for the search model." % \
              inverse_freq
            return seq_freq
        # too much model
        # XXX is this actually possible the way I've written the function?
        elif (inverse_freq < model_copies) :
          error = "The number of copies of the model expected exceeds the "+\
            "number indicated by the sequence file (%d versus %d)." % \
            (model_copies, int(inverse_freq)) + \
            "This may mean that either the sequence file or the predicted "+\
            "number of copies is wrong."
        # model_copies times model contents matches sequence
        else :
          print >> out, "Assuming only one copy of sequence file contents."
          return 1
      else :
        error = "The sequence file does not appear to contain a round "+\
          "number of copies of the model (%g)." % inverse_freq
  if (error is not None) :
    if (force_accept_composition) :
      print >> out, "WARNING: " + error
      print >> out, ""
      print >> out, "Ambiguous input, defaulting to 1 copy of sequence file"
      return 1
    else :
      raise_sorry(error)
  else :
    raise RuntimeError("No error but frequency not determined")

def get_sequence_n_copies_from_files (seq_file, pdb_file, **kwds) :
  from iotbx import file_reader
  seq_in = file_reader.any_file(seq_file,
    raise_sorry_if_errors=True,
    raise_sorry_if_not_expected_format=True)
  if (seq_in.file_type != "seq") :
    raise Sorry("Can't parse %s as a sequence file.")
  pdb_in = file_reader.any_file(pdb_file,
    raise_sorry_if_errors=True,
    raise_sorry_if_not_expected_format=True)
  if (pdb_in.file_type != "pdb") :
    raise Sorry("Can't parse %s as a PDB or mmCIF file.")
  kwds['pdb_hierarchy'] = pdb_in.file_object.construct_hierarchy()
  kwds['sequences'] = seq_in.file_object
  return get_sequence_n_copies(**kwds)

# XXX test needed
def group_chains_and_sequences (seq_file, pdb_file, **kwds) :
  from iotbx import file_reader
  seq_in = file_reader.any_file(seq_file,
    raise_sorry_if_errors=True,
    raise_sorry_if_not_expected_format=True)
  if (seq_in.file_type != "seq") :
    raise Sorry("Can't parse %s as a sequence file.")
  pdb_in = file_reader.any_file(pdb_file,
    raise_sorry_if_errors=True,
    raise_sorry_if_not_expected_format=True)
  if (pdb_in.file_type != "pdb") :
    raise Sorry("Can't parse %s as a PDB or mmCIF file.")
  kwds['pdb_hierarchy'] = pdb_in.file_object.construct_hierarchy()
  kwds['sequences'] = seq_in.file_object
  v = validation(**kwds)
  chain_to_sequence_mappings = {}
  sequence_to_chain_mappings = {}
  for chain in v.chains :
    seq_id = chain.sequence_id
    chain_id = chain.chain_id
    if (seq_id is None) :
      raise Sorry("Can't map chain %s to a sequence in %s." % (chain_id,
        seq_file))
    sequence = seq_in.file_object[seq_id].sequence
    if (chain_id in chain_to_sequence_mappings) :
      if (chain_to_sequence_mappings[chain_id] != sequence) :
        raise Sorry("Multiple unique chains named '%s'" % chain_id)
    else :
      chain_to_sequence_mappings[chain_id] = sequence
    if (not chain.sequence in sequence_to_chain_mappings) :
      sequence_to_chain_mappings[sequence] = []
    sequence_to_chain_mappings[sequence].append(chain_id)
  return sequence_to_chain_mappings

########################################################################
# REGRESSION TESTING
def exercise () :
  import libtbx.utils
  if (libtbx.utils.detect_multiprocessing_problem() is not None) :
    print "multiprocessing not available, skipping this test"
    return
  if (os.name == "nt"):
    print "easy_mp fixed_func not supported under Windows, skipping this test"
    return
  import iotbx.bioinformatics
  import iotbx.pdb
  from iotbx import file_reader
  import libtbx.load_env # import dependency
  from libtbx.test_utils import Exception_expected, contains_lines, approx_equal
  from cStringIO import StringIO
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
ATOM      2  CA  ARG A  10      -6.299  36.344   7.806  1.00 55.20           C
ATOM     25  CA  TYR A  11      -3.391  33.962   7.211  1.00 40.56           C
ATOM     46  CA  ALA A  12      -0.693  34.802   4.693  1.00 67.95           C
ATOM     56  CA  ALA A  13       0.811  31.422   3.858  1.00 57.97           C
ATOM     66  CA  GLY A  14       4.466  31.094   2.905  1.00 49.24           C
ATOM     73  CA  ALA A  15       7.163  28.421   2.671  1.00 54.70           C
ATOM     83  CA  ILE A  16       6.554  24.685   2.957  1.00 51.79           C
ATOM    102  CA  LEU A  17       7.691  23.612   6.406  1.00 42.30           C
ATOM    121  CA  PTY A  18       7.292  19.882   5.861  1.00 36.68           C
ATOM    128  CA  PHE A  19       5.417  16.968   4.327  1.00 44.99           C
ATOM    148  CA  GLY A  20       3.466  14.289   6.150  1.00 41.99           C
ATOM    155  CA  GLY A  21       1.756  11.130   4.965  1.00 35.77           C
ATOM    190  CA  ALA A  24       1.294  19.658   3.683  1.00 47.02           C
ATOM    200  CA  VAL A  24A      2.361  22.009   6.464  1.00 37.13           C
ATOM    216  CA  HIS A  25       2.980  25.633   5.535  1.00 42.52           C
ATOM    234  CA  LEU A  26       4.518  28.425   7.577  1.00 47.63           C
ATOM    253  CA  ALA A  27       2.095  31.320   7.634  1.00 38.61           C
ATOM    263  CA  ARG A  28       1.589  34.719   9.165  1.00 37.04           C
END""")
  seq1 = iotbx.bioinformatics.sequence("MTTPSHLSDRYELGEILGFGGMSEVHLARD".lower())
  v = validation(
    pdb_hierarchy=pdb_in.construct_hierarchy(),
    sequences=[seq1],
    log=null_out(),
    nproc=1)
  out = StringIO()
  v.show(out=out)
  assert contains_lines(out.getvalue(), """\
  sequence identity: 76.47%
  12 residue(s) missing from PDB chain (9 at start, 1 at end)
  2 gap(s) in chain
  4 mismatches to sequence
    residue IDs:  12 13 15 24""")
  cif_block = v.as_cif_block()
  assert list(cif_block['_struct_ref.pdbx_seq_one_letter_code']) == [
    'MTTPSHLSDRYELGEILGFGGMSEVHLARD']
  assert approx_equal(cif_block['_struct_ref_seq.pdbx_auth_seq_align_beg'],
                      ['10', '14', '16', '19', '24'])
  assert approx_equal(cif_block['_struct_ref_seq.pdbx_auth_seq_align_end'],
                      ['11', '14', '17', '21', '28'])
  assert approx_equal(cif_block['_struct_ref_seq.db_align_beg'],
                      ['10', '14', '16', '19', '25'])
  assert approx_equal(cif_block['_struct_ref_seq.db_align_end'],
                      ['11', '14', '17', '21', '29'])
  assert cif_block['_struct_ref_seq.pdbx_seq_align_beg_ins_code'][4] == 'A'
  seq2 = iotbx.bioinformatics.sequence("MTTPSHLSDRYELGEILGFGGMSEVHLA")
  v = validation(
    pdb_hierarchy=pdb_in.construct_hierarchy(),
    sequences=[seq2],
    log=null_out(),
    nproc=1)
  out = StringIO()
  v.show(out=out)
  assert contains_lines(out.getvalue(), """\
  1 residues not found in sequence
    residue IDs:  28""")
  try :
    v = validation(
      pdb_hierarchy=pdb_in.construct_hierarchy(),
      sequences=[],
      log=null_out(),
      nproc=1)
  except AssertionError :
    pass
  else :
    raise Exception_expected
  cif_block = v.as_cif_block()
  assert list(cif_block['_struct_ref.pdbx_seq_one_letter_code']) == [
    'MTTPSHLSDRYELGEILGFGGMSEVHLA-']
  assert approx_equal(cif_block['_struct_ref_seq.pdbx_auth_seq_align_end'],
                      ['11', '14', '17', '21', '27'])
  assert approx_equal(cif_block['_struct_ref_seq.db_align_end'],
                      ['11', '14', '17', '21', '28'])
  #
  pdb_in2 = iotbx.pdb.input(source_info=None, lines="""\
ATOM      2  CA  ARG A  10      -6.299  36.344   7.806  1.00 55.20           C
ATOM     25  CA  TYR A  11      -3.391  33.962   7.211  1.00 40.56           C
ATOM     46  CA  ALA A  12      -0.693  34.802   4.693  1.00 67.95           C
ATOM     56  CA  ALA A  13       0.811  31.422   3.858  1.00 57.97           C
ATOM     66  CA  GLY A  14       4.466  31.094   2.905  1.00 49.24           C
ATOM     73  CA  ALA A  15       7.163  28.421   2.671  1.00 54.70           C
ATOM     83  CA  ILE A  16       6.554  24.685   2.957  1.00 51.79           C
ATOM    102  CA  LEU A  17       7.691  23.612   6.406  1.00 42.30           C
TER
ATOM   1936  P     G B   2     -22.947 -23.615  15.323  1.00123.20           P
ATOM   1959  P     C B   3     -26.398 -26.111  19.062  1.00110.06           P
ATOM   1979  P     U B   4     -29.512 -30.638  21.164  1.00101.06           P
ATOM   1999  P     C B   5     -30.524 -36.109  21.527  1.00 92.76           P
ATOM   2019  P     U B   6     -28.684 -41.458  21.223  1.00 87.42           P
ATOM   2062  P     G B   8     -18.396 -45.415  21.903  1.00 80.35           P
ATOM   2085  P     A B   9     -13.852 -43.272  24.156  1.00 77.76           P
ATOM   2107  P     G B  10      -8.285 -44.242  26.815  1.00 79.86           P
END
""")
  seq3 = iotbx.bioinformatics.sequence("AGCUUUGGAG")
  v = validation(
    pdb_hierarchy=pdb_in2.construct_hierarchy(),
    sequences=[seq2,seq3],
    log=null_out(),
    nproc=1,
    extract_coordinates=True)
  out = StringIO()
  v.show(out=out)
  cif_block = v.as_cif_block()
  assert approx_equal(cif_block['_struct_ref.pdbx_seq_one_letter_code'],
                      ['MTTPSHLSDRYELGEILGFGGMSEVHLA', 'AGCUUUGGAG'])
  assert approx_equal(cif_block['_struct_ref_seq.pdbx_auth_seq_align_beg'],
                      ['10', '14', '16', '2', '6', '8'])
  assert approx_equal(cif_block['_struct_ref_seq.pdbx_auth_seq_align_end'],
                      ['11', '14', '17', '4', '6', '10'])
  assert (len(v.chains[0].get_outliers_table()) == 3)
  assert (len(v.get_table_data()) == 4)
  assert approx_equal(
    v.chains[0].get_mean_coordinate_for_alignment_range(11,11),
    (-0.693, 34.802, 4.693))
  assert approx_equal(
    v.chains[0].get_mean_coordinate_for_alignment_range(11,14),
    (2.93675, 31.43475, 3.53175))
  assert (v.chains[0].get_highlighted_residues() == [11,12,14])
  assert contains_lines(out.getvalue(), """\
  3 mismatches to sequence
    residue IDs:  12 13 15""")
  assert contains_lines(out.getvalue(), """\
  sequence identity: 87.50%
  2 residue(s) missing from PDB chain (1 at start, 0 at end)
  1 gap(s) in chain
  1 mismatches to sequence
    residue IDs:  5""")
  s = easy_pickle.dumps(v)
  seq4 = iotbx.bioinformatics.sequence("")
  try :
    v = validation(
      pdb_hierarchy=pdb_in2.construct_hierarchy(),
      sequences=[seq4],
      log=null_out(),
      nproc=1,
      extract_coordinates=True)
  except AssertionError :
    pass
  else :
    raise Exception_expected
  # check that nucleic acid chain doesn't get aligned against protein sequence
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
ATOM  18932  P  B DG D   1     -12.183  60.531  25.090  0.50364.79           P
ATOM  18963  P  B DG D   2      -9.738  55.258  20.689  0.50278.77           P
ATOM  18994  P  B DA D   3     -10.119  47.855  19.481  0.50355.17           P
ATOM  19025  P  B DT D   4     -13.664  42.707  21.119  0.50237.06           P
ATOM  19056  P  B DG D   5     -19.510  39.821  21.770  0.50255.45           P
ATOM  19088  P  B DA D   6     -26.096  40.001  21.038  0.50437.49           P
ATOM  19120  P  B DC D   7     -31.790  41.189  18.413  0.50210.00           P
ATOM  19149  P  B DG D   8     -34.639  41.306  12.582  0.50313.99           P
ATOM  19179  P  B DA D   9     -34.987  38.244   6.813  0.50158.92           P
ATOM  19210  P  B DT D  10     -32.560  35.160   1.082  0.50181.38           P
HETATM19241  P  BTSP D  11     -27.614  30.137   0.455  0.50508.17           P
""")
  sequences, _ = iotbx.bioinformatics.fasta_sequence_parse.parse(
    """>4GFH:A|PDBID|CHAIN|SEQUENCE
MSTEPVSASDKYQKISQLEHILKRPDTYIGSVETQEQLQWIYDEETDCMIEKNVTIVPGLFKIFDEILVNAADNKVRDPS
MKRIDVNIHAEEHTIEVKNDGKGIPIEIHNKENIYIPEMIFGHLLTSSNYDDDEKKVTGGRNGYGAKLCNIFSTEFILET
ADLNVGQKYVQKWENNMSICHPPKITSYKKGPSYTKVTFKPDLTRFGMKELDNDILGVMRRRVYDINGSVRDINVYLNGK
SLKIRNFKNYVELYLKSLEKKRQLDNGEDGAAKSDIPTILYERINNRWEVAFAVSDISFQQISFVNSIATTMGGTHVNYI
TDQIVKKISEILKKKKKKSVKSFQIKNNMFIFINCLIENPAFTSQTKEQLTTRVKDFGSRCEIPLEYINKIMKTDLATRM
FEIADANEENALKKSDGTRKSRITNYPKLEDANKAGTKEGYKCTLVLTEGDSALSLAVAGLAVVGRDYYGCYPLRGKMLN
VREASADQILKNAEIQAIKKIMGLQHRKKYEDTKSLRYGHLMIMTDQDHDGSHIKGLIINFLESSFPGLLDIQGFLLEFI
TPIIKVSITKPTKNTIAFYNMPDYEKWREEESHKFTWKQKYYKGLGTSLAQEVREYFSNLDRHLKIFHSLQGNDKDYIDL
AFSKKKADDRKEWLRQYEPGTVLDPTLKEIPISDFINKELILFSLADNIRSIPNVLDGFKPGQRKVLYGCFKKNLKSELK
VAQLAPYVSECTAYHHGEQSLAQTIIGLAQNFVGSNNIYLLLPNGAFGTRATGGKDAAAARYIYTELNKLTRKIFHPADD
PLYKYIQEDEKTVEPEWYLPILPMILVNGAEGIGTGWSTYIPPFNPLEIIKNIRHLMNDEELEQMHPWFRGWTGTIEEIE
PLRYRMYGRIEQIGDNVLEITELPARTWTSTIKEYLLLGLSGNDKIKPWIKDMEEQHDDNIKFIITLSPEEMAKTRKIGF
YERFKLISPISLMNMVAFDPHGKIKKYNSVNEILSEFYYVRLEYYQKRKDHMSERLQWEVEKYSFQVKFIKMIIEKELTV
TNKPRNAIIQELENLGFPRFNKEGKPYYGSPNDEIAEQINDVKGATSDEEDEESSHEDTENVINGPEELYGTYEYLLGMR
IWSLTKERYQKLLKQKQEKETELENLLKLSAKDIWNTDLKAFEVGYQEFLQRDAEAR
>4GFH:D|PDBID|CHAIN|SEQUENCE
GGATGACGATX
""")
  v = validation(
    pdb_hierarchy=pdb_in.construct_hierarchy(),
    sequences=sequences,
    log=null_out(),
    nproc=1,)
  out = StringIO()
  v.show(out=out)
  assert v.chains[0].n_missing == 0
  assert v.chains[0].n_missing_end == 0
  assert v.chains[0].n_missing_start == 0
  assert len(v.chains[0].alignment.matches()) == 11
  #
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
ATOM      2  CA  GLY A   1       1.367   0.551   0.300  1.00  7.71           C
ATOM      6  CA  CYS A   2       2.782   3.785   1.683  1.00  5.18           C
ATOM     12  CA  CYS A   3      -0.375   5.128   3.282  1.00  5.21           C
ATOM     18  CA  SER A   4      -0.870   2.048   5.492  1.00  7.19           C
ATOM     25  CA  LEU A   5       2.786   2.056   6.642  1.00  6.78           C
ATOM     33  CA  PRO A   6       3.212   4.746   9.312  1.00  7.03           C
ATOM     40  CA  PRO A   7       6.870   5.690   8.552  1.00  7.97           C
ATOM     47  CA  CYS A   8       6.021   6.070   4.855  1.00  6.48           C
ATOM     53  CA  ALA A   9       2.812   8.041   5.452  1.00  7.15           C
ATOM     58  CA  LEU A  10       4.739  10.382   7.748  1.00  8.36           C
ATOM     66  CA  SER A  11       7.292  11.200   5.016  1.00  7.00           C
ATOM     73  CA  ASN A  12       4.649  11.435   2.264  1.00  5.40           C
ATOM     81  CA  PRO A  13       1.879  13.433   3.968  1.00  5.97           C
ATOM     88  CA  ASP A  14       0.485  15.371   0.986  1.00  7.70           C
ATOM     96  CA  TYR A  15       0.565  12.245  -1.180  1.00  6.55           C
ATOM    108  CA  CYS A  16      -1.466  10.260   1.363  1.00  7.32           C
ATOM    113  N   NH2 A  17      -2.612  12.308   2.058  1.00  8.11           N
""")
  seq = iotbx.bioinformatics.sequence("GCCSLPPCALSNPDYCX")
  v = validation(
    pdb_hierarchy=pdb_in.construct_hierarchy(),
    sequences=[seq],
    log=null_out(),
    nproc=1,)
  out = StringIO()
  v.show(out=out)
  assert v.chains[0].n_missing == 0
  assert v.chains[0].n_missing_end == 0
  assert v.chains[0].n_missing_start == 0
  assert len(v.chains[0].alignment.matches()) == 17
  #
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
ATOM   2518  CA  PRO C   3      23.450  -5.848  45.723  1.00 85.24           C
ATOM   2525  CA  GLY C   4      20.066  -4.416  44.815  1.00 79.25           C
ATOM   2529  CA  PHE C   5      19.408  -0.913  46.032  1.00 77.13           C
ATOM   2540  CA  GLY C   6      17.384  -1.466  49.208  1.00 83.44           C
ATOM   2544  CA  GLN C   7      17.316  -5.259  49.606  1.00 89.25           C
ATOM   2553  CA  GLY C   8      19.061  -6.829  52.657  1.00 90.67           C
""")
  sequences, _ = iotbx.bioinformatics.fasta_sequence_parse.parse(
    """>1JN5:A|PDBID|CHAIN|SEQUENCE
MASVDFKTYVDQACRAAEEFVNVYYTTMDKRRRLLSRLYMGTATLVWNGNAVSGQESLSEFFEMLPSSEFQISVVDCQPV
HDEATPSQTTVLVVICGSVKFEGNKQRDFNQNFILTAQASPSNTVWKIASDCFRFQDWAS
>1JN5:B|PDBID|CHAIN|SEQUENCE
APPCKGSYFGTENLKSLVLHFLQQYYAIYDSGDRQGLLDAYHDGACCSLSIPFIPQNPARSSLAEYFKDSRNVKKLKDPT
LRFRLLKHTRLNVVAFLNELPKTQHDVNSFVVDISAQTSTLLCFSVNGVFKEVDGKSRDSLRAFTRTFIAVPASNSGLCI
VNDELFVRNASSEEIQRAFAMPAPTPSSSPVPTLSPEQQEMLQAFSTQSGMNLEWSQKCLQDNNWDYTRSAQAFTHLKAK
GEIPEVAFMK
>1JN5:C|PDBID|CHAIN|SEQUENCE
GQSPGFGQGGSV
""")
  v = validation(
    pdb_hierarchy=pdb_in.construct_hierarchy(),
    sequences=sequences,
    log=null_out(),
    nproc=1,)
  out = StringIO()
  v.show(out=out)
  assert v.chains[0].n_missing_start == 3
  assert v.chains[0].n_missing_end == 3
  assert v.chains[0].identity == 1.0
  assert v.chains[0].alignment.match_codes == 'iiimmmmmmiii'
  #
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
ATOM      2  CA  ALA A   2      -8.453  57.214 -12.754  1.00 52.95           C
ATOM      7  CA  LEU A   3      -8.574  59.274  -9.471  1.00 24.33           C
ATOM     15  CA  ARG A   4     -12.178  60.092  -8.575  1.00 28.40           C
ATOM     26  CA  GLY A   5     -14.170  61.485  -5.667  1.00 26.54           C
ATOM     30  CA  THR A   6     -17.784  60.743  -4.783  1.00 31.78           C
ATOM     37  CA  VAL A   7     -19.080  64.405  -4.464  1.00 21.31           C
""")
  seq = iotbx.bioinformatics.sequence("XALRGTV")
  v = validation(
    pdb_hierarchy=pdb_in.construct_hierarchy(),
    sequences=[seq],
    log=null_out(),
    nproc=1,)
  out = StringIO()
  v.show(out=out)
  assert v.chains[0].n_missing_start == 1
  assert v.chains[0].n_missing_end == 0
  assert v.chains[0].identity == 1.0
  assert v.chains[0].alignment.match_codes == 'immmmmm'
  #
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
ATOM   2171  CA  ASP I 355       5.591 -11.903   1.133  1.00 41.60           C
ATOM   2175  CA  PHE I 356       7.082  -8.454   0.828  1.00 39.82           C
ATOM   2186  CA  GLU I 357       5.814  -6.112  -1.877  1.00 41.12           C
ATOM   2195  CA  GLU I 358       8.623  -5.111  -4.219  1.00 42.70           C
ATOM   2199  CA  ILE I 359      10.346  -1.867  -3.363  1.00 43.32           C
ATOM   2207  CA  PRO I 360      11.658   0.659  -5.880  1.00 44.86           C
ATOM   2214  CA  GLU I 361      14.921  -0.125  -7.592  1.00 44.32           C
ATOM   2219  CA  GLU I 362      15.848   3.489  -6.866  1.00 44.27           C
HETATM 2224  CA  TYS I 363      16.482   2.005  -3.448  1.00 44.52           C
""")
  seq = iotbx.bioinformatics.sequence("NGDFEEIPEEYL")
  v = validation(
    pdb_hierarchy=pdb_in.construct_hierarchy(),
    sequences=[seq],
    log=null_out(),
    nproc=1,)
  out = StringIO()
  v.show(out=out)
  assert v.chains[0].n_missing_start == 2
  assert v.chains[0].n_missing_end == 1
  assert v.chains[0].identity == 1.0
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
ATOM    450  CA  ASN A   1      37.242  41.665  44.160  1.00 35.89           C
ATOM    458  CA  GLY A   2      37.796  38.269  42.523  1.00 30.13           C
HETATM  463  CA AMSE A   3      35.878  39.005  39.326  0.54 22.83           C
HETATM  464  CA BMSE A   3      35.892  39.018  39.323  0.46 22.96           C
ATOM    478  CA  ILE A   4      37.580  38.048  36.061  1.00 22.00           C
ATOM    486  CA  SER A   5      37.593  40.843  33.476  1.00 18.73           C
ATOM    819  CA  ALA A   8      25.982  34.781  27.220  1.00 18.43           C
ATOM    824  CA  ALA A   9      23.292  32.475  28.614  1.00 19.60           C
HETATM  830  CA BMSE A  10      22.793  30.814  25.223  0.41 22.60           C
HETATM  831  CA CMSE A  10      22.801  30.850  25.208  0.59 22.54           C
ATOM    845  CA  GLU A  11      26.504  30.054  24.966  1.00 25.19           C
ATOM    854  CA  GLY A  12      25.907  28.394  28.320  1.00 38.88           C
""")
  seq = iotbx.bioinformatics.sequence("NGMISAAAAMEG")
  v = validation(
    pdb_hierarchy=pdb_in.construct_hierarchy(),
    sequences=[seq],
    log=null_out(),
    nproc=1,)
  out = StringIO()
  v.show(out=out)
  assert v.chains[0].alignment.a == 'NGMISXXAAMEG'
  assert v.chains[0].alignment.b == 'NGMISAAAAMEG'
  # all tests below here have additional dependencies
  if (not libtbx.env.has_module("ksdssp")) :
    print "Skipping advanced tests (require ksdssp module)"
    return
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1ywf.pdb",
    test=os.path.isfile)
  if (pdb_file is not None) :
    seq = iotbx.bioinformatics.sequence("MGSSHHHHHHSSGLVPRGSHMAVRELPGAWNFRDVADTATALRPGRLFRSSELSRLDDAGRATLRRLGITDVADLRSSREVARRGPGRVPDGIDVHLLPFPDLADDDADDSAPHETAFKRLLTNDGSNGESGESSQSINDAATRYMTDEYRQFPTRNGAQRALHRVVTLLAAGRPVLTHCFAGKDRTGFVVALVLEAVGLDRDVIVADYLRSNDSVPQLRARISEMIQQRFDTELAPEVVTFTKARLSDGVLGVRAEYLAAARQTIDETYGSLGGYLRDAGISQATVNRMRGVLLG")
    pdb_in = file_reader.any_file(pdb_file, force_type="pdb")
    hierarchy = pdb_in.file_object.construct_hierarchy()
    v = validation(
      pdb_hierarchy=hierarchy,
      sequences=[seq],
      log=null_out(),
      nproc=1,
      include_secondary_structure=True,
      extract_coordinates=True)
    out = StringIO()
    v.show(out=out)
    aln1, aln2, ss = v.chains[0].get_alignment(include_sec_str=True)
    assert ("HHH" in ss) and ("LLL" in ss) and ("---" in ss)
    cif_block = v.as_cif_block()
    assert cif_block['_struct_ref.pdbx_seq_one_letter_code'] == seq.sequence
    assert list(
      cif_block['_struct_ref_seq.pdbx_auth_seq_align_beg']) == ['4', '117']
    assert list(
      cif_block['_struct_ref_seq.pdbx_auth_seq_align_end']) == ['85', '275']
    assert list(cif_block['_struct_ref_seq.seq_align_beg']) == ['1', '114']
    assert list(cif_block['_struct_ref_seq.seq_align_end']) == ['82', '272']
    # determine relative counts of sequences and chains
    n_seq = get_sequence_n_copies(
      pdb_hierarchy=hierarchy,
      sequences=[seq] * 4,
      copies_from_xtriage=4,
      out=null_out())
    assert (n_seq == 1)
    hierarchy = hierarchy.deep_copy()
    chain2 = hierarchy.only_model().chains()[0].detached_copy()
    hierarchy.only_model().append_chain(chain2)
    n_seq = get_sequence_n_copies(
      pdb_hierarchy=hierarchy,
      sequences=[seq] * 4,
      copies_from_xtriage=2,
      out=null_out())
    assert (n_seq == 1)
    n_seq = get_sequence_n_copies(
      pdb_hierarchy=hierarchy,
      sequences=[seq],
      copies_from_xtriage=2,
      out=null_out())
    assert (n_seq == 4)
    try :
      n_seq = get_sequence_n_copies(
        pdb_hierarchy=hierarchy,
        sequences=[seq] * 3,
        copies_from_xtriage=2,
        out=null_out())
    except Sorry, s :
      assert ("round number" in str(s))
    else :
      raise Exception_expected
    n_seq = get_sequence_n_copies(
      pdb_hierarchy=hierarchy,
      sequences=[seq] * 3,
      copies_from_xtriage=2,
      force_accept_composition=True,
      out=null_out())
    assert (n_seq == 1)
    try :
      n_seq = get_sequence_n_copies(
        pdb_hierarchy=hierarchy,
        sequences=[seq] * 4,
        copies_from_xtriage=1,
        out=null_out())
    except Sorry, s :
      assert ("less than" in str(s))
    else :
      raise Exception_expected
    n_seq = get_sequence_n_copies(
      pdb_hierarchy=hierarchy,
      sequences=[seq] * 4,
      copies_from_xtriage=1,
      assume_xtriage_copies_from_sequence_file=True,
      out=null_out())
    assert (n_seq == 0.5)
    hierarchy = hierarchy.deep_copy()
    chain2 = hierarchy.only_model().chains()[0].detached_copy()
    hierarchy.only_model().append_chain(chain2)
    try :
      n_seq = get_sequence_n_copies(
        pdb_hierarchy=hierarchy,
        sequences=[seq] * 2,
        copies_from_xtriage=2,
        out=null_out())
    except Sorry, s :
      assert ("round number" in str(s))
    else :
      raise Exception_expected
    n_seq = get_sequence_n_copies(
      pdb_hierarchy=hierarchy,
      sequences=[seq],
      copies_from_xtriage=1,
      out=null_out())
    assert (n_seq == 3)
    hierarchy = hierarchy.deep_copy()
    chain2 = hierarchy.only_model().chains()[0].detached_copy()
    hierarchy.only_model().append_chain(chain2)
    n_seq = get_sequence_n_copies(
      pdb_hierarchy=hierarchy,
      sequences=[seq] * 2,
      copies_from_xtriage=2,
      out=null_out())
    assert (n_seq == 4)
    # now with files as input
    seq_file = "tmp_mmtbx_validation_sequence.fa"
    open(seq_file, "w").write(">1ywf\n%s" % seq.sequence)
    n_seq = get_sequence_n_copies_from_files(
      pdb_file=pdb_file,
      seq_file=seq_file,
      copies_from_xtriage=4,
      out=null_out())
    try :
      assert (n_seq == 4)
    finally :
      os.remove(seq_file)

if (__name__ == "__main__") :
  exercise()
  print "OK"
