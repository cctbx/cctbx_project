
from libtbx import easy_mp
from libtbx import str_utils
import libtbx.phil
from libtbx.utils import Sorry, null_out
from libtbx import adopt_init_args, Auto
import sys

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

class chain (object) :
  def __init__ (self, chain_id, sequence, resids, chain_type) :
    adopt_init_args(self, locals())
    assert (chain_type in [PROTEIN, NUCLEIC_ACID])
    assert (len(sequence) == len(resids))
    self.alignment = None
    self.sequence_name = None
    self.identity = 0.
    self.n_missing = 0
    self.n_missing_start = 0
    self.n_missing_end = 0
    self.n_gaps = 0
    self.extra = []
    self.unknown = []
    self.mismatch = []

  def set_alignment (self, alignment, sequence_name) :
    assert (len(alignment.a) == len(alignment.b)) and (len(alignment.a) > 0)
    self.alignment = alignment
    self.sequence_name = sequence_name
    if (self.sequence_name in [None, ""]) :
      self.sequence_name = "(unnamed)"
    if (alignment is not None) :
      self.identity = alignment.calculate_sequence_identity(skip_chars=['X'])
    i = j = 0
    chain_started = False
    prev_char = None
    while (i < len(alignment.a)) :
      a = alignment.a[i]
      b = alignment.b[i]
      resid = None
      if (a != '-') :
        resid = self.resids[j]
        j += 1
      if (a == 'X') :
        if (not b in ["X","-"]) :
          self.n_missing += 1
          if (not chain_started) :
            self.n_missing_start += 1
          elif (prev_char != 'X') :
            self.n_gaps += 1
      if (not a in ["X","-"]) :
        chain_started = True
        if (b == "-") :
          self.extra.append(resid)
        elif (not b in ["X","-"]) and (b != a) :
          if (((a == UNK_AA) and (self.chain_type == PROTEIN)) or
              ((a == UNK_NA) and (self.chain_type == NUCLEIC_ACID))) :
            self.unknown.append(resid)
          else :
            self.mismatch.append(resid)
      prev_char = a
      i += 1
    i = len(alignment.a) - 1
    while (alignment.a[i] in ["X", "-"]) :
      self.n_missing_end += 1
      i -= 1

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
      nproc=Auto) :
    assert (len(sequences) > 0)
    if (log is None) :
      log = sys.stdout
    if (params is None) :
      params = master_phil.extract()
    self.n_protein = 0
    self.n_rna_dna = 0
    self.n_other = 0
    self.chains = []
    self.sequences = sequences
    if (len(pdb_hierarchy.models()) > 1) :
      raise Sorry("Multi-model PDB files not supported.")
    for pdb_chain in pdb_hierarchy.models()[0].chains() :
      unk = UNK_AA
      chain_id = pdb_chain.id
      main_conf = pdb_chain.conformers()[0]
      if (main_conf.is_na()) :
        self.n_rna_dna += 1
        unk = UNK_NA
      elif (main_conf.is_protein()) :
        self.n_protein += 1
      else :
        self.n_other += 1
        print >> log, "Skipping non-polymer chain '%s'" % chain_id
        continue
      seq = main_conf.as_padded_sequence(substitute_unknown=unk)
      resids = main_conf.get_residue_ids()
      assert (len(seq) == len(resids))
      c = chain(chain_id=chain_id,
        sequence=seq,
        resids=resids,
        chain_type=PROTEIN)
      self.chains.append(c)
    alignments_and_names = easy_mp.pool_map(
      fixed_func=self.align_chain,
      args=range(len(self.chains)),
      processes=nproc)
    assert (len(alignments_and_names) == len(self.chains))
    for c, (alignment,seq_name) in zip(self.chains, alignments_and_names) :
      try :
        c.set_alignment(alignment, seq_name)
      except Exception, e :
        print "Error processing chain %s" % c.chain_id
        print e

  def align_chain (self, i) :
    import mmtbx.alignment
    chain = self.chains[i]
    seq1 = chain.sequence
    best_alignment = None
    best_sequence = None
    best_identity = 0.
    for seq_object in self.sequences :
      seq2 = seq_object.sequence
      alignment = mmtbx.alignment.align(
        seq_a=chain.sequence,
        seq_b=seq_object.sequence).extract_alignment()
      identity = alignment.calculate_sequence_identity(skip_chars=['X'])
      if (identity > best_identity) :
        best_identity = identity
        best_alignment = alignment
        best_sequence = seq_object.name
    return best_alignment, best_sequence

  def show (self, out=None) :
    if (out is None) :
      out = sys.stdout
    for chain in self.chains :
      chain.show_summary(out)

########################################################################
# REGRESSION TESTING
def exercise () :
  import iotbx.bioinformatics
  import iotbx.pdb
  from libtbx.test_utils import Exception_expected, contains_lines
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
  seq1 = iotbx.bioinformatics.sequence("MTTPSHLSDRYELGEILGFGGMSEVHLARD")
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
    log=None,#null_out(),
    nproc=1)
  out = StringIO()
  v.show(out=out)
  assert contains_lines(out.getvalue(), """\
  3 mismatches to sequence
    residue IDs:  12 13 15""")
  assert contains_lines(out.getvalue(), """\
  sequence identity: 87.50%
  2 residue(s) missing from PDB chain (1 at start, 0 at end)
  1 gap(s) in chain
  1 mismatches to sequence
    residue IDs:  5""")

if (__name__ == "__main__") :
  exercise()
  print "OK"
