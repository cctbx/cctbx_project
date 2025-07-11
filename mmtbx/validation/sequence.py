from __future__ import absolute_import, division, print_function
from libtbx import easy_mp
from libtbx import str_utils
import libtbx.phil
from libtbx.utils import Sorry
from libtbx.test_utils import approx_equal_core
from libtbx import adopt_init_args, Auto
from mmtbx.alignment import align
import sys

import iotbx.cif.model
from iotbx.pdb import modified_aa_names, modified_rna_dna_names
from iotbx.pdb.amino_acid_codes import one_letter_given_three_letter, \
    three_letter_l_given_three_letter_d, three_letter_given_one_letter
from six.moves import zip
from six.moves import range

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

def get_mean_coordinate(sites):
  if (len(sites) == 0):
    return None
  elif (len(sites) == 1):
    return sites[0]
  else :
    from scitbx.array_family import flex
    v = flex.vec3_double(sites)
    return v.mean()

class chain(object):
  """
  Stores information on a protein or nucleic acid chain, its alignment to
  the target sequence, and the coordinates of each residue.  For command-line
  use, much of this is irrelevant.  In the PHENIX GUI, much of this data is
  fed to the sequence/alignment viewer (wxtbx.sequence_view), which controls
  the graphics window(s).
  """
  def __init__(self, chain_id, sequence, resids, chain_type, sec_str=None,
                resnames=None):
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
    self.actual_code = []
    self._xyz = []
    self._table = None
    self._flag_indices = []

  def set_alignment(self, alignment, sequence_name, sequence_id):
    assert (len(alignment.a) == len(alignment.b)) and (len(alignment.a) > 0)
    self.alignment = alignment
    self.sequence_name = sequence_name
    self.sequence_id = sequence_id
    if (self.sec_str is not None):
      raw_sec_str = list(self.sec_str)
    else :
      raw_sec_str = list("-" * len(self.sequence))
    self.sec_str = ""
    if (self.sequence_name in [None, ""]):
      self.sequence_name = "(unnamed)"
    if (alignment is not None):
      self.identity = alignment.calculate_sequence_identity(skip_chars=['X'])
    i_aln = i_resid = 0
    chain_started = False
    prev_char = None
    alignment_length = len(alignment.a)
    while (i_aln < alignment_length):
      symbol_pdb = alignment.a[i_aln]
      symbol_seq = alignment.b[i_aln]
      resid = None
      if (symbol_pdb != '-'):
        resid = self.resids[i_resid]
        self.sec_str += raw_sec_str[i_resid]
        i_resid += 1
      else :
        self.sec_str += "-"
      j_resid = i_aln - self.n_missing_start
      if (symbol_pdb == 'X'):
        if (not symbol_seq in ["X","-"]):
          self.n_missing += 1
          if (not chain_started and
              (self.resnames is None or self.resnames[j_resid] is None)):
            self.n_missing_start += 1
          elif (prev_char != 'X'):
            self.n_gaps += 1
      elif (symbol_pdb == '-'):
        if chain_started:
          if (j_resid < len(self.resids)):
            self.resids.insert(j_resid, None)
            raw_sec_str.insert(j_resid, "-")
            if self.resnames is not None:
              self.resnames.insert(j_resid, None)
            i_resid += 1
        elif (not symbol_seq in ['-']):
          self.n_missing += 1
          self.n_missing_start += 1
      if (not symbol_pdb in ["X","-"]):
        chain_started = True
        if (symbol_seq == "-"):
          self.extra.append(resid)
          self._flag_indices.append(i_aln)
        elif (not symbol_seq in ["X","-"]) and (symbol_seq != symbol_pdb):
          if (((symbol_pdb == UNK_AA) and (self.chain_type == PROTEIN)) or
              ((symbol_pdb == UNK_NA) and (self.chain_type == NUCLEIC_ACID))):
            self.unknown.append(resid)
          else :
            self.mismatch.append(resid)
            self.actual_code.append(symbol_seq)
          self._flag_indices.append(i_aln)
      prev_char = symbol_pdb
      i_aln += 1
    i_aln = alignment_length - 1
    matches = alignment.matches()
    while (alignment.a[i_aln] in ["X", "-"]):
      i_resid = i_aln - self.n_missing_start
      if (matches[i_aln] == " " and
          (self.resnames is None or
           (i_resid >= len(self.resnames) or
           self.resnames[i_resid] is None))):
        self.n_missing_end += 1
        self.n_missing += 1
      i_aln -= 1
      if i_aln < 0:
        break
    assert (len(self.sec_str) == len(alignment.a))

  def extract_coordinates(self, pdb_chain):
    """
    Collect the coordinate of the central atom (CA or P) in each residue,
    padding the array with None so it matches the sequence and resid arrays.
    """
    assert (self.chain_id == pdb_chain.id)
    self._table = []
    k = 0
    for residue_group in pdb_chain.residue_groups():
      resid = residue_group.resid()
      while (self.resids[k] != resid):
        self._xyz.append(None)
        k += 1
      res_class = None
      if (resid in self.extra):
        res_class = "not in sequence"
      elif (resid in self.mismatch):
        res_class = "mismatch to sequence"
      elif (resid in self.unknown):
        res_class = "special residue"
      xyz = None
      for atom in residue_group.atoms():
        if (atom.name == " CA ") or (atom.name == " P  "):
          xyz = atom.xyz
          break
      else :
        print("WARNING: can't find center of residue %s" % resid)
        xyz = residue_group.atoms()[0].xyz
      self._xyz.append(xyz)
      if (res_class is not None):
        self._table.append([self.chain_id, resid, res_class,
          "chain '%s' and resid %s" % (self.chain_id, resid), xyz])
      k += 1
    while (k < len(self.resids)):
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

  def get_outliers_table(self):
    """Used in PHENIX validation GUI"""
    return self._table

  def get_highlighted_residues(self):
    """Used for wxtbx.sequence_view to highlight mismatches, etc."""
    return self._flag_indices

  def get_coordinates_for_alignment_range(self, i1, i2):
    assert (len(self._xyz) > 0)
    k = 0
    sites = []
    for j, a in enumerate(self.alignment.a):
      if (j > i2):
        break
      elif (j >= i1):
        sites.append(self._xyz[k])
      if (a != '-'):
        k += 1
    return sites

  def get_mean_coordinate_for_alignment_range(self, *args, **kwds):
    sites = self.get_coordinates_for_alignment_range(*args, **kwds)
    return get_mean_coordinate(sites)

  def get_coordinates_for_alignment_ranges(self, ranges):
    sites = []
    for i1, i2 in ranges :
      sites.extend(self.get_coordinates_for_alignment_range(i1,i2))
    return sites

  def get_mean_coordinate_for_alignment_ranges(self, *args, **kwds):
    sites = self.get_coordinates_for_alignment_ranges(*args, **kwds)
    return get_mean_coordinate(sites)

  def get_alignment(self, include_sec_str=False):
    if (include_sec_str):
      return [self.alignment.a, self.alignment.b, self.sec_str]
    else :
      return [self.alignment.a, self.alignment.b]

  def show_summary(self, out, verbose=True):
    def print_resids(resids):
      assert (len(resids) > 0)
      w = 0
      line = " ".join([ resid.strip() for resid in resids ])
      lines = str_utils.wordwrap(line, 60).split("\n")
      print("    residue IDs: %s" % lines[0], file=out)
      for line in lines[1:] :
        print("                 %s" % line, file=out)
    print("Chain '%s':" % self.chain_id, file=out)
    if (self.alignment is None):
      print("  No appropriate sequence match found!", file=out)
    else :
      print("  best matching sequence: %s" % self.sequence_name, file=out)
      print("  sequence identity: %.2f%%" % (self.identity*100), file=out)
      if (self.n_missing > 0):
        print("  %d residue(s) missing from PDB chain (%d at start, %d at end)" % (self.n_missing, self.n_missing_start, self.n_missing_end), file=out)
      if (self.n_gaps > 0):
        print("  %d gap(s) in chain" % self.n_gaps, file=out)
      if (len(self.mismatch) > 0):
        print("  %d mismatches to sequence" % len(self.mismatch), file=out)
        print_resids(self.mismatch)
      if (len(self.extra) > 0):
        print("  %d residues not found in sequence" % len(self.extra), file=out)
        print_resids(self.extra)
      if (len(self.unknown) > 0):
        print("  %d residues of unknown type" % len(self.unknown), file=out)
        print_resids(self.unknown)
      if (verbose):
        self.alignment.pretty_print(out=out,
          block_size=60,
          top_name="PDB file",
          bottom_name="sequence",
          show_ruler=False)

class validation(object):
  def __init__(self, pdb_hierarchy, sequences, params=None, log=None,
      nproc=Auto, include_secondary_structure=False,
      extract_coordinates=False, extract_residue_groups=False,
      minimum_identity=0, custom_residues=[], ignore_hetatm=False):
    #
    # XXX This is to stop assertion crash that checks lengths of provided
    # XXX sequence with the length of sequence from model (which can happen to
    # XXX be different due to presence of altlocs!
    # XXX Also this assumes altlocs within residue groups are the same resname's
    # XXX And of course making a copy is a bad idea, obviously!
    # XXX No test.
    #
    pdb_hierarchy = pdb_hierarchy.deep_copy()
    pdb_hierarchy.remove_alt_confs(always_keep_one_conformer=True)
    pdb_hierarchy.atoms().reset_i_seq()
    #
    assert (len(sequences) > 0)
    for seq_object in sequences :
      assert (seq_object.sequence != "")
    if (log is None):
      log = sys.stdout
    if (params is None):
      params = master_phil.extract()
    self.n_protein = 0
    self.n_rna_dna = 0
    self.n_other = 0
    self.chains = []
    self.minimum_identity = minimum_identity
    self.sequences = sequences
    self.custom_residues = custom_residues
    if self.custom_residues is None:
      self.custom_residues = list()
    self.sequence_mappings = [ None ] * len(sequences)
    for i_seq in range(1, len(sequences)):
      seq_obj1 = sequences[i_seq]
      for j_seq in range(0, len(sequences)):
        if (j_seq == i_seq):
          break
        else :
          seq_obj2 = sequences[j_seq]
          if (seq_obj1.sequence == seq_obj2.sequence):
            self.sequence_mappings[i_seq ] = j_seq
            break
    if (len(pdb_hierarchy.models()) > 1):
      raise Sorry("Multi-model PDB files not supported.")
    helix_selection = sheet_selection = None
    if (include_secondary_structure):
      import mmtbx.secondary_structure
      ssm = mmtbx.secondary_structure.manager(
        pdb_hierarchy=pdb_hierarchy,
        sec_str_from_pdb_file=None)
      helix_selection = ssm.helix_selection()
      sheet_selection = ssm.beta_selection()
    pdb_chains = []
    chain_seq = {}
    for pdb_chain in pdb_hierarchy.models()[0].chains():
      unk = UNK_AA
      chain_id = pdb_chain.id
      main_conf = pdb_chain.conformers()[0]
      if (main_conf.is_na()):
        self.n_rna_dna += 1
        unk = UNK_NA
        chain_type = NUCLEIC_ACID
      elif (main_conf.is_protein()):
        self.n_protein += 1
        chain_type = PROTEIN
      else :
        self.n_other += 1
        print("Skipping non-polymer chain '%s'" % chain_id, file=log)
        continue
      pad = True
      pad_at_start = False
      seq = pdb_chain.as_padded_sequence(
        substitute_unknown=unk, pad=pad, pad_at_start=pad_at_start, ignore_hetatm=ignore_hetatm)
      chain_seq[chain_id] = seq
      resids = pdb_chain.get_residue_ids(pad=pad, pad_at_start=pad_at_start,
        ignore_hetatm=ignore_hetatm)
      resnames = pdb_chain.get_residue_names_padded(
        pad=pad, pad_at_start=pad_at_start, ignore_hetatm=ignore_hetatm)
      assert (len(seq) == len(resids) == len(resnames))
      sec_str = None
      if (helix_selection is not None) and (main_conf.is_protein()):
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
      try:
        alignments_and_names = easy_mp.pool_map(
          fixed_func=self.align_chain,
          args=range(len(self.chains)),
          processes=nproc)
      except Exception as e:
        # retry without mp
        print("Failed to get alignments with multiprocessing, retrying")
        alignments_and_names = []
        for i in range(len(self.chains)):
          alignments_and_names.append(self.align_chain(i))
    assert (len(alignments_and_names) == len(self.chains) == len(pdb_chains))
    for i, c in enumerate(self.chains):
      alignment, seq_name, seq_id = alignments_and_names[i]
      # if no alignment was found, just use the sequence from the model
      if alignment is None:
        alignment = align(chain_seq[c.chain_id], chain_seq[c.chain_id]).extract_alignment()
        seq_name = 'model'
        seq_id = 0
      pdb_chain = pdb_chains[i]
      try :
        c.set_alignment(alignment, seq_name, seq_id)
      except Exception as e :
        print("Error processing chain %s" % c.chain_id)
        raise
        print(e)
      else :
        if (extract_coordinates):
          c.extract_coordinates(pdb_chain)
        if extract_residue_groups:
          c.extract_residue_groups(pdb_chain)
    self.sequences = None

  def align_chain(self, i):
    import mmtbx.alignment
    chain = self.chains[i]
    best_alignment = None
    best_sequence = None
    best_seq_id = None
    best_identity = self.minimum_identity
    best_width = sys.maxsize
    best_length = sys.maxsize
    for i_seq, seq_object in enumerate(self.sequences):
      alignment = mmtbx.alignment.align(
        seq_a=chain.sequence,
        seq_b=seq_object.sequence).extract_alignment()
      identity = alignment.calculate_sequence_identity(skip_chars=['X'])
      # if the identities of two alignments are equal, then we prefer the
      # alignment that has the narrowest range for the match and the
      # shortest sequence
      width = alignment.match_codes.rfind('m') - alignment.match_codes.find('m')
      length = len(seq_object.sequence)
      if ((identity > best_identity) or
          (approx_equal_core(identity, best_identity, 1.e-6, 1.e10, None, "")
           and width <= best_width and length < best_length)):
        best_identity = identity
        best_alignment = alignment
        best_sequence = seq_object.name
        best_seq_id = i_seq
        best_width = width
        best_length = length
    return best_alignment, best_sequence, best_seq_id

  def get_table_data(self):
    table = []
    for c in self.chains :
      outliers = c.get_outliers_table()
      if (outliers is not None):
        table.extend(outliers)
    return table

  def get_missing_chains(self):
    missing = []
    for c in self.chains :
      if (c.alignment is None):
        missing.append((c.chain_id, c.sequence))
    return missing

  def show(self, out=None):
    if (out is None):
      out = sys.stdout
    for chain in self.chains :
      chain.show_summary(out)

  def get_relative_sequence_copy_number(self):
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
      if (c.sequence_id is not None):
        counts[c.sequence_id] += 1
    redundancies = [ 1 ] * n_seq
    for i_seq in range(n_seq):
      j_seq = self.sequence_mappings[i_seq]
      if (j_seq is not None):
        redundancies[j_seq] += 1
        redundancies[i_seq] = 0
    counts_relative = [ 0 ] * n_seq
    for i_seq in range(n_seq):
      if (redundancies[i_seq] == 0):
        counts_relative[i_seq] = -1
      else :
        counts_relative[i_seq] = counts[i_seq] / redundancies[i_seq]
    return counts_relative

  def sequence_as_cif_block(self, custom_residues=None):
    """
    Export sequence information as mmCIF block
    Version 5.0 of mmCIF/PDBx dictionary

    Parameters
    ----------
    custom_residues: list of str
      List of custom 3-letter residues to keep in pdbx_one_letter_sequence
      The 3-letter residue must exist in the model. If None, the value
      from self.custom_residues is used.

    Returns
    -------
    cif_block: iotbx.cif.model.block
    """

    if custom_residues is None:
      custom_residues = self.custom_residues

    dna = set(['DA', 'DT', 'DC', 'DG', 'DI'])
    rna = set(['A', 'U', 'C', 'G'])
    rna_to_dna = {'A': 'DA', 'U': 'DT', 'T': 'DT', 'C': 'DC', 'G': 'DG',
                  'I': 'DI'}
    modified_dna = set()
    modified_rna = set()
    for key in modified_rna_dna_names.lookup.keys():
      value = modified_rna_dna_names.lookup[key]
      if value in dna:
        modified_dna.add(key)
      elif value in rna:
        modified_rna.add(key)

    # http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/entity.html
    entity_loop = iotbx.cif.model.loop(header=(
      '_entity.id',
      '_entity.pdbx_description'
    ))

    # http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/entity_poly.html
    entity_poly_loop = iotbx.cif.model.loop(header=(
      '_entity_poly.entity_id',
      '_entity_poly.nstd_linkage',
      '_entity_poly.nstd_monomer',
      '_entity_poly.pdbx_seq_one_letter_code',
      '_entity_poly.pdbx_seq_one_letter_code_can',
      '_entity_poly.pdbx_strand_id',
      '_entity_poly.pdbx_target_identifier',
      '_entity_poly.type',
    ))

    # http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/entity_poly_seq.html
    entity_poly_seq_loop = iotbx.cif.model.loop(header=(
      '_entity_poly_seq.entity_id',
      '_entity_poly_seq.num',
      '_entity_poly_seq.mon_id',
      '_entity_poly_seq.hetero',
    ))

    # http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/struct_ref.html
    struct_ref_loop = iotbx.cif.model.loop(header=(
      '_struct_ref.id',
      '_struct_ref.db_code',
      '_struct_ref.db_name',
      '_struct_ref.entity_id',
      '_struct_ref.pdbx_align_begin',
      '_struct_ref.pdbx_db_accession',
      '_struct_ref.pdbx_db_isoform',
      '_struct_ref.pdbx_seq_one_letter_code',
    ))

    # http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/struct_ref_seq.html
    struct_ref_seq_loop = iotbx.cif.model.loop(header=(
      '_struct_ref_seq.align_id',
      '_struct_ref_seq.db_align_beg',
      '_struct_ref_seq.db_align_end',
      '_struct_ref_seq.pdbx_PDB_id_code',
      '_struct_ref_seq.pdbx_auth_seq_align_beg',
      '_struct_ref_seq.pdbx_auth_seq_align_end',
      '_struct_ref_seq.pdbx_db_accession',
      '_struct_ref_seq.pdbx_db_align_beg_ins_code',
      '_struct_ref_seq.pdbx_db_align_end_ins_code',
      '_struct_ref_seq.pdbx_seq_align_beg_ins_code',
      '_struct_ref_seq.pdbx_seq_align_end_ins_code',
      '_struct_ref_seq.pdbx_strand_id',
      '_struct_ref_seq.ref_id',
      '_struct_ref_seq.seq_align_beg',
      '_struct_ref_seq.seq_align_end',
    ))

    entity_id = 0
    # entity_poly
    sequence_to_entity_id = dict()
    nstd_linkage = dict()
    nstd_monomer = dict()
    seq_one_letter_code = dict()
    seq_one_letter_code_can = dict()
    strand_id = dict()
    target_identifier = dict()
    sequence_type = dict()
    # entity_poly_seq
    num = dict()
    mon_id = dict()
    hetero = dict()
    # struct_ref (work in progress)
    chain_id = dict()
    db_code = '?'
    db_name = '?'
    align_begin = '?'
    db_accession = '?'
    db_isoform = '?'
    # struct_ref_seq (work in progress)
    db_align_beg = '?'
    db_align_end = '?'
    PDB_id_code = '?'
    align_beg_ins_code = '?'
    align_end_ins_code = '?'

    for i_chain, chain in enumerate(self.chains):
      seq_can = chain.alignment.b
      # entity_id
      if seq_can not in sequence_to_entity_id:
        entity_id += 1
        sequence_to_entity_id[seq_can] = entity_id
      else:
        # subsequent matches just add strand_id
        entity_id = sequence_to_entity_id[seq_can]
        strand_id[entity_id].append(chain.chain_id)
        continue

      # entity_poly items
      # nstd_linkage (work in progress)
      if entity_id not in nstd_linkage:
        nstd_linkage[entity_id] = 'no'
      # nstd_monomer
      if entity_id not in nstd_monomer:
        nstd_monomer[entity_id] = 'no'
      # pdbx_seq_one_letter_code
      if entity_id not in seq_one_letter_code:
        seq_one_letter_code[entity_id] = list()
      # type (work in progress)
      if entity_id not in sequence_type:
        sequence_type[entity_id] = '?'
      has_protein = False
      has_rna = False
      has_dna = False
      has_sugar = False
      is_d = False
      # chain.alignment.a is the model
      # chain.alignment.b is the sequence
      for i_a, i_b in zip(chain.alignment.i_seqs_a, chain.alignment.i_seqs_b):
        # sequence does not have residue in model
        if i_b is None:
          continue
        # model does not have residue in sequence
        if i_a is None or chain.resnames[i_a] is None:
          letter = seq_can[i_b]
        else:
          resname = chain.resnames[i_a].strip()
          # check for modified residues
          if (resname in modified_aa_names.lookup or
              resname in modified_rna_dna_names.lookup or
              resname in custom_residues):
            letter = '({resname})'.format(resname=resname)
            nstd_monomer[entity_id] = 'yes'
          elif resname in three_letter_l_given_three_letter_d:
            letter = '({resname})'.format(resname=resname)
            nstd_monomer[entity_id] = 'yes'
          # check for nucleic acid
          elif resname in dna:
            letter = '({resname})'.format(resname=resname)
          elif resname in rna:
            letter = resname
          # regular protein
          else:
            letter = one_letter_given_three_letter.get(resname)
            if letter is None:
              letter = 'X'

          # check for protein
          if (resname in one_letter_given_three_letter or
              resname in modified_aa_names.lookup):
            has_protein = True
          # check for DNA
          # hybrid protein/DNA/RNA chains are not allowed
          if resname in dna or resname in modified_dna:
            has_dna = True
            has_protein = False
          # check for RNA
          # does not handle hybrid DNA/RNA chains
          if resname in rna or resname in modified_rna:
            has_rna = True
            has_dna = False
            has_protein = False
          # check chirality
          # hybrid D/L handed chains are not allowed
          if resname in three_letter_l_given_three_letter_d:
            is_d = True
        # pdbx_seq_one_letter_code
        seq_one_letter_code[entity_id].append(letter)

      # pdbx_seq_one_letter_code_can
      seq_one_letter_code_can[entity_id] = seq_can.replace('-', '')
      # strand_id
      if entity_id not in strand_id:
        strand_id[entity_id] = list()
      strand_id[entity_id].append(chain.chain_id)
      # target_identifier (work in progress)
      if entity_id not in target_identifier:
        target_identifier[entity_id] = '?'
      # type
      #   polypeptide(L)
      #   polypeptide(D)
      #   polydeoxyribonucleotide,
      #   polyribonucleotide
      # missing
      #   cyclic-psuedo-peptide
      #   other
      #   peptide nucleic acid
      #   polydeoxyribonucleotide/polyribonucleotide
      #   polysaccharide(D)
      #   polysaccahride(L)
      if has_protein:
        choice = 'polypeptide'
        if is_d:
          choice += '(D)'
        else:
          choice += '(L)'
      if has_dna:
        choice = 'polydeoxyribonucleotide'
      if has_rna:
        choice = 'polyribonucleotide'
      sequence_type[entity_id] = choice

      # entity_poly_seq items
      if entity_id not in mon_id:
        mon_id[entity_id] = list()
      if entity_id not in num:
        num[entity_id] = list()
      if entity_id not in hetero:
        hetero[entity_id] = list()

      # struct_ref items
      if entity_id not in chain_id:
        chain_id[entity_id] = i_chain + 1

      for i_a, i_b in zip(chain.alignment.i_seqs_a, chain.alignment.i_seqs_b):
        # sequence does not have residue in model
        if i_b is None:
          continue
        seq_resname = None
        if has_protein:
          seq_resname = three_letter_given_one_letter.get(seq_can[i_b])
        if has_dna:
          seq_resname = rna_to_dna.get(seq_can[i_b])
        if has_rna:
          seq_resname = seq_can[i_b]
        if seq_resname is None:
          seq_resname = 'UNK'
        # model does not have residue in sequence
        if i_a is None or chain.resnames[i_a] is None:
          resname = seq_resname
        else:
          resname = chain.resnames[i_a]
        mon_id[entity_id].append(resname.strip())
        if len(num[entity_id]) == 0:
          num[entity_id].append(1)
        else:
          num[entity_id].append(num[entity_id][-1] + 1)
        hetero[entity_id].append('no')

    # build loops
    ids = list(sequence_to_entity_id.values())
    ids.sort()
    align_id = 1
    for entity_id in ids:
      # construct entity_poly loop
      if len(strand_id[entity_id]) == 1:
        chains = strand_id[entity_id][0]
      else:
        chains = strand_id[entity_id]
        #chains.sort()
        chains = ','.join(chains)
      entity_poly_loop.add_row((
        entity_id,
        nstd_linkage[entity_id],
        nstd_monomer[entity_id],
        ';' + ''.join(seq_one_letter_code[entity_id]) + '\n;',
        ';' + seq_one_letter_code_can[entity_id] + '\n;',
        chains,
        target_identifier[entity_id],
        sequence_type[entity_id]
      ))

      # construct entity loop
      entity_loop.add_row((
        entity_id,
        'Chains: ' + chains
      ))

      # construct entity_poly_seq loop
      chain_length = len(mon_id[entity_id])
      for i in range(chain_length):
        entity_poly_seq_loop.add_row((
          entity_id,
          num[entity_id][i],
          mon_id[entity_id][i],
          hetero[entity_id][i]
        ))

      # construct struct_ref loop
      struct_ref_loop.add_row((
        chain_id[entity_id],
        db_code,
        db_name,
        entity_id,
        align_begin,
        db_accession,
        db_isoform,
        ';' + seq_one_letter_code_can[entity_id] + '\n;'
      ))

      # construct struct_ref_seq loop
      for chain in strand_id[entity_id]:
        struct_ref_seq_loop.add_row((
          align_id,
          db_align_beg,
          db_align_end,
          PDB_id_code,
          '1',
          len(seq_one_letter_code_can[entity_id]) - 1,
          db_accession,
          align_beg_ins_code,
          align_end_ins_code,
          align_beg_ins_code,
          align_end_ins_code,
          chain,
          chain_id[entity_id],
          '1',
          len(seq_one_letter_code_can[entity_id]) - 1
        ))
        align_id += 1

    # construct block
    cif_block = iotbx.cif.model.block()
    cif_block.add_loop(entity_loop)
    cif_block.add_loop(entity_poly_loop)
    cif_block.add_loop(entity_poly_seq_loop)
    cif_block.add_loop(struct_ref_loop)
    cif_block.add_loop(struct_ref_seq_loop)

    return cif_block

# XXX I am not particularly proud of this code
def get_sequence_n_copies(
    sequences,
    pdb_hierarchy,
    force_accept_composition=False,
    copies_from_xtriage=None,
    copies_from_user=Auto,
    minimum_identity=0.3, # okay for MR, may not be suitable in all cases
    assume_xtriage_copies_from_sequence_file=None,
    out=sys.stdout,
    nproc=1):
  """
  Utility function for reconciling the contents of the sequence file, the
  chains in the model, and the ASU.  Returns the number of copies of the
  sequence file to tell Phaser are present, or raises an error if this is
  either ambiguous or in conflict with the search model multiplicity.
  This is intended to allow the user to specify any combination of inputs -
  for instance, given a tetrameric search model, 2 copies in the ASU, and a
  monomer sequence file, the ASU contains 8 copies of the sequence(s).
  """
  print("Guessing the number of copies of the sequence file in the ASU", file=out)
  v = validation(
    pdb_hierarchy=pdb_hierarchy,
    sequences=sequences,
    log=out,
    nproc=1,
    include_secondary_structure=True,
    extract_coordinates=False,
    minimum_identity=0.3)
  missing = v.get_missing_chains()
  def raise_sorry(msg):
    raise Sorry(msg + " (Add the parameter force_accept_composition=True to "+
      "disable this error message and continue assuming 1 copy of sequence "+
      "file.)")
  if (len(missing) > 0):
    error = [
      "%d chain(s) not found in sequence file.  If the sequence does " % \
        len(missing),
      "not accurately represent the contents of the crystal after adjusting ",
      "for copy number, molecular replacement and building may fail.", ]
    if (force_accept_composition):
      print("WARNING: %s" % error, file=out)
    else :
      raise_sorry(error)
  counts = v.get_relative_sequence_copy_number()
  unique_counts = set(counts)
  if (-1 in unique_counts) : unique_counts.remove(-1)
  if (0 in unique_counts):
    unique_counts.remove(0)
    print("WARNING: %d sequence(s) not found in model.  This is not" %\
      counts.count(0), file=out)
    print("usually a problem, but it may indicate an incomplete model.", file=out)
  error = freq = None
  model_copies = copies_from_user
  if (model_copies is Auto):
    model_copies = copies_from_xtriage
  if (model_copies is None):
    model_copies = 1
    print("WARNING: assuming 1 copy of model", file=out)
  if (model_copies < 1):
    raise Sorry("Must have at least one copy of model.")
  if (len(unique_counts) == 0):
    error = "The sequence file provided does not appear to match the " +\
            "search model."
  elif (len(unique_counts) > 1):
    error = "The sequence file provided does not map evenly to the "+\
      "search model; this usually means an error in the sequence file (or "+\
      "model) such as accidental duplication or deletion of a chain. "+\
      "[counts: %s]" % list(unique_counts)
  else :
    # XXX float is_integer introduced in Python 2.6
    def is_integer(x) : return (int(x) == x)
    seq_freq = list(unique_counts)[0]
    assert seq_freq > 0
    # 1-1 mapping between PDB hierarchy and sequence
    if (seq_freq == 1.0):
      print("Assuming %d copies of sequence file (as well as model)" %\
        model_copies, file=out)
      return model_copies
    elif (seq_freq > 1):
      # hierarchy contains N copies of sequence
      if is_integer(seq_freq):
        freq = int(seq_freq) * model_copies
        print("Assuming %d copies of sequence file present" % freq, file=out)
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
      if is_integer(inverse_freq):
        # too much sequence
        if (inverse_freq > model_copies):
          error = "The number of copies of the model expected is less than "+\
            "the number indicated by the sequence file (%d versus %d).  " % \
            (model_copies, int(inverse_freq)) + \
            "This may mean that either the sequence file or the predicted "+\
            "number of copies is wrong."
          if ((copies_from_user is Auto) and
              (copies_from_xtriage is not None) and
              (assume_xtriage_copies_from_sequence_file)):
            print(error, file=out)
            print("", file=out)
            print("Since the number of copies was guessed by Xtriage", file=out)
            print("based on the sequence file, it will be scaled", file=out)
            print("by %g to be appropriate for the search model." % \
              inverse_freq, file=out)
            return seq_freq
        # too much model
        # XXX is this actually possible the way I've written the function?
        elif (inverse_freq < model_copies):
          error = "The number of copies of the model expected exceeds the "+\
            "number indicated by the sequence file (%d versus %d)." % \
            (model_copies, int(inverse_freq)) + \
            "This may mean that either the sequence file or the predicted "+\
            "number of copies is wrong."
        # model_copies times model contents matches sequence
        else :
          print("Assuming only one copy of sequence file contents.", file=out)
          return 1
      else :
        error = "The sequence file does not appear to contain a round "+\
          "number of copies of the model (%g)." % inverse_freq
  if (error is not None):
    if (force_accept_composition):
      print("WARNING: " + error, file=out)
      print("", file=out)
      print("Ambiguous input, defaulting to 1 copy of sequence file", file=out)
      return 1
    else :
      raise_sorry(error)
  else :
    raise RuntimeError("No error but frequency not determined")

def get_sequence_n_copies_from_files(seq_file, pdb_file, **kwds):
  from iotbx import file_reader
  import iotbx.pdb
  seq_in = file_reader.any_file(seq_file,
    raise_sorry_if_errors=True,
    raise_sorry_if_not_expected_format=True)
  if (seq_in.file_type != "seq"):
    raise Sorry("Can't parse %s as a sequence file.")
  try:
    pdb_in = iotbx.pdb.input(pdb_file)
  except Exception:
    raise Sorry("Can't parse %s as a PDB or mmCIF file.")
  kwds['pdb_hierarchy'] = pdb_in.construct_hierarchy()
  kwds['sequences'] = seq_in.file_object
  return get_sequence_n_copies(**kwds)

# XXX test needed
def group_chains_and_sequences(seq_file, pdb_file, **kwds):
  from iotbx import file_reader
  import iotbx.pdb
  seq_in = file_reader.any_file(seq_file,
    raise_sorry_if_errors=True,
    raise_sorry_if_not_expected_format=True)
  if (seq_in.file_type != "seq"):
    raise Sorry("Can't parse %s as a sequence file.")
  try:
    pdb_in = iotbx.pdb.input(pdb_file)
  except Exception:
    raise Sorry("Can't parse %s as a PDB or mmCIF file.")
  kwds['pdb_hierarchy'] = pdb_in.construct_hierarchy()
  kwds['sequences'] = seq_in.file_object
  v = validation(**kwds)
  chain_to_sequence_mappings = {}
  sequence_to_chain_mappings = {}
  for chain in v.chains :
    seq_id = chain.sequence_id
    chain_id = chain.chain_id
    if (seq_id is None):
      raise Sorry("Can't map chain %s to a sequence in %s." % (chain_id,
        seq_file))
    sequence = seq_in.file_object[seq_id].sequence
    if (chain_id in chain_to_sequence_mappings):
      if (chain_to_sequence_mappings[chain_id] != sequence):
        raise Sorry("Multiple unique chains named '%s'" % chain_id)
    else :
      chain_to_sequence_mappings[chain_id] = sequence
    if (not chain.sequence in sequence_to_chain_mappings):
      sequence_to_chain_mappings[sequence] = []
    sequence_to_chain_mappings[sequence].append(chain_id)
  return sequence_to_chain_mappings
