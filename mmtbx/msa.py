
from libtbx import easy_run
import libtbx.phil
import libtbx.load_env
from libtbx.utils import Sorry
from libtbx import adopt_init_args
import os
import sys

# XXX this is more complex than ideal, in order to support long and
# potentially problematic "sequence names", which in some applications will
# actually be file paths + chain IDs.  to ensure that MUSCLE doesn't choke
# on these, the option is given to substitute numerical sequence names
# internally, while still allowing retrieval using the original names.
class align_pdb_residues (object) :
  """
  Provides mapping between original residue numbers or IDs and a reference
  numbering, via multiple sequence alignment.  Depending on whether or not
  the insertion code was used in the original PDB files, it can use either
  the resseq or the resid as the lookup key.
  """
  def __init__ (self,
                pdb_sequences,
                pdb_names,
                pdb_offsets=None,
                pdb_resids=None,
                reference_sequence=None,
                reference_sequence_name="reference_seq",
                reference_sequence_offset=0,
                reference_sequence_resids=None,
                reference_index=None,
                substitute_names=False,
                out=None) :
    adopt_init_args(self, locals())
    n_models = len(pdb_sequences)
    assert (pdb_resids is not None) or (pdb_offsets is not None)
    assert ((n_models >= 1) and (n_models==len(pdb_names)))
    if (pdb_offsets is not None) :
      assert (pdb_resids is None)
      assert (len(pdb_names) == len(pdb_offsets))
      assert (reference_sequence is None) or (reference_sequence_resids is None)
    else :
      assert (len(pdb_names) == len(pdb_resids))
      assert ((reference_sequence is None) or
              (reference_sequence_resids is not None))
    self.fasta_names = []
    self._name_lookup = None
    self.run_alignment(out=out)
    self.build_lookup_table()
    self.out = None # XXX required for pickling

  def run_alignment (self, out=None) :
    use_pdb_sequence = False
    # build index of sequence names given to MUSCLE
    if self.substitute_names :
      self._name_lookup = {}
      for i, name in enumerate(self.pdb_names) :
        self._name_lookup[name] = str(i)
        self.fasta_names.append(str(i))
    else :
      self.fasta_names = self.pdb_names
    # determine sequence for reference numbering
    if (self.reference_sequence is not None) :
      assert (self.reference_index is None)
      assert ((self.reference_sequence_resids is not None) or
              (self.reference_sequence_offset is not None))
    else :
      i_ref = self.reference_index
      use_pdb_sequence = True
      if (i_ref is None) :
        i_ref = self.reference_index = 0
      self.reference_sequence = self.pdb_sequences[i_ref]
      if (self.pdb_offsets is not None) :
        self.reference_sequence_offset = self.pdb_offsets[i_ref]
      else :
        self.reference_sequence_resids = self.pdb_resids[i_ref]
      if self.substitute_names :
        self.reference_sequence_name = self.fasta_names[i_ref]
      else :
        self.reference_sequence_name = self.pdb_names[i_ref]
    fasta = "\n".join([ ">%s\n%s" % (n,s) for n,s in zip(self.fasta_names,
                        self.pdb_sequences) ])
    if (not use_pdb_sequence) :
      ref_seq_fasta = ">%s\n%s" % (self.reference_sequence_name,
        self.reference_sequence)
      fasta = ref_seq_fasta + "\n" + fasta
    self.muscle_aln = get_muscle_alignment(fasta, out=out)
    assert (self.muscle_aln is not None)

  # I am not proud of this.
  def build_lookup_table (self) :
    # find sequences and determine equivalent numbering
    self._lookup_table = {}
    self._indices = {}
    self._resids_padded = {}
    i_ref = self.muscle_aln.names.index(self.reference_sequence_name)
    self._reference_alignment = self.muscle_aln.alignments[i_ref]
    assert (i_ref is not None)
    for i, name in enumerate(self.muscle_aln.names) :
      if (i == i_ref) :
        self._lookup_table[name] = None
        if (self.pdb_offsets is None) :
          indices = {}
          resids_padded = [None] * self.get_alignment_size()
          h = k = 0
          for resi in self._reference_alignment :
            if (resi != '-') :
              resid = self.reference_sequence_resids[k]
              if (resid is not None) :
                resid = resid.strip()
                indices[resid] = h
                resids_padded[h] = resid
              k += 1
            h += 1
          self._indices[name] = indices
          self._resids_padded[name] = resids_padded
      else :
        i_pdb = self.fasta_names.index(name)
        aln = self.muscle_aln.alignments[i]
        h = j = k = 0
        # case 1: index by resseq
        if (self.pdb_offsets is not None) :
          new_resseqs = []
          for res1, res2 in zip(aln, self._reference_alignment) :
            if (res2 != '-') :
              k += 1
            if (res1 == '-') :
              continue
            else :
              if (res2 == '-') :
                new_resseqs.append(None)
              else :
                new_resseqs.append(k - self.reference_sequence_offset)
              j += 1
          assert (j == len(new_resseqs))
          self._lookup_table[name] = new_resseqs
        # case 2: index by resid
        else :
          new_resids = {}
          indices = {}
          resids_padded = [None] * self.get_alignment_size()
          for res1, res2 in zip(aln, self._reference_alignment) :
            resid1 = resid2 = None
            if (res2 != '-') :
              resid2 = self.reference_sequence_resids[k]
              if (resid2 is not None) :
                resid2 = resid2.strip()
              k += 1
            if (res1 != '-') :
              resid1 = self.pdb_resids[i_pdb][j]
              if (resid1 is not None) :
                resid1 = resid1.strip()
                new_resids[resid1] = resid2
                indices[resid1] = h
                resids_padded[h] = resid1
              j += 1
            h += 1
          assert (len(new_resids) == (len(self.pdb_resids[i_pdb]) -
                                      self.pdb_resids[i_pdb].count(None)))
          self._lookup_table[name] = new_resids
          self._indices[name] = indices
          self._resids_padded[name] = resids_padded

  def get_alignment_size (self) :
    return len(getattr(self, "_reference_alignment", []))

  def write_file (self, file_name) :
    f = open(file_name, "w")
    f.write(str(self.muscle_aln))
    f.close()

  def get_name_index (self, pdb_name) :
    if (self._name_lookup is not None) :
      return self._name_lookup[pdb_name]
    else :
      return pdb_name
      #return self.pdb_names.index(pdb_name)

  def get_residue_position (self, pdb_name, resseq=None, resid=None) :
    if (resseq is not None) :
      assert (resid is None)
      return self.convert_residue_number(pdb_name, resseq)
    else :
      assert (resid is not None)
      return self.convert_resid(pdb_name, resid)

  def get_resid_array_index (self, pdb_name, resid) :
    assert (resid is not None)
    seq_name = self.get_name_index(pdb_name)
    indices = self._indices[seq_name]
    return indices[resid]

  def get_all_resids_at_index (self, index) :
    resids = []
    for name in self.fasta_names : #self._resids_padded.keys() :
      resids.append(self._resids_padded[name][index])
    return resids

  def get_reference_resid (self, index) :
    if (self.pdb_offsets is None) :
      resids = self._resids_padded[self.reference_sequence_name]
      return resids[index]
    else :
      return index + 1

  def convert_resid (self, pdb_name, resid) :
    assert (self.pdb_resids is not None)
    assert (isinstance(resid, str))
    seq_name = self.get_name_index(pdb_name)
    if (self._lookup_table[seq_name] is None) :
      return resid
    try :
      return self._lookup_table[seq_name][resid.strip()]
    except KeyError, e :
      raise RuntimeError("""\
Encountered IndexError attempting to convert residue ID!
Values:
  pdb_name = %s
  resid = %s
  seq_name = %s

Dump of full alignment:
  %s
""" % (pdb_name, resseq, seq_name, self.muscle_aln))

  def convert_residue_number (self, pdb_name, resseq) :
    assert (self.pdb_offsets is not None)
    seq_name = self.get_name_index(pdb_name)
    assert isinstance(resseq, int)
    if (self._lookup_table[seq_name] is None) :
      return resseq
    i_pdb = self.pdb_names.index(pdb_name)
    offset = self.pdb_offsets[i_pdb]
    i_res = resseq + offset - 1
    try :
      return self._lookup_table[seq_name][i_res]
    except IndexError, e :
      raise RuntimeError("""\
Encountered IndexError attempting to convert residue number!
Values:
  pdb_name = %s
  resseq = %s
  seq_name = %s
  i_res = %s

Dump of full alignment:
  %s
""" % (pdb_name, resseq, seq_name, i_res, self.muscle_aln))

def align_pdb_hierarchies (hierarchies,
                           hierarchy_names,
                           reference_hierarchy=None,
                           substitute_names=True,
                           log=None) :
  """
  Convenience function: takes a collection of *single-model, single-chain* PDB
  hierarchies, plus optional reference hierarchy, extracts sequences and
  resseq offsets or resids, and generates alignment object.  By default, if
  any of the hierarchies contain atoms with insertion codes, the resid mapping
  will be used automatically.
  """
  assert (reference_hierarchy is None)
  if (log is None) :
    log = sys.stdout
  pdb_resids = []
  i = 0
  reference_index = None
  pdb_sequences = []
  def strip (string_or_none) :
    if (string_or_none is None) :
      return None
    return string_or_none.strip()
  for hierarchy, name in zip(hierarchies, hierarchy_names) :
    assert (hierarchy.overall_counts().n_chains == 1)
    main_conf = hierarchy.models()[0].chains()[0].conformers()[0]
    chain_seq = main_conf.as_padded_sequence(skip_insertions=False)
    pdb_sequences.append(chain_seq)
    resids = main_conf.get_residue_ids(skip_insertions=False)
    pdb_resids.append([ strip(resid) for resid in resids ])
    if (hierarchy is reference_hierarchy) :
      reference_index = i
      print >> log, "  Using %s for sequence numbering" % name
    elif (reference_hierarchy is None) and (reference_index is None) :
      reference_index = i
      print >> log, "  Using %s for sequence numbering" % name
  if (reference_index is not None) :
    pdb_names = hierarchy_names
    msa_manager = align_pdb_residues(
      pdb_sequences=pdb_sequences,
      pdb_names=pdb_names,
      pdb_resids=pdb_resids,
      reference_index=reference_index,
      substitute_names=substitute_names,
      out=log)
  else :
    assert (reference_hierarchy.overall_counts().n_chains == 1)
    main_conf = reference_hierarchy.models()[0].chains()[0].conformers()[0]
    reference_sequence = main_conf.as_padded_sequence(skip_insertions=False)
    reference_sequence_offset = reference_sequence_resids = None
    resids = main_conf.get_residue_ids(skip_insertions=False)
    reference_sequence_resids = [ strip(resid) for resid in resids ]
    msa_manager = align_pdb_residues(
      pdb_sequences=pdb_sequences,
      pdb_names=pdb_names,
      pdb_resids=pdb_resids,
      reference_sequence=reference_sequence,
      reference_sequence_resids=reference_sequence_resids,
      substitute_names=substitute_names,
      out=self.log)
  return msa_manager

def run_muscle (fasta_sequences, group_sequences=True) :
  assert group_sequences # XXX this isn't actually optional!
  if not libtbx.env.has_module(name="muscle") :
    raise RuntimeError("MUSCLE not available or not configured.")
  exe_path = libtbx.env.under_build("muscle/exe/muscle")
  if (os.name == "nt") :
    exe_path += ".exe"
  if (not os.path.isfile(exe_path)) :
    raise RuntimeError("muscle executable is not available.")
  assert isinstance(fasta_sequences, str)
  cmd = "%s -quiet -clw" % exe_path
  if group_sequences :
    cmd += " -group"
  #else :
  #  cmd += " -stable"
  muscle_out = easy_run.fully_buffered(cmd,
    stdin_lines=fasta_sequences)
  return muscle_out.stdout_lines

def get_muscle_alignment (fasta_sequences, group_sequences=True, out=None) :
  muscle_out = run_muscle(fasta_sequences, group_sequences)
  from iotbx.bioinformatics import clustal_alignment_parse
  alignment, null = clustal_alignment_parse("\n".join(muscle_out))
  if (out is not None) :
    print >> out, "\n".join(muscle_out)
  return alignment

########################################################################
# PHENIX GUI ADAPTOR
master_phil = libtbx.phil.parse("""
muscle
  .caption = PHENIX includes the open-source multiple sequence alignment \
    program MUSCLE, written by Bob Edgar.  It can produce output identical in \
    format to CLUSTALW, which is suitable for input to the Sculptor model \
    preparation program.  You may provide all sequences in a single file, \
    or in as many different files as desired.  Alternately, you may provide \
    one or more PDB files from which the chain sequence(s) will be extracted \
    (note however that this is only useful for homogenous models).
  .style = caption_img:icons/custom/msa64.png caption_width:480
{
  seq_file = None
    .type = path
    .multiple = True
    .short_caption = Sequence or PDB file
    .style = file_type:seq,pdb use_list input_file
  output_file = None
    .type = path
    .style = bold file_type:aln new_file
  load_in_text_editor = True
    .type = bool
    .short_caption = Open alignment in text editor when complete
}""")

def run (args=(), params=None, out=sys.stdout) :
  assert (params is not None)
  seq_files = params.muscle.seq_file
  output_file = params.muscle.output_file
  if (output_file is None) or (output_file == "") :
    output_file = os.path.join(os.getcwd(), "muscle.aln")
  from iotbx import file_reader
  from iotbx.bioinformatics import any_sequence_format, sequence
  seqs = []
  for file_name in seq_files :
    if (file_name.endswith(".pdb") or file_name.endswith(".ent") or
        file_name.endswith(".pdb.gz") or file_name.endswith(".ent.gz")) :
      pdb_in = file_reader.any_file(file_name, force_type="pdb").file_object
      hierarchy = pdb_in.construct_hierarchy()
      first_model = hierarchy.models()[0]
      found_protein = False
      for chain in first_model.chains() :
        first_conf = chain.conformers()[0]
        if first_conf.is_protein() :
          chain_seq = first_conf.as_padded_sequence()
          base_name = os.path.basename(file_name)
          seq_name = "%s_%s" % (os.path.splitext(base_name)[0], chain.id)
          seqs.append(sequence(chain_seq, seq_name))
          found_protein = True
      if (not found_protein) :
        raise Sorry(("The PDB file %s does not contain any recognizable "+
          "protein chains.") % file_name)
    else :
      try :
        seq_objects, non_compliant = any_sequence_format(file_name,
          assign_name_if_not_defined=True)
        seqs.extend(seq_objects)
      except Exception, e :
        raise Sorry(("Error parsing '%s' - not a recognizable sequence "+
          "format.  (Original message: %s)") % (file_name, str(e)))
  if (len(seqs) < 2) :
    raise Sorry("Need at least two valid sequences to run MUSCLE.")
  combined = "\n".join([ seq.format(80) for seq in seqs ])
  muscle_out = run_muscle(combined)
  open(output_file, "w").write("\n".join(muscle_out))
  return (output_file, "\n".join(muscle_out))

def validate_params (params) :
  if (len(params.muscle.seq_file) == 0) :
    raise Sorry("No sequence files provided!")
