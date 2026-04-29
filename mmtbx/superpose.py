"""Model superposition tools"""

from __future__ import absolute_import, division, print_function
import sys
import os
import itertools

import cctbx.array_family.flex
import iotbx.pdb
import iotbx.pdb.amino_acid_codes
import iotbx.phil
import libtbx
import libtbx.math_utils
import libtbx.phil
import libtbx.runtime_utils
import mmtbx.alignment
import scitbx.math
import scitbx.math.superpose

from libtbx.utils import Sorry
from six.moves import zip
from six.moves import map

# Ideas:
#   preset options: Ca, backbone, all atom
#   sieve fit -- scitbx.math.superpose
#   SCEDS fit -- http://www.ncbi.nlm.nih.gov/pubmed/24189233
#   re-alignment of models using a different selection
#   automatic determination of subunit equivalences for multimers
#   chain independent superposition
# TODO: Merge back into a single pdb.
# TODO: Provide different alignment methods? Create subclasses override lsq?
# TODO: Port GUI

##############################################
# PHIL and Phenix Launcher infrastructure
##############################################

# PHIL Arguments.
PHIL_PARAMS = """
input {
  pdb_file_name_fixed = None
    .type = path
    .short_caption = Fixed model
    .help = Name of PDB file with model to fit to
    .style = bold file_type:pdb input_file
  pdb_file_name_moving = None
    .type = path
    .short_caption = Moving model
    .help = Name of PDB file with model that will be fit to pdb_file_name_fixed
    .style = bold file_type:pdb input_file
}

output {
  file_name = None
    .help = Name of PDB file with model that best fits to pdb_file_name_fixed
    .short_caption = Output model
    .type = path
    .style = bold new_file file_type:pdb
  include scope libtbx.phil.interface.tracking_params
}

selection_fixed = None
  .type = str
  .short_caption = Selection for fixed model
  .input_size = 400
  .help = Selection of the target atoms to fit to (optional)
selection_moving = None
  .type = str
  .short_caption = Selection for moving model
  .input_size = 400
  .help = Selection of the atoms that will be fit to selection_fixed (optional)
selection_fixed_preset = * ca backbone all
  .type = choice
  .help = Selection preset for fixed model.
selection_moving_preset = * ca backbone all
  .type = choice
  .help = Selection preset for moving model.
selection_fixed_chain = None
  .type = str
  .help = Selection chain for fixed model.
selection_moving_chain = None
  .type = str
  .help = Selection chain for moving model.

alignment
  .help = Set of parameters for sequence alignment. Defaults are good for most \
          of cases
  .short_caption = Sequence alignment
  .style = box auto_align
{
  alignment_style = local *global
    .type = choice
  gap_opening_penalty = 1
    .type = float
  gap_extension_penalty = 1
    .type = float
  similarity_matrix = blosum50  dayhoff *identity
    .type = choice
}"""

# Preset selections.
PRESETS = {
    'backbone': """pepnames and (name ca or name n or name c) and altloc \" \" chain %(chain)s""",
    'ca': """pepnames and (name ca) and altloc \" \" chain %(chain)s""",
    'all': """all chain %(chain)s""",
}
PRESET_ORDER = ['backbone', 'ca', 'all']

class SuperposePDB(object):
  """Superimpose PDB files.

  Instantiate with a filename or pdb instance. You may also provide a default
  selection string, selection preset, or chain to use for atom selection.
  Otherwise, the select method will attempt to find a reasonable default.

  Use the superpose(target) method to perform superposition. This will return
  the RMSD, and the transformation used. To write output, use output(lsq,
  filename).

  The selectomatic(target) method performs a pairwise comparison of chain
  sequence alignments, and updates the selections to use the chains with the
  highest similarity.

  You can also update the selection using the select_update() method, or return
  a given selection using select().

  Multiple models are supported by creating multiple instances. The class method
  open_models() helper can be used to open a file with multiple models.

  To create a subclass with a modified fitting routine, override fit().
  """

  # Alignment cache.
  _seqalign_cache = {}

  def __init__(self, filename=None, pdb=None, selection=None, preset=None,
               chain=None, quiet=False, log=None, desc=None):
    """A filename or pdb instance is required. You may also provide arguments
    for the atom selection. The quiet option will suppress informational
    output. The desc option provides a label for the log output.
    """
    self._quiet = quiet
    self._log = log
    self.desc = desc or os.path.basename(filename)

    # Validate input.
    if (filename is None) and (pdb is None):
      raise Sorry("Filename or PDB instance required.")

    # Read PDB or used specified pdb_hierarchy.
    self.pdb_input = pdb or iotbx.pdb.input(filename)
    self.pdb = self.pdb_input.construct_hierarchy()

    if len(self.pdb.models()) > 1:
      raise Sorry("Only one model supported; use SuperposePDB.open_models() to open multiple models.")

    # Select atoms.
    self.select_update(selection=selection, chain=chain, preset=preset)

  @classmethod
  def open_models(cls, filename=None, pdb=None, **kwargs):
    """Class method to return instances for each model in a file."""
    # pdb = pdb or iotbx.pdb.input(filename).construct_hierarchy()
    pdb = iotbx.pdb.input(filename)
    models = iotbx.pdb.input(filename).construct_hierarchy().models()
    desc = kwargs.pop('desc', 'moving')
    # print "Model check:"
    # for count, i in enumerate(models):
    #   print "\tmodel: %s, atoms: %s, chains: %s"%(count, len(i.atoms()), len(i.chains()))
    ret = []
    if len(models) == 1:
      ret.append(cls(pdb=pdb, desc=desc, **kwargs))
    elif len(models) > 1:
      for count, model in enumerate(models):
        pdb_new = iotbx.pdb.hierarchy.root()
        pdb_new.append_model(model.detached_copy())
        ret.append(cls(pdb=pdb_new, desc='%s-%s'%(desc, count), **kwargs))
    else:
      raise Sorry("No models found!")
    return ret

  def get_transformed(self):
    """Return a transformed model."""
    raise NotImplemented

  def selectomatic(self, target):
    """Perform pairwise sequence alignments and find the best aligned chain pair.

    This method performs pairwise sequence alignments between all chains in
    itself and the target. The selections are updated to the pair with the
    highest sequence similarity, with chain IDs in the return value.
    """
    alignments = []
    for moving_chain in self.pdb.chains():
      _, moving_selection, _ = self.select(selection="all chain %s"%moving_chain.id)
      for target_chain in target.pdb.chains():
        # For each pair,
        #   ... get the selection array
        _, target_selection, _ = target.select(selection="all chain %s"%target_chain.id)
        #   ... extract sequence from selection
        moving_seq, _ = self._extract_sequence_chain(moving_chain, moving_selection)
        target_seq, _ = target._extract_sequence_chain(target_chain, target_selection)
        #   ... perform alignment
        alignment = self._seqalign(moving_seq, target_seq)
        #   ... calculate score
        score = self._seqalign_score(alignment)
        alignments.append((score, moving_chain.id, target_chain.id, alignment))

    if not alignments:
      raise Sorry("No alignments?")

    best = sorted(alignments)[-1]
    self.log("Select-o-matic: Aligning chain %s to target chain %s"%(best[1], best[2]))
    self._print_seqalign(best[3])
    self.select_update(chain=best[1])
    target.select_update(chain=best[2])
    return best[1], best[2]

  def select(self, selection=None, chain=None, preset=None):
    """Update the selection.

    Specify either an explicit selection, or one of the handy preset selections
    (ca, backbone, all). Presets will use the longest chain by default, but
    this can also be specified using chain. If no selection or preset is
    specified, the various presets will be tested before falling back to all
    atom selection.
    """
    # Check for various Nones
    if selection in [None, "None", "none"]:
      selection = None

    # Use the specified chain, or find the longest chain.
    chains = sorted(self.pdb.chains(), key=lambda x:x.atoms_size())
    chain = chain or chains[-1].id

    if selection:
      # If a selection was specified, use it directly.
      selected_atoms = self.pdb.atom_selection_cache().selection(string=selection)
    elif PRESETS.get(preset):
      # Use a preset if it was specified.
      selection = PRESETS.get(preset)
      selection = selection%{'chain':chain}
      selected_atoms = self.pdb.atom_selection_cache().selection(string=selection)
    else:
      # Automagical selection; try the various presets,
      #   eventually falling back to 'all'
      for preset in PRESET_ORDER:
        selection = PRESETS[preset]
        selection = selection%{'chain':chain}
        selected_atoms = self.pdb.atom_selection_cache().selection(string=selection)
        if selected_atoms.count(True):
          break

    # Resequence the atoms and get the xyz coords.
    atoms = self.pdb.atoms()
    atoms.reset_i_seq()
    selected_xyz = atoms.extract_xyz().select(selected_atoms)
    return selection, selected_atoms, selected_xyz

  def select_update(self, selection=None, chain=None, preset=None):
    """select() and update state."""
    self.selection, self.selected_atoms, self.selected_xyz = self.select(selection, chain, preset)
    self.log("Using selection: %s, atoms total: %s, selected: %s"%(self.selection, len(self.selected_atoms), len(self.selected_xyz)))

  def log(self, *msg):
    """Log."""
    if self._quiet:
      return
    print("%s: "%self.desc, " ".join(map(str, msg)), file=self._log)

  def _print_rmsd(self, sites_moving, sites_fixed, desc=''):
    """Print statistics on the RMSD between sites_moving and sites_fixed."""
    rmsd = sites_fixed.rms_difference(sites_moving)
    deltas = (sites_fixed - sites_moving).norms()
    try:
      pbs = libtbx.math_utils.percentile_based_spread(deltas)
      self.log("Percentile-based spread (%s): %-.3f"%(desc, pbs))
    except Exception:
      self.log("Could not calculate percentile-based spread... Please fix!")
    self.log("RMSD between fixed and moving atoms (%s): %-.3f"%(desc, rmsd))

  def _print_lsq(self, lsq=None):
    """Print the transformation described by a lsq."""
    self.log("Rotation:")
    self.log("\n"+lsq.r.mathematica_form(label="r", one_row_per_line=True, format="%8.5f"))
    self.log("Translation:")
    self.log(lsq.t.mathematica_form(label="t", format="%8.5f"))

  def _print_seqalign(self, alignment, quiet=False):
    """Print a sequence alignment details."""
    matches = alignment.matches()
    self.log("Alignment details:")
    self.log("\tmatches after alignment: %s"%(matches.count("|") + matches.count("*")))
    self.log("\tsequence alignment:")
    # Since this prints directly, check if quiet.
    if not (quiet or self._quiet):
      # Change the labels to target.desc
      alignment.pretty_print(
        matches   = matches,
        block_size  = 50,
        n_block   = 1,
        top_name  = "moving",
        bottom_name = "fixed")

  def fit(self, sites_moving, sites_fixed):
    """Perform a least squares fit between sites_moving and sites_fixed. Return
    the rmsd, lsq, and updated sites_moving, sites_fixed."""
    # Check moving and fixed are the same size.
    if len(sites_moving) != len(sites_fixed):
      raise Sorry("Cannot align different number of atoms: %s and %s"%(len(sites_moving), len(sites_fixed)))
    if len(sites_moving) == 0:
      raise Sorry("Cannot align an empty set of atoms.")
    if len(sites_moving) < 2:
      raise Sorry("Cannot align a single atom.")

    # Perform the least squares fit
    # print "len:", len(sites_moving), len(sites_fixed)
    # for i,j in zip(sites_moving, sites_fixed):
    #    print "moving/fixed:", i, j
    lsq = scitbx.math.superpose.least_squares_fit(reference_sites=sites_fixed, other_sites=sites_moving)
    sites_moving2 = lsq.other_sites_best_fit()
    rmsd = sites_fixed.rms_difference(sites_moving2)
    return rmsd, lsq, sites_moving, sites_fixed

  def superpose(self, target):
    """Superpose to target. Return rmsd and lsq.

    This method will perform a sequence alignment first if the target is not
    sequence identical.
    """
    seq_moving, str_moving = self._extract_sequence_rgs()
    seq_fixed, str_fixed = target._extract_sequence_rgs()
    self.log("Let's go!")

    # Perform the sequence alignment
    alignment = self._seqalign(seq_moving, seq_fixed)
    self._print_seqalign(alignment)

    # Too many arguments; should it just extract seq/str in _align?
    rmsds = self._superpose_align(
      alignment,
      seq_moving,
      str_moving,
      seq_fixed,
      str_fixed,
      self.selected_atoms,
      target.selected_atoms
      )
    rmsds = sorted(rmsds, key=lambda x:x[0], reverse=True)
    if not rmsds:
      raise Sorry("Failed to get a good sequence alignment, or failed to fit models.")
    rmsd, lsq, xyz_moving, xyz_fixed = rmsds.pop()

    self.log("Initial stats:")
    self._print_rmsd(xyz_moving, xyz_fixed, desc='start')
    self.log("Best fit:")
    self.log("Number of atoms for LS fitting: %s"%xyz_moving.size())
    self._print_rmsd(lsq.other_sites_best_fit(), xyz_fixed, desc='final')
    self._print_lsq(lsq)
    return rmsd, lsq

  def _superpose_align(self, alignment, seq_moving, str_moving, seq_fixed, str_fixed, sel_moving, sel_fixed, windows=None):
    """Perform a sequence alignment, try to find a sequence alignment window
    that produces the lowest lsq fit. Note that this is a generator."""
    # Find the alignment window with the lowest RMSD.
    rmsds = []
    windows = (2,1,0) # windows or (0, 1)
    for count, window in enumerate(windows):
      # Create new coordinates based on the alignment.
      xyz_moving = cctbx.array_family.flex.vec3_double()
      xyz_fixed  = cctbx.array_family.flex.vec3_double()
      aligned = self._seqalign_pickmatches(alignment, window=window)
      for i,j,k in aligned:
        # alignment index, fixed index, moving index
        # assert seq_moving[j] == seq_fixed[k]
        # Add the atoms for this residue.
        for atom1 in str_moving[j].rg.atoms():
          for atom2 in str_fixed[k].rg.atoms():
            if (atom1.name == atom2.name) and sel_moving[atom1.i_seq] and sel_fixed[atom2.i_seq]:
              # print "checked:", atom1.name, atom2.name, atom1.i_seq, atom2.i_seq
              xyz_moving.append(atom1.xyz)
              xyz_fixed.append(atom2.xyz)

      # Yield the RMSD for this window.
      try:
        yield self.fit(xyz_moving, xyz_fixed)
      except Exception as e:
        pass
        # self.log("Fitting error: %s"%e)

  def _seqalign(self, seq_moving, seq_fixed, quiet=False, cache=True):
    """Perform sequence alignment using mmtbx.alignment. Return alignment.

    Sequence alignments are cached; use cache=False to disable.
    """
    key = (seq_moving, seq_fixed)
    if self._seqalign_cache.get(key) and cache:
      return self._seqalign_cache[key]

    opts = {
      'gap_opening_penalty': 12,
      'gap_extension_penalty': 1,
      'similarity_function': 'blosum50', # identity
      'style': 'global'
    }
    opts['seq_a'] = seq_moving
    opts['seq_b'] = seq_fixed
    align_obj = mmtbx.alignment.align(**opts)
    alignment = align_obj.extract_alignment()
    if cache:
      self._seqalign_cache[key] = alignment
    return alignment

  def _seqalign_pickmatches(self, alignment, window=0):
    """Find the indexes of aligned residues.

    Returns a list of lists with the indexes:
      [match index, sequence a index, sequence b index]
    """
    # Find sequence alignment matches -- "|" or "*"
    # Return the alignment index, the source index, and target index
    #   of each match.
    matches = alignment.matches()
    aligned = []
    for i, ia, ib, im, a, b in zip(
        itertools.count(0),
        alignment.i_seqs_a,
        alignment.i_seqs_b,
        matches,
        alignment.a,
        alignment.b):
      # print i, ia, ib, im, a, b
      left = i-window
      if left < 0:
        left = 0
      right = i+window+1
      if right > len(matches):
        right = len(matches)
      w = matches[left:right]
      append = all((i == '*' or i == '|') for i in w)
      if append:
        # print "window: %s -- left, right: %s %s -- i: %s -- append? %s"%(w, left, right, i, append)
        aligned.append((i, ia, ib))
    return aligned

  def _seqalign_score(self, alignment):
    matches = alignment.matches()
    total = len(alignment.a) - alignment.a.count("-")
    equal = matches.count("|")
    similar = matches.count("*")
    score = 100.*(equal+similar) / max(1,total)
    return score

  def _extract_sequence_rgs(self):
    """Extract the sequence and residue groups of the selected atoms."""
    sequence = []
    rgs = []
    for chain in self.pdb.chains():
      a, b = self._extract_sequence_chain(chain, self.selected_atoms)
      sequence.append(a)
      rgs.extend(b)
    return "".join(sequence), rgs

  def _extract_sequence_chain(self, chain, selection):
    """Extract the sequence and residue groups of a chain."""
    sequence = []
    rgs = []
    counter = 0
    for rg in chain.residue_groups():
      good = False
      for ai in rg.atoms().extract_i_seq():
        # print ai, selection[ai]
        if selection[ai]:
          good = True
          break
      if good:
        if len(rg.unique_resnames()) == 1:
          resname = rg.unique_resnames()[0]
          olc = iotbx.pdb.amino_acid_codes.one_letter_given_three_letter.get(resname, "X")
          if olc != "X":
            # proteins only
            sequence.append(olc)
            rgs.append(libtbx.group_args(i_seq=counter, rg=rg))
            counter += 1
          elif iotbx.pdb.common_residue_names_get_class(name=resname.strip()) == "common_rna_dna":
            # rna/dna
            sequence.append(resname.strip())
            rgs.append(libtbx.group_args(i_seq=counter, rg=rg))
            counter += 1
    return "".join(sequence), rgs

  def output(self, lsq, filename):
    """Output PDB model to filename given transformation lsq."""
    output = self.pdb.atoms()
    output.set_xyz(lsq.r.elems * output.extract_xyz() + lsq.t.elems)
    output.set_uij(cctbx.array_family.flex.sym_mat3_double(output.size(), [-1,-1,-1,-1,-1,-1]))
    self.pdb.write_pdb_file(
      file_name=filename,
      crystal_symmetry=self.pdb_input.crystal_symmetry()
    )
    return 0.0 # Return RMSD to output # TODO

  @classmethod
  def output_merge(cls, instances, filename):
    """Merge multiple instances back into a single PDB file."""
    pdb = iotbx.pdb.hierarchy.root()
    for i in instances:
      pdb.append_model(i.get_transformed().detached_copy())
    pdb.write_pdb_file(file_name=filename)

class SuperposePDBSieve(SuperposePDB):
  def fit(self, sites_moving, sites_fixed):
    lsq = scitbx.math.superpose.sieve_fit(sites_fixed=sites_fixed, sites_moving=sites_moving)
    sites_moving2 = lsq.other_sites_best_fit()
    rmsd = sites_fixed.rms_difference(sites_moving2)
    return rmsd, lsq, sites_moving, sites_fixed

class SuperposePDBSCEDS(SuperposePDB):
  def fit(self, sites_moving, sites_fixed):
    raise NotImplemented

def run(args, command_name="phenix.superpose_pdbs", log=None):
  import mmtbx.utils
  quiet = False
  #if log is None:
  #  log = mmtbx.utils.set_log(args)
  # Process the parameters.
  params = iotbx.phil.parse(PHIL_PARAMS, process_includes=True)
  args = mmtbx.utils.process_command_line_args(args=args, master_params=params, suppress_symmetry_related_errors=True)
  params = args.params

  # Copy command line args back into PHIL.
  fixed, moving = None, None
  if len(args.pdb_file_names) >= 2:
    fixed, moving = args.pdb_file_names[0], args.pdb_file_names[1]
    sources = []
    parameter_interpreter = libtbx.phil.command_line.argument_interpreter(master_phil=params, home_scope=None)
    sources.append(parameter_interpreter.process(arg="pdb_file_name_fixed=%s"%fixed))
    sources.append(parameter_interpreter.process(arg="pdb_file_name_moving=%s"%moving))
    params, _ = params.fetch(sources=sources, track_unused_definitions=True)

  # Get input and output filenames.
  pe = params.extract()
  fixed, moving, output = pe.input.pdb_file_name_fixed, pe.input.pdb_file_name_moving, pe.output.file_name
  if fixed is None and moving is None:
    raise Sorry("Need two input PDB files.")

  # Information.
  params.show(out=log)

  print("\n===== Init =====")
  # The fixed model can only contain a single model.
  # It will raise an Exception if there is more than one!
  fixed = SuperposePDB(
    fixed,
    selection=pe.selection_fixed,
    preset=pe.selection_fixed_preset,
    log=log,
    quiet=quiet,
    desc=fixed
  )

  # The moving pdb can contain many models. These will each be aligned to the
  # fixed model and output as a separate file...
  moving_args = dict(
      selection=pe.selection_moving,
      preset=pe.selection_moving_preset,
      desc=moving,
      log=log,
      quiet=quiet
  )
  for count, moving in enumerate(SuperposePDB.open_models(moving, **moving_args)):
    print("\n===== Aligning %s to %s ====="%(moving.desc, fixed.desc))
    if not pe.selection_moving:
      moving.selectomatic(fixed)
    rmsd, lsq = moving.superpose(fixed)
    if output:
      filename = '%s-%s.pdb'%(output, count)
      print("\n===== Writing %s output to %s ====="%(moving.desc, filename))
      moving.output(lsq, filename=filename)

class launcher(libtbx.runtime_utils.target_with_save_result):
  def run(self):
    return run(args=list(self.args), log=sys.stdout)

def finish_job(result):
  if (result is not None):
    (output_file, rmsd) = result
    output_files = [(output_file, "Superposed model")]
    statistics = [("RMSD", format_value("%.3f", rmsd))]
    return (output_files, statistics)
  return ([], [])

def validate_params(params):
  if (params.input.pdb_file_name_fixed is None or
      params.input.pdb_file_name_moving is None):
    raise Sorry("One or both PDB files missing.")
  if params.output.file_name is None :
    raise Sorry("Please specify a name for the output file.")

if __name__ == "__main__":
  run(args=sys.argv[1:])
