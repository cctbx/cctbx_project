"""Replace values in B-factor field with estimated B values.
Optionally remove low-confidence residues and split into domains."""
# -*- coding: utf-8 -*-

from __future__ import division, print_function
import os
from libtbx.utils import Sorry
from libtbx import group_args
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate

# =============================================================================
class Program(ProgramTemplate):
  description = '''
Replace values in B-factor field with estimated B values.
Optionally remove low-confidence residues and split into domains.

Inputs: Model file (PDB, mmCIF)
'''

  datatypes = ['phil', 'model', 'sequence']

  master_phil_str = """

  job_title = None
    .type = str
    .input_size = 400
    .help = Job title in PHENIX GUI, not used on command line
    .style = noauto bold

  input_files {
    chain_id = None
      .type = str
      .short_caption = chain_id
      .help = If specified, find domains in this chain only. NOTE: only one \
                chain can be used for finding domains at a time.
      .input_size = 400

    selection = None
      .type = str
      .short_caption = Selection
      .help = If specified, use only selected part of model
       .input_size = 400

    pae_file = None
      .type = path
      .help = Optional input json file with matrix of inter-residue estimated\
               errors (pae file)
      .short_caption = PAE file

    distance_model_file = None
      .type = path
      .help = Distance_model_file. A PDB or mmCIF file containing the model \
               corresponding to the PAE  \
               matrix. Only needed if weight_by_ca_ca_distances is True and \
               you want to specify a file other than your model file. (Default\
               is to use the model file)
      .short_caption = Distance model file

    model = None
      .type = path
      .help = Input predicted model (e.g., AlphaFold model).  Assumed to \
                have LDDT values in B-value field (or RMSD values).
      .style = file_type:pdb input_file
      .short_caption = Predicted model
      .expert_level = 3

  }
  output_files {
     target_output_format = *None pdb mmcif
       .type = choice
       .help = Desired output format (if possible). Choices are None (\
                try to use input format), pdb, mmcif.  If output model\
                 does not fit in pdb format, mmcif will be used. \
                 Default is pdb.
       .short_caption = Desired output format

    processed_model_prefix = None
      .type = str
      .help = Output file with processed models will begin with this prefix.\
               If not specified, the input model file name will be used with\
               the suffix _processed.
      .short_caption = Output file prefix (optional)

    remainder_seq_file_prefix = None
      .type = str
      .help = Output file with sequences of deleted parts of model \
          will begin with this prefix
      .short_caption = Output remainder seq file prefix

    mark_atoms_to_keep_with_occ_one = False
       .type = bool
       .help = Mark atoms to keep with occupancy of 1 and those to remove \
                  with zero (default is to remove those that are not desired)
       .short_caption = Mark atoms to keep with occupancy of 1


     maximum_output_b = 999.
       .type = float
       .help = Limit output B values (so that they fit in old-style PDB \
              format). Note that this limit is applied just before writing \
              output model files, it does not affect anything else. Also \
              if output model is trimmed normally, high-B value residues are\
              already removed, so in most cases this keyword has no effect.
       .short_caption = Maximum output B

     remove_hydrogen = True
       .type = bool
       .help = Remove hydrogen atoms from model on input
       .short_caption = Remove hydrogen

     single_letter_chain_ids = False
       .type = bool
       .help = Write output files with all chain IDS as single characters.\
                Default is to use original chain ID and to add digits (1-9)\
                for domains.
       .short_caption = Use only single-letter chain ID

  }

  include scope mmtbx.process_predicted_model.master_phil_str
  control {

    write_files = True
      .type = bool
      .help = Write output files
      .short_caption = Write output files
   }


  gui
    .help = "GUI-specific parameter required for output directory"
  {
    output_dir = None
    .type = path
    .style = output_dir
  }


  """

  def run(self):

    # print version and date
    # self.print_version_info()
    self.data_manager.set_overwrite(True)

    #
    self.get_data_inputs()  # get any file-based information
    from iotbx.pdb.utils import set_target_output_format_in_params
    set_target_output_format_in_params(self.params)

    # self.print_params()
    self.starting_model = self.model
    self.model_list = []
    self.processed_model = None
    self.processed_model_text = None
    self.processed_model_file_name = None
    self.processed_model_file_name_list = []

    # Split by chain ID. Normally only one chain ID but run separately if more
    working_model_list = self.model.as_model_manager_each_chain()
    from mmtbx.process_predicted_model import process_predicted_model

    processed_model_list = []
    remainder_sequence_str_list = []
    for model in working_model_list:
      # Run on each unique chain ID
      info = process_predicted_model(model = model,
       distance_model = self.distance_model,
       params = self.params,
       pae_matrix = self.pae_matrix,
       mark_atoms_to_keep_with_occ_one = \
          self.params.output_files.mark_atoms_to_keep_with_occ_one,
       log = self.logger,
       )

      if not info:
        print("Unable to process predicted model", file = self.logger)
        return

      self.model_list += info.model_list
      processed_model_list.append(info.model)
      remainder_sequence_str_list.append(info.remainder_sequence_str)


    if len(processed_model_list) > 1:
      # rename chains to match orig and merge. Limited functionality here
      for m, orig_m in zip(processed_model_list, working_model_list):
        chain_id = orig_m.first_chain_id()
        for model_m in m.get_hierarchy().models()[:1]:
          for chain in model_m.chains():
            chain.id = chain_id

      from iotbx.pdb.utils import add_model
      self.processed_model = processed_model_list[0]
      for m in processed_model_list[1:]:
        self.processed_model = add_model(self.processed_model, m)
    else: # usual
      self.processed_model = processed_model_list[0]

    self.dock_and_rebuild = group_args(
      group_args_type = 'dummy dock_and_rebuild for summary',
      processed_model = self.processed_model)

    if not self.params.control.write_files:
      return  # done


    starting_residues = self.model.get_hierarchy().overall_counts().n_residues
    print("\nStarting residues: %s" %(starting_residues), file = self.logger)

    if self.params.output_files.processed_model_prefix:
      prefix = self.params.output_files.processed_model_prefix
    else:
      prefix, ext  = os.path.splitext(self.params.input_files.model)
      prefix = "%s_processed" %(prefix)
      prefix = os.path.basename(prefix)
    self.processed_model_file_name = "%s.pdb" %(prefix) # PDB OK updated below
    if not self.processed_model or \
         self.processed_model.overall_counts().n_residues < 1:
      print("No residues obtained after processing...", file = self.logger)
      return None

    # Special case: write out split models marking residues as missing with
    #   occ = 0
    if self.params.output_files.mark_atoms_to_keep_with_occ_one:
      ii = 0
      for m in self.model_list:
        ii += 1
        fn = "%s_%s.pdb" %(prefix,ii)
        fn = self.data_manager.write_model_file(
           m, fn,
           format=self.params.output_files.target_output_format)
        print("Wrote model with all residues present to %s" %(fn)
          +"\n marking domain %s with" %(
           ii) + " occupancy=1 and rest with occupancy=0")
      return
    if (self.params.output_files.maximum_output_b is not None) and (
       self.processed_model.get_b_iso().min_max_mean().max >
       self.params.output_files.maximum_output_b):
      print("Limiting output B values to %.0f" %(
        self.params.output_files.maximum_output_b), file = self.logger)
    mm = limit_output_b(self.processed_model,
         maximum_output_b = self.params.output_files.maximum_output_b)

    if self.params.output_files.single_letter_chain_ids:
      # convert all chain_ids to single character (A)
      mm_to_split = convert_chain_ids_to_single_character(mm, chain_id = "A")
    else:
      mm_to_split = mm

    # original (multi-char chain IDs)

    self.processed_model_file_name = self.data_manager.write_model_file(
      mm, self.processed_model_file_name,
      format=self.params.output_files.target_output_format)

    final_residues = mm.get_hierarchy().overall_counts().n_residues
    print("Final residues: %s\n" %(final_residues), file = self.logger)

    # Split up processed model and write each chain as well
    if len(mm_to_split.chain_ids()) > 1 or \
         self.params.output_files.single_letter_chain_ids:
      model_list = mm_to_split.as_model_manager_each_chain()
    else:
      model_list = [mm_to_split]
    count = 0
    for m in model_list:
      count += 1
      chain_id = m.first_chain_id().strip()
      if not chain_id:
        if len(model_list) > 1:
          raise Sorry(
           "Input model cannot have a blank chain ID and non-blank chain IDS")
        chain_id = "A"
      fn = "%s_%s_%s.pdb" %(prefix,chain_id, count)      # PDB OK (updated below)
      if not m or not m.overall_counts().n_residues:
        print("Skipping #%s (no residues)" %(count), file = self.logger)
        continue
      mm = limit_output_b(m,
           maximum_output_b = self.params.output_files.maximum_output_b)
      fn = self.data_manager.write_model_file(mm,fn,
        format=self.params.output_files.target_output_format)
      print("Copying predicted model chain %s (#%s) to %s" %(
           chain_id,count,fn), file = self.logger)

      self.processed_model_file_name_list.append(fn)


    # Write out seq file for remainder (unused part) of model
    if remainder_sequence_str_list:
      remainder_sequence_str = "\n".join(remainder_sequence_str_list)
      if self.params.output_files.remainder_seq_file_prefix:
        prefix = self.params.output_files.remainder_seq_file_prefix
      else:
        prefix, ext  = os.path.splitext(self.params.input_files.model)
        prefix = "%s_remainder" %(prefix)
        prefix = os.path.basename(prefix)
      self.remainder_sequence_file_name = os.path.join(
        os.getcwd(), "%s.seq" %(prefix))
      sequence_str = remainder_sequence_str
      self.data_manager.write_sequence_file(sequence_str,
        filename = self.remainder_sequence_file_name)

    # Summarize each step

    self.summarize_predicted_model()

    print ('\nFinished with process_predicted_model', file=self.logger)
  # ---------------------------------------------------------------------------
  def get_results(self):

    return group_args(
      group_args_type = 'results of process_predicted_model',
      starting_model = self.starting_model,
      processed_model = self.processed_model,
      processed_model_file_name = self.processed_model_file_name,
      processed_model_file_name_list = self.processed_model_file_name_list,
      model_list = self.model_list,
      processed_model_text = self.processed_model_text,
      params = self.params)

# =============================================================================
#    Custom operations
# =============================================================================
#

  def summarize_predicted_model(self):
    #  Process predicted model
    dr = self.dock_and_rebuild
    processed_model = dr.processed_model
    if processed_model:
      chain_id_list = processed_model.chain_ids()
      text = ""
      text += "Processed model with %s domains and %s residues " %(
           len(chain_id_list),
           processed_model.get_hierarchy().overall_counts().n_residues,)
      text += "\n\nResidues by domain (as chains):"

      self.processed_model_text = ""
      for chain_id in chain_id_list:
        m = processed_model.apply_selection_string("chain '%s'" %(chain_id))
        n_found = m.get_hierarchy().overall_counts().n_residues
        text += "\nCHAIN: %s   Residues: %s " %(
           chain_id,n_found)
        print(text, file = self.logger)
        self.processed_model_text += text
    else:
        self.processed_model_text = "No processed model information available"


  def set_defaults(self):
    # set params for files identified automatically and vice versa
    params=self.params
    self.titles=[]

    # Read in default model as pdb_in
    file_name = getattr(params.input_files,'model',None)
    if file_name:
      if not os.path.isfile(file_name):
        raise Sorry("The file %s is missing?" %(file_name))
      else:
        pass # ok so far
    else: # guess it
      try:
        file_name=self.data_manager.get_default_model_name()
        self.params.input_files.model=file_name
      except Exception as e:
        pass # did not work

  def get_data_inputs(self):  # get any file-based information
    self.set_defaults()
    file_name=self.params.input_files.model
    if not file_name:
      raise Sorry("Unable to guess model file name...please specify")
    if not os.path.isfile(file_name):
      raise Sorry("Missing the model file '%s'" %(file_name))
    try:
      self.model=self.data_manager.get_model(filename=file_name)
    except Exception as e:
      raise Sorry("Failed to read model file '%s'" %(file_name))

    print("Read model from %s" %(file_name), file = self.logger)

    if not self.model:
      raise Sorry("Missing model")

    self.model.add_crystal_symmetry_if_necessary()

    # Remove waters and hetero atoms (ligands)
    self.model = self.model.apply_selection_string(
       "(not hetero) and (not water)")

    # Remove hydrogens and apply user selection
    selections = []
    if self.params.output_files.remove_hydrogen:
      selections.append("(not (element H))")
    if self.params.input_files.selection:
      selections.append("(%s)" %(self.params.input_files.selection))
    if selections:
      self.model = self.model.apply_selection_string(" and ".join(selections))

    if self.params.process_predicted_model.weight_by_ca_ca_distance:
      if not self.params.input_files.distance_model_file:
        self.params.input_files.distance_model_file = \
           self.params.input_files.model
      file_name=self.params.input_files.distance_model_file
      self.distance_model = self.data_manager.get_model(filename=file_name)
      print("Read distance model from %s" %(file_name), file = self.logger)
      self.distance_model.add_crystal_symmetry_if_necessary()
    else:
      self.distance_model = None

    if self.params.input_files.pae_file and \
       os.path.isfile(self.params.input_files.pae_file):
      from mmtbx.domains_from_pae import parse_pae_file
      self.pae_matrix = parse_pae_file(self.params.input_files.pae_file)
    else:
      self.pae_matrix = None

    if len(self.model.chain_ids()) > 1:
      if self.distance_model:
        raise Sorry(
          "The distance model cannot currently be used for multiple chains")
      if self.pae_matrix:
        raise Sorry(
          "The PAE matrix cannot currently be used for multiple chains")

  def validate(self):  # make sure we have files
    return True

  def print_params(self):
    import iotbx.phil
    master_phil = iotbx.phil.parse(master_phil_str)
    print ("\nInput parameters for process_predicted_model:\n", file = self.logger)
    master_phil.format(python_object = self.params).show(out = self.logger)

  def print_version_info(self):

    # Print version info
    import time
    print ("\n"+60*"*"+"\n"+10*" "+"PHENIX process_predicted_model" +\
      "  "+str(time.asctime())+"\n"+60*"*"+"\n",file=self.logger)
    print ("Working directory: ",os.getcwd(),"\n",file=self.logger)
    print ("PHENIX VERSION: ",os.environ.get('PHENIX_VERSION','svn'),"\n",
     file=self.logger)

def convert_chain_ids_to_single_character(m, chain_id = "A"):

  from mmtbx.secondary_structure.find_ss_from_ca import set_chain_id
  # rename all chains to "A" but keep separate chains for different segments
  mm = m.deep_copy()
  set_chain_id(mm.get_hierarchy(), chain_id)
  mm.reset_after_changing_hierarchy()
  return mm

def get_available_letter(c, all_letters, used_ids):

  if len(c) == 1 and not c in used_ids:
    return c
  if c and (not (c[0] in used_ids)):
    return c[0]

  for a in all_letters:
    if not (a in used_ids):
      return a
  return None



def limit_output_b(m, maximum_output_b = None):
  """ create deep copy of model m in which all isotropic
      values > maximum_output_b
      are set to maximum_output_b. If maximum_output_b is None or there are no
      b-values > maximum_output_b, return original model (not deep copy)"""

  if (maximum_output_b is not None) and (
       m.get_b_iso().min_max_mean().max > maximum_output_b):
    b_values = m.get_b_iso()
    b_values.set_selected((b_values > maximum_output_b), maximum_output_b)
    mm = m.deep_copy() # REQUIRED so we do not modify b-values in m itself
    mm.set_b_iso(b_values)
    return mm
  else:
    return m
# =============================================================================
# for reference documentation keywords
master_phil_str = Program.master_phil_str
