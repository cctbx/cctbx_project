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
Process Predicted Model

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
    # self.print_params()

    self.starting_model = self.model
    self.model_list = []
    self.processed_model = None
    self.processed_model_file_name = None

    from mmtbx.process_predicted_model import process_predicted_model
    info = process_predicted_model(model = self.model,
     distance_model = self.distance_model,
     params = self.params,
     pae_matrix = self.pae_matrix,
     log = self.logger,
       )

    if not info:
      print("Unable to process predicted model", file = self.logger)
      return

    self.model_list = info.model_list
    self.processed_model = info.model
    self.dock_and_rebuild = group_args(
      group_args_type = 'dummy dock_and_rebuild for summary',
      processed_model = self.processed_model)

    if not self.params.control.write_files:
      return  # done

    starting_residues = self.model.get_hierarchy().overall_counts().n_residues
    print("Starting residues: %s" %(starting_residues), file = self.logger)

    if self.params.output_files.processed_model_prefix:
      prefix = self.params.output_files.processed_model_prefix
    else:
      prefix, ext  = os.path.splitext(self.params.input_files.model)
      prefix = "%s_processed" %(prefix)
      prefix = os.path.basename(prefix)
    self.processed_model_file_name = "%s.pdb" %(prefix)
    self.data_manager.write_model_file(
      info.model, self.processed_model_file_name)

    # Write out seq file for remainder (unused part) of model
    if info.remainder_sequence_str:
      if self.params.output_files.remainder_seq_file_prefix:
        prefix = self.params.output_files.remainder_seq_file_prefix
      else:
        prefix, ext  = os.path.splitext(self.params.input_files.model)
        prefix = "%s_remainder" %(prefix)
      self.processed_model_file_name = os.path.join(
        os.getcwd(), "%s.seq" %(prefix))
      sequence_str = info.remainder_sequence_str
      self.data_manager.write_sequence_file(sequence_str,
        filename = self.processed_model_file_name)

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
      for chain_id  in chain_id_list:
        m = processed_model.apply_selection_string("chain %s" %(chain_id))
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
    file_name=self.data_manager.get_default_model_name()
    if not file_name in [
          getattr(params.input_files,'model',None),]:
      # set it by default
      self.params.input_files.model=file_name
      print ("Set model=%s " %(file_name),file=self.logger)


  def get_data_inputs(self):  # get any file-based information
    self.set_defaults()
    file_name=self.params.input_files.model
    self.model=self.data_manager.get_model(filename=file_name)
    print("Read model from %s" %(file_name), file = self.logger)
    if not self.model:
      raise Sorry("Missing model")

    self.model.add_crystal_symmetry_if_necessary()
    if self.params.input_files.selection:
      self.model = self.model.apply_selection_string(
        self.params.input_files.selection)

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

# =============================================================================
# for reference documentation keywords
master_phil_str = Program.master_phil_str
