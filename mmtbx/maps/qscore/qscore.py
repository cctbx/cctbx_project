from __future__ import absolute_import, division, print_function
from phenix.program_template import ProgramTemplate
from mmtbx.maps.qscore.qscore_lib import qscore_np


class Program(ProgramTemplate):

  description = """
  Program to identify calculate Q-score, a map-model validation metric
  Usage: phenix.qscore model.pdb map.mrc

  """

  # Define the data types that will be used
  datatypes = ['model', "real_map",'phil']

  # Input parameters

  # Note: include of map_model_phil_str allows specification of
  #   full_map, half_map, another half_map, and model and Program Template
  #   produces input_files.map_model_manager containing these

  master_phil_str = """

  input_files {
    include scope iotbx.map_model_manager.map_model_phil_str
    }
  

    n_probes = 16
      .type = int
      .help = number of probes to use
      .short_caption = "This many probes will be used when sampling the density in radial shells around each atom"

    algorithm_version = 2
      .type = int
      .help = "Q-score algorithm version"
      .short_caption = "Version 1 will guarantee n_probes per each atom. Version 2 is much faster, but does not.
    

  """

  # Define how to determine if inputs are ok
  def validate(self):

    # Expect exactly one map and one model. Stop if not the case
    self.data_manager.has_models(
      raise_sorry = True,
      expected_n  = 1,
      exact_count = True)
    mmm = self.data_manager.get_map_model_manager()
    q = qscore_np(mmm,n_probes=self.params.n_probes,version=self.params.algorithm_version)

  # Set any defaults
  def set_defaults(self):
    pass


  # Run the program
  def run(self):
    print("Running Hotspots...",file=self.logger)
    model = self._model
    
    
  # Get the results
  def get_results(self):
      return self.get_results_JSON()

  def get_results_JSON(self):
    json_str = self.hotspots.to_json()
    return json.dumps(json.loads(json_str),cls=NpEncoder)