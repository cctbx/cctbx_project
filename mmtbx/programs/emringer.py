"""Side-chain-directed model and map validation for 3D Electron Cryomicroscopy"""
from __future__ import absolute_import, division, print_function
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate
import os
import libtbx.phil
from libtbx.utils import Sorry
from libtbx import easy_pickle
import mmtbx.ringer.emringer

# =============================================================================

program_citations = libtbx.phil.parse('''
citation {
  article_id = emringer1
  authors = Barad BA, Echols N, Wang RY, Cheng Y, DiMaio F, Adams PD, Fraser JS
  title = Side-chain-directed model and map validation for 3D Electron Cryomicroscopy.
  journal = Nature Methods
  volume = 10
  pages = 943-46
  year = 2015
  doi_id = "10.1038/nmeth.3541"
  pmid = 26280328
  external = True
}

citation {
  article_id = emringer2
  authors = Lang PT, Ng HL, Fraser JS, Corn JE, Echols N, Sales M, Holton JM, Alber T
  title = Automated electron-density sampling reveals widespread conformational polymorphism in proteins.
  journal = Protein Sci.
  volume = 7
  pages = 1420-31
  year = 2010
  doi_id = ""
  pmid = 20499387
  external = True
}
''')

# =============================================================================

master_phil_str = '''
include scope libtbx.phil.interface.tracking_params
include scope mmtbx.ringer.emringer.master_params
map_label = 2FOFCWT,PH2FOFCWT
  .type = str
  .input_size = 200
  .short_caption = 2Fo-FC map labels
  .help = Labels for 2Fo-Fc map coefficients
show_gui = False
  .type = bool
output_base = None
  .type = str
output_dir = None
  .type = path
  .short_caption = Output directory
quiet = False
  .type = bool
  .short_caption = no graphs
  .help = Don't output files or graphs
'''

# =============================================================================

class Program(ProgramTemplate):

  description = '''
Program for calculating the EMRinger score.\n

Minimum required inputs:
  Model file
  Map file (or file with map coefficients)

How to run:
  phenix.emringer model.pdb map.ccp4
'''

  datatypes = ['model', 'real_map', 'phil', 'map_coefficients']

  citations = program_citations
  master_phil_str = master_phil_str


  # ---------------------------------------------------------------------------
  def validate(self):
    print('Validating inputs', file=self.logger)
    self.data_manager.has_models(raise_sorry=True)
    if not (self.data_manager.has_real_maps() or
        self.data_manager.has_map_coefficients()):
      raise Sorry("Supply a map file or a file with map coefficients.")
    elif (self.data_manager.has_real_maps() and
        self.data_manager.has_map_coefficients()):
      raise Sorry("Supply either a map file or a file with map coefficients.")

  # ---------------------------------------------------------------------------
  def run(self):
    map_inp = None
    miller_array = None

    print('Using model: %s' % self.data_manager.get_default_model_name(),
      file=self.logger)
    model = self.data_manager.get_model()

    if self.data_manager.has_map_coefficients():
      miller_arrays = self.data_manager.get_miller_arrays()
      miller_array = self.find_label(miller_arrays = miller_arrays)
      print('Using miller array: %s' % miller_array.info().label_string(),
        file=self.logger)
    elif self.data_manager.has_real_maps():
      print('Using map: %s' % self.data_manager.get_default_real_map_name(),
        file=self.logger)
      map_inp = self.data_manager.get_real_map()
      print("CCP4 map statistics:", file=self.logger)
      map_inp.show_summary(out=self.logger, prefix="  ")

    if (self.params.output_base is None):
      pdb_base = os.path.basename(self.data_manager.get_default_model_name())
      self.params.output_base = os.path.splitext(pdb_base)[0] + "_emringer"

    if not self.params.quiet:
      plots_dir = self.params.output_base + "_plots"
      if (not os.path.isdir(plots_dir)):
        os.makedirs(plots_dir)

    task_obj = mmtbx.ringer.emringer.emringer(
      model        = model,
      miller_array = miller_array,
      map_inp      = map_inp,
      params       = self.params,
      out          = self.logger)
    task_obj.validate()
    task_obj.run()
    self.results = task_obj.get_results()

    ringer_result = self.results.ringer_result

    if not self.params.quiet:
      # save as pickle
      easy_pickle.dump("%s.pkl" % self.params.output_base, ringer_result)
      print ('Wrote %s.pkl' % self.params.output_base, file=self.logger)
      # save as CSV
      csv = "\n".join([ r.format_csv() for r in ringer_result])
      open("%s.csv" % self.params.output_base, "w").write(csv)
      print ('Wrote %s.csv' % self.params.output_base, file=self.logger)

    scoring_result = self.results.scoring_result
    scoring_result.show_summary(out = self.logger)

    #rolling_result = self.results.rolling_result

  # It would be good to have central code for this
  # ---------------------------------------------------------------------------
  def find_label(self, miller_arrays):
    best_guess = None
    best_labels = []
    all_labels = []
    miller_array = None
    for array in miller_arrays:
      label = array.info().label_string().replace(" ", "")
      if (self.params.map_label is not None):
        if (label == self.params.map_label.replace(" ", "")):
          miller_array = array
          return miller_array
      elif (self.params.map_label is None):
        if (array.is_complex_array()):
          all_labels.append(label)
          if (label.startswith("2FOFCWT") or label.startswith("2mFoDFc") or
              label.startswith("FWT")):
            best_guess = array
            best_labels.append(label)
    if (miller_array is None):
      if (len(all_labels) == 0):
        raise Sorry("No valid (pre-weighted) map coefficients found in file.")
      elif (len(best_labels) == 0):
        raise Sorry("Couldn't automatically determine appropriate map labels. "+
          "Choices:\n  %s" % "  \n".join(all_labels))
      elif (len(best_labels) > 1):
        raise Sorry("Multiple appropriate map coefficients found in file. "+
          "Choices:\n  %s" % "\n  ".join(best_labels))
      elif (len(best_labels) == 1):
        miller_array = best_guess
        print("  Guessing %s for input map coefficients"% best_labels[0], file=self.logger)
        return miller_array

  # ---------------------------------------------------------------------------
  def get_results(self):
    return self.results

