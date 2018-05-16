from __future__ import division, print_function
from phenix.program_template import ProgramTemplate
import time, os
#import iotbx.phil
import libtbx.phil
from libtbx import easy_pickle
import mmtbx.ringer.emringer

#from libtbx.utils import null_out
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
map_coeffs = None
  .type = path
  .short_caption = Map coefficients
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
  Map file

How to run:
  phenix.emringer model.pdb map.ccp4
'''

  datatypes = ['model', 'real_map', 'phil']

  citations = program_citations
  master_phil_str = master_phil_str


  # ---------------------------------------------------------------------------
  def validate(self):
    print('Validating inputs', file=self.logger)
    self.data_manager.has_models(raise_sorry=True)
    #TODO here should be test if map or mtz, one of them, available
    self.data_manager.has_real_maps(raise_sorry=True)

  # ---------------------------------------------------------------------------
  def run(self):
    t0 = time.time()
    print('Using model: %s' % self.data_manager.get_default_model_name())
    print('Using map: %s' % self.data_manager.get_default_real_map_name())

    model = self.data_manager.get_model()
    map_inp = self.data_manager.get_real_map()

    if (self.params.output_base is None) :
      pdb_base = os.path.basename(self.data_manager.get_default_model_name())
      self.params.output_base = os.path.splitext(pdb_base)[0] + "_emringer"

    if not self.params.quiet:
      plots_dir = self.params.output_base + "_plots"
      if (not os.path.isdir(plots_dir)) :
        os.makedirs(plots_dir)

    task_obj = mmtbx.ringer.emringer.emringer(
      model      = model,
      map_coeffs = None,
      ccp4_map   = map_inp,
      params     = self.params,
      out        = self.logger,
      quiet      = self.params.quiet)
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

    rolling_result = self.results.rolling_result

#    if (self.params.show_gui) :
#      run_app(results)
#    else :
#      return (ringer_result, scoring_result, rolling_result)

  # ---------------------------------------------------------------------------
  def get_results(self):
    return self.results

