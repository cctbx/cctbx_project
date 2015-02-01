from __future__ import division
import iotbx.phil
from libtbx.utils import Sorry
import sys, os

master_phil = iotbx.phil.parse("""
project_name = project
  .type = str
  .help = Project name specified as a path to store different runs.
  .optional = False
run_name = run
  .type = str
  .help = Run name is used as folder name that stores output files.
  .optional = False
title = None
  .type = str
  .help = Title of the run.
flag_log_verbose_on = False
    .type = bool
    .help = Turn this flag on to output more information.
hkl
  .help = Parameters related to hkl input file.
{
  phibin = None
    .type = path
    .help = Mtz file containing FP SIGFP PHIB FOM HLA HLB HLC HLD columns.
  refin = None
    .type = path
    .help = Mtz file used as a reference to calculate phase errors.
  column_names = FP,SIGFP,PHIB,FOM,HLA,HLB,HLC,HLD
    .type = str
    .help = List of column names if not FP, SIGFP, PHIB, FOM, HLA, HLB, HLC, HLD.
  column_phic = PHIC
    .type = str
    .help = Column name for phases in refin mtz file.
}
autobuild
  .help = Parameters for autobuild.
{
  flag_on = False
    .type = bool
    .help = Turn this flag on to run phenix.autobuil automatically.
  seq_file = None
    .type = path
    .help = Path to sequence file.
  ha_file = None
    .type = path
    .help = Path to heavy atom pdb file (preliminary results for previous searches).
  model_file = None
    .type = path
    .help = Path to pre-built model.
  flag_semet_on = False
    .type = bool
    .help = Turn this flag on to specify that the heavy atoms are semet.
}
ga
  .help = Genetic algorithm parameters.
{
  fom_min = 0.0
    .type = float
    .help = Minimum cut-off figure of merits.
  d_min = 0.1
    .type = float
    .help = Minimum resolution.
  d_max = 99.0
    .type = float
    .help = Maximum resolution.
  sigma_min = 1.5
    .type = float
    .help = Minimum I/sigI.
  n_stacks = 5
    .type = int
    .help = No. of stacks used in optimization (each stack contribute to specified percent_f_squared).
  n_micro_cycles = 10
    .type = int
    .help = No. of microcycles.
  n_macro_cycles = 3
    .type = int
    .help = No. of macrocycles.
  n_processors = 32
    .type = int
    .help = No. processors.
  percent_f_squared = 10
    .type = float
    .help = Percent of F^2 used in selecting no. of reflections in each stack.
  flag_excl_centric = False
    .type = bool
    .help = Turn this flag on to exclude centric reflections from the search.
  fit_params
    .help = search parameters
  {
    solvent_content = 0.5
      .type = float
      .help = Solvent content (fractions).
    wang_radius = 4.0
      .type = float
      .help = Wang radius
    ed_sigma_thres = 5.0
      .type = float
      .help = Electrion density sigma cut-off.
    w_llmap = 0.0
      .type = float
      .help = Weight of log-likelihood map score.
    w_skew = 0.0
      .type = float
      .help = Weight of map skew score.
  }
  ga_params
    .help = genetic algorithm parameters.
  {
    pop_size = 400
      .type = int
      .help = No. of chromosomes in the population.
    max_gen = 100
      .type = int
      .help = Maximum number of generations.
    prob_of_cross = 0.95
      .type = float
      .help = Probability of crossover.
    prob_of_mut = 0.01
      .type = float
      .help = Probability of mutation.
    ratio_cross = 0.2
      .type = float
      .help = Fractions used in crossing over.
    num_point_mut = 1
      .type = int
      .help = No. of mutation points.
    xmap_radius = 3
      .type = int
      .help = No. of cells used during map local search.
    num_sel_mate = 15
      .type = int
      .help = NA
    crossover_start_rate = 0.2
      .type = float
      .help = Starting crossover rate.
    crossover_end_rate = 0.2
      .type = float
      .help = Final crossover rate.
    crossover_slope = 1.0
      .type = float
      .help = Slope of the crossover rate function.
    skew_sigma_sel_lo = 0.0
      .type = float
      .help = NA
    skew_sigma_sel_hi = 1.0
      .type = float
      .help = NA
  }
}
""")


def process_input(argv=None):
  if argv == None:
    argv = sys.argv[1:]

  user_phil = []
  for arg in sys.argv[1:]:
    if os.path.isfile(arg):
      user_phil.append(iotbx.phil.parse(open(arg).read()))
    else :
      try :
        user_phil.append(iotbx.phil.parse(arg))
      except RuntimeError, e :
        raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))

  working_phil = master_phil.fetch(sources=user_phil)
  params = working_phil.extract()

  #generate run_no folder
  project_run_path = params.project_name + '/' + params.run_name
  if os.path.exists(project_run_path):
    import shutil
    shutil.rmtree(project_run_path)
    print 'Folder :', project_run_path, ' exists. This folder will be removed'

  os.makedirs(project_run_path)

  #capture input read out by phil
  from cStringIO import StringIO
  class Capturing(list):
    def __enter__(self):
      self._stdout = sys.stdout
      sys.stdout = self._stringio = StringIO()
      return self
    def __exit__(self, *args):
      self.extend(self._stringio.getvalue().splitlines())
      sys.stdout = self._stdout

  with Capturing() as output:
    working_phil.show()

  txt_out = 'sisa.optimize input:\n'
  for one_output in output:
    txt_out += one_output + '\n'

  return params, txt_out
