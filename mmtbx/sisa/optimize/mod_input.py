from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import iotbx.phil
from libtbx.utils import Sorry
import sys, os

master_phil = iotbx.phil.parse("""
data = None
  .type = path
  .help = Path to an input mtz file.
project_name = .
  .type = str
  .help = Project name specified as a path to store different runs.
run_name = None
  .type = str
  .help = Run name is used as folder name that stores output files.
title = None
  .type = str
  .help = Title of the run.
flag_log_verbose_on = False
  .type = bool
  .help = Turn this flag on to output more information.
hklrefin = None
  .type = path
  .help = Mtz file used as a reference to calculate phase errors.
column_names = FP,SIGFP,PHIB,FOM,HLA,HLB,HLC,HLD
  .type = str
  .help = List of column names if not FP, SIGFP, PHIB, FOM, HLA, HLB, HLC, HLD.
column_phic = PHIC
  .type = str
  .help = Column name for phases in refin mtz file.
autodm = True
  .type = bool
  .help = Set this flag to False to disable automatic density modification.
seq_file = None
  .type = path
  .help = Path to sequence file.
ha_file = None
  .type = path
  .help = Path to heavy atom pdb file (preliminary results for previous searches).
model_file = None
  .type = path
  .help = Path to pre-built model.
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
flag_apply_b_factor = False
  .type = bool
  .help = Turn this flag on to apply Wilson B-factor (wavelength is needed).
wavelength = 1.0
  .type = float
  .help = Wavelength.
n_stacks = 1
  .type = int
  .help = No. of stacks used in optimization (each stack contribute to specified percent_f_squared).
n_micro_cycles = 10
  .type = int
  .help = No. of microcycles.
n_macro_cycles = 1
  .type = int
  .help = No. of macrocycles.
n_processors = 32
  .type = int
  .help = No. processors.
percent_f_squared = 25
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
""")

txt_help = """**************************************************************************************************

Sisa helps resampling a new combination of experimental phases for a few strongest reflections.
The output mtz file can be used to perform density modification.

Usage: phenix.sisa parameter.phil

With this command, you can specify all parameters required by sisa in your parameter.phil file.
To obtain the template of these parameters, you can perform a dry run (simply run phenix.sisa).
You can then change the values of the parameters.

For feedback, please contact monarin@stanford.edu.

**************************************************************************************************

List of available parameters:
"""
def process_input(argv=None):
  if argv == None:
    argv = sys.argv[1:]

  user_phil = []
  for arg in sys.argv[1:]:
    if os.path.isfile(arg):
      if arg.lower().find('.mtz') > 0:
        user_phil.append(iotbx.phil.parse("""data=\"%s\"""" % arg))
      else:
        try:
          user_phil.append(iotbx.phil.parse(open(arg).read()))
        except RuntimeError as e :
          print('Error reading input: run phenix.sisa -h for help')
          raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))
    else :
      print(arg)
      if arg == '--help' or arg == '-h':
        print(txt_help)
        master_phil.show(attributes_level=1)
        exit()
      else:
        try :
          user_phil.append(iotbx.phil.parse(arg))
        except RuntimeError as e :
          print('Error reading input: run phenix.sisa -h for help')
          raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))

  working_phil = master_phil.fetch(sources=user_phil)
  params = working_phil.extract()

  #check dry-run
  if len(user_phil) == 0:
    master_phil.show()
    print('Use the above list of parameters to generate your input file (.phil). For more information, run phenix.sisa -h.')
    exit()

  if params.data is None:
    print('MTZ file with amplitudes, HL coefficients, and PHIB is required. For more information, run phenix.sisa -h.')
    exit()

  #capture input read out by phil
  from six.moves import cStringIO as StringIO
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

  txt_out = 'sisa.optimize (v.0 150326a) input:\n'
  for one_output in output:
    txt_out += one_output + '\n'

  print(txt_out)


  if params.autodm:
    if params.seq_file is None:
      raise Sorry("Sequence file is needed for automatic density modification. You can set autodm=False to disable this step.")

  #determine automatic run_no if run_name is not given
  run_seq_list = flex.int()
  if params.run_name is None:
    for dirs in os.walk(params.project_name):
      current_dir = dirs[0]
      if current_dir.find('Sisa_run_') > 0:
        try:
          current_dir_arr = current_dir.split('_')
          run_seq_list.append(int(current_dir_arr[len(current_dir_arr)-1]))
        except Exception:
          dummy = 1
    if len(run_seq_list) == 0:
      run_seq = 1
    else:
      run_seq = flex.max(run_seq_list) + 1
    params.run_name = 'Sisa_run_'+str(run_seq)

  #generate run_no folder only if the folder does not exist.
  project_run_path = params.project_name + '/' + params.run_name
  if os.path.exists(project_run_path):
    raise Sorry("The run number %s already exists."%project_run_path)

  os.makedirs(project_run_path)
  return params, txt_out
