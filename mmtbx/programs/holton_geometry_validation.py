"""
Holton geometry validation score calculation.  Based on a
script by James Holton (untangle_score.csh).  See details in
mmtbx/validation/holton_geometry_validation.py
"""

from __future__ import absolute_import, division, print_function

import os
from mmtbx.validation.holton_geometry_validation import \
    holton_geometry_validation
from libtbx.program_template import ProgramTemplate

try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  pass

class Program(ProgramTemplate):
  prog = os.getenv('LIBTBX_DISPATCHER_NAME')
  description="""
%(prog)s file.pdb [params.eff] [options ...]

Options:

  model=input_file          input PDB file

Example:

  %(prog)s model=1ubq.pdb
""" % locals()


  master_phil_str = """
    model = None
      .type = path
      .optional = False
      .help = '''input PDB file'''

    other_models = None
      .type = path
      .help = "Models to compare with model"
      .multiple = True

    get_individual_residue_scores = None
      .type = bool
      .short_caption = Residue scores
      .help = Calculate individual residue scores in addition to overall score

    round_numbers = True
      .type = bool
      .short_caption = Round numbers
      .help = Round numbers before calculation

    worst_clash_is_one_plus_n_times_worst_clash = True
      .type = bool
      .short_caption = Scale clashes
      .help = Scale worst clash score by (1 + n) where n is total clashes

    clash_energy_add_n = True
      .type = bool
      .short_caption = Add N to clash energy
      .help = Add number of clashes to clash energy score

    minimum_nonbond_score_to_be_worst = -0.1
      .type = float
      .short_caption = Minimum worst nonbond_score
      .help = Only include worst nonbond score if it is at least this value

    minimum_nonbond_score_to_be_included_in_average = 0
      .type = float
      .short_caption = Minimum nonbond_score
      .help = Only include nonbond score in average if it is at least this value

    include_full_nonbond_score = True
      .type = bool
      .short_caption = Include full nonbond score
      .help = If set, add additional scoring term in which all nonbond \
              values (even those less than the minimums) are included

    keep_hydrogens = True
      .type = bool
      .short_caption = Keep hydrogens
      .help = If set, keep input hydrogens, but add any necessary riding H.

    ignore_cis_peptides = False
      .type = bool
      .short_caption = Ignore cis peptides
      .help = If set, ignore cis peptides. Otherwise (if False), penalize them.

    ignore_h_except_in_nonbond = True
      .type = bool
      .short_caption = Ignore H except in nonbond
      .help = If set, ignore H atoms except in nonbond term

    ignore_arg_h_nonbond = True
      .type = bool
      .short_caption = Ignore H in Arg
      .help = If set, ignore H atoms in Arginine

    ignore_bond_lengths_with_h = False
      .type = bool
      .short_caption = Ignore bond lengths with H
      .help = If set, ignore bond lengths involving H atoms

    ignore_water_h_bonds = False
      .type = bool
      .short_caption = Ignore water H bonds
      .help = If set, ignore nonbonded contacts with water hydrogens

    rotalyze_ramalyze_max_energy = 99
      .type = float
      .short_caption = Maximum ramalyze/rotalyze energy
      .help = Maximum value of energy for a bad rotamer or phi-psi combination

    overall_max_energy = None
      .type = float
      .short_caption = Maximum overall energy
      .help = Maximum value of energy

    omega_angle_sigma = 4
      .type = float
      .short_caption = Omega angle sigma
      .help = Sigma for omega angle
    cbetadev_sigma = 0.05
      .type = float
      .short_caption = CB position sigma
      .help = Sigma for CB position

    clashscore_ideal_dist = 3
      .type = float
      .short_caption = Clashscore ideal distance
      .help = Ideal distance in Lennard-Jones potential for clashscore result

    lj_dist_that_yields_zero = 6 # Distance for modified LJ to cross zero
      .type = float
      .short_caption = LJ distance yielding zero
      .help = Distance at which modified Lennard-Jones potential crosses zero

    const_shrink_donor_acceptor = 0
      .type = float
      .help = Allow contacts closer by const_shrink_donor_acceptor from \
                normal target for H-bonding atoms. NOTE: matches \
                behavior of Phenix pre-2024.
      .short_caption = Shrink donor acceptor distance

    remove_waters = None
      .type = bool
      .help = Remove waters before analysis
      .short_caption = Remove waters

    score_this_altloc_only = None
      .type = str
      .help = Score only this altloc if specified
      .short_caption = Score this altloc only

    softPnna_params {
      y0 = 1
      .type = float
      .help = Parameter y0 for softPnna calculation
      a2 = -0.0192266
      .type = float
      .help = Parameter a2 for softPnna calculation
      a1 = 0.751694
      .type = float
      .help = Parameter a1 for softPnna calculation
      a0 = 1.12482
      .type = float
      .help = Parameter a0 for softPnna calculation
      mx = 0.21805
      .type = float
      .help = Parameter mx for softPnna calculation
      my = 0.736621
      .type = float
      .help = Parameter my for softPnna calculation
     }


    include_random = True
      .type = bool
      .help = Estimate expected value for structure with normal errors
      .short_caption = Estimate expected value

    n_random = 20
      .type = int
      .help = Number of samples for estimation of expected values
      .short_caption = Number of samples for estimates

    random_seed = 171927
      .type = int
      .help = Random seed
      .short_caption = Random seed

    variable_params {
      clashscore_ideal_dist_values = 2.5 3 3 3 3.5
         .type = floats
         .help = Values to try for variables that are uncertain in comparison \
               between models
      lj_dist_that_yields_zero_values = 5 6 6 6 7
         .type = floats
      omega_angle_sigma_values = 4 5 6 7 8 9 10
         .type = floats
      cbetadev_sigma_values = .04  .05  .06
         .type = floats
      worst_clash_is_one_plus_n_times_worst_clash_values = 0 1
          .type = floats
      clash_energy_add_n_values = 0 1
          .type = floats
      minimum_nonbond_score_to_be_worst_values = -0.1  0.0
         .type = floats
      minimum_nonbond_score_to_be_included_in_average_values = 0 0 0 0  -1.1
         .type = floats
      include_full_nonbond_score_values = 0 1
         .type = floats
      ignore_arg_h_nonbond_values = 0 1
         .type = floats
      rotalyze_ramalyze_max_energy_values = 20  99 99 99  1000
         .type = floats
    }

"""


  datatypes = ['model','phil']

  def validate(self):
    self.set_defaults()
    self.data_manager.has_models(raise_sorry=True)

  def set_defaults(self):

    if not self.params.model:
      self.data_manager.has_models(raise_sorry=True)
      self.params.model = self.data_manager.get_default_model_name()
    self.model = self.data_manager.get_model(self.params.model)
    print("\nModel read from %s" %(self.params.model), file = self.logger)

    # decide if there are other models
    if not self.params.other_models:
      self.params.other_models = []
      for fn in self.data_manager.get_model_names():
        if fn != self.params.model:
          self.params.other_models.append(fn)
    if self.params.other_models:
      print("Other models to compare: %s" %(" ".join(self.params.other_models)),
        file = self.logger)

  def run(self):
    if self.params.other_models:  # Run comparison analysis
      return self.run_comparison()

    self.results = holton_geometry_validation(
      dm = self.data_manager,
     model = self.model,
     get_individual_residue_scores = self.params.get_individual_residue_scores,
     round_numbers = self.params.round_numbers,
     worst_clash_is_one_plus_n_times_worst_clash =
       self.params.worst_clash_is_one_plus_n_times_worst_clash,
     clash_energy_add_n = self.params.clash_energy_add_n,
     minimum_nonbond_score_to_be_worst =
         self.params.minimum_nonbond_score_to_be_worst,
     minimum_nonbond_score_to_be_included_in_average =
        self.params.minimum_nonbond_score_to_be_included_in_average,
     include_full_nonbond_score = self.params.include_full_nonbond_score,
     keep_hydrogens = self.params.keep_hydrogens,
     ignore_cis_peptides = self.params.ignore_cis_peptides,
     ignore_h_except_in_nonbond = self.params.ignore_h_except_in_nonbond,
     ignore_arg_h_nonbond = self.params.ignore_arg_h_nonbond,
     ignore_bond_lengths_with_h = self.params.ignore_bond_lengths_with_h,
     ignore_water_h_bonds = self.params.ignore_water_h_bonds,
     rotalyze_ramalyze_max_energy = self.params.rotalyze_ramalyze_max_energy,
     overall_max_energy = self.params.overall_max_energy,
     omega_angle_sigma = self.params.omega_angle_sigma,
     cbetadev_sigma = self.params.cbetadev_sigma,
     clashscore_ideal_dist = self.params.clashscore_ideal_dist,
     lj_dist_that_yields_zero = self.params.lj_dist_that_yields_zero,
     const_shrink_donor_acceptor = self.params.const_shrink_donor_acceptor,
     remove_waters = self.params.remove_waters,
     score_this_altloc_only = self.params.score_this_altloc_only,
     softPnna_params = self.params.softPnna_params,
     include_random = self.params.include_random,
     n_random = self.params.n_random,
      log =self.logger)

  def get_results(self):
    return self.results

  def get_results_as_JSON(self):
    return self.results.as_JSON(self.info_json)

  def run_comparison(self):
    """ Analyze model and other models with various values of parameters and
     determine which are convincingly better"""

    all_model_list = [self.model]
    for fn in self.params.other_models:
      all_model_list.append(self.data_manager.get_model(fn))
    self.params.other_models = []
    from copy import deepcopy
    base_params = deepcopy(self.params)

    base_params.include_random = False # skip this

    variable_params = {}
    for x in dir(self.params.variable_params):
      if x.startswith("__"): continue
      xx = x.replace("_values","")
      current_value = getattr(self.params,xx, None)
      value_list = getattr(self.params.variable_params, x)
      if current_value is None: continue
      elif current_value in [True, False]: # bool
        new_value_list = []
        for v in value_list:
          new_value_list.append(True if v==1 else False)
        value_list = new_value_list
      variable_params[xx] = value_list

    import random
    random.seed(self.params.random_seed)
    next_seed = random.randint(0,10000000)
    all_result_dict = {}
    for i in range(self.params.n_random):
      random.seed(next_seed)
      self.params = deepcopy(base_params)
      for key in list(variable_params.keys()):
        value_list = variable_params[key]
        n = len(value_list)
        k = random.randint(0,n-1)
        value = value_list[k]
        setattr(self.params,key,value)
      next_seed = random.randint(0,10000000)
      self.params.random_seed = next_seed
      next_seed = random.randint(0,10000000)

      # Ready to run on all the models
      all_result_dict[i] = []
      for model in all_model_list:
        self.model = model
        self.run()
        r = self.get_results()
        all_result_dict[i].append(r)
        print("RESULT:",self.params.omega_angle_sigma,r.sum_energy,model.info().file_name)
    value_list_dict = {}
    from scitbx.array_family import flex
    for k in range(len(all_model_list)):
      value_list = flex.double()
      for i in range(self.params.n_random):
        value_list.append(all_result_dict[i][k].sum_energy)
      value_list_dict[k] = value_list
      delta_list = value_list - value_list_dict[0]
      mean_delta = delta_list.min_max_mean().mean
      sd_delta = delta_list.standard_deviation_of_the_sample()
      print("Model %s:  delta: %.2f  SD %.2f" %(
        all_model_list[k].info().file_name, mean_delta, sd_delta),
         file = self.logger)

# =============================================================================
# for reference documentation keywords
master_phil_str = Program.master_phil_str
