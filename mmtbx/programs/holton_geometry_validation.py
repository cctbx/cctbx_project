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

  def run(self):

    self.results = holton_geometry_validation(
      dm = self.data_manager,
      filename = self.params.model,
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
     softPnna_params = self.params.softPnna_params,
      log =self.logger)

  def get_results(self):
    return self.results

  def get_results_as_JSON(self):
    return self.results.as_JSON(self.info_json)

# =============================================================================
# for reference documentation keywords
master_phil_str = Program.master_phil_str
