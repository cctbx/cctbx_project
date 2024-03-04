'''
File to store PHIL migrations for the CCTBXParser

The original full path of the migrated PHIL is necessary for matching
'''
from __future__ import absolute_import, division, print_function

from libtbx.utils import null_out, Sorry

# =============================================================================
class PhilMigrator():

  name = 'PHIL migrator'
  description = None
  original_phil_str = None
  end_date = None

  def __init__(self, parser_class=None):
    # recursively use CCTXParser for finding deprecated PHIL
    if parser_class is None:
      from iotbx.cli_parser import CCTBXParser
      parser_class = CCTBXParser
    from libtbx.program_template import ProgramTemplate

    class program_class(ProgramTemplate):
      master_phil_str = self.original_phil_str

      def validate(self):
        return True

      def run(self):
        pass

    self.parser = parser_class(program_class=program_class,
                               custom_process_arguments=None,
                               unused_phil_raises_sorry=True,
                               logger=null_out())

  def matches(self, unused_phil):
    '''
    _summary_

    Parameters
    ----------
    unused_phil : str
      Unused PHIL string as returned by CCTBXParser

    Returns
    -------
      Returns True if the unused PHIL matches the original_phil in this
      migrator
    '''
    try:
      namespace = self.parser.parse_args([unused_phil])
    except Sorry as s:
      if 'Some PHIL parameters are not recognized' in str(s):
        return False
      else:
        raise
    return True

  def update(self, working_phil, master_phil):
    raise NotImplementedError('PHIL migrator is not implemented')

# =============================================================================
class ramachandran_restraints_phil(PhilMigrator):

  name = 'Ramachandran plot restraints'
  original_phil_str = '''pdb_interpretation.peptide_link {
  rama_weight = 1.0
    .type = float
    .short_caption = Ramachandran gradients weight
    .expert_level = 1
    .style = hidden
  scale_allowed = 1.0
    .type = float
    .short_caption = Rescale allowed region pseudo-energy by
    .style = hidden
  rama_potential = *oldfield emsley
    .type = choice(multi=False)
    .short_caption = Ramachandran potential
    .caption = Oldfield Coot
    .style = hidden
  oldfield
    .short_caption = Oldfield potential parameters
    .style = box auto_align hidden
  {
    esd = 10.0
      .type = float
      .expert_level = 2
      .short_caption = E.S.D.
    weight_scale = 1.0
      .type = float
      .expert_level = 2
    dist_weight_max = 10.0
      .type = float
      .expert_level = 2
    weight = None
      .type = float
      .expert_level = 2
    plot_cutoff = 0.027
      .type = float
      .expert_level = 2
  }
  rama_selection = None
    .type = atom_selection
    .short_caption = Atom selection for Ramachandran restraints
    .help = Selection of part of the model for which \
        Ramachandran restraints will be set up.
    .expert_level = 1
    .style = hidden
  restrain_rama_outliers = True
    .type = bool
    .help = Apply restraints to Ramachandran outliers
    .style = hidden
  restrain_rama_allowed = True
    .type = bool
    .help = Apply restraints to residues in allowed region on Ramachandran plot
    .style = hidden
  restrain_allowed_outliers_with_emsley = False
    .type = bool
    .help = In case of restrain_rama_outliers=True and/or restrain_rama_allowed=True \
      still restrain these residues with emsley. Make sense only in case of \
      using oldfield potential.
    .style = hidden
}
'''

  def update(self, working_phil, master_phil):
    # new PHIL
    new_params = working_phil.extract()
    w_params = new_params.pdb_interpretation.ramachandran_plot_restraints
    # old PHIL
    old_params = self.parser.working_phil.extract()
    params = old_params.pdb_interpretation.peptide_link
    # old params, make transfer
    w_params.selection = params.rama_selection
    # oldfield
    w_params.enabled = True
    w_params.oldfield.weight = \
        params.oldfield.weight if (params.oldfield.weight is None or params.oldfield.weight > 0) else 0
    w_params.oldfield.weight_scale = \
        1/(params.oldfield.esd**2) * params.oldfield.weight_scale
    w_params.oldfield.distance_weight_min = 2.0
    w_params.oldfield.distance_weight_max = params.oldfield.dist_weight_max

    # emsley
    w_params.emsley.weight = params.rama_weight
    w_params.emsley.scale_allowed = params.scale_allowed
    # strategy
    if params.rama_potential == 'oldfield':
      pass
    elif params.rama_potential == 'emsley':
      w_params.favored = 'emsley'
      w_params.allowed = 'emsley'
      w_params.outlier = 'emsley'
    if params.restrain_rama_outliers:
      w_params.outlier = params.rama_potential
    else:
      w_params.outlier = None
    if params.restrain_rama_allowed:
      w_params.allowed = params.rama_potential
    else:
      w_params.allowed = None
    if params.restrain_allowed_outliers_with_emsley:
      if not params.restrain_rama_allowed:
        w_params.allowed = 'emsley'
      if not params.restrain_rama_outliers:
        w_params.outlier = 'emsley'

    # update working_phil
    new_params.pdb_interpretation.ramachandran_plot_restraints = w_params
    return master_phil.format(python_object=new_params)

# ============================================================================
# end
