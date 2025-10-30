from __future__ import absolute_import, division, print_function
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate
from mmtbx import real_space_correlation
from libtbx.utils import Sorry
from cctbx.array_family import flex
import mmtbx.utils

def broadcast(m, log):
  print("-"*79, file=log)
  print(m, file=log)
  print("*"*len(m), file=log)

class Program(ProgramTemplate):

  description = '''
Compute map correlation coefficient given input PDB model and reflection data.

Examples:

  phenix.real_space_correlation m.pdb d.mtz
  phenix.real_space_correlation m.pdb d.mtz detail=atom
  phenix.real_space_correlation m.pdb d.mtz detail=residue
  phenix.real_space_correlation m.pdb d.mtz data_labels=FOBS
  phenix.real_space_correlation m.pdb d.mtz scattering_table=neutron
  phenix.real_space_correlation m.pdb d.mtz detail=atom use_hydrogens=true
  phenix.real_space_correlation m.pdb d.mtz map_1.type=Fc map_2.type="2mFo-DFc"

  phenix.real_space_correlation m.pdb d.mtz map_coefficients_label="2FOFCWT,PH2FOFCWT"
  phenix.real_space_correlation m.pdb d.ccp4
'''

  datatypes = ['model', 'phil', 'miller_array', 'real_map', 'map_coefficients']

  master_phil_str = real_space_correlation.master_params_str
  #data_manager_options = ['model_skip_expand_with_mtrix',
  #                        'model_skip_ss_annotations']

  # ---------------------------------------------------------------------------

  def validate(self):
    print('Validating inputs...\n', file=self.logger)
    if self.params.pdb_file_name is None:
      self.data_manager.has_models(
        raise_sorry = True,
        expected_n  = 1,
        exact_count = True)
      self.params.pdb_file_name = self.data_manager.get_default_model_name()

    has_miller = self.data_manager.has_miller_arrays()
    has_map = self.data_manager.has_real_maps()

    if self.params.reflection_file_name is None and self.params.map_file_name is None:
      if (has_miller and has_map) or not (has_miller or has_map):
        raise Sorry("Specify either data, map_coefficients or a map.")

    if self.data_manager.has_miller_arrays():
      if len(self.data_manager.get_miller_array_names())>1:
        raise Sorry('Supply only one data file or map_coefficients file.')
      miller_fn = self.data_manager.get_default_miller_array_name()
      if self.data_manager.get_miller_array_array_type(miller_fn) != 'complex':
        self.params.reflection_file_name = miller_fn
      else:
        self.params.map_coefficients_file_name = miller_fn
      if ( (self.params.map_coefficients_file_name is not None) and
         (self.params.map_coefficients_label is None) ):
        raise Sorry('Please specify map coefficient labels for %s.' %
                  self.params.map_coefficients_file_name)

    if self.data_manager.has_real_maps():
      if len(self.data_manager.get_real_map_names())>1:
        raise Sorry('Supply only one map file.')
      self.params.map_file_name = self.data_manager.get_default_real_map_name()

  # ---------------------------------------------------------------------------

  def run(self):

    broadcast(m="Input model file name: %s"%self.params.pdb_file_name,
      log=self.logger)

    m = self.data_manager.get_model(self.params.pdb_file_name)
    xray_structure = m.get_xray_structure()
    pdb_hierarchy = m.get_hierarchy()
    #pdb_hierarchy.atoms().reset_i_seq() # VERY important to do.
    mmtbx.utils.setup_scattering_dictionaries(
      scattering_table = self.params.scattering_table,
      xray_structure = xray_structure,
      d_min = None)
    xray_structure.show_summary(f=self.logger, prefix="  ")

    data_and_flags = None
    if (self.params.reflection_file_name is not None):
      broadcast(
        m="Input reflection file name: %s"%self.params.reflection_file_name,
        log=self.logger)
      data_and_flags = real_space_correlation.extract_data_and_flags(params = self.params)
      data_and_flags.f_obs.show_comprehensive_summary(f=self.logger, prefix="  ")

    if self.params.map_coefficients_file_name is not None:
      if (self.params.map_coefficients_label not in
        self.data_manager.get_map_coefficients_all_labels()):
        raise Sorry('%s labels for map_coefficients were not found in %s' %
                    (self.params.map_coefficients_label,
                     self.params.map_coefficients_file_name))
      print('map coefficients', file=self.logger)
      print('  Map labels:', self.params.map_coefficients_label, file=self.logger)
    elif self.params.map_file_name is not None:
      print('Using CCP4-format', self.params.map_file_name, file=self.logger)

    # check crystal symmetry
    cs1 = xray_structure.crystal_symmetry()
    if self.params.reflection_file_name:
      miller_fn = self.params.reflection_file_name
    elif self.params.map_coefficients_file_name:
      miller_fn = self.params.map_coefficients_file_name
    if miller_fn:
      mas = self.data_manager.get_miller_array(miller_fn).as_miller_arrays()
      for m in mas:
        cs2 = m.crystal_symmetry()
        break
    if self.params.map_file_name:
      real_map = self.data_manager.get_real_map(self.params.map_file_name)
      cs2 = real_map.crystal_symmetry()

    if (cs1.is_similar_symmetry(cs2) is False):
      print(cs1, cs2)
      raise Sorry('The symmetry of the model and map/data files is not similar')

    # create fmodel with f_obs (if available)
    fmodel = None
    if (self.params.reflection_file_name is not None):
      r_free_flags = data_and_flags.f_obs.array(
        data = flex.bool(data_and_flags.f_obs.size(), False))
      fmodel = mmtbx.utils.fmodel_simple(
        xray_structures     = [xray_structure],
        scattering_table    = self.params.scattering_table,
        f_obs               = data_and_flags.f_obs,
        r_free_flags        = r_free_flags)
      broadcast(m="R-factors, reflection counts and scales", log=self.logger)
      fmodel.show(log=self.logger, show_header=False)

    # compute cc
    results = real_space_correlation.simple(
      fmodel        = fmodel,
      pdb_hierarchy = pdb_hierarchy,
      params        = self.params,
      show_results  = True,
      log           = self.logger)

