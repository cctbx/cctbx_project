import boost.python
ext = boost.python.import_ext("smtbx_structure_factors_direct_ext")

class _linearisation(object):

  def __init__(self, xray_structure):
    xs = xray_structure
    super(_linearisation, self).__init__(xs.n_parameters_XXX(),
                                         xs.unit_cell(),
                                         xs.space_group(),
                                         xs.scatterers(),
                                         xs.scattering_type_registry())
    self.xray_structure = xray_structure


class linearisation_of_f_calc_modulus_squared(
  _linearisation, ext.linearisation_of_f_calc_modulus_squared): pass

class linearisation_of_f_calc_modulus(
  _linearisation, ext.linearisation_of_f_calc_modulus): pass
