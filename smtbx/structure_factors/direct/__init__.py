import boost.python
ext = boost.python.import_ext("smtbx_structure_factors_direct_ext")

class constructed_with_xray_structure(object):

  def __init__(self, xray_structure, *args, **kwds):
    xs = xray_structure
    args = (xs.unit_cell(),
            xs.space_group(),
            xs.scatterers(),
            xs.scattering_type_registry()) + args
    super(constructed_with_xray_structure, self).__init__(*args, **kwds)
    self.xray_structure = xray_structure


class f_calc_modulus_squared_with_std_trigonometry(
  constructed_with_xray_structure,
  ext.f_calc_modulus_squared_with_std_trigonometry):
  pass

class f_calc_modulus_squared_with_custom_trigonometry(
  constructed_with_xray_structure,
  ext.f_calc_modulus_squared_with_custom_trigonometry):
  pass

class f_calc_modulus_with_std_trigonometry(
  constructed_with_xray_structure,
  ext.f_calc_modulus_with_std_trigonometry):
  pass

class f_calc_modulus_with_custom_trigonometry(
  constructed_with_xray_structure,
  ext.f_calc_modulus_with_custom_trigonometry):
  pass


def f_calc_modulus_squared(xray_structure,
                           exp_i_2pi_functor=None):
  if exp_i_2pi_functor is None:
    return f_calc_modulus_squared_with_std_trigonometry(xray_structure)
  else:
    return f_calc_modulus_squared_with_custom_trigonometry(xray_structure,
                                                           exp_i_2pi_functor)

def f_calc_modulus(xray_structure,
                   exp_i_2pi_functor=None):
  if exp_i_2pi_functor is None:
    return f_calc_modulus_with_std_trigonometry(xray_structure)
  else:
    return f_calc_modulus_with_custom_trigonometry(xray_structure,
                                                   exp_i_2pi_functor)
