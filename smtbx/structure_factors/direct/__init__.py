from __future__ import division
from __future__ import absolute_import, print_function
import boost.python
ext = boost.python.import_ext("smtbx_structure_factors_direct_ext")

class constructed_with_xray_structure(object):

  def __init__(self, xray_structure, table_file_name=None, reflections=None,
               *args, **kwds):
    xs = xray_structure
    if not table_file_name:
      if reflections:
          self.scatterer_contribution = ext.isotropic_scatterer_contribution(
            xs.scatterers(),
            xs.scattering_type_registry(),
            unit_cell=xs.unit_cell(),
            reflections=reflections)
      else:
          self.scatterer_contribution = ext.isotropic_scatterer_contribution(
            xs.scatterers(),
            xs.scattering_type_registry())
    else:
      self.scatterer_contribution = ext.table_based_scatterer_contribution.build(
        xs.scatterers(),
        table_file_name,
        xs.space_group(),
        xs.space_group().is_origin_centric())

    args = (xs.unit_cell(),
            xs.space_group(),
            xs.scatterers(),
            self.scatterer_contribution) + args
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
                           table_file_name=None,
                           reflections=None,
                           exp_i_2pi_functor=None):
  if exp_i_2pi_functor is None:
    return f_calc_modulus_squared_with_std_trigonometry(xray_structure,
                                                        table_file_name=table_file_name,
                                                        reflections=reflections)
  else:
    return f_calc_modulus_squared_with_custom_trigonometry(xray_structure,
                                                           exp_i_2pi_functor,
                                                           table_file_name=table_file_name,
                                                           reflections=reflections)
def f_calc_modulus(xray_structure,
                   exp_i_2pi_functor=None):
  if exp_i_2pi_functor is None:
    return f_calc_modulus_with_std_trigonometry(xray_structure)
  else:
    return f_calc_modulus_with_custom_trigonometry(xray_structure,
                                                   exp_i_2pi_functor)
#for tests
def generate_isc_table_file(file_name,
                            xray_structure,
                            indices):
  xs = xray_structure
  isc = ext.isotropic_scatterer_contribution(
    xs.scatterers(),
    xs.scattering_type_registry())
  with open(file_name, "w") as out:
    out.write("Title: generated from isotropic AFF")
    out.write("\nScatterers:")
    for sc in xs.scatterers():
      out.write(" %s" %sc.label)
    out.write("\nAD accounted: true")
    out.write("\nSymm:")
    for idx in indices:
      d_star_sq = xs.unit_cell().d_star_sq(idx)
      isc.at_d_star_sq(d_star_sq)
      out.write("\n%s %s %s" %(idx[0], idx[1], idx[2]))
      for sci in xrange(xs.scatterers().size()):
        val = isc.get(sci, idx)
        out.write(" %.6f,%.6f" %(val.real, val.imag))
