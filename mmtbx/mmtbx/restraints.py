import cctbx.adp_restraints
import math
from cctbx.array_family import flex
import scitbx.restraints
from cctbx import adptbx

class manager(object):

  def __init__(self,
        geometry=None,
        ncs_groups=None,
        normalization=False):
    self.geometry = geometry
    self.ncs_groups = ncs_groups
    self.normalization = normalization

  def energies_sites(self,
        sites_cart,
        geometry_flags=None,
        compute_gradients=False,
        gradients=None,
        disable_asu_cache=False,
        lock_for_line_search=False):
    result = scitbx.restraints.energies(
      compute_gradients=compute_gradients,
      gradients=gradients,
      gradients_size=sites_cart.size(),
      gradients_factory=flex.vec3_double,
      normalization=self.normalization)
    if (self.geometry is None):
      result.geometry = None
    else:
      result.geometry = self.geometry.energies_sites(
        sites_cart=sites_cart,
        flags=geometry_flags,
        compute_gradients=compute_gradients,
        gradients=result.gradients,
        disable_asu_cache=disable_asu_cache,
        lock_pair_proxies=lock_for_line_search,
        normalization=False)
      result += result.geometry
    if (self.ncs_groups is None):
      result.ncs_groups = None
    else:
      result.ncs_groups = self.ncs_groups.energies_sites(
        sites_cart=sites_cart,
        compute_gradients=compute_gradients,
        gradients=result.gradients,
        lock_operators=lock_for_line_search,
        normalization=False)
      result += result.ncs_groups
    result.finalize_target_and_gradients()
    return result

  def energies_adp_iso(self,
        xray_structure,
        parameters,
        wilson_b=None,
        compute_gradients=False,
        tan_b_iso_max=None,
        mean_displacements=None,
        gradients=None):
    result = scitbx.restraints.energies(
      compute_gradients=compute_gradients,
      gradients=gradients,
      gradients_size=xray_structure.scatterers().size(),
      gradients_factory=flex.double,
      normalization=self.normalization)
    if (self.geometry is None):
      result.geometry = None
    else:
      result.geometry = cctbx.adp_restraints.energies_iso(
        geometry_restraints_manager=self.geometry,
        xray_structure=xray_structure,
        parameters=parameters,
        wilson_b=wilson_b,
        compute_gradients=compute_gradients,
        gradients=result.gradients)
      result += result.geometry
    if (self.ncs_groups is None):
      result.ncs_groups = None
    else:
      result.ncs_groups = self.ncs_groups.energies_adp_iso(
        u_isos=xray_structure.extract_u_iso_or_u_equiv(),
        average_power=parameters.average_power,
        compute_gradients=compute_gradients,
        gradients=result.gradients)
      result += result.ncs_groups
    result.finalize_target_and_gradients()
    if(compute_gradients):
       if(tan_b_iso_max != 0):
          u_iso_max = adptbx.b_as_u(tan_b_iso_max)
          if(mean_displacements is not None):
             scale = u_iso_max / math.pi / (flex.pow2(mean_displacements)+1.0)
          else:
             mean_displacements = flex.tan(math.pi*(xray_structure.scatterers().extract_u_iso()/u_iso_max-1./2.))
             scale = u_iso_max / math.pi / (flex.pow2(mean_displacements)+1.0)

          #a =u_iso_max / math.pi / (flex.pow2(mean_displacements)+1.0)
          #b = u_iso_max / math.pi / (flex.pow2(flex.tan(math.pi*(xray_structure.scatterers().extract_u_iso()/u_iso_max-1./2.)))+1.0)
          #print a[0],b[0]
       else:
          scale = 1.0
       result.gradients = result.gradients * scale
    return result
