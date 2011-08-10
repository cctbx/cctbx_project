import cctbx.adp_restraints
import math
from cctbx.array_family import flex
import scitbx.restraints
from cctbx import adptbx
from libtbx.utils import Sorry

class manager(object):

  def __init__(self,
        geometry=None,
        ncs_groups=None,
        torsion_ncs_groups=None,
        normalization=False):
    self.geometry = geometry
    self.ncs_groups = ncs_groups
    self.torsion_ncs_groups = torsion_ncs_groups
    self.normalization = normalization

  def select(self, selection):
    if (self.geometry is None):
      geometry = None
    else:
      if(isinstance(selection, flex.size_t)):
        geometry = self.geometry.select(iselection=selection)
      else:
        geometry = self.geometry.select(selection=selection)
    if (self.ncs_groups is None):
      ncs_groups = None
    else:
      ncs_groups = self.ncs_groups.select(iselection=selection.iselection())
    if (self.torsion_ncs_groups is None):
      torsion_ncs_groups = None
    else:
      torsion_ncs_groups = \
        self.torsion_ncs_groups.select(iselection=selection.iselection())
    return manager(
      geometry=geometry,
      ncs_groups=ncs_groups,
      torsion_ncs_groups=torsion_ncs_groups,
      normalization=self.normalization)

  def energies_sites(self,
        sites_cart,
        geometry_flags=None,
        external_energy_function=None,
        custom_nonbonded_function=None,
        compute_gradients=False,
        gradients=None,
        disable_asu_cache=False):
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
        external_energy_function=external_energy_function,
        custom_nonbonded_function=custom_nonbonded_function,
        compute_gradients=compute_gradients,
        gradients=result.gradients,
        disable_asu_cache=disable_asu_cache,
        normalization=False)
      result += result.geometry
    if (self.ncs_groups is None):
      result.ncs_groups = None
    else:
      result.ncs_groups = self.ncs_groups.energies_sites(
        sites_cart=sites_cart,
        compute_gradients=compute_gradients,
        gradients=result.gradients,
        normalization=False)
      result += result.ncs_groups
    result.finalize_target_and_gradients()
    return result

  def energies_adp_iso(self,
        xray_structure,
        parameters,
        use_u_local_only,
        use_hd,
        wilson_b=None,
        compute_gradients=False,
        tan_b_iso_max=None,
        u_iso_refinable_params=None,
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
        use_hd = use_hd,
        use_u_local_only=use_u_local_only,
        compute_gradients=compute_gradients,
        gradients=result.gradients)
      result += result.geometry
    if (self.ncs_groups is not None and \
        self.torsion_ncs_groups is not None):
      raise Sorry("Cannot have both Cartesian and torsion NCS restraints"+\
                  " at the same time.")
    if (self.ncs_groups is None):
      result.ncs_groups = None
    else:
      result.ncs_groups = self.ncs_groups.energies_adp_iso(
        u_isos=xray_structure.extract_u_iso_or_u_equiv(),
        average_power=parameters.average_power,
        compute_gradients=compute_gradients,
        gradients=result.gradients)
      result += result.ncs_groups
    if (self.torsion_ncs_groups is None):
      result.torsion_ncs_groups = None
    else:
      result.torsion_ncs_groups = self.torsion_ncs_groups.energies_adp_iso(
        u_isos=xray_structure.extract_u_iso_or_u_equiv(),
        average_power=parameters.average_power,
        compute_gradients=compute_gradients,
        gradients=result.gradients)
      result += result.torsion_ncs_groups
    result.finalize_target_and_gradients()
    if(compute_gradients):
       #XXX highly inefficient code: do something asap by adopting new scatters flags
       if(tan_b_iso_max is not None and tan_b_iso_max != 0):
          u_iso_max = adptbx.b_as_u(tan_b_iso_max)
          if(u_iso_refinable_params is not None):
             chain_rule_scale = u_iso_max / math.pi / (flex.pow2(u_iso_refinable_params)+1.0)
          else:
             u_iso_refinable_params = flex.tan(math.pi*(xray_structure.scatterers().extract_u_iso()/u_iso_max-1./2.))
             chain_rule_scale = u_iso_max / math.pi / (flex.pow2(u_iso_refinable_params)+1.0)
       else:
          chain_rule_scale = 1.0
       result.gradients = result.gradients * chain_rule_scale
    return result

  def energies_adp_aniso(self,
        xray_structure,
        use_hd = None,
        selection = None,
        compute_gradients=False,
        gradients=None):
    result = cctbx.adp_restraints.adp_aniso_restraints(
        restraints_manager = self.geometry,
        xray_structure = xray_structure,
        selection = selection,
        use_hd = use_hd
        #compute_gradients=compute_gradients,
        #gradients=result.gradients
        )
    if(self.normalization):
       normalization_scale = 1.0
       if(result.number_of_restraints > 0):
          normalization_scale /= result.number_of_restraints
       result.target *= normalization_scale
       result.gradients_aniso_cart *= normalization_scale
       result.gradients_aniso_star *= normalization_scale
       if(result.gradients_iso is not None):
          result.gradients_iso *= normalization_scale
    return result
