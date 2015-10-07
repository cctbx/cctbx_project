
"""
Wrapper module for computing targets and gradients for restraints (or other
energy functions) on coordinates and B-factors; used in phenix.refine.
"""

from __future__ import division
import cctbx.adp_restraints
import math
from cctbx.array_family import flex
import scitbx.restraints
from cctbx import adptbx
from libtbx.utils import Sorry
import sys

class manager(object):
  """
  Central management of restraints on molecular geometry and ADPs.  This
  includes the standard stereochemistry restraints (i.e. Engh & Huber),
  non-crystallographic symmetry restraints, reference model restraints, and
  optionally.
  """

  def __init__(self,
        geometry=None,
        ncs_groups=None,
        torsion_ncs_groups=None,
        normalization=False,
        use_amber=False,
        use_sander=False,
        amber_structs=None,
        use_afitt=False, #afitt
        afitt_object=None,
        use_rosetta=False,
        rosetta_manager=None) :
    self.geometry = geometry
    self.ncs_groups = ncs_groups
    self.torsion_ncs_groups = torsion_ncs_groups
    self.normalization = normalization
    # amber
    self.use_amber = use_amber
    self.amber_structs = amber_structs
    #afitt
    self.use_afitt = use_afitt
    self.afitt_object = afitt_object
    #rosetta
    self.use_rosetta = use_rosetta
    self.rosetta_manager = rosetta_manager

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
      if(not isinstance(selection, flex.size_t)):
        selection = selection.iselection()
      ncs_groups = self.ncs_groups.select(iselection=selection)
    if (self.torsion_ncs_groups is None):
      torsion_ncs_groups = None
    else:
      torsion_ncs_groups = \
        self.torsion_ncs_groups.select(iselection=selection)

    return manager(
      geometry=geometry,
      ncs_groups=ncs_groups,
      torsion_ncs_groups=torsion_ncs_groups,
      normalization=self.normalization,
      use_amber=self.use_amber,
      amber_structs=self.amber_structs,
      )

  def energies_sites(self,
        sites_cart,
        geometry_flags=None,
        external_energy_function=None,
        custom_nonbonded_function=None,
        compute_gradients=False,
        gradients=None,
        force_restraints_model=False,
        disable_asu_cache=False,
        hd_selection=None,
        ):
    """
    Compute energies for coordinates.  Originally this just used the standard
    geometry restraints from the monomer library, but it has since been
    extended to optionally incorporate a variety of external energy functions.

    :returns: scitbx.restraints.energies object
    """
    result = scitbx.restraints.energies(
      compute_gradients=compute_gradients,
      gradients=gradients,
      gradients_size=sites_cart.size(),
      gradients_factory=flex.vec3_double,
      normalization=self.normalization)
    if (self.geometry is None):
      result.geometry = None
    else:
      if (self.use_amber and not force_restraints_model) :
        geometry_energy = self.geometry.energies_sites(
          sites_cart=sites_cart,
          flags=geometry_flags,
          external_energy_function=external_energy_function,
          custom_nonbonded_function=custom_nonbonded_function,
          compute_gradients=False,
          gradients=None,
          disable_asu_cache=disable_asu_cache,
          normalization=False)
        ##################################################################
        #                                                                #
        # AMBER CALL - Amber force field gradients and target            #
        #                                                                #
        ##################################################################
        import amber_adaptbx
        amber_geometry_manager = amber_adaptbx.geometry_manager(
          sites_cart=sites_cart,
          number_of_restraints=geometry_energy.number_of_restraints,
          gradients_factory=flex.vec3_double,
          amber_structs=self.amber_structs)
        result.geometry = amber_geometry_manager.energies_sites(
          crystal_symmetry = self.geometry.crystal_symmetry,
          compute_gradients = compute_gradients)

      elif (self.use_afitt and
            len(sites_cart)==self.afitt_object.total_model_atoms
            ):
        ##################################################################
        #                                                                #
        # AFITT CALL - OpenEye AFITT gradients and target                #
        #                                                                #
        ##################################################################
        from mmtbx.geometry_restraints import afitt
        result.geometry = self.geometry.energies_sites(
          sites_cart=sites_cart,
          flags=geometry_flags,
          external_energy_function=external_energy_function,
          custom_nonbonded_function=custom_nonbonded_function,
          compute_gradients=compute_gradients,
          gradients=result.gradients,
          disable_asu_cache=disable_asu_cache,
          normalization=False,
          )
        result = afitt.apply(result, self.afitt_object, sites_cart)
        result = afitt.adjust_energy_and_gradients(
          result,
          self,
          sites_cart,
          hd_selection,
          self.afitt_object,
        )
        result.target = result.residual_sum
        result.afitt_energy=result.residual_sum
      elif (self.use_rosetta) :
        ##################################################################
        #                                                                #
        # ROSETTA energy function - requires rosetta_adaptbx             #
        #                                                                #
        ##################################################################
        assert (self.rosetta_manager is not None)
        geometry_energy = self.geometry.energies_sites(
          sites_cart=sites_cart,
          flags=geometry_flags,
          external_energy_function=external_energy_function,
          custom_nonbonded_function=custom_nonbonded_function,
          compute_gradients=False,
          gradients=None,
          disable_asu_cache=disable_asu_cache,
          normalization=False)
        gradients_rosetta = flex.vec3_double(sites_cart.size(), (0,0,0))
        rosetta_energy = scitbx.restraints.energies(
          compute_gradients=compute_gradients,
          gradients=gradients_rosetta,
          gradients_size=sites_cart.size(),
          gradients_factory=flex.vec3_double,
          normalization=False)
        rosetta_energy.residual_sum = self.rosetta_manager.compute_target(
          sites_cart=sites_cart,
          gradients=gradients_rosetta)
        rosetta_energy.number_of_restraints = \
          geometry_energy.number_of_restraints
        result.geometry = rosetta_energy
      # default restraints manager
      else :
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
    """
    Compute target and gradients for isotropic ADPs/B-factors relative to
    restraints.

    :returns: scitbx.restraints.energies object
    """
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
    """
    Compute target and gradients for isotropic ADPs/B-factors relative to
    restraints.

    :returns: scitbx.restraints.energies object
    """
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

  def write_geo_file(self,
      sites_cart=None,
      site_labels=None,
      file_name=None,
      file_descriptor=sys.stdout,
      header="# Geometry restraints\n",
      # Stuff for outputting ncs_groups
      excessive_distance_limit = 1.5,
      xray_structure=None,
      processed_pdb_file=None):
    """
    This should make complete .geo file with geometry and NCS if present.
    Instead of sites_cart and site_labels one may pass xray_structure and they
    will be extracted from it.
    Content of header will be outputted before restraints in the file.
    If file_name is not specified, everything will be outputted to
    file_descriptor which is stdout by default. Instead of file_name one may
    directly pass file_descriptor and it will be used for output. The caller
    will have take care of saving, closing or flushing it separately.
    """

    outf_descriptor = None
    if file_name is None:
      outf_descriptor = file_descriptor
    else:
      outf_descriptor = open(file_name, "w")
    if xray_structure is not None:
      if sites_cart is None:
        sites_cart = xray_structure.sites_cart()
      if site_labels is None:
        site_labels = xray_structure.scatterers().extract_labels()
    self.geometry.write_geo_file(
      sites_cart=sites_cart,
      site_labels=site_labels,
      header=header,
      file_descriptor=outf_descriptor)
    n_excessive = 0
    if [self.ncs_groups, xray_structure, processed_pdb_file].count(None) == 0:
      n_excessive = self.ncs_groups.show_sites_distances_to_average(
          sites_cart=sites_cart,
          site_labels=site_labels,
          excessive_distance_limit=excessive_distance_limit,
          out=outf_descriptor)
      print >> outf_descriptor
      self.ncs_groups.show_adp_iso_differences_to_average(
          u_isos=xray_structure.extract_u_iso_or_u_equiv(),
          site_labels=site_labels,
          out=outf_descriptor)
      print >> outf_descriptor
      processed_pdb_file.show_atoms_without_ncs_restraints(
        ncs_restraints_groups=self.ncs_groups,
        out=outf_descriptor)
      print >> outf_descriptor
    if file_name is not None:
      outf_descriptor.close()
    if (n_excessive != 0):
      raise Sorry("Excessive distances to NCS averages:\n"
        + "  Please inspect the resulting .geo file\n"
        + "  for a full listing of the distances to the NCS averages.\n"
        + '  Look for the word "EXCESSIVE".\n'
        + "  The current limit is defined by parameter:\n"
        + "    refinement.ncs.excessive_distance_limit\n"
        + "  The number of distances exceeding this limit is: %d\n"
            % n_excessive
        + "  Please correct your model or redefine the limit.\n"
        + "  To disable this message completely define:\n"
        + "    refinement.ncs.excessive_distance_limit=None")
