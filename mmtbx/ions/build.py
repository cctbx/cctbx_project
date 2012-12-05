"""
Deals with modifying a structure to include unbuilt and misidentified ions.
"""

from __future__ import division
from libtbx.str_utils import make_sub_header
from libtbx import Auto
import sys

ion_building_params_str = """
debug = False
  .type = bool
elements = Auto
  .type = strings
ion_chain_id = X
  .type = str
initial_occupancy = 1.0
  .type = float
  .help = Occupancy for newly placed ions - if less than 1.0, the occupancy \
    may be refined automatically in future runs of phenix.refine.
refine_ion_occupancies = True
  .type = bool
  .help = Toggles refinement of occupancies for newly placed ions.  This \
    will only happen if the occupancy refinement strategy is selected.
refine_ion_adp = *Auto isotropic anisotropic none
  .type = choice
  .help = B-factor refinement type for newly placed ions.  At medium-to-high \
    resolution, anisotropic refinement may be preferrable for the heavier \
    elements.
"""

def find_and_build_ions (
      manager,
      fmodel,
      model,
      wavelength,
      params,
      nproc=1,
      out=None,
      run_ordered_solvent=False) :
  import mmtbx.ions
  from cctbx.eltbx import sasaki
  assert (1.0 >= params.initial_occupancy >= 0)
  if (out is None) : out = sys.stdout
  if (manager is None) :
    manager = mmtbx.ions.create_manager(
      pdb_hierarchy=model.pdb_hierarchy(sync_with_xray_structure=True),
      geometry_restraints_manager=model.restraints_manager.geometry,
      fmodel=fmodel,
      wavelength=wavelength,
      params=params,
      nproc=nproc,
      verbose=params.debug,
      log=out)
  else :
    grm = model.restraints_manager.geometry
    connectivity = grm.shell_sym_tables[0].full_simple_connectivity()
    manager.update_structure(
      pdb_hierarchy=model.pdb_hierarchy(sync_with_xray_structure=True),
      xray_structure=fmodel.xray_structure,
      connectivity=connectivity,
      log=out)
    manager.update_maps()
  make_sub_header("Analyzing water molecules", out=out)
  elements = params.elements
  # XXX somehow comma-separation of phil strings fields doesn't work
  if (isinstance(elements, list)) and (len(elements) == 1) :
    elements = elements[0].split(",")
  water_ion_candidates = manager.analyze_waters(
    out=out,
    candidates=elements)
  structure_was_modified = False
  # Build in the identified ions
  for i_seq, final_choices in water_ion_candidates :
    atom = manager.pdb_atoms[i_seq]
    if len(final_choices) > 1:
      # Ambiguous results
      pass
    elif (len(final_choices) == 1) :
      final_choice = final_choices[0]
      print >> out, "  %s becomes %s%+d" % \
          (atom.id_str(), final_choice.element, final_choice.charge)
      refine_adp = params.refine_ion_adp
      if (refine_adp == "Auto") :
        if (fmodel.f_obs().d_min() <= 1.5) :
          refine_adp = "anisotropic"
        elif (fmodel.f_obs().d_min() < 2.5) :
          atomic_number = sasaki.table(final_choice.element).atomic_number()
          if (atomic_number >= 19) :
            refine_adp = "anisotropic"
      # Modify the atom object
      model.convert_atom(
        i_seq=i_seq,
        scattering_type=final_choice.scattering_type(),
        atom_name=final_choice.element,
        element=final_choice.element,
        charge=final_choice.charge,
        residue_name=final_choice.element,
        initial_occupancy=params.initial_occupancy,
        chain_id=params.ion_chain_id,
        segid="ION",
        refine_adp=refine_adp,
        refine_occupancies=params.refine_ion_occupancies)
      structure_was_modified = True
  if (structure_was_modified) :
    fmodel.update_xray_structure(update_f_calc=True)
  return manager
