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
initial_b_iso = 20.0
  .type = float
refine_ion_occupancies = True
  .type = bool
  .help = Toggles refinement of occupancies for newly placed ions.  This \
    will only happen if the occupancy refinement strategy is selected.
refine_ion_adp = *Auto isotropic anisotropic none
  .type = choice
  .help = B-factor refinement type for newly placed ions.  At medium-to-high \
    resolution, anisotropic refinement may be preferrable for the heavier \
    elements.
anomalous = True
  .type = bool
  .help = If True and the wavelength is specified, any newly placed ions will \
    have anomalous scattering factors set using theoretical values.  This is \
    unlikely to affect R-factors but should flatten the anomalous LLG map.
refine_anomalous = False
  .type = bool
"""

def find_and_build_ions (
      manager,
      fmodels,
      model,
      wavelength,
      params,
      nproc=1,
      out=None,
      run_ordered_solvent=False,
      group_anomalous_strategy_enabled=False) :
  import mmtbx.refinement.minimization
  from mmtbx.refinement.anomalous_scatterer_groups import \
    get_single_atom_selection_string
  from mmtbx.refinement import anomalous_scatterer_groups
  import mmtbx.ions
  from cctbx.eltbx import sasaki
  from cctbx import xray
  from scitbx.array_family import flex
  import scitbx.lbfgs
  assert (1.0 >= params.initial_occupancy >= 0)
  fmodel = fmodels.fmodel_xray()
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
  manager.show_current_scattering_statistics(out=out)
  elements = params.elements
  anomalous_groups = []
  # XXX somehow comma-separation of phil strings fields doesn't work
  if (isinstance(elements, list)) and (len(elements) == 1) :
    elements = elements[0].split(",")
  water_ion_candidates = manager.analyze_waters(
    out=out,
    candidates=elements)
  modified_iselection = flex.size_t()
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
      # FIXME this is really insufficient - I need to also group them into
      # a chain
      modified_atom = model.convert_atom(
        i_seq=i_seq,
        scattering_type=final_choice.scattering_type(),
        atom_name=final_choice.element,
        element=final_choice.element,
        charge=final_choice.charge,
        residue_name=final_choice.element,
        initial_occupancy=params.initial_occupancy,
        initial_b_iso=params.initial_b_iso,
        chain_id=params.ion_chain_id,
        segid="ION",
        refine_adp=refine_adp,
        refine_occupancies=False) #params.refine_ion_occupancies)
      if (params.anomalous) :
        scatterer = fmodel.xray_structure.scatterers()[i_seq]
        if (wavelength is not None) :
          fp_fdp_info = sasaki.table(final_choice.element).at_angstrom(
            wavelength)
          scatterer.fp = fp_fdp_info.fp()
          scatterer.fdp = fp_fdp_info.fdp()
          print >> out, "    setting f'=%g, f''=%g" % (scatterer.fp,
            scatterer.fdp)
        if (params.refine_anomalous) :
          group = xray.anomalous_scatterer_group(
            iselection=flex.size_t([i_seq]),
            f_prime=scatterer.fp,
            f_double_prime=scatterer.fdp,
            refine=["f_prime","f_double_prime"],
            selection_string=get_single_atom_selection_string(modified_atom))
          anomalous_groups.append(group)
      modified_iselection.append(i_seq)
  if (len(modified_iselection) > 0) :
    def show_r_factors () :
       return "r_work=%6.4f r_free=%6.4f" % (fmodel.r_work(), fmodel.r_free())
    fmodel.update_xray_structure(
      update_f_calc=True,
      update_f_mask=True)
    if (params.refine_ion_occupancies) :
      # XXX not ideal - need to determine whether the occupancy refinement
      # strategy will be used, and only refine here if it won't happen
      # otherwise.
      print >> out, "occupancy refinement (new ions only): start %s" % \
        show_r_factors()
      fmodel.xray_structure.scatterers().flags_set_grads(state = False)
      fmodel.xray_structure.scatterers().flags_set_grad_occupancy(
        iselection = modified_iselection)
      lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
        max_iterations = 25)
      minimized = mmtbx.refinement.minimization.lbfgs(
        restraints_manager       = None,
        fmodels                  = fmodels,
        model                    = model,
        is_neutron_scat_table    = False,
        lbfgs_termination_params = lbfgs_termination_params)
      fmodel.xray_structure.adjust_occupancy(
        occ_max   = 1.0,
        occ_min   = 0,
        selection = modified_iselection)
      fmodel.update_xray_structure(
        update_f_calc=True,
        update_f_mask=True)
      print >> out, "occupancy refinement (new ions only): final %s" % \
        show_r_factors()
    if (params.refine_anomalous) and (len(anomalous_groups) > 0) :
      if ((model.anomalous_scatterer_groups is not None) and
          (group_anomalous_strategy_enabled)) :
        model.anomalous_scatterer_groups.extend(anomalous_groups)
      else :
        print >> out, "anomalous refinement (new ions only): start %s" % \
          show_r_factors()
        anomalous_scatterer_groups.minimizer(
          fmodel=fmodel,
          groups=anomalous_groups)
        print >> out, "anomalous refinement (new ions only): final %s" % \
          show_r_factors()
  return manager

# XXX is there any circumstance in which this could reorder atoms?
def clean_up_ions (model, params) :
  atoms = model.pdb_hierarchy().atoms()
  resseq = 1
  ion_selection = model.pdb_hierarchy().atom_selection_cache().selection(
    "segid ION").iselection()
  if (len(ion_selection) == 0) :
    return
  for chain in model.pdb_hierarchy().only_model().chains() :
    if (chain.id == params.ion_chain_id) :
      for residue in chain.residue_groups() :
        residue.resseq = "%4d" % resseq
        resseq += 1
      for i_seq in ion_selection :
        atom = atoms[i_seq]
        if (atom.segid == "ION") :
          residue_group = atom.parent().parent()
          assert (len(residue_group.atoms()) == 1)
          rg_chain = residue_group.parent()
          if (rg_chain != chain) :
            rg_chain.remove_residue_group(residue_group)
            residue_group.resseq = "%4d" % resseq
            resseq += 1
            chain.append_residue_group(residue_group)
      break
