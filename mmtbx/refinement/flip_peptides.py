
# XXX liberally copied from mmtbx.refinement.fit_rotamers - it would probably
# be a good idea to consolidate these modules into one framework.

from mmtbx.refinement import fit_rotamers
from mmtbx.refinement import print_statistics
import mmtbx.refinement.real_space
import mmtbx.model
import mmtbx.restraints
import iotbx.phil
import iotbx.pdb
from scitbx.array_family import flex
from libtbx.utils import Sorry, Usage
import sys
import os

torsion_search_params_str = fit_rotamers.torsion_search_params_str

master_params_str = """\
validate_change = True
  .type = bool
filter_residual_map_value = 2.0
  .type = float
filter_2fofc_map = None
  .type = float
target_map = 2mFo-DFc
  .type = str
residual_map = mFo-DFc
  .type = str
model_map = Fc
  .type = str
ignore_water_when_flipping = True
  .type = bool
number_of_macro_cycles = 1
  .type = int
real_space_refine_overall = False
  .type = bool
residue_iteration
  .style = box auto_align
{
  poor_cc_threshold = 0.9
    .type = float
    .short_caption = Poor CC threshold
  real_space_refine_peptide = True
    .type = bool
  real_space_refine_window = 1
    .type = int
  real_space_refine_max_iterations = 25
    .type = int
  real_space_refine_target_weight = 100.
    .type = float
  torsion_grid_search = True
    .type = bool
    .expert_level = 2
  ignore_alt_conformers = True
    .type = bool
    .short_caption = Ignore alternate conformers
  %s
}
""" % torsion_search_params_str

def torsion_search_params () :
  return iotbx.phil.parse(input_string = torsion_search_params_str)

def master_params():
  return iotbx.phil.parse(input_string = master_params_str)

def extract_calpha_atoms (residue1_isel, residue2_isel, atom_names) :
  (ca1, ca2) = (None, None)
  for atom_i_seq in residue1_isel :
    atom_name = atom_names[atom_i_seq].strip()
    if (atom_name == "CA") :
      ca1 = atom_i_seq
  for atom_i_seq in residue2_isel :
    atom_name = atom_names[atom_i_seq].strip()
    if (atom_name == "CA") :
      ca2 = atom_i_seq
  return (ca1, ca2)

def extract_peptide_atoms (residue1_isel, residue2_isel, atom_names) :
  (c1, o1, n2, h2) = (None, None, None, None)
  isel = flex.size_t()
  for atom_i_seq in residue1_isel :
    atom_name = atom_names[atom_i_seq].strip()
    if (atom_name == "C") :
      c1 = atom_i_seq
      isel.append(atom_i_seq)
    elif (atom_name == "O") :
      o1 = atom_i_seq
      isel.append(atom_i_seq)
  for atom_i_seq in residue2_isel :
    atom_name = atom_names[atom_i_seq].strip()
    if (atom_name == "N") :
      n2 = atom_i_seq
      isel.append(atom_i_seq)
    elif (atom_name in ["H", "D"]) :
      h2 = atom_i_seq
      isel.append(atom_i_seq)
  if None in [c1, o1, n2] :
    return None
  return isel

def rotate_peptide (ca1, ca2, sites_cart, selected, angle_deg) :
  axis_point_1 = sites_cart[ca1]
  axis_point_2 = sites_cart[ca2]
  new_sites = flex.vec3_double()
  for i_seq in selected :
    xyz_new = fit_rotamers.rotate_point_around_axis(
      axis_point_1=sites_cart[ca1],
      axis_point_2=sites_cart[ca2],
      point=sites_cart[i_seq],
      angle_deg=angle_deg)
    new_sites.append(xyz_new)
  return new_sites

class select_map (fit_rotamers.select_map) :
  def is_refinement_needed(self, peptide_isel, sites_cart, atoms, cc_limit) :
    # TODO: add Ramachandran check
    for i_seq in peptide_isel :
      if (not atoms[i_seq].element.strip().lower() in ["h","d"]) :
        xyz = sites_cart[i_seq]
        sel_map = self.select(sites_cart = flex.vec3_double([xyz]))
        m1 = self.target_map_data.select(sel_map)
        m2 = self.model_map_data.select(sel_map)
        cc = flex.linear_correlation(x = m1, y = m2).coefficient()
        if(cc < cc_limit): return True
    return False

def rotation_search (residue_evaluator, ca1, ca2, sites_cart, peptide_isel,
    params) :
  angle_deg_best = None
  sites_cart_best = None
  angles = fit_rotamers.generate_range(
    start=180+params.range_start,
    stop=180+params.range_stop,
    step=params.step)
  for angle_deg in angles :
    sites_cart_flipped = rotate_peptide(
      ca1=ca1,
      ca2=ca2,
      sites_cart=sites_cart,
      selected=peptide_isel,
      angle_deg=angle_deg)
    if residue_evaluator.is_better(sites_cart=sites_cart_flipped) :
      angle_deg_best = angle_deg
      sites_cart_best = sites_cart_flipped
  return (sites_cart_best, angle_deg_best)

def residue_iteration (pdb_hierarchy,
                       xray_structure,
                       selection,
                       target_map_data,
                       model_map_data,
                       residual_map_data,
                       rsr_manager,
                       params,
                       log) :
  assert target_map_data.focus() == model_map_data.focus()
  assert target_map_data.all() == model_map_data.all()
  fmt1 = "               |--------START--------| |--------FINAL--------|"
  fmt2 = "     residue   map_cc 2mFo-DFc mFo-DFc map_cc 2mFo-DFc mFo-DFc " \
      "rotation"
  fmt3 = "  %12s%7.4f %8.2f %7.2f %6.4f  %7.2f %7.2f %4.1fdeg"
  print >> log, fmt1
  print >> log, fmt2
  unit_cell = xray_structure.unit_cell()
  map_selector = select_map(
    unit_cell  = xray_structure.unit_cell(),
    target_map_data = target_map_data,
    model_map_data = model_map_data)
  get_class = iotbx.pdb.common_residue_names_get_class
  n_other_residues = 0
  n_amino_acids_ignored = 0
  n_amino_acids_scored = 0
  sites_cart_start = xray_structure.sites_cart()
  atoms = pdb_hierarchy.atoms()
  atom_names = atoms.extract_name()
  result = []
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for conformer in chain.conformers() :
        residues = conformer.residues()
        for i_res in range(len(residues) - 1) :
          residue1 = residues[i_res]
          residue2 = residues[i_res+1]
          if (residue2.resseq_as_int() != (residue1.resseq_as_int() + 1)) :
            continue
          if(get_class(residue1.resname) == "common_amino_acid"):
            residue_id_str = residue1.id_str(suppress_segid=1)[-12:]
            residue1_iselection = residue1.atoms().extract_i_seq()
            residue2_iselection = residue2.atoms().extract_i_seq()
            (ca1, ca2) = extract_calpha_atoms(
              residue1_isel=residue1_iselection,
              residue2_isel=residue2_iselection,
              atom_names=atom_names)
            if (None in [ca1, ca2]) :
              continue # missing flanking C-alphas
            peptide_isel = extract_peptide_atoms(
              residue1_isel=residue1_iselection,
              residue2_isel=residue2_iselection,
              atom_names=atom_names)
            if (peptide_isel is None) :
              continue # don't waste time on incomplete peptides
            residues_iselection = flex.size_t()
            residues_iselection.extend(residue1_iselection)
            residues_iselection.extend(residue2_iselection)
            if(map_selector.is_refinement_needed(
                peptide_isel=peptide_isel,
                sites_cart=sites_cart_start,
                atoms=atoms,
                cc_limit=params.poor_cc_threshold)) :
              sites_cart_peptide = sites_cart_start.select(peptide_isel)
              rsel, rs = fit_rotamers.include_residue_selection(
                selection=selection,
                residue_iselection=residues_iselection)
              rev = fit_rotamers.rotamer_evaluator(
                sites_cart_start=sites_cart_peptide,
                unit_cell=unit_cell,
                two_mfo_dfc_map=target_map_data,
                mfo_dfc_map=residual_map_data)
              cc_start = map_selector.get_cc(
                sites_cart=sites_cart_peptide,
                residue_iselection=residues_iselection)
              rm = fit_rotamers.residue_rsr_monitor(
                residue_id_str=residue_id_str,
                selection=peptide_isel.deep_copy(),
                sites_cart=sites_cart_peptide.deep_copy(),
                twomfodfc=rev.t1_start,
                mfodfc=rev.t2_start,
                cc=cc_start)
              if params.torsion_grid_search :
                torsion_params = params.torsion_search
              else :
                torsion_params = torsion_search_params().extract()
                torsion_params.start = 0
                torsion_params.stop = 0
              sites_cart_best, angle_best = rotation_search(
                residue_evaluator=rev,
                ca1=ca1,
                ca2=ca2,
                sites_cart=sites_cart_start,
                peptide_isel=peptide_isel,
                params=params.torsion_search)
              if (angle_best is None) :
                continue
              cc_flipped = map_selector.get_cc(
                sites_cart=sites_cart_best,
                residue_iselection=residues_iselection)
              cc_end = None
              if (not params.real_space_refine_peptide) :
                sites_cart_start = sites_cart_start.set_selected(
                  peptide_isel, sites_cart_best)
                cc_end = cc_flipped
              else :
                tmp = sites_cart_start.set_selected(
                  peptide_isel, sites_cart_best)
                sites_cart_refined = rsr_manager.refine_restrained(
                  tmp.select(rsel), rsel, rs)
                if rev.is_better(sites_cart_refined) :
                  sites_cart_start = sites_cart_start.set_selected(
                    peptide_isel, sites_cart_best)
                  cc_end = map_selector.get_cc(
                    sites_cart=sites_cart_refined,
                    residue_iselection=residues_iselection)
                else : # undo changes
                  sites_cart_start = sites_cart_start.set_selected(
                    peptide_isel, sites_cart_peptide)
              if (cc_end is not None) :
                result.append(rm)
                print >> log, fmt3 % (
                  residue_id_str, cc_start, rev.t1_start, rev.t2_start,
                  cc_end, rev.t1_best, rev.t2_best, angle_best)
  xray_structure.set_sites_cart(sites_cart_start)
  return result

def run(fmodel,
        geometry_restraints_manager,
        pdb_hierarchy,
        solvent_selection,
        log = None,
        params = None):
  if(log is None): log = sys.stdout
  if(params is None): params = master_params().extract()
  print_statistics.make_sub_header(text="Peptide bond flips", out=log)
  restraints_manager = mmtbx.restraints.manager(
    geometry      = geometry_restraints_manager,
    normalization = True)
  model = mmtbx.model.manager(
    restraints_manager = restraints_manager,
    xray_structure = fmodel.xray_structure,
    pdb_hierarchy = pdb_hierarchy)
  backbone_selections = model.backbone_selections()
  if(params.ignore_water_when_flipping):
    selection = ~model.solvent_selection()
  else:
    selection = flex.bool(model.xray_structure.scatterers().size(), True)
  selection &= backbone_selections
  fmt = "Macro-cycle %2d: r_work=%6.4f r_free=%6.4f"
  print >> log, fmt%(0, fmodel.r_work(), fmodel.r_free())
  for macro_cycle in range(1,params.number_of_macro_cycles+1):
    target_map_data, fft_map_1 = fit_rotamers.get_map_data(
      fmodel=fmodel,
      map_type=params.target_map,
      kick=False)
    model_map_data,fft_map_2 = fit_rotamers.get_map_data(
      fmodel=fmodel,
      map_type=params.model_map)
    residual_map_data,fft_map_3 = fit_rotamers.get_map_data(
      fmodel=fmodel,
      map_type=params.residual_map,
      kick=False)
    if(params.filter_residual_map_value is not None): #XXX use filtering....
      map_sel = flex.abs(residual_map_data) < params.filter_residual_map_value
      residual_map_data = residual_map_data.set_selected(map_sel, 0)
    if(params.filter_2fofc_map is not None):
      map_sel = flex.abs(target_map_data) < params.filter_2fofc_map
      target_map_data = target_map_data.set_selected(map_sel, 0)
    rsr_manager = fit_rotamers.refiner(
      pdb_hierarchy=pdb_hierarchy,
      target_map=target_map_data,
      geometry_restraints_manager=geometry_restraints_manager,
      real_space_target_weight=params.residue_iteration.real_space_refine_target_weight,
      real_space_gradients_delta=fmodel.f_obs.d_min()/4,
      max_iterations=params.residue_iteration.real_space_refine_max_iterations)
    residue_rsr_monitor = residue_iteration(
      pdb_hierarchy     = pdb_hierarchy,
      xray_structure    = fmodel.xray_structure,
      selection         = selection,
      target_map_data   = target_map_data,
      model_map_data    = model_map_data,
      residual_map_data = residual_map_data,
      rsr_manager       = rsr_manager,
      params            = params.residue_iteration,
      log               = log)
    fmodel.update_xray_structure(update_f_calc=True, update_f_mask=True)
    print >> log, "1:", fmt%(macro_cycle, fmodel.r_work(), fmodel.r_free())
    del target_map_data, model_map_data, residual_map_data, fft_map_1, fft_map_2, fft_map_3
    if(params.real_space_refine_overall):
      assert model.xray_structure is fmodel.xray_structure
      rsr_params = mmtbx.refinement.real_space.master_params().extract()
      rsr_params.real_space_refinement.mode="diff_map"
      if(params.exclude_hydrogens and optimize_hd):
        hd_selection = fmodel.xray_structure.hd_selection()
        occupancies_cache= fmodel.xray_structure.scatterers().extract_occupancies()
        fmodel.xray_structure.set_occupancies(value=0, selection=hd_selection)
      mmtbx.refinement.real_space.run(
        fmodel = fmodel, # XXX neutron ?
        model  = model,
        params = rsr_params.real_space_refinement,
        log    = log)
      if(params.exclude_hydrogens and optimize_hd):
        fmodel.xray_structure.set_occupancies(value = occupancies_cache)
      assert model.xray_structure is fmodel.xray_structure
      fmodel.update_xray_structure(update_f_calc=True, update_f_mask=True)
      print >> log, "2:",fmt%(macro_cycle, fmodel.r_work(), fmodel.r_free())
    if params.validate_change :
      fit_rotamers.validate_changes(
        fmodel=fmodel,
        residue_rsr_monitor=residue_rsr_monitor,
        log=log)
    pdb_hierarchy.atoms().set_xyz(fmodel.xray_structure.sites_cart())
