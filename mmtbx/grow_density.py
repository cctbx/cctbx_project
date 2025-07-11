"""Local improvement of electron density map by growing density"""
from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import mmtbx.f_model
from mmtbx import utils
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from six.moves import cStringIO as StringIO
import iotbx.phil
from iotbx import crystal_symmetry_from_any
from cctbx import adptbx
from libtbx.utils import Sorry
import os, math
from libtbx import group_args
from iotbx import pdb
from cctbx import xray
import mmtbx.refinement.minimization
import scitbx.lbfgs
from mmtbx import maps
import mmtbx.masks
from six.moves import zip
from six.moves import range
from iotbx import extract_xtal_data


master_params_str = """\
pdb_file_name = None
  .type = str
  .multiple = False
  .help = PDB file name (model only).
external_da_pdb_file_name = None
  .type = str
  .multiple = False
  .help = PDB file name (dummy atoms only).
reflection_file_name = None
  .type = str
  .help = File with experimental data (most of formats: CNS, SHELX, MTZ, etc).
data_labels = None
  .type = str
  .help = Labels for experimental data.
refine = *occupancies *adp *sites
  .type = choice(multi=True)
  .help = Parameters of DA to refine.
stop_reset_occupancies_at_macro_cycle = 5
  .type = int
  .help = Reset occupancies during first stop_reset_occupancies_at_macro_cycle \
          macro-cycles.
stop_reset_adp_at_macro_cycle = 5
  .type = int
  .help = Reset ADPs during first stop_reset_adp_at_macro_cycle macro-cycles.
start_filtering_at_macro_cycle=10
  .type = int
  .help = Macro-cycle number at which filtering takes off.
sphere
  .multiple=True
  {
    center = None
      .type = floats(size=3)
      .help = Approximate coordinates of the center of problem density.
    radius = None
      .type = float
      .help = Sphere radius of where the dummy atoms will be placed.
  }
mode = build *build_and_refine
  .type = choice(multi=False)
high_resolution = None
  .type = float
low_resolution = None
  .type = float
output_file_name_prefix = None
  .type = str
initial_occupancy = 0.1
  .type = float
  .help = Starting occupancy value for a newly placed DA
initial_b_factor = None
  .type = float
  .help = Starting B-factor value for a newly palced DA. Of None, then \
          determined automatically as mean B.
atom_gap = 0.7
  .type = float
  .help = Thje gap between the atoms in the grid.
  .expert_level = 2
overlap_interval = 0.25
  .type = float
  .help = Grid interval size - the gap between two grids.
  .expert_level = 2
atom_type = H
  .type = str
  .help = Atom to use to make grid. Nitrogen is a good choice
atom_name = " DA "
  .type = str
residue_name = DUM
  .type = str
  .help = Residue name for a DA
scattering_table = *n_gaussian wk1995 it1992 neutron
  .type = choice
  .help = Scattering table for structure factors calculations
number_of_refinement_cycles = 20
  .type = int
  .help = Number of refinement cycles
number_of_minimization_iterations = 25
  .type = int
  .help = Number of refinement iterations
filter
  {
    b_iso_min = 1.0
      .type = float
      .help = Min B-factor value to remove DA.
    b_iso_max = 100.0
      .type = float
      .help = Max B-factor value to remove DA.
    occupancy_min = 0.1
      .type = float
      .help = Min occupancy value to remove DA.
    occupancy_max = 10.0
      .type = float
      .help = Max occupancy value to remove DA.
    inside_sphere_scale = 1.2
      .type = float
      .help = DA with distance larger than sphere.radius*inside_sphere_scale \
              from sphere.center will be removed.
    da_model_min_dist = 1.
      .type = float
      .help = Max allowable distance between DA and any model atom.
  }
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=False)

def create_da_xray_structures(xray_structure, params):
  def grid(sphere, gap, overlap):
    c = flex.double(sphere.center)
    x_start, y_start, z_start = c - float(sphere.radius)
    x_end, y_end, z_end       = c + float(sphere.radius)
    x_range = frange(c[0], c[0]+gap, overlap)
    y_range = frange(c[1], c[1]+gap, overlap)
    z_range = frange(c[2], c[2]+gap, overlap)
    return group_args(x_range = x_range, y_range = y_range, z_range = z_range)
  grids = []
  for sphere in params.sphere:
    grids.append(grid(sphere = sphere, gap = params.atom_gap,
      overlap = params.overlap_interval))
  initial_b_factor = params.initial_b_factor
  if(initial_b_factor is None):
    initial_b_factor = flex.mean(
      xray_structure.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.))
  da_xray_structures = []
  counter = 0
  for grid, sphere in zip(grids, params.sphere):
    cntr_g = 0 # XXX
    for x_start in grid.x_range:
      for y_start in grid.y_range:
        for z_start in grid.z_range:
          cntr_g += 1
          if(cntr_g>1): continue # XXX
          counter += 1
          new_center = flex.double([x_start, y_start, z_start])
          atom_grid = make_grid(
            center          = new_center,
            radius          = sphere.radius,
            gap             = params.atom_gap,
            occupancy       = params.initial_occupancy,
            b_factor        = initial_b_factor,
            atom_name       = params.atom_name,
            scattering_type = params.atom_type,
            resname         = params.residue_name)
          da_xray_structure = pdb_atoms_as_xray_structure(pdb_atoms = atom_grid,
            crystal_symmetry = xray_structure.crystal_symmetry())
          closest_distances_result = xray_structure.closest_distances(
            sites_frac      = da_xray_structure.sites_frac(),
            distance_cutoff = 5)
          selection = closest_distances_result.smallest_distances > 0
          selection &= closest_distances_result.smallest_distances < 1
          da_xray_structure = da_xray_structure.select(~selection)
          #print counter, da_xray_structure.scatterers().size()
          if(cntr_g==1): # XXX
            da_xray_structures.append(da_xray_structure)
  ###
  result = []
  for i, x1 in enumerate(da_xray_structures):
    for j, x2 in enumerate(da_xray_structures):
      if(x1 is not x2):
        closest_distances_result = x1.closest_distances(
          sites_frac      = x2.sites_frac(),
          distance_cutoff = 5) # XXX ???
        selection = closest_distances_result.smallest_distances > 0
        selection &= closest_distances_result.smallest_distances < params.atom_gap
        da_xray_structures[j] = x2.select(~selection)
  return da_xray_structures

def grow_density(f_obs,
                 r_free_flags,
                 xray_structure,
                 xray_structure_da,
                 params):
    if(xray_structure_da is not None):
      print("Using external DA...")
      da_xray_structures = [xray_structure_da]
    else:
      print("Start creating DAs...")
      da_xray_structures = create_da_xray_structures(xray_structure =
        xray_structure, params = params)
    n_da = 0
    for daxrs in da_xray_structures:
      n_da += daxrs.scatterers().size()
    print("Total number of dummy atoms:", n_da)
    #
    xray_structure_start = xray_structure.deep_copy_scatterers()
    #
    if(params.mode == "build_and_refine"):
      mask_params = mmtbx.masks.mask_master_params.extract()
      mask_params.ignore_hydrogens=True
      fmodel = mmtbx.f_model.manager(
        xray_structure = xray_structure_start,
        r_free_flags   = r_free_flags,
        target_name    = "ml",
        mask_params    = mask_params,
        f_obs          = f_obs)
      fmodel.update_all_scales()
      print("START R-work and R-free: %6.4f %6.4f"%(fmodel.r_work(),
        fmodel.r_free()))
    xray_structure_current = xray_structure_start.deep_copy_scatterers()
    da_sel_all = flex.bool(xray_structure_start.scatterers().size(), False)
    for i_model, da_xray_structure in enumerate(da_xray_structures):
      print("Model %d, adding %d dummy atoms" % (i_model,
        da_xray_structure.scatterers().size()))
      da_sel_refinable = flex.bool(xray_structure_current.scatterers().size(), False)
      xray_structure_current = xray_structure_current.concatenate(da_xray_structure)
      da_sel_all.extend(flex.bool(da_xray_structure.scatterers().size(), True))
      da_sel_refinable.extend(flex.bool(da_xray_structure.scatterers().size(), True))
      if(params.mode == "build_and_refine"):
        fmodel.update_xray_structure(update_f_calc=True, update_f_mask=False,
          xray_structure = xray_structure_current)
        #XXX print fmodel.r_work()
        #XXX assert 0 # update_f_mask=True would not work
        da_sel_all = refine_da(
          fmodel           = fmodel,
          selection        = da_sel_all,
          da_sel_refinable = da_sel_all,#XXXda_sel_refinable,
          params           = params)
        xray_structure_current = fmodel.xray_structure
    da_sel_all = flex.bool(xray_structure_start.scatterers().size(), False)
    da_sel_all.extend(flex.bool(xray_structure_current.scatterers().size()-
                            xray_structure_start.scatterers().size(), True))
    all_da_xray_structure = xray_structure_current.select(da_sel_all)
    ofn = params.output_file_name_prefix
    if(ofn is None):
      ofn = "DA.pdb"
    else:
      ofn = ofn+"_DA.pdb"
    ofn = open(ofn,"w")
    pdb_file_str = all_da_xray_structure.as_pdb_file()
    print(pdb_file_str, file=ofn)
    #
    if(params.mode == "build_and_refine"):
      map_type_obj = mmtbx.map_names(map_name_string = "2mFo-DFc")
      map_params = maps.map_and_map_coeff_master_params().fetch(
        maps.cast_map_coeff_params(map_type_obj)).extract()
      coeffs = maps.map_coefficients_from_fmodel(fmodel = fmodel,
        params = map_params.map_coefficients[0])
      lbl_mgr = maps.map_coeffs_mtz_label_manager(map_params =
        map_params.map_coefficients[0])
      mtz_dataset = coeffs.as_mtz_dataset(
        column_root_label = lbl_mgr.amplitudes(),
        label_decorator   = lbl_mgr)
      #
      f_calc_da = coeffs.structure_factors_from_scatterers(xray_structure =
        all_da_xray_structure).f_calc()
      mtz_dataset.add_miller_array(
        miller_array      = f_calc_da,
        column_root_label = "Fcalc_da")
      #
      mtz_object = mtz_dataset.mtz_object()
      output_map_file_name = "map_coeffs.mtz"
      if(params.output_file_name_prefix is not None):
        output_map_file_name = params.output_file_name_prefix+"_map_coeffs.mtz"
      mtz_object.write(file_name = output_map_file_name)
    print("Finished")

def refinery(fmodels, number_of_iterations, iselection, parameter):
  if(iselection.size()==0): return
  fmodels.fmodel_xray().xray_structure.scatterers().flags_set_grads(state=False)
  lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
    max_iterations = number_of_iterations)
  if(parameter == "occupancies"):
    fmodels.fmodel_xray().xray_structure.scatterers().flags_set_grad_occupancy(
      iselection = iselection)
  if(parameter == "adp"):
    fmodels.fmodel_xray().xray_structure.scatterers().flags_set_grad_u_iso(
      iselection = iselection)
  if(parameter == "sites"):
    fmodels.fmodel_xray().xray_structure.scatterers().flags_set_grad_site(
      iselection = iselection)
  try:
    minimized = mmtbx.refinement.minimization.lbfgs(
      fmodels                  = fmodels,
      lbfgs_termination_params = lbfgs_termination_params,
      collect_monitor          = False)
  except Exception as e:
    print("ERROR: Refinement failed... ")
    print(str(e))
    print("... carry on")
    return
  xrs = fmodels.fmodel_xray().xray_structure
  fmodels.update_xray_structure(xray_structure = xrs, update_f_calc=True,
    update_f_mask=False)
  assert minimized.xray_structure is fmodels.fmodel_xray().xray_structure
  fmodels.fmodel_xray().xray_structure.scatterers().flags_set_grads(state=False)

def reset_occupancies(fmodels, selection, params):
  xrs = fmodels.fmodel_xray().xray_structure
  occ = xrs.scatterers().extract_occupancies()
  sel = occ < params.filter.occupancy_min
  sel &= selection
  occ = occ.set_selected(sel, params.filter.occupancy_min)
  sel = occ > params.filter.occupancy_max
  sel &= selection
  occ = occ.set_selected(sel, params.filter.occupancy_max)
  xrs.set_occupancies(occ)
  fmodels.update_xray_structure(xray_structure = xrs, update_f_calc=True,
    update_f_mask=False)

def reset_adps(fmodels, selection, params):
  xrs = fmodels.fmodel_xray().xray_structure
  b = xrs.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.)
  sel = b > params.filter.b_iso_max
  sel &= selection
  b = b.set_selected(sel, params.filter.b_iso_max)
  sel = b < params.filter.b_iso_min
  sel &= selection
  b = b.set_selected(sel, params.filter.b_iso_min)
  xrs = xrs.set_b_iso(values=b)
  fmodels.update_xray_structure(xray_structure = xrs, update_f_calc=True,
    update_f_mask=False)

def show_refinement_update(fmodels, selection, da_sel_refinable, prefix):
  fmt1 = "%s Rwork= %8.6f Rfree= %8.6f Number of: non-DA= %d DA= %d all= %d"
  print(fmt1%(prefix, fmodels.fmodel_xray().r_work(),
    fmodels.fmodel_xray().r_free(),
    selection.count(False),selection.count(True),
    fmodels.fmodel_xray().xray_structure.scatterers().size()))
  occ = fmodels.fmodel_xray().xray_structure.scatterers().extract_occupancies()
  occ_da = occ.select(selection)
  if(occ_da.size()>0):
    occ_ma = occ.select(~selection)
    print("         non-da: occ(min,max,mean)= %6.3f %6.3f %6.3f"%(
      flex.min(occ_ma),flex.max(occ_ma),flex.mean(occ_ma)))
    print("             da: occ(min,max,mean)= %6.3f %6.3f %6.3f"%(
      flex.min(occ_da),flex.max(occ_da),flex.mean(occ_da)))
    b = fmodels.fmodel_xray().xray_structure.extract_u_iso_or_u_equiv()*\
      adptbx.u_as_b(1.)
    b_da = b.select(selection)
    b_ma = b.select(~selection)
    print("         non-da: ADP(min,max,mean)= %7.2f %7.2f %7.2f"%(
      flex.min(b_ma),flex.max(b_ma),flex.mean(b_ma)))
    print("             da: ADP(min,max,mean)= %7.2f %7.2f %7.2f"%(
      flex.min(b_da),flex.max(b_da),flex.mean(b_da)))
    print("da_sel_refinable:", da_sel_refinable.size(), da_sel_refinable.count(True))

def refine_occupancies(fmodels, selection, da_sel_refinable, params, macro_cycle):
  refinery(
    fmodels              = fmodels,
    number_of_iterations = params.number_of_minimization_iterations,
    iselection           = da_sel_refinable.iselection(),
    parameter            = "occupancies")
  if(params.stop_reset_occupancies_at_macro_cycle > macro_cycle):
    reset_occupancies(fmodels = fmodels, selection = selection, params = params)
  if(params.start_filtering_at_macro_cycle <= macro_cycle):
    fmodels, selection, da_sel_refinable = filter_da(fmodels, selection, da_sel_refinable, params)
  show_refinement_update(fmodels, selection, da_sel_refinable, "occ(%2d):"%macro_cycle)
  return fmodels, selection, da_sel_refinable

def refine_adp(fmodels, selection, da_sel_refinable, params, macro_cycle):
  refinery(
    fmodels              = fmodels,
    number_of_iterations = params.number_of_minimization_iterations,
    iselection           = da_sel_refinable.iselection(),
    parameter            = "adp")
  if(params.stop_reset_adp_at_macro_cycle > macro_cycle):
    reset_adps(fmodels = fmodels, selection = selection, params = params)
  if(params.start_filtering_at_macro_cycle <= macro_cycle):
    fmodels, selection, da_sel_refinable = filter_da(fmodels, selection, da_sel_refinable, params)
  show_refinement_update(fmodels, selection, da_sel_refinable, "adp(%2d):"%macro_cycle)
  return fmodels, selection, da_sel_refinable

def refine_sites(fmodels, selection, da_sel_refinable, params, macro_cycle):
  refinery(
    fmodels              = fmodels,
    number_of_iterations = params.number_of_minimization_iterations,
    iselection           = da_sel_refinable.iselection(),
    parameter            = "sites")
  if(params.start_filtering_at_macro_cycle <= macro_cycle):
    fmodels, selection, da_sel_refinable = filter_da(fmodels, selection, da_sel_refinable, params)
  show_refinement_update(fmodels, selection, da_sel_refinable, "xyz(%2d):"%macro_cycle)
  return fmodels, selection, da_sel_refinable

def filter_da(fmodels, selection, da_sel_refinable, params):
  xrs = fmodels.fmodel_xray().xray_structure
  xrs_d = xrs.select(selection)
  xrs_m = xrs.select(~selection)
  #
  occ = xrs_d.scatterers().extract_occupancies()
  adp = xrs_d.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.)
  sel  = occ < params.filter.occupancy_max
  sel &= occ > params.filter.occupancy_min
  sel &= adp < params.filter.b_iso_max
  sel &= adp > params.filter.b_iso_min
  uc = xrs_d.unit_cell()
  sphere_defined = True
  if(len(params.sphere)==0): sphere_defined = False
  for sphere in params.sphere:
    if([sphere.center, sphere.radius].count(None) != 0): sphere_defined = False
  for i_seq, site_frac_d in enumerate(xrs_d.sites_frac()):
    if(sphere_defined):
      inside = False
      for sphere in params.sphere:
        if([sphere.center, sphere.radius].count(None)==0):
          dist = uc.distance(site_frac_d, uc.fractionalize(sphere.center))
          if(dist <= sphere.radius*params.filter.inside_sphere_scale):
            inside = True
            break
      if(not inside): sel[i_seq] = False
    for site_frac_m in xrs_m.sites_frac():
      dist = uc.distance(site_frac_d, site_frac_m)
      if(dist < params.filter.da_model_min_dist):
        sel[i_seq] = False
        break
  #
  sel_all = flex.bool(xrs_m.scatterers().size(),True)
  sel_all.extend(sel)
  da_sel_refinable = da_sel_refinable.select(sel_all)
  #
  xrs_d = xrs_d.select(sel)
  xrs = xrs_m.concatenate(xrs_d)
  selection = flex.bool(xrs_m.scatterers().size(), False)
  selection.extend(flex.bool(xrs_d.scatterers().size(), True))
  fmodels.update_xray_structure(xray_structure = xrs, update_f_calc=True,
    update_f_mask=False)
  assert selection.size() == da_sel_refinable.size()
  #XXXX
  da_sel_refinable = selection
  #XXXX
  return fmodels, selection, da_sel_refinable

def refine_da(fmodel, selection, params, da_sel_refinable):
  fmodels = mmtbx.fmodels(fmodel_xray = fmodel)
  show_refinement_update(fmodels, selection, da_sel_refinable, "  START:")
  nrm = params.number_of_refinement_cycles
  assert nrm >= params.stop_reset_occupancies_at_macro_cycle
  assert nrm >= params.stop_reset_adp_at_macro_cycle
  assert nrm >= params.start_filtering_at_macro_cycle
  for macro_cycle in range(nrm):
    if("occupancies" in params.refine):
      fmodels, selection, da_sel_refinable = refine_occupancies(
        fmodels     = fmodels,
        selection   = selection,
        da_sel_refinable = da_sel_refinable,
        params      = params,
        macro_cycle = macro_cycle)
    if("adp" in params.refine):
      fmodels, selection, da_sel_refinable = refine_adp(
        fmodels     = fmodels,
        selection   = selection,
        da_sel_refinable = da_sel_refinable,
        params      = params,
        macro_cycle = macro_cycle)
    if("sites" in params.refine):
      fmodels, selection, da_sel_refinable = refine_sites(
        fmodels     = fmodels,
        selection   = selection,
        da_sel_refinable = da_sel_refinable,
        params      = params,
        macro_cycle = macro_cycle)
    assert fmodel.xray_structure is fmodels.fmodel_xray().xray_structure
  return selection

def pdb_atoms_as_xray_structure(pdb_atoms, crystal_symmetry):
  xray_structure = xray.structure(crystal_symmetry = crystal_symmetry)
  unit_cell = xray_structure.unit_cell()
  for atom in pdb_atoms:
    scatterer = xray.scatterer(
      label           = atom.name,
      site            = unit_cell.fractionalize(atom.xyz),
      b               = atom.b,
      occupancy       = atom.occ,
      scattering_type = atom.element)
    xray_structure.add_scatterer(scatterer)
  return xray_structure

def make_grid(center, radius, gap, occupancy, b_factor, atom_name,
              scattering_type, resname):
  x_start, y_start, z_start = center - float(radius)
  x_end, y_end, z_end       = center + float(radius)
  x_range = frange(x_start, x_end, gap)
  y_range = frange(y_start, y_end, gap)
  z_range = frange(z_start, z_end, gap)
  result = []
  counter = 1
  for x in x_range:
    for y in y_range:
      for z in z_range:
        d = math.sqrt((x-center[0])**2+(y-center[1])**2+(z-center[2])**2)
        if(d < radius):
          atom = iotbx.pdb.hierarchy.atom_with_labels()
          atom.serial  = counter
          atom.name    = atom_name
          atom.resname = resname
          atom.resseq  = counter
          atom.xyz     = (x,y,z)
          atom.occ     = occupancy
          atom.b       = b_factor
          atom.element = scattering_type
          result.append(atom)
          counter += 1
  return result

def frange(start, end=None, inc=None):
    if end == None:
        end = start + 0.0
        start = 0.0
    if inc == None:
        inc = 1.0
    L = []
    while 1:
        next = start + len(L) * inc
        if inc > 0 and next >= end:
            break
        elif inc < 0 and next <= end:
            break
        L.append(next)
    return L

def reflection_file_server(crystal_symmetry, reflection_files):
  return reflection_file_utils.reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files,
    err=StringIO())

def cmd_run(args, command_name):
  msg = """\

Tool for local improvement of electron density map.

How to use:
1: Run this command: phenix.grow_density
2: Copy, save into a file and edit the parameters shown between the
   lines *** below. Do not include *** lines.
3: Run the command with this parameters file:
   phenix.grow_density parameters.txt
"""
  if(len(args) == 0):
    print(msg)
    print("*"*79)
    master_params().show()
    print("*"*79)
    return
  else:
    if(not os.path.isfile(args[0]) or len(args)>1):
      print("Parameter file is expected at input. This is not a parameter file:\n", \
        args)
      print("Run phenix.grow_density without arguments for running instructions.")
      return
    processed_args = utils.process_command_line_args(args = args,
      master_params = master_params(), log = None)
    params = processed_args.params.extract()
    if(params.pdb_file_name is None):
      assert len(processed_args.pdb_file_names)==1
      params.pdb_file_name = processed_args.pdb_file_names[0]
    run(processed_args = processed_args, params = params)

def run(processed_args, params):
  if(params.scattering_table not in ["n_gaussian","wk1995",
     "it1992","neutron"]):
    raise Sorry("Incorrect scattering_table.")
  crystal_symmetry = None
  crystal_symmetries = []
  for f in [str(params.pdb_file_name), str(params.reflection_file_name)]:
    cs = crystal_symmetry_from_any.extract_from(f)
    if(cs is not None): crystal_symmetries.append(cs)
  if(len(crystal_symmetries) == 1): crystal_symmetry = crystal_symmetries[0]
  elif(len(crystal_symmetries) == 0):
    raise Sorry("No crystal symmetry found.")
  else:
    if(not crystal_symmetries[0].is_similar_symmetry(crystal_symmetries[1])):
      raise Sorry("Crystal symmetry mismatch between different files.")
    crystal_symmetry = crystal_symmetries[0]
  f_obs = None
  r_free_flags = None
  if(params.reflection_file_name is not None):
    reflection_file = reflection_file_reader.any_reflection_file(
      file_name = params.reflection_file_name, ensure_read_access = True)
    rfs = reflection_file_server(
      crystal_symmetry = crystal_symmetry,
      reflection_files = [reflection_file])
    parameters = extract_xtal_data.data_and_flags_master_params().extract()
    if(params.data_labels is not None):
      parameters.labels = [processed_args.data_labels]
    determine_data_and_flags_result = extract_xtal_data.run(
      reflection_file_server = rfs,
      parameters             = parameters,
      keep_going             = True)
    f_obs = determine_data_and_flags_result.f_obs
    r_free_flags = determine_data_and_flags_result.r_free_flags
    if(r_free_flags is None):
      r_free_flags=f_obs.array(data=flex.bool(f_obs.data().size(), False))
      test_flag_value=None
  xray_structure = iotbx.pdb.input(file_name =
    params.pdb_file_name).xray_structure_simple()
  xray_structure_da = None
  if(params.external_da_pdb_file_name is not None):
    xray_structure_da = iotbx.pdb.input(file_name =
      params.external_da_pdb_file_name).xray_structure_simple()
  if(f_obs is not None):
    f_obs = f_obs.resolution_filter(d_min = params.high_resolution,
      d_max = params.low_resolution)
    r_free_flags = r_free_flags.resolution_filter(d_min = params.high_resolution,
      d_max = params.low_resolution)
  #
  assert params.mode in ["build", "build_and_refine"]
  grow_density(f_obs             = f_obs,
               r_free_flags      = r_free_flags,
               xray_structure    = xray_structure,
               xray_structure_da = xray_structure_da,
               params            = params)
