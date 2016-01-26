from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.polder

import time
import sys
from cctbx.array_family import flex
import mmtbx.f_model
from mmtbx import utils
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from cStringIO import StringIO
import iotbx.phil
from iotbx import crystal_symmetry_from_any
from libtbx.utils import Sorry
import os, math
from iotbx import pdb
from cctbx import xray
from mmtbx import maps
from mmtbx import map_tools
import mmtbx.masks

#*******************
# questions
#
# - where should hierarchy be read in?
# run() or in the function where it is used?
#
# - careful with ligands which have 2 letters for ligand code (PI)
#
# - add possibility of two ligands --> atom selection!
#
# - add possibility of selecting ligand of one chain --> atom selection!
#
# - add feature enhanced map (FEM) function?
#
# - for large ligands, DA procedure is super slow
#
# - error if there are no Rfree flags


master_params_str = """\
pdb_file_name = None
  .type = str
  .multiple = False
  .help = PDB file name (model only).
external_da_pdb_file_name = None
  .type = str
  .multiple = False
  .help = PDB file name (dummy atoms only).
ligand_code = None
  .type = str
  .help = 3letter code for the ligand to be replaced by dummy atoms
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
sphere_radius = 5
  .type = float
  .help = radius of sphere around ligand atoms in which dummy atoms will be placed
high_resolution = None
  .type = float
low_resolution = None
  .type = float
output_file_name_prefix = None
  .type = str
atom_gap = 0.7
  .type = float
  .help = The gap between the atoms in the grid.
  .expert_level = 2
overlap_interval = 0.25
  .type = float
  .help = Grid interval size - the gap between two grids.
  .expert_level = 2
scattering_table = *n_gaussian wk1995 it1992 neutron
  .type = choice
  .help = Scattering table for structure factors calculations
"""


def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=False)

def remove_ligand(pdb_hierarchy, params):
  print "Start deleting ligand atoms..."
  ligand_coord=[]
  for ag in pdb_hierarchy.atom_groups():
    if (ag.resname == params.ligand_code):
      print "The following atoms of ligand %s will be deleted" % params.ligand_code
      for atom in ag.atoms():
        print atom.quote(), atom.xyz
        ligand_coord.append(atom.xyz)
        ag.remove_atom(atom=atom)
  return ligand_coord
  print '...done removing ligand atoms'

def create_da_xray_structures(xray_structure, params, ligand_coord):
  print 'Looping through atoms...'
  da_xray_structures = []
  counter = 0
  for coord in ligand_coord:
          counter += 1
          new_center = flex.double(coord)
          atom_grid = make_grid(
            center          = new_center,
            radius          = params.sphere_radius,
            gap             = params.atom_gap)
          da_xray_structure = pdb_atoms_as_xray_structure(pdb_atoms = atom_grid,
            crystal_symmetry = xray_structure.crystal_symmetry())
          closest_distances_result = xray_structure.closest_distances(
            sites_frac      = da_xray_structure.sites_frac(),
            distance_cutoff = 5)
          print counter, len(atom_grid)
          selection = closest_distances_result.smallest_distances > 0
          selection &= closest_distances_result.smallest_distances < 1
          da_xray_structure = da_xray_structure.select(~selection)
          print counter, da_xray_structure.scatterers().size() #delete
          da_xray_structures.append(da_xray_structure)
  ###
  result = []
  for i, x1 in enumerate(da_xray_structures):
    print 'now deleting superflouos from DA', i
    for j, x2 in enumerate(da_xray_structures):
      if(x1 is not x2):
        closest_distances_result = x1.closest_distances(
          sites_frac      = x2.sites_frac(),
          distance_cutoff = 5) # XXX ???
        selection = closest_distances_result.smallest_distances > 0
        selection &= closest_distances_result.smallest_distances < params.atom_gap
        da_xray_structures[j] = x2.select(~selection)
  return da_xray_structures



def polder(f_obs,
           r_free_flags,
           xray_structure,
           xray_structure_da,
           pdb_hierarchy,
           params):
  crystal_symmetry = xray_structure.crystal_symmetry()
  pdb_hierarchy.write_pdb_file(file_name="before.pdb" )
  ligand_coord = remove_ligand(pdb_hierarchy, params)
  pdb_hierarchy.write_pdb_file(file_name="after.pdb" )
  #print dir(pdb_hierarchy)
  xray_structure_noligand = pdb_hierarchy.extract_xray_structure(crystal_symmetry = crystal_symmetry)
  #print dir(xray_structure_noligand)
  print xray_structure_noligand.crystal_symmetry().show_summary()
  if(xray_structure_da is not None):
    print "Using external DA..."
    da_xray_structures = [xray_structure_da]
  else:
    print "Start creating DAs..."
    da_xray_structures = create_da_xray_structures(xray_structure =
        xray_structure_noligand, params = params, ligand_coord = ligand_coord)
  n_da = 0
  for daxrs in da_xray_structures:
    n_da += daxrs.scatterers().size()
  print "Total number of dummy atoms:", n_da
  xray_structure_start = xray_structure_noligand.deep_copy_scatterers() #XXXXXXXXXXXXXXX
  xray_structure_current = xray_structure_start.deep_copy_scatterers()
   # not so clear
  #da_sel_all = flex.bool(xray_structure_start.scatterers().size(), False)
  print "xray_structure", xray_structure.scatterers().size()
  print "start", xray_structure_start.scatterers().size()
  print "current", xray_structure_current.scatterers().size()
  for i_model, da_xray_structure in enumerate(da_xray_structures):
      print "Atom %d, added %d dummy atoms" % (i_model,
         da_xray_structure.scatterers().size())
      #da_sel_refinable = flex.bool(xray_structure_current.scatterers().size(), False)
      xray_structure_current = xray_structure_current.concatenate(da_xray_structure)
      print "now it is longer", xray_structure_current.scatterers().size()
      #da_sel_all.extend(flex.bool(da_xray_structure.scatterers().size(), True))
      #da_sel_refinable.extend(flex.bool(da_xray_structure.scatterers().size(), True))
  da_sel_all = flex.bool(xray_structure_start.scatterers().size(), False)
  da_sel_all.extend(flex.bool(xray_structure_current.scatterers().size()-
                            xray_structure_start.scatterers().size(), True))
  sel_all = flex.bool(xray_structure_current.scatterers().size(), True) # DL
  all_da_xray_structure = xray_structure_current.select(da_sel_all)
  sel_all_xray_structure = xray_structure_current.select(sel_all) # DL
  # not so clear
  #print dir(xray_structure_current)
  ofn = params.output_file_name_prefix
  if(ofn is None):
    ofn = "DA.pdb"
  else:
    ofn = ofn+"_DA.pdb"
  ofn = open(ofn,"w")
  pdb_file_str = all_da_xray_structure.as_pdb_file()
  print >> ofn, pdb_file_str
  pdbroot = params.pdb_file_name.split(".")[0]
  ofn_all = open(pdbroot+"_DA.pdb","w")
  pdb_file_str = sel_all_xray_structure.as_pdb_file()
  print >> ofn_all, pdb_file_str
  #
  print 'now calculating maps...'
  mask_params = mmtbx.masks.mask_master_params.extract()
  mask_params.ignore_zero_occupancy_atoms = False # This is critical to do!
  fmodel = mmtbx.f_model.manager(
    f_obs          = f_obs,
    r_free_flags   = r_free_flags,
    mask_params    = mask_params,
    xray_structure = xray_structure_current)
  fmodel.update_all_scales(update_f_part1=False) # XXX
  print fmodel.r_work(), fmodel.r_free()
  fmodel.show(show_header=False, show_approx=False)

  mc_diff = map_tools.electron_density_map(
    fmodel=fmodel).map_coefficients(
      map_type         = "mFo-DFc",
      isotropize       = True,
      fill_missing     = False)
  #
  mtz_dataset = mc_diff.as_mtz_dataset(column_root_label="mFo-DFc")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "polder_map_coeffs.mtz")

  print 'now calculating maps without dummy atoms'
  mask_params2 = mmtbx.masks.mask_master_params.extract()
  mask_params2.ignore_zero_occupancy_atoms = True # This time ignore DA atoms
  fmodel2 = mmtbx.f_model.manager(
    f_obs          = f_obs,
    r_free_flags   = r_free_flags,
    mask_params    = mask_params2,
    xray_structure = xray_structure_current)
  fmodel2.update_all_scales(update_f_part1=False) # XXX
  print fmodel2.r_work(), fmodel2.r_free()
  fmodel2.show(show_header=False, show_approx=False)

  mc_diff2 = map_tools.electron_density_map(
    fmodel=fmodel2).map_coefficients(
      map_type         = "mFo-DFc",
      isotropize       = True,
      fill_missing     = False)
  #
  mtz_dataset2 = mc_diff2.as_mtz_dataset(column_root_label="mFo-DFc")
  mtz_object2 = mtz_dataset2.mtz_object()
  mtz_object2.write(file_name = "no-polder_map_coeffs.mtz")


  ######
  # add something to calculate the map from initial model, in order to compare


def reflection_file_server(crystal_symmetry, reflection_files):
  return reflection_file_utils.reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files,
    err=StringIO())


def make_grid(center, radius, gap):
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
          atom.name    = " DA "
          atom.resname = "DUM"
          atom.resseq  = counter
          atom.xyz     = (x,y,z)
          atom.occ     = 0.0
          atom.b       = 10.0
          atom.element = "O"
          result.append(atom)
          counter += 1
  return result


# needs from cctbx import xray
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

# returns a list L = [start, start + 1 x inc, start + 2 x inc, ... ]
# until start + n x inc < end
def frange(start, end=None, inc=None):
    if end == None:
        end = start + 0.0
        start = 0.0
    if inc == None:
        inc = 1.0
    L = []
 # while 1 is faster than while True
    while 1:
        next = start + len(L) * inc
        if inc > 0 and next >= end:
            break
        elif inc < 0 and next <= end:
            break
        L.append(next)
    return L


def cmd_run(args, command_name):
  msg = """\

Tool for improvement of ligand omit map.

How to use:
1: Run this command: phenix.polder
2: Copy, save into a file (e.g. parameters.txt) and edit the parameters
   shown between the lines *** below. Do not include *** lines.
3: Run the command with this parameters file:
   phenix.polder parameters.txt
"""
  if(len(args) == 0):
    print msg
    print "*"*79
    master_params().show()
    print "*"*79
    return
  else:
    if(not os.path.isfile(args[0]) or len(args)>1):
      print "Parameter file is expected at input. This is not a parameter file:\n", \
        args
      print "Run phenix.polder without argumets for running instructions."
      return
    processed_args = utils.process_command_line_args(args = args,
      master_params = master_params(), log = None)
    params = processed_args.params.extract()
    if(params.pdb_file_name is None):
      assert len(processed_args.pdb_file_names)==1
      params.pdb_file_name = processed_args.pdb_file_names[0]
    run(processed_args = processed_args, params = params)

def run(processed_args, params):
  print "phenix.polder is running..."
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
    parameters = utils.data_and_flags_master_params().extract()
    if(params.data_labels is not None):
      parameters.labels = [processed_args.data_labels]
    determine_data_and_flags_result = utils.determine_data_and_flags(
      reflection_file_server  = rfs,
      parameters              = parameters,
      keep_going              = True,
      log                     = StringIO())
    f_obs = determine_data_and_flags_result.f_obs
    r_free_flags = determine_data_and_flags_result.r_free_flags
    if(r_free_flags is None):
      r_free_flags=f_obs.array(data=flex.bool(f_obs.data().size(), False))
      test_flag_value=None
    f_obs, r_free_flags = f_obs.common_sets(r_free_flags) #DL
  pdb_input = iotbx.pdb.input(file_name = params.pdb_file_name)
  pdb_hierarchy = pdb_input.construct_hierarchy()
  #print dir(pdb_hierarchy)
  xray_structure = pdb_input.xray_structure_simple()
  print xray_structure.crystal_symmetry().show_summary()
  xray_structure_da = None
  #pdb_inp = iotbx.pdb.hierarchy.input(file_name=params.pdb_file_name) #not read file twice
  if(params.external_da_pdb_file_name is not None):
    xray_structure_da = iotbx.pdb.input(file_name =
      params.external_da_pdb_file_name).xray_structure_simple()
  if(f_obs is not None):
    f_obs = f_obs.resolution_filter(d_min = params.high_resolution,
      d_max = params.low_resolution)
    r_free_flags = r_free_flags.resolution_filter(d_min = params.high_resolution,
      d_max = params.low_resolution)
  #assert params.mode in ["build", "build_and_refine"]
  polder(f_obs             = f_obs,
         r_free_flags      = r_free_flags,
         xray_structure    = xray_structure,
         xray_structure_da = xray_structure_da,
         pdb_hierarchy     = pdb_hierarchy,
         params            = params)





if(__name__ == "__main__"):
  t0 = time.time()
  cmd_run(
    args         = sys.argv[1:],
    command_name = "mmtbx.dl_polder")
  print "Time:", time.time()-t0
