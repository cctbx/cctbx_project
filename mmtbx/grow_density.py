from cctbx.array_family import flex
import mmtbx.f_model
from cctbx.array_family import flex
from mmtbx import utils
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx.pdb import crystal_symmetry_from_pdb
from cStringIO import StringIO
from cctbx.xray import ext
from mmtbx import utils
from scitbx.array_family import shared
from cctbx import maptbx
import iotbx.phil
from iotbx import crystal_symmetry_from_any
from cctbx import adptbx
from libtbx.utils import Sorry
import os, math, time
from cctbx import miller
from mmtbx import map_tools
from libtbx import adopt_init_args
from libtbx import Auto, group_args
import mmtbx.bulk_solvent.bulk_solvent_and_scaling as bss
from iotbx import pdb
from libtbx.str_utils import format_value
from libtbx import smart_open
from cctbx import xray
from cctbx.array_family import flex
import time
from cctbx import adptbx
import mmtbx.tls.tools
from cctbx.development import random_structure
from cctbx import sgtbx
import mmtbx.refinement
import scitbx.lbfgs
import random, sys
from libtbx import group_args
import iotbx.pdb


master_params_str = """\
pdb_file_name = None
  .type = str
  .multiple = False
  .help = PDB file name.
reflection_file_name = None
  .type = str
  .help = File with experimental data (most of formats: CNS, SHELX, MTZ, etc).
data_labels = None
  .type = str
  .help = Labels for experimental data.
high_resolution = None
  .type = float
low_resolution = None
  .type = float

x_center = 0.5
  .type = float
  .help = X coordinate center of problem density. Determine in Coot by placing a temporary atom in area.
  .expert_level = 2

y_center = 0.5
  .type = float
  .help = X coordinate center of problem density. Determine in Coot by placing a temporary atom in area.
  .expert_level = 2

z_center = 0.5
  .type = float
  .help = X coordinate center of problem density. Determine in Coot by placing a temporary atom in area.
  .expert_level = 2

radius = 10
  .type = float
  .help = X coordinate center of problem density. Determine in Coot by placing a temporary atom in area.
  .expert_level = 2

atom_gap = 1.2
  .type = float
  .help = Thje gap between the atoms in the grid.
  .expert_level = 2

overlap_interval = 0.4
  .type = float
  .help = Grid interval size - the gap between two grids.
  .expert_level = 2

atom_type = N
    .type = str
    .help = Atom to use to make grid.  Nitrogen is a good choice


scattering_table = *n_gaussian wk1995 it1992 neutron
  .type = choice
  .help = Scattering table for structure factors calculations

cycles = 5
    .type = int
    .help = Number of refinement cycles

iterations = 20
    .type = int
    .help = Number of refinement iterations

bfac_dummy_atoms = 25
  .type = float
  .help = The initial B factor for the dummy atoms. Should be around protein average.
  .expert_level = 2

bfac_cutoff = 40.0
  .type = float
  .help = The value above which to remove atoms for the final model.
  .expert_level = 2

occ_cutoff = 1.0
  .type = float
  .help = The value below which to remove atoms for the final model.
  .expert_level = 2



"""
master_params = iotbx.phil.parse(master_params_str, process_includes=False)

def grow_density(f_obs,
                 r_free_flags,
                 file_name,
                 xray_structure,
                 x_center,
                 y_center,
                 z_center,
                 radius,

                overlap_interval=0.2,
                step_size=1,
                atom_type = "N",
                number_of_iterations = 20,
                number_of_cycles     = 5,

                bfac_dummy_atoms = 25,
                bfac_cutoff = 60,
                occ_cutoff = 1.0,

                sf_algorithm         = "fft",
                free_r_fraction      = 0.05,
                target_name          = "ls_wunit_kunit",

                 ):
    """ method to improve local density using dummy atoms """
    """TODO: trying to work out how to add atom correctly to pdb, so can create models /
    need to make new atom group?  Currently adds atoms to last residue (as shown in junk.pdb) /
    which is obviously not correct.  Code below will be able to produce atom grids /
    """
    step = overlap_interval
    atom_gap = step_size

    kept_atoms = []
    x_start = float(x_center) - float(radius)
    x_end   = float(x_center) + float(radius)
    y_start = float(y_center) - float(radius)
    y_end   = float(y_center) + float(radius)
    z_start = float(z_center) - float(radius)
    z_end   = float(z_center) + float(radius)
    #
    step_size = step_size
    overlap_interval = overlap_interval
    overlap_list = [0.0]

    kept_atoms = []

    center = [x_center,y_center,z_center]
    atom_count = xray_structure.scatterers().size()
    """still not sure this is right for making grid of grids"""
    x_range = frange(center[0],center[0]+atom_gap,step)
    y_range = frange(center[1],center[1]+atom_gap,step)
    z_range = frange(center[2],center[2]+atom_gap,step)
    all_atoms = []
    number_of_grids = len(x_range) * len(y_range) * len(z_range)
    print "Creating %s grids with atom spacing %s, each grid is %s apart" %(str(number_of_grids),str(atom_gap), str(step))
    print "Criteria for removing final atoms is b factor above %s, or occupancy below %s" %(str(bfac_cutoff), str(occ_cutoff))
    xray_structure_start = xray_structure
    fmodel = mmtbx.f_model.manager(
      xray_structure = xray_structure_start,
      r_free_flags   = r_free_flags,
      target_name    = "ml",
      f_obs          = f_obs)
    fmodel.update_solvent_and_scale()
    print "Start R-wok and R-free: %6.4f %6.4f"%(fmodel.r_work(), fmodel.r_free())
    counter = 0
    all_da_xray_structure = None
    for x_start in x_range:
      for y_start in y_range:
        for z_start in z_range:
          da_file_name = "tmp_%s.pdb"%str(counter)
          counter += 1
          new_center = [x_start, y_start, z_start]
          print "Number of grids left: ", number_of_grids
          number_of_grids = number_of_grids - 1
          atom_grid = make_grid(center = new_center, radius = radius,
            gap = atom_gap, occupancy = 0.1, b_factor = 10.)
          da_xray_structure = pdb_atoms_as_xray_structure(pdb_atoms = atom_grid,
            crystal_symmetry = xray_structure_start.crystal_symmetry())
          closest_distances_result = xray_structure_start.closest_distances(
            sites_frac      = da_xray_structure.sites_frac(),
            distance_cutoff = 3)
          selection = closest_distances_result.smallest_distances > 0
          selection &= closest_distances_result.smallest_distances < 1.2
          da_xray_structure = da_xray_structure.select(~selection)
          #print da_xray_structure.scatterers().size()
          xray_structure_current = xray_structure_start.concatenate(da_xray_structure)
          fmodel.update_xray_structure(update_f_calc=True, update_f_mask=True,
            xray_structure = xray_structure_current)
          da_sel = flex.bool(xray_structure_start.scatterers().size(), False)
          da_sel.extend(flex.bool(da_xray_structure.scatterers().size(), True))
          if 1:
            refine_atoms(
              fmodel               = fmodel,
              number_of_iterations = number_of_iterations,
              number_of_cycles     = number_of_cycles,
              selection            = da_sel)
            xray_structure_current = fmodel.xray_structure
          assert da_sel.size() == xray_structure_current.scatterers().size()
          if(all_da_xray_structure is None):
            all_da_xray_structure = xray_structure_current.select(da_sel)
          else:
            all_da_xray_structure = all_da_xray_structure.concatenate(
              xray_structure_current.select(da_sel))

    #occupancies = all_da_xray_structure.scatterers().extract_occupancies()
    #b_isos = all_da_xray_structure.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.)
    #selection = occupancies > 0.3
    #selection &= occupancies < 1.5
    #selection &= b_isos < 60.
    #selection &= b_isos > 1.0
    #all_da_xray_structure = all_da_xray_structure.select(selection)
    #print selection.size(), selection.count(True)

    ###
    ofn = open("zzz.pdb","w")
    #xrs_tmp = iotbx.pdb.input(file_name = "1akg.pdb").xray_structure_simple()
    #print >> ofn, all_da_xray_structure.concatenate(xrs_tmp).as_pdb_file()
    print >> ofn, all_da_xray_structure.as_pdb_file()
    print "Finished"


def refine_atoms(fmodel, number_of_iterations, number_of_cycles, selection):
  fmodels = mmtbx.fmodels(fmodel_xray = fmodel)
  selection = selection.iselection()
  lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
    max_iterations = number_of_iterations)
  for i in xrange(number_of_cycles):
    fmodel.xray_structure.scatterers().flags_set_grad_u_iso(
      iselection = selection)
    minimized = mmtbx.refinement.minimization.lbfgs(
      fmodels                  = fmodels,
      lbfgs_termination_params = lbfgs_termination_params,
      collect_monitor          = False)
    print "Refined B-factors     Rwork = %8.6f Rfree = %8.6f"%(fmodel.r_work(),
      fmodel.r_free())
    fmodel.xray_structure.scatterers().flags_set_grad_occupancy(
      iselection = selection)
    minimized = mmtbx.refinement.minimization.lbfgs(
      fmodels                  = fmodels,
      lbfgs_termination_params = lbfgs_termination_params,
      collect_monitor          = False)
    print "Refined occupancies   Rwork = %8.6f Rfree = %8.6f"%(fmodel.r_work(),
      fmodel.r_free())
  assert minimized.xray_structure is fmodels.fmodel_xray().xray_structure
  assert minimized.xray_structure is fmodel.xray_structure

def pdb_atoms_as_xray_structure(pdb_atoms, crystal_symmetry):
  xray_structure = xray.structure(crystal_symmetry = crystal_symmetry)
  unit_cell = xray_structure.unit_cell()
  for atom in pdb_atoms:
    scatterer = xray.scatterer(
      label           = atom.name,
      site            = unit_cell.fractionalize(atom.xyz),
      b               = atom.b,
      scattering_type = atom.element)
    xray_structure.add_scatterer(scatterer)
  return xray_structure

def make_grid(center, radius, gap, occupancy, b_factor, atom_name = " DA ",
              scattering_type = "N", resname = "DUM"):
  x_start = float(center[0]) - (float(radius))
  x_end   = float(center[0]) + (float(radius))
  y_start = float(center[1]) - (float(radius))
  y_end   = float(center[1]) + (float(radius))
  z_start = float(center[2]) - (float(radius))
  z_end   = float(center[2]) + (float(radius))
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
    """A range function, that does accept float increments..."""
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

Description:
Tool to compute local map correlation coefficient.

How to use:
1: Run this command: phenix.real_space_correlation;
2: Copy, save into a file and edit the parameters shown between the lines *** below;
3: Run the command with this parameters file:
   phenix.grow_density parameters.txt
"""
  if(len(args) == 0):
    print msg
    print "*"*79
    master_params.show()
    print "*"*79
    return
  else :
    arg = args[-1]
    if(not os.path.isfile(arg)):
      raise Sorry("%s is not a file."%arg)
    parsed_params = None
    try: parsed_params = iotbx.phil.parse(file_name=arg)
    except KeyboardInterrupt: raise
    except RuntimeError, e:
      print e
    params, unused_definitions = master_params.fetch(
        sources = [parsed_params],
        track_unused_definitions = True)
    if(len(unused_definitions)):
      print "*"*79
      print "ERROR:",
      print "Unused parameter definitions:"
      for obj_loc in unused_definitions:
        print " ", str(obj_loc)
      print "*"*79
      raise Sorry("Fix parameters file and run again.")
      print
    run(params = params.extract())

def run(params, d_min_default=1.5, d_max_default=999.9) :
  if(params.scattering_table not in ["n_gaussian","wk1995",
     "it1992","neutron"]):
    raise Sorry("Incorrect scattering_table.")
  if(params.x_center is None):
    raise Sorry("Need to specify x center.")
  if(params.y_center is None):
    raise Sorry("Need to specify y center.")
  if(params.z_center is None):
    raise Sorry("Need to specify z center.")
  if(params.radius is None):
    raise Sorry("Need to specify radius.")
  if(params.atom_gap is None):
    raise Sorry("Need to specify atom_gap.")
  if(params.overlap_interval is None):
    raise Sorry("Need to specify overlap_interval.")

  if(params.atom_type is None):
    raise Sorry("Need to specify atom_type.")
  if(params.cycles is None):
    raise Sorry("Need to specify cycles.")
  if(params.iterations is None):
    raise Sorry("Need to specify iterations.")
  if(params.bfac_dummy_atoms is None):
    raise Sorry("Need to specify bfac_dummy_atoms.")
  if(params.bfac_cutoff is None):
    raise Sorry("Need to specify bfac_cutoff.")
  if(params.occ_cutoff is None):
    raise Sorry("Need to specify iterations.")

  crystal_symmetry = None
  crystal_symmetry_data = None
  crystal_symmetry_model = None
  hkl_file_name = None
  pdb_file_name = None
  reflection_file = None
  files = [str(params.pdb_file_name), str(params.reflection_file_name)]

  print params.pdb_file_name, params.reflection_file_name
  for arg in files:
    arg_is_processed = False
    if(not os.path.isfile(arg)):
      raise Sorry(" %s is not a file."%arg)
    else:
      if(pdb.is_pdb_file(file_name=arg)):
        pdb_file_name = arg
        arg_is_processed = True
        try:
          crystal_symmetry_model = crystal_symmetry_from_pdb.extract_from(
            file_name=arg)
        except RuntimeError, e:
          if(str(e) == "No CRYST1 record."): pass
      if(not arg_is_processed):
        reflection_file = reflection_file_reader.any_reflection_file(
          file_name=arg, ensure_read_access=False)
        if(reflection_file.file_type() is not None):
          hkl_file_name = arg
          arg_is_processed = True
          crystal_symmetry_data = crystal_symmetry_from_any.extract_from(arg)
      if(not arg_is_processed):
        raise Sorry(
          "The command line argument %s is not a valid PDB or data file."%arg)
  if(hkl_file_name is None): raise Sorry("No file with Fobs is given.")
  if(pdb_file_name is None): raise Sorry("No PDB file is given.")
  print "Model and data files: %s %s"%(
    format_value("%5s",os.path.basename(pdb_file_name)),
    format_value("%5s",os.path.basename(hkl_file_name)))
  if([crystal_symmetry_model, crystal_symmetry_data].count(None)==0):
    if([crystal_symmetry_model.space_group_info(),
        crystal_symmetry_data.space_group_info()].count(None)==0):
      if(not crystal_symmetry_model.is_similar_symmetry(crystal_symmetry_data)):
        raise Sorry("Crystal symmetry mismatch between data and PDB files.")
  if(crystal_symmetry_model is not None):
    crystal_symmetry = crystal_symmetry_model
  if(crystal_symmetry_data is not None):
    crystal_symmetry = crystal_symmetry_data
  if(crystal_symmetry is None): raise Sorry("Crystal symmetry is not defined.")
  if(pdb_file_name is None): raise Sorry("No PDB file given.")
  if(reflection_file is None): raise Sorry("No reflection file given.")
  rfs = reflection_file_server(
    crystal_symmetry = crystal_symmetry,
    reflection_files = [reflection_file])
  parameters = utils.data_and_flags.extract()
  determine_data_and_flags_result = utils.determine_data_and_flags(
    reflection_file_server  = rfs,
    parameters              = parameters,
    data_parameter_scope    = "refinement.input.xray_data",
    flags_parameter_scope   = "refinement.input.xray_data.r_free_flags",
    data_description        = "X-ray data",
    keep_going              = True,
    log                     = StringIO())
  f_obs = determine_data_and_flags_result.f_obs
  number_of_reflections = f_obs.indices().size()
  #
  r_free_flags = determine_data_and_flags_result.r_free_flags
  test_flag_value = determine_data_and_flags_result.test_flag_value
  if(r_free_flags is None):
    r_free_flags=f_obs.array(data=flex.bool(f_obs.data().size(), False))
    test_flag_value=None
  #
  mmtbx_pdb_file = mmtbx.utils.pdb_file(
    pdb_file_names = [pdb_file_name],
    cryst1         = pdb.format_cryst1_record(crystal_symmetry = crystal_symmetry),
    log            = sys.stdout)
  mmtbx_pdb_file.set_ppf()
  processed_pdb_file = mmtbx_pdb_file.processed_pdb_file
  pdb_raw_records = mmtbx_pdb_file.pdb_raw_records
  pdb_inp = mmtbx_pdb_file.pdb_inp
  #
  xsfppf = mmtbx.utils.xray_structures_from_processed_pdb_file(
    processed_pdb_file = processed_pdb_file,
    scattering_table   =params.scattering_table,
    d_min              = f_obs.d_min())
  xray_structure = xsfppf.xray_structures[0]
  assert len(xsfppf.xray_structures) == 1
  #
  print "  Unit cell:       ", f_obs.unit_cell()
  print "  Space group:     ", f_obs.crystal_symmetry().space_group_info().\
    symbol_and_number(), "number of symmetry operators:",\
    f_obs.crystal_symmetry().space_group_info().type().group().order_z()
  print "  Unit cell volume: %-15.4f" % f_obs.unit_cell().volume()
  #
  grow_density(f_obs,
               r_free_flags,
               pdb_file_name,
               xray_structure,
               x_center=params.x_center,
               y_center=params.y_center,
               z_center=params.z_center,
               radius=params.radius,
               step_size=params.atom_gap,
               overlap_interval=params.overlap_interval,
               atom_type= params.atom_type,
               number_of_cycles = params.cycles,
               number_of_iterations = params.iterations,
               bfac_dummy_atoms = params.bfac_dummy_atoms,
               bfac_cutoff = params.bfac_cutoff,
               occ_cutoff = params.occ_cutoff)
