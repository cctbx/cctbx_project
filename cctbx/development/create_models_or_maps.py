from __future__ import absolute_import, division, print_function
import sys,os
import iotbx.phil
from scitbx.array_family import flex


'''
 Utilities to generate models and maps from standard or user-selected PDB files

 Use map_model_manager to access these utilities

'''

def get_map_manager(map_data,wrapping,
     unit_cell_dimensions=None,crystal_symmetry=None,):
  '''
      Get a minimal map_manager in SG p1 from any map_data and
      if unit_cell_dimensions.  Assume cell angles are 90,90,90
      Shift origin to (0,0,0)
  '''
  assert unit_cell_dimensions or crystal_symmetry
  assert isinstance(wrapping, bool)
  from iotbx.map_manager import map_manager
  map_data = map_data.shift_origin()

  if not crystal_symmetry:
    from cctbx import crystal
    crystal_symmetry=crystal.symmetry(
        tuple(list(unit_cell_dimensions[:3])+[90,90,90]), 1)
  mm=map_manager(map_data=map_data,
    unit_cell_grid = map_data.all(),
    unit_cell_crystal_symmetry = crystal_symmetry,
    origin_shift_grid_units = (0,0,0),
    wrapping = wrapping)
  return mm

def read_map_and_model(file_name_1,file_name_2, regression_directory = None,
    prefix = None):
  '''
    Identify which file is map and which is model, read in and
    create map_model_manager
    If regression_directory is specified, look there for these files, assuming
    prefix of $PHENIX/modules/phenix_regression/
  '''

  if regression_directory and not prefix:
    import libtbx.load_env
    prefix = libtbx.env.under_dist(
      module_name="phenix_regression",
      path=regression_directory,
      test=os.path.isdir)

  if prefix:
    file_name_1 = os.path.join(prefix,file_name_1)
    file_name_2 = os.path.join(prefix,file_name_2)

  map_file_name = None
  model_file_name = None
  for f in [file_name_1,file_name_2]:
    for ending in ['.ccp4','.mrc','.map']:
      if f.endswith(ending):
        map_file_name = f
    for ending in ['.pdb','.cif']:
      if f.endswith(ending):
        model_file_name = f
  if not map_file_name or not model_file_name:
    raise Sorry("Unable to guess map and model from %s and %s" %(
     file_name_1, file_name_2))

  from iotbx.data_manager import DataManager
  from iotbx.map_model_manager import map_model_manager
  dm = DataManager()
  dm.process_real_map_file(map_file_name)
  mm = dm.get_real_map(map_file_name)

  dm.process_model_file(model_file_name)
  model = dm.get_model(model_file_name)
  mam = map_model_manager(model=model,map_manager=mm)
  return mam

def generate_model(
      file_name=None,
      n_residues=None,
      start_res=None,
      b_iso=30,
      box_cushion=5,
      space_group_number=1,
      output_model_file_name=None,
      shake=None,
      random_seed=None,
      log=sys.stdout):

  '''
    generate_model: Simple utility for generating a model for testing purposes.

    This function typically accessed and tested through map_model_manager

     Summary
    -------

    Generate a model from a user-specified file or from some examples available
    in the cctbx.  Cut out specified number of residues, shift to place on
    positive side of origin, optionally set b values to b_iso,
    place in box with buffering of box_cushion on all
    edges, optionally randomly shift (shake) atomic positions by rms of shake A,
    and write out to output_model_file_name and return model object.

    Parameters:

      file_name (string, None):  File containing model (PDB, CIF format)
      n_residues (int, 10):      Number of residues to include
      start_res (int, None):     Starting residue number
      b_iso (float, 30):         B-value (ADP) to use for all atoms
      box_cushion (float, 5):     Buffer (A) around model
      space_group_number (int, 1):  Space group to use
      output_model_file_name (string, None):  File for output model
      shake (float, None):       RMS variation to add (A) in shake
      random_seed (int, None):    Random seed for shake

    Returns:
      model.manager object (model) in a box defined by a crystal_symmetry object

  '''

  # Get the parameters

  space_group_number=int(space_group_number)
  if n_residues is not None:
    n_residues=int(n_residues)
  box_cushion=float(box_cushion)
  if start_res:
    start_res=int(start_res)
  if shake:
    shake=float(shake)
  if random_seed:
    random_seed=int(random_seed)
    import random
    random.seed(random_seed)
    random_seed=random.randint(1,714717)
    flex.set_random_seed(random_seed)

  # Choose file with coordinates

  if not file_name:
    if not n_residues:
      n_residues = 10 # default
    import libtbx.load_env
    iotbx_regression = os.path.join(libtbx.env.find_in_repositories("iotbx"),
      'regression')
    if n_residues < 25:
      file_name=os.path.join(iotbx_regression,'secondary_structure', # PDB OK
         '5a63_chainBp.pdb')  # starts at 219
      if not start_res: start_res=219
    elif n_residues < 167:
      file_name=os.path.join(iotbx_regression,'secondary_structure', # PDB OK
         '3jd6_noh.pdb') # starts at 58
      if not start_res:start_res=58
    else:
      file_name=os.path.join(iotbx_regression,'secondary_structure', # PDB OK
         '4a7h_chainC.pdb') # starts at 9
      if not start_res:start_res=9
  else: # have file_name                              # PDB OK
    if start_res is None:
      start_res=1
    if not n_residues:
      n_residues = 100000 #  a big number

  # Read in coordinates and cut out the part of the model we want

  from iotbx.data_manager import DataManager

  dm = DataManager(['model'])
  dm.process_model_file(file_name)
  model = dm.get_model(file_name)

  if not model.crystal_symmetry() or not model.crystal_symmetry().unit_cell():
    from cctbx.maptbx.box import shift_and_box_model
    model = shift_and_box_model(model = model,
      box_cushion = box_cushion)

  selection=model.selection('resseq %s:%s' %(start_res,start_res+n_residues-1))
  model=model.select(selection)

  # shift the model and return it with new crystal_symmetry
  from cctbx.maptbx.box import shift_and_box_model
  model = shift_and_box_model(model = model,
    box_cushion = box_cushion)

  if b_iso is not None:
    b_values=flex.double(model.get_sites_cart().size(), b_iso)
    ph = model.get_hierarchy()
    ph.atoms().set_b(b_values)

  # Optionally shake model
  if shake:
    model=shake_model(model,shake=shake)

  if output_model_file_name:
    f=open(output_model_file_name,'w')
    print ("%s" %(model.model_as_pdb()),file=f)  # PDB OK
    f.close()
    print ("Writing model with %s residues and b_iso=%s from %s to %s" %(
      n_residues,b_iso,file_name,output_model_file_name),file=log)
  else:
    print ("Generated model with %s residues and b_iso=%s from %s " %(
      n_residues,b_iso,file_name),file=log)
  return model

def shake_model(model,shake=None,log=sys.stdout):
  # Shake model with rmsd of shake
  from mmtbx.pdbtools import master_params_str
  params=iotbx.phil.parse(master_params_str).extract()
  params.modify.sites[0].shake=shake
  from mmtbx.pdbtools import modify
  new_model = modify(
      model  = model,
      params = params.modify,
      log    = log).get_results().model
  return new_model

def generate_map_coefficients(
      model=None,  # Required model
      output_map_coeffs_file_name=None,
      d_min=3,
      k_sol = None,
      b_sol = None,
      scattering_table='electron',
      f_obs_array = None,
      log=sys.stdout):
  '''
    Convenience function to create map coefficients from a model.

    This function typically accessed and tested through map_model_manager

    Summary:  Supply model.manager object or parameters for generate_model,
    high_resolution limit of d_min and optional scattering_table ("electron"
    for cryo-EM, "n_gaussian" for x-ray) and optional
    output_map_coeffs_file_name.

    Writes map coefficients and returns miller_array object
    containing the map coefficients

    Parameters:
    -----------

      model (model.manager object, None):    model to use.
         contains crystal_symmetry
      output_map_coeffs_file_name (string, None): output model file name
      d_min (float, 3):   High-resolution limit for map coeffs (A)
      scattering_table (choice, 'electron'): choice of scattering table
           All choices: wk1995 it1992 n_gaussian neutron electron
      f_obs_array:  array used to match indices of fcalc
      crystal_symmetry:  used for crystal symmetry (overrides model value but
        not f_obs)

  '''

  assert model is not None

  # get map coefficients
  from mmtbx.utils import fmodel_from_xray_structure
  import iotbx.phil
  from mmtbx.programs.fmodel import master_phil
  fmodel_params=master_phil.extract()

  if f_obs_array:
    assert f_obs_array.crystal_symmetry().is_similar_symmetry(
      model.crystal_symmetry())

  assert model.crystal_symmetry() is not None # Need crystal_symmetry in model

  xrs=model.get_xray_structure()

  fmodel_params.high_resolution=float(d_min)
  fmodel_params.scattering_table=scattering_table
  if k_sol is not None and b_sol is not None:
    fmodel_params.fmodel.k_sol = k_sol
    fmodel_params.fmodel.b_sol = b_sol
  if f_obs_array:
     fmodel_params.high_resolution=f_obs_array.d_min()-0.0001 # different cut
  fmodel=fmodel_from_xray_structure(
    xray_structure = xrs,
    params         = fmodel_params,
    f_obs          = f_obs_array,
    out            = log)
  f_model=fmodel.f_model
  if output_map_coeffs_file_name:
    f_model.as_mtz_dataset(column_root_label='FWT').mtz_object().write(
      file_name=output_map_coeffs_file_name)
    print("Writing map coefficients to resolution of %s A to %s" %(
       d_min,output_map_coeffs_file_name),file=log)
  else:
    print("Generated map coefficients to resolution of %s A " %(
      d_min),file=log)
  return f_model

def generate_map(
      output_map_file_name=None,
      map_coeffs=None,  # Required
      d_min=None,
      map_manager = None, # source of info, not map
      gridding=None,
      wrapping=False,
      resolution_factor=None,
      origin_shift_grid_units=None,
      low_resolution_fourier_noise_fraction=0,
      high_resolution_fourier_noise_fraction=0,
      low_resolution_real_space_noise_fraction=0,
      high_resolution_real_space_noise_fraction=0,
      low_resolution_noise_cutoff=None,
      random_seed = None,
      log=sys.stdout):


  '''
      Generate map from map_coefficients and add noise in Fourier or real space

      This function typically accessed and tested through map_model_manager

      Summary:
      --------

      Calculate a map and optionally add noise to it.  Supply map
      coefficients (miller_array object) and types of noise to add,
      along with optional gridding (nx,ny,nz), and origin_shift_grid_units.
      Optionally create map coefficients from a model and optionally
      generate a model.

      Unique aspect of this noise generation is that it can be specified
      whether the noise is local in real space (every point in a map
      gets a random value before Fourier filtering), or local in Fourier
      space (every Fourier coefficient gets a complex random offset).
      Also the relative contribution of each type of noise vs resolution
      can be controlled.

      Parameters:
      -----------


      output_map_file_name (string, None):  Output map file (MRC/CCP4 format)
      map_coeffs (miller.array object, None) : map coefficients
      d_min(float):      high_resolution limit (A)
      gridding (tuple (nx,ny,nz), None):  Gridding of map (optional)
      origin_shift_grid_units (tuple (ix,iy,iz), None):  Move location of
          origin of resulting map to (ix,iy,iz) before writing out
      low_resolution_fourier_noise_fraction (float, 0): Low-res Fourier noise
      high_resolution_fourier_noise_fraction (float, 0): High-res Fourier noise
      low_resolution_real_space_noise_fraction(float, 0): Low-res
          real-space noise
      high_resolution_real_space_noise_fraction (float, 0): High-res
          real-space noise
      low_resolution_noise_cutoff (float, None):  Low resolution where noise
          starts to be added

  '''

  if random_seed:
    random_seed=int(random_seed)
    import random
    random.seed(random_seed)
    random_seed=random.randint(1,714717)
    flex.set_random_seed(random_seed)

  if map_manager:
    origin_shift_grid_units=map_manager.origin_shift_grid_units
    gridding=map_manager.map_data().all()
    wrapping=map_manager.wrapping()
  if gridding:
    if type(gridding) in [type((1,2,3)), type([1,2,3])] and type(gridding[0])==type(1):
      pass # already fine
    else:
      new_gridding=[]
      for x in str(gridding).replace("(","").replace(")","").replace("[","").replace("]","").replace(",","").split():
        new_gridding.append(int(x))
      gridding = new_gridding
  low_resolution_fourier_noise_fraction=float(
     low_resolution_fourier_noise_fraction)
  high_resolution_fourier_noise_fraction=float(
     high_resolution_fourier_noise_fraction)
  low_resolution_real_space_noise_fraction=float(
     low_resolution_real_space_noise_fraction)
  high_resolution_real_space_noise_fraction=float(
     high_resolution_real_space_noise_fraction)
  if low_resolution_noise_cutoff:
     low_resolution_noise_cutoff=float(low_resolution_noise_cutoff)

  if d_min:
    d_min=float(d_min)
    map_coeffs=map_coeffs.resolution_filter(d_min=d_min)

  # Calculate a map from Fourier coefficients:
  map_data=get_map_from_map_coeffs(
     map_coeffs=map_coeffs,
     crystal_symmetry=map_coeffs.crystal_symmetry(),
     n_real=gridding,
     resolution_factor=resolution_factor,
     apply_sigma_scaling=False)

  # Optionally add noise to this map as an additive noise map
  # Noise can be added in Fourier space (leads to correlated errors
  #    in real space)  or in real space (leads to correlated errors
  #    in Fourier space).
  # Noise is Fourier-weighted as function of resolution.
  # RMS noise to add at low-resolution (Fourier based noise) as fraction of RMS
  #   value in map at low-resolution is: low_resolution_fourier_noise_fraction
  # RMS Fourier high-res noise:is high_resolution_fourier_noise_fraction
  # RMS real-space low-res noise:is low_resolution_real_space_noise_fraction
  # RMS real-space high-res noise:is high_resolution_real_space_noise_fraction
  # Low-resolution where noise begins to be added is low_resolution_noise_cutoff

  if (low_resolution_fourier_noise_fraction or
      high_resolution_fourier_noise_fraction):
    fourier_noise_map=get_fourier_noise_map(n_real=map_data.all(),
     map_coeffs=map_coeffs,
     low_resolution_fourier_noise_fraction=
        low_resolution_fourier_noise_fraction,
     high_resolution_fourier_noise_fraction=
        high_resolution_fourier_noise_fraction,
     d_min=d_min,
     low_resolution_noise_cutoff=low_resolution_noise_cutoff,
     log=log)
  else:
    fourier_noise_map=None

  if (low_resolution_real_space_noise_fraction or
      high_resolution_real_space_noise_fraction):
    real_space_noise_map=get_real_space_noise_map(map_data=map_data,
     map_coeffs=map_coeffs,
     low_resolution_real_space_noise_fraction=
        low_resolution_real_space_noise_fraction,
     high_resolution_real_space_noise_fraction=
        high_resolution_real_space_noise_fraction,
     d_min=d_min,
     low_resolution_noise_cutoff=low_resolution_noise_cutoff,
     log=log)
  else:
    real_space_noise_map=None

  if fourier_noise_map:
    map_data+=fourier_noise_map
  if real_space_noise_map:
    map_data+=real_space_noise_map

  if map_manager:
    mm=map_manager.customized_copy(map_data=map_data)
  else:
    # Create a map_manager object directly (unusual use of map_manager)
    from iotbx.map_manager import map_manager
    mm=map_manager(map_data=map_data,
      unit_cell_grid=map_data.all(),
      unit_cell_crystal_symmetry=map_coeffs.crystal_symmetry(),
      origin_shift_grid_units=origin_shift_grid_units,
      wrapping=wrapping)

  if output_map_file_name:
    mm.write_map(output_map_file_name)
  else:
    print("Generated map with origin at %s and size of %s" %(
      mm.map_data().origin(),mm.map_data().all()),file=log)

  return mm

def squares_of_complex(m1):
  a1=flex.pow2(m1.parts()[0])
  a2=flex.pow2(m1.parts()[1])
  a3=a1+a2
  return a3

def norm(m1):
  a3=squares_of_complex(m1)
  return a3.min_max_mean().mean**0.5


def map_coeffs_as_fp_phi(map_coeffs):
  amplitudes=map_coeffs.amplitudes()
  amplitudes.set_observation_type_xray_amplitude()
  assert amplitudes.is_real_array()
  phases=map_coeffs.phases(deg=True)
  return amplitudes,phases

def fp_phi_as_map_coeffs(fp,phi):
  return fp.phase_transfer(phase_source=phi,deg=True)

def get_real_space_noise_map(map_data=None,
     map_coeffs=None,
     low_resolution_real_space_noise_fraction=None,
     high_resolution_real_space_noise_fraction=None,
     d_min=None,
     low_resolution_noise_cutoff=None,
     log=sys.stdout):

  '''
    Procedure for generating a map with noise in real space

    NOTE: only applies for space group P1

    Parameters:
    -----------

    map_data (map_data object, None): map_data with current values in map

    map_coeffs (miller.array object, None): map coefficients assumed to
      match map_data. Used to get RMS values in Fourier representation of
      map vs resolution and used to specify which Fourier coefficients to
      calculate
    low_resolution_real_space_noise_fraction (float, None): ratio of
      RMS value in output map to input map at low resolution
    high_resolution_real_space_noise_fraction (float, None): ratio of
      RMS value in output map to input map at high resolution
    low_resolution_noise_cutoff (float, None): resolution where noise is added
  '''

  # Get random values at grid points in map.  Then obtain Fourier
  #  coefficients by back-transformation.  Then scale Fourier coefficients to
  #  yield low_resolution_real_space_noise_fraction at low resolution and
  #    high_resolution_real_space_noise_fraction at high resolution and
  #    linear in 1/d in between.

  assert map_coeffs.crystal_symmetry().space_group_number() in [
    0,1]

  print ("\nGenerating random map in real space, then Fourier filtering",
     file=log)

  # Create a map with mean zero and rms 1 with random values
  random_values=flex.random_double(map_data.size())
  acc=map_data.accessor()
  random_values.reshape(acc)  # map with random values rms 1 mean zero

  # Get map as fourier coefficients
  from cctbx import miller
  randomized_map_coeffs=map_coeffs.structure_factors_from_map(random_values)

  return scale_map_coeffs(
    n_real=map_data.all(),
    randomized_map_coeffs=randomized_map_coeffs,
    map_coeffs=map_coeffs,
    high_resolution_noise_fraction=high_resolution_real_space_noise_fraction,
    low_resolution_noise_fraction=low_resolution_real_space_noise_fraction,
    random_selection_within_bins=False,
    low_resolution_noise_cutoff=low_resolution_noise_cutoff,
    log=log)

def get_fourier_noise_map(n_real=None,
     map_coeffs=None,
     low_resolution_fourier_noise_fraction=None,
     high_resolution_fourier_noise_fraction=None,
     d_min=None,
     low_resolution_noise_cutoff=None,
     log=sys.stdout):

  '''
    Procedure for generating a map with noise in Fourier space

    NOTE: only applies for space group P1

    Parameters:
    -----------

    n_real (tuple (nx,ny,nz), None):  Gridding of output map. Usually
      obtained from map_data.all()
    map_coeffs (miller.array object, None): map coefficients assumed to
      match n_real used to get RMS values in Fourier representation of
      map vs resolution and used to specify which Fourier coefficients to
      calculate
    low_resolution_fourier_noise_fraction (float, None): ratio of
      RMS value in output map to input map at low resolution
    high_resolution_fourier_noise_fraction (float, None): ratio of
      RMS value in output map to input map at high resolution
    low_resolution_noise_cutoff (float, None): resolution where noise is added
  '''

  # Get random values of Fourier coefficients with rms values scaled to
  #  yield low_resolution_fourier_noise_fraction at low resolution and
  #    high_resolution_fourier_noise_fraction at high resolution and
  #    linear in 1/d in between.

  assert map_coeffs.crystal_symmetry().space_group_number() in [
    0,1]

  return scale_map_coeffs(
    n_real=n_real,
    map_coeffs=map_coeffs,
    high_resolution_noise_fraction=high_resolution_fourier_noise_fraction,
    low_resolution_noise_fraction=low_resolution_fourier_noise_fraction,
    random_selection_within_bins=True,
    low_resolution_noise_cutoff=low_resolution_noise_cutoff,
    log=log)

def scale_map_coeffs(
   n_real=None,
   randomized_map_coeffs=None,
   map_coeffs=None,
   high_resolution_noise_fraction=None,
   low_resolution_noise_fraction=None,
   random_selection_within_bins=False,
   low_resolution_noise_cutoff=None,
   log=sys.stdout):

  '''
    Scale map coefficients to target value vs resolution, optionally randomize

    Scales map coefficients in resolution bins, scale factor adjusted
    to yield high_resolution_noise_fraction at high-resolution limit
    and low_resolution_noise_fraction at low_resolution_noise_cutoff and
    linearly between in 1/resolution.

    Optionally randomizes amplitudes and phases by shuffling within bins

  '''

  assert random_selection_within_bins or randomized_map_coeffs

  if not hasattr(map_coeffs,'binner') or not map_coeffs.binner():
    map_coeffs.setup_binner(auto_binning=True)

  d_max,d_min=map_coeffs.d_max_min()
  if d_max < 0: d_max = 1.e+10

  if random_selection_within_bins:
    new_map_coeffs=map_coeffs.customized_copy(
      data=flex.complex_double(map_coeffs.size(),(0+0.j)))
    print ("\nGenerating map randomized in Fourier space",file=log)
  else:
    new_map_coeffs=randomized_map_coeffs
  print ("Relative error added at high-resolution: %.3f" %(
     high_resolution_noise_fraction),file=log)
  print ("Relative error added at low-resolution: %.3f" %(
     low_resolution_noise_fraction),file=log)

  print("Resolution  Noise ratio   RMS original            RMS error ",file=log)

  for i_bin in map_coeffs.binner().range_used():
    sel=map_coeffs.binner().selection(i_bin)
    dd=map_coeffs.d_spacings().data().select(sel)
    local_d_mean     = dd.min_max_mean().mean
    local_s_mean=1/local_d_mean
    s_max=1/max(1.e-10,d_min)
    if low_resolution_noise_cutoff:
      s_min=1/max(1.e-10,min(d_max,low_resolution_noise_cutoff))
    else:
      s_min=1/max(1.e-10,d_max)
    fraction_high= max(0.,min(1.,
       (local_s_mean-s_min)/max(1.e-10,s_max-s_min)))

    noise_ratio=low_resolution_noise_fraction+\
      fraction_high * (
       high_resolution_noise_fraction-
        low_resolution_noise_fraction)

    mc=map_coeffs.select(sel)
    rms_original=norm(mc.data())
    if random_selection_within_bins: # randomize, normalize, scale
      fp,phi=map_coeffs_as_fp_phi(mc)
      sel_fp=flex.random_permutation(fp.size())
      new_fp=fp.select(sel_fp)
      data=new_fp.data()
      data*=noise_ratio
      sel_phi=flex.random_permutation(phi.size())
      new_phi=phi.select(sel_phi)
      new_mc=fp_phi_as_map_coeffs(new_fp,new_phi)
    else:  # just normalize and scale
      randomized_mc=randomized_map_coeffs.select(sel)
      rms_new=norm(randomized_mc.data())
      scale=rms_original/max(1.e-10,rms_new)
      new_fp,new_phi=map_coeffs_as_fp_phi(randomized_mc)
      data=new_fp.data()
      data*=noise_ratio*scale
      new_mc=fp_phi_as_map_coeffs(new_fp,new_phi)

    rms_new=norm(new_mc.data())
    new_map_coeffs.data().set_selected(sel,new_mc.data())
    print ("  %.3f         %.3f         %.3f               %.3f " %(
      local_d_mean, noise_ratio,rms_original,rms_new),file=log)

  # Convert to map
  from cctbx import maptbx
  return maptbx.map_coefficients_to_map(
      map_coeffs       = new_map_coeffs,
      crystal_symmetry = new_map_coeffs.crystal_symmetry(),
      n_real           = n_real)

def get_map_from_map_coeffs(map_coeffs = None, crystal_symmetry = None,
     n_real = None, resolution_factor = None, apply_sigma_scaling = True):
    if resolution_factor is None:
      resolution_factor = 0.25
    from cctbx import maptbx
    from cctbx.maptbx import crystal_gridding
    if not crystal_symmetry:
      crystal_symmetry = map_coeffs.crystal_symmetry()
    if map_coeffs.crystal_symmetry().space_group_info()!=  \
       crystal_symmetry.space_group_info():
      assert str(map_coeffs.crystal_symmetry().space_group_info()
          ).replace(" ", "").lower() == 'p1'
      # use map_coeffs.crystal_symmetry
      crystal_symmetry = map_coeffs.crystal_symmetry()
    if n_real:
      cg = crystal_gridding(
        unit_cell = crystal_symmetry.unit_cell(),
        space_group_info = crystal_symmetry.space_group_info(),
        pre_determined_n_real = n_real)
    else:
      cg = None
    fft_map = map_coeffs.fft_map( resolution_factor = resolution_factor,
       crystal_gridding = cg,
       symmetry_flags = maptbx.use_space_group_symmetry)
    if apply_sigma_scaling:
      fft_map.apply_sigma_scaling()
    else:
      fft_map.apply_volume_scaling()
    map_data = fft_map.real_map_unpadded()
    return map_data
