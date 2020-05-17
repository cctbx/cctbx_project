from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry
import sys,os
from libtbx import group_args
import iotbx.phil


'''
 Utility to generate models and maps from standard or user-selected PDB files
'''

def generate_model(
      file_name=None,
      n_residues=10,
      b_iso=30,
      box_buffer=5,
      start_res=None,
      space_group_number=1,
      output_model_file_name=None,
      random_seed=None,
      shake=None,
      log=sys.stdout,**ignored_kw):

  '''
    generate_model: Simple utility for generating a model for testing purposes.

    Generate a model from a user-specified file or from some examples available
    in the cctbx.  Cut out specified number of residues, shift to place on
    positive side of origin, optionally set b values to b_iso,
    place in box with buffering of box_buffer on all
    edges, optionally randomly shift (shake) atomic positions by rms of shake A,
    and write out to output_model_file_name and return model object.
  '''

  # Get the parameters

  space_group_number=int(space_group_number)
  n_residues=int(n_residues)
  box_buffer=float(box_buffer)
  if start_res:
    start_res=int(start_res)
  if shake:
    shake=float(shake)
  if random_seed:
    random_seed=int(random_seed)
    import random
    random.seed(random_seed)
    random_seed=random.randint(1,714717)
    from scitbx.array_family import flex
    flex.set_random_seed(random_seed)

  # Choose file with coordinates

  if not file_name:
    import libtbx.load_env
    iotbx_regression = os.path.join(libtbx.env.find_in_repositories("iotbx"),
      'regression')
    if n_residues < 25:
      file_name=os.path.join(iotbx_regression,'secondary_structure',
         '5a63_chainBp.pdb')  # starts at 219
      if not start_res: start_res=219
    elif n_residues < 167:
      file_name=os.path.join(iotbx_regression,'secondary_structure',
         '3jd6_noh.pdb') # starts at 58
      if not start_res:start_res=58
    else:
      file_name=os.path.join(iotbx_regression,'secondary_structure',
         '4a7h_chainC.pdb') # starts at 9
      if not start_res:start_res=9

  # Read in coordinates and cut out the part of the model we want

  from mmtbx.model import manager as model_manager
  model=model_manager(file_name)
  selection=model.selection('resseq %s:%s' %(start_res,start_res+n_residues-1))
  model=model.select(selection)

  # shift the model and return it with new crystal_symmetry

  import iotbx.pdb
  from scitbx.matrix import col
  from cctbx import crystal
  ph=model.get_hierarchy()
  sites_cart=ph.atoms().extract_xyz()
  sites_cart=sites_cart-col(sites_cart.min())+col(
      (box_buffer,box_buffer,box_buffer))
  box_end=col(sites_cart.max())+col((box_buffer,box_buffer,box_buffer))
  a,b,c=box_end
  crystal_symmetry=crystal.symmetry((a,b,c, 90,90,90),1)
  ph.atoms().set_xyz(sites_cart)

  if b_iso is not None:
    from scitbx.array_family import flex
    b_values=flex.double(sites_cart.size(), b_iso)
    ph.atoms().set_b(b_values)
  model=model_manager(
     ph.as_pdb_input(),
     crystal_symmetry = crystal_symmetry,
     log = log)

  # Optionally shake model
  if shake:
    model=shake_model(model,shake=shake)

  if output_model_file_name:
    f=open(output_model_file_name,'w')
    print ("%s" %(model.model_as_pdb()),file=f)
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

def generate_map_coefficients(model=None,output_map_coeffs_file_name=None,
        high_resolution=3,scattering_table='electron',
       log=sys.stdout,**ignored_kw):

  '''
    generate_map_coefficients  Convenience function to create
    map coefficients from a model.  Supply model.manager object,
    high_resolution limit and optional scattering_table ("electron"
    for cryo-EM, "n_gaussian" for x-ray) and optional
    output_map_coeffs_file_name.

    Writes map coefficients and returns miller_array object
    containing the map
  '''

  if not model:
    model=generate_model(log=log,**ignored_kw)

  # get map coefficients
  from mmtbx.utils import fmodel_from_xray_structure
  import iotbx.phil
  from mmtbx.command_line.fmodel import master_phil
  fmodel_params=master_phil.extract()

  xrs=model.get_xray_structure()

  fmodel_params.high_resolution=float(high_resolution)
  fmodel_params.scattering_table=scattering_table

  fmodel=fmodel_from_xray_structure(
    xray_structure = xrs,
    params         = fmodel_params,
    out            = log)
  f_model=fmodel.f_model
  if output_map_coeffs_file_name:
    f_model.as_mtz_dataset(column_root_label='FWT').mtz_object().write(
      file_name=output_map_coeffs_file_name)
    print("Writing map coefficients to resolution of %s A to %s" %(
       high_resolution,output_map_coeffs_file_name),file=log)
  else:
    print("Generated map coefficients to resolution of %s A " %(
      high_resolution),file=log)
  return f_model

def generate_map(map_coeffs=None,
   high_resolution=3,
   gridding=None,
   low_resolution_fourier_noise_fraction=0,
   high_resolution_fourier_noise_fraction=0,
   low_resolution_real_space_noise_fraction=0,
   high_resolution_real_space_noise_fraction=0,
   output_map_file_name=None,
   log=sys.stdout, **ignored_kw):

  '''
    generate_map

    Convenience method to calculate a map and
    optionally add noise to it.  Supply map coefficients (miller_array
    object) and types of noise to add, along with optional gridding (nx,ny,nz).

    Not implemented:
    Unique aspect of this noise generation is that it can be specified
    whether the noise is local in real space (every point in a map
    gets a random value before Fourier filtering), or local in Fourier
    space (every Fourier coefficient gets a complex random offset).
    Also the relative contribution of each type of noise vs resolution
    can be controlled.

  '''

  if gridding:
    if type(gridding)==type([1,2,3]) and type(gridding[0])==type(1):
      pass # already fine
    else:
      gridding=[]
      for x in str(gridding).replace(",","").split():
        gridding.append(int(x))


  if not map_coeffs:  # get map coefficients
    map_coeffs=generate_map_coefficients(
     high_resolution=high_resolution,log=log,**ignored_kw)

  if high_resolution:
    map_coeffs=map_coeffs.resolution_filter(d_min=high_resolution)

  # Calculate a map from Fourier coefficients:
  from cctbx.maptbx.segment_and_split_map import get_map_from_map_coeffs
  map_data=get_map_from_map_coeffs(
     map_coeffs=map_coeffs,
     crystal_symmetry=map_coeffs.crystal_symmetry(),
     n_real=gridding,
     apply_sigma_scaling=False)

  from iotbx.map_manager import map_manager
  mm=map_manager(map_data=map_data,
    unit_cell_grid=map_data.all(),
    unit_cell_parameters=map_coeffs.crystal_symmetry().unit_cell().parameters(),
    space_group_number=map_coeffs.crystal_symmetry().space_group().info(
         ).symbol_and_number().split('(')[0],
    origin_shift_grid_units=(0,0,0),
    working_map_n_xyz=map_data.all())

  if output_map_file_name:
    mm.write_map(output_map_file_name)
  else:
    print("Generated map with origin at %s and size of %s" %(
      mm.map_data().origin(),mm.map_data().all()),file=log)

  return mm


def run(**kw):

  if kw.get('generate_model','None').lower()=='true':
    model=generate_model(**kw)

  elif kw.get('generate_map_coefficients','None').lower()=='true':
    map_coeffs=generate_map_coefficients(**kw)

  elif kw.get('generate_map','None').lower()=='true':
    mm=generate_map(**kw)


if __name__=="__main__":
  kw={}
  args=[]
  for x in sys.argv[1:]:
    args.append(x.lower())
  if args.count(['generate_model=true','generate_map_coefficients=true','generate_map=true'])==0:
    extra=['generate_map=true']
  else:
    extra=[]

  if  kw.get('output_model_file_name','None').lower()=='none':
      kw['output_model_file_name']='dummy_model.pdb'
  if  kw.get('output_map_coeffs_file_name','None').lower()=='none':
      kw['output_map_coeffs_file_name']='dummy_model.mtz'
  if  kw.get('output_map_file_name','None').lower()=='none':
      kw['output_map_file_name']='dummy_model.ccp4'

  for x in args+extra:
    spl=x.split("=")
    if len(spl)==2:
      kw[spl[0]]=spl[1]
  run(**kw)
