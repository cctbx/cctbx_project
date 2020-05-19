from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry
import sys
from iotbx.map_manager import map_manager

class map_model_manager:

  '''
   map_model_manager

   Use: Container to hold maps, models, ncs objects, with high-level tools to
   generate them, manipulate them, and keep track of origin shifts and
   cell dimensions of boxed maps.

   Main tools:
     Read and write maps and models and NCS objects in their original
       coordinate systems
     Box maps and models

   See notes in iotbx.map_manager for information about MRC/CCP4 maps and
     how they represent part or all of a "unit cell".


   Normal usage:

     Initialize empty, then read in or add a group of model.manager,
     map_manager, and ncs objects

     Optional: read in the models, maps, ncs objects

     Optional: box the maps, models, ncs objects and save boxed versions

     Shift origin to (0,0,0) and save position of this (0,0,0) point in the
        original coordinate system so that everything can be written out
        superimposed on the original locations. This is origin_shift_grid_units
        in grid units

     Can add modified maps/models later if they have the same value of
        origin_shift_grid_units and crystal_symmetry() matches existing
        crystal_symmetry() or unit_cell_crystal_symmetry()
        (dimensions of box of data that is present or dimensions of
        full unit cell)

     Implemented:  Only one map, model at present

     NOTE: modifies the model, map_manager, and ncs objects. Call with
     deep_copy() of these if originals need to be preserved.


  '''

  def __init__(self,):
    self._map_manager=None
    self._model=None
    self._ncs_object=None

  def show_summary(self,log=sys.stdout):
    print ("Summary of maps and models",file=log)
    if self._map_manager:
      print("Map summary:",file=log)
      self._map_manager.show_summary(out=log)
    if self._model:
      print("Model summary:",file=log)
      print("Residues: %s" %(
       self._model.get_hierarchy().overall_counts().n_residues),file=log)

  def generate_map(self,
     log=sys.stdout,
     **kw):  # Accepts keywords for generate_model, generate_map_coefficients,
             #  and generate_map and all are passed to generate_model etc.


    '''
      Simple interface to iotbx.create_models_or_maps to generate a map
      Simple use:
       mm=map_manager.generate_map()  # generates a small map (and model)
       map is in mm.map_data()
       model is in mm.model()
      Advanced use...see keywords in iotbx.create_models_or_maps.generate_map
    '''

    print("\nGenerating new map data\n",file=log)
    if self._map_manager:
      print("NOTE: replacing existing map data\n",file=log)
    if self._model and  not kw.get('file_name','None').lower()=='none':
      print("NOTE: using existing model to generate map data\n",file=log)
      model=self._model
    else:
      model=None

    from iotbx.create_models_or_maps import generate_model,\
       generate_map_coefficients,generate_map

    if not model:
      model=generate_model(**kw)
    map_coeffs=generate_map_coefficients(model=model,**kw)
    mm=generate_map(map_coeffs=map_coeffs,**kw)
    mm.show_summary()
    self._map_manager=mm
    self._model=model

  def expected_crystal_symmetries(self):
    expected=[]
    if self.crystal_symmetry():
      expected.append(self.crystal_symmetry())
    if self.unit_cell_crystal_symmetry():
      expected.append(self.unit_cell_crystal_symmetry())
    return expected

  def crystal_symmetry(self):
    # Return crystal symmetry of first map, or if not present, of first model
    if self._map_manager:
      return self._map_manager.crystal_symmetry()
    elif self._model:
      return self._model.crystal_symmetry()
    else:
      return None

  def unit_cell_crystal_symmetry(self):
    # Return unit_cell crystal symmetry of first map
    if self._map_manager:
      return self._map_manager.unit_cell_crystal_symmetry()
    else:
      return None

  def map_manager(self):
    return self._map_manager

  def model(self):
    return self._model

  def ncs_object(self):
    return self._ncs_object

  def write_map(self,file_name=None,log=sys.stdout):
    if not self._map_manager:
      print ("No map to write out",file=log)
    elif not file_name:
      print ("Need file name to write map",file=log)
    else:
      self._map_manager.write_map(file_name=file_name,log=log)

  def write_model(self,
     file_name=None,
     log=sys.stdout):
    if not self._model:
      print ("No model to write out",file=log)
    elif not file_name:
      print ("Need file name to write model",file=log)
    else:
      # See if we need to shift model
      if self._map_manager and \
           self._map_manager.origin_shift_grid_units != (0,0,0):  # shift it
        model=self._model.deep_copy()
        print("Shifting model to match original map",file=log)
        model=self.shift_model_to_match_original_map(
          model=model)
      else:
        model=self._model

      f=open(file_name,'w')
      print(model.model_as_pdb(),file=f)
      f.close()
      print("Wrote model with %s residues to %s" %(
         model.get_hierarchy().overall_counts().n_residues,
         file_name),file=log)

  def add_map(self,map_manager=None,log=sys.stdout):
    # Add a map and make sure its symmetry is similar to others
    self._map_manager=map_manager
    self.check_crystal_symmetry(map_manager.crystal_symmetry(),log=log)
    # if existing:  check that map_manager.is_similar(other_map_manager)

  def add_model(self,model=None,log=sys.stdout):
    # Add a model and make sure its symmetry is similar to others
    # Check that model crystal_symmetry matches either full or working
    # crystal_symmetry of map
    self._model=model
    self.check_crystal_symmetry(model.crystal_symmetry(),log=log)

  def check_crystal_symmetry(self,crystal_symmetry,
     text_on_failure='Symmetry of model',log=sys.stdout):
    ok=False
    for cs in self.expected_crystal_symmetries():
      if crystal_symmetry.is_similar_symmetry(cs):
        ok=True
        break
    if 1: # ZZZnot ok:
      raise Sorry("Crystal symmetry of %s" %(text_on_failure)+
      "\n%s\n" %(crystal_symmetry) +
      "does not match " %(
       "that of the map that is present: "+
        " \n%s\n or that of the full map:\n%s\n" %(
       self.crystal_symmetry(),self.unit_cell_crystal_symmetry())))

  def add_ncs_object(self,ncs_object=None,log=sys.stdout):
    # Add an NCS object
    self._ncs_object=ncs_object

  def read_map(self,file_name=None,log=sys.stdout):
    # Read in a map and make sure its symmetry is similar to others
    mm=map_manager(file_name)
    self.add_map(mm,log=log)

  def read_model(self,file_name=None,log=sys.stdout):
    print("Reading model from %s " %(file_name),file=log)
    from iotbx.pdb import input
    inp=input(file_name=file_name)
    from mmtbx.model import manager as model_manager
    model = model_manager(model_input=inp)
    self.add_model(model,log=log)


  def read_ncs_object(self,file_name=None,log=sys.stdout):
    # Read in an NCS object and make sure its symmetry is similar to others
    from mmtbx.ncs.ncs import ncs
    ncs_object=ncs()
    ncs_object.read_ncs(file_name=file_name, log=log)
    if ncs_object.max_operators()<2:
       self.ncs_object.set_unit_ncs()
    self.add_ncs_object(ncs_object)

  def shift_origin(self,update_shift=None,log=sys.stdout):
    # shift the origin of all maps/models to (0,0,0)
    self._map_manager.shift_origin()
    self._model=self.shift_model_to_match_working_map(
       model=self._model,log=log)
    self._ncs_object=self.shift_ncs_to_match_working_map(
       ncs_object=self._ncs_object, log=log)

  def shift_ncs_to_match_working_map(self,ncs_object=None,reverse=False,
    log=sys.stdout):
    # Shift an ncs object to match the working map (based
    #    on self._map_manager.origin_shift_grid_units)
    coordinate_shift=self.get_coordinate_shift(reverse=reverse)
    ncs_object=ncs_object.coordinate_shift(coordinate_shift) # ZZZ
    return ncs_object

  def shift_ncs_to_match_original_map(self,ncs_object=None,log=sys.stdout):
    return shift_ncs_to_match_working_map(ncs_object=ncs_object,
      reverse=True,log=log)

  def get_coordinate_shift(self,reverse=False):
    if not reverse: # usual
      origin_shift=self._map_manager.origin_shift_grid_units
    else:  # go backwards
      a=self._map_manager.origin_shift_grid_units
      origin_shift=[-a[0],-a[1],-a[2]]

    coordinate_shift=[]
    for shift_grid_units,spacing in zip(origin_shift,self._map_manager.pixel_sizes()):
      coordinate_shift.append(shift_grid_units*spacing)
    return coordinate_shift

  def shift_model_to_match_working_map(self,model=None,reverse=False,
     log=sys.stdout):

    # Shift a model object to match the working map (based
    #    on self._map_manager.origin_shift_grid_units)
    coordinate_shift=self.get_coordinate_shift(
       reverse=reverse)
    if(coordinate_shift !=(0,0,0)):
      print ("\nMap origin is not at (0,0,0): shifting the model by %s" %(
         str(coordinate_shift)), file=log)
      sites_cart=model.get_sites_cart()
      sites_cart+=coordinate_shift
      model.set_sites_cart(sites_cart)
      model._process_input_model()
    return model

  def shift_model_to_match_original_map(self,model=None):
    # Shift a model object to match the original map (based
    #    on -self._map_manager.origin_shift_grid_units)
    return self.shift_model_to_match_working_map(model=model,reverse=True)


  def map_as_fourier_coefficients(self,
    high_resolution=None,
    log=sys.stdout):

    '''
       Convert a map to Fourier coefficients to a resolution of high_resolution.
       Generates error if the resolution is not possible given the gridding
       of the map.
    '''


