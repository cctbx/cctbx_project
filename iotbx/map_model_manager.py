from __future__ import absolute_import, division, print_function
from libtbx.utils import null_out
import sys
from iotbx.map_manager import map_manager
from libtbx.test_utils import approx_equal

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

     Shift origin to (0, 0, 0) and save position of this (0, 0, 0) point in the
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

     Input models, maps, and ncs_object must all match in crystal_symmetry,
     original (unit_cell) crystal_symmetry, and shift_cart (-origin_shift_cart
     for maps)

  '''

  def __init__(self, ):
    self._map_manager = None
    self._model = None
    self._ncs_object = None

  def deep_copy(self):
    from copy import deepcopy
    new_mmm = map_model_manager()
    if self._model:
      new_mmm.add_model(self._model.deep_copy())
    if self._map_manager:
      new_mmm.add_map_manager(self._map_manager.deep_copy())
    if self._ncs_object:
      new_mmm.add_model(self._model.deep_copy())
    return new_mmm

  def show_summary(self, log = sys.stdout):
    print ("Summary of maps and models", file = log)
    if self._map_manager:
      print("Map summary:", file = log)
      self._map_manager.show_summary(out = log)
    if self._model:
      print("Model summary:", file = log)
      print("Residues: %s" %(
       self._model.get_hierarchy().overall_counts().n_residues), file = log)

    if self._ncs_object:
      print("NCS summary:", file = log)
      print("Operators: %s" %(
       self._ncs_object.max_operators()), file = log)

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

  def write_map(self, file_name = None, log = sys.stdout):
    if not self._map_manager:
      print ("No map to write out", file = log)
    elif not file_name:
      print ("Need file name to write map", file = log)
    else:
      self._map_manager.write_map(file_name = file_name)

  def write_model(self,
     file_name = None,
     log = sys.stdout):
    if not self._model:
      print ("No model to write out", file = log)
    elif not file_name:
      print ("Need file name to write model", file = log)
    else:
      # See if we need to shift model
      if self._map_manager and \
           self._map_manager.origin_shift_grid_units !=  (0, 0, 0):  # shift it
        model = self._model.deep_copy()
        print("Shifting model to match original map", file = log)
        model = self.shift_model_to_match_original_map(
          model = model, log = log)
      else:
        model = self._model

      f = open(file_name, 'w')
      print(model.model_as_pdb(), file = f)
      f.close()
      print("Wrote model with %s residues to %s" %(
         model.get_hierarchy().overall_counts().n_residues,
         file_name), file = log)

  def add_map_manager(self, map_manager = None, log = sys.stdout):
    # Add a map and make sure its symmetry is similar to others
    self._map_manager = map_manager
    if self.model():
      self.check_model_and_set_to_match_map()

  def check_model_and_set_to_match_map(self):
    # Map, model and ncs_object all must have same symmetry and shifts at end

    if self.map_manager() and self.model():
      # Must be compatible...then set model symmetry if not set
      ok=self.map_manager().is_compatible_model(self.model(),
        require_similar=False)
      if ok:
        self.map_manager().set_model_symmetries_and_shift_cart_to_match_map(
          self.model())  # modifies self.model() in place
      else:
         raise AssertionError(self.map_manager().warning_message())

    if self.map_manager() and self.ncs_object():
      # Must be similar...
      if not self.map_manager().is_similar_ncs_object(self.ncs_object()):
        raise AssertionError(self.map_manager().warning_message())

  def add_model(self, model = None, set_model_log_to_null = True,
     log = sys.stdout):
    # Add a model and make sure its symmetry is similar to others
    # Check that model original crystal_symmetry matches full
    #    crystal_symmetry of map
    if set_model_log_to_null:
      model.set_log(null_out())
    self._model = model
    if self.map_manager():
      self.check_model_and_set_to_match_map()

  def add_ncs_object(self, ncs_object = None, log = sys.stdout):
    # Add an NCS object
    self._ncs_object = ncs_object
    # Check to make sure its shift_cart matches

  def read_map(self, file_name = None, log = sys.stdout):
    # Read in a map and make sure its symmetry is similar to others
    mm = map_manager(file_name)
    self.add_map_manager(mm, log = log)

  def read_model(self, file_name = None, log = sys.stdout):
    print("Reading model from %s " %(file_name), file = log)
    from iotbx.pdb import input
    inp = input(file_name = file_name)
    from mmtbx.model import manager as model_manager
    model = model_manager(model_input = inp)
    self.add_model(model, log = log)


  def read_ncs_file(self, file_name = None, log = sys.stdout):
    # Read in an NCS file and make sure its symmetry is similar to others
    from mmtbx.ncs.ncs import ncs
    ncs_object = ncs()
    ncs_object.read_ncs(file_name = file_name, log = log)
    if ncs_object.max_operators()<2:
       self.ncs_object.set_unit_ncs()
    self.add_ncs_object(ncs_object)

  def set_original_origin_and_gridding(self,
      original_origin = None,
      gridding = None):
    '''
     Use map_manager to reset (redefine) the original origin and gridding
     of the map.
     You can supply just the original origin in grid units, or just the
     gridding of the full unit_cell map, or both.

     Update shift_cart for model and ncs object if present.

    '''

    assert self._map_manager is not None

    self._map_manager.set_original_origin_and_gridding(
         original_origin = original_origin,
         gridding = gridding)

    # Get the current origin shift based on this new original origin
    shift_cart = tuple([-x for x in self._map_manager.origin_shift_cart()])
    if self._model:
      if self._model.shift_cart() is None:
        self._model.set_unit_cell_crystal_symmetry_and_shift_cart(
          unit_cell_crystal_symmetry = \
           self._map_manager.unit_cell_crystal_symmetry())
      self._model.set_shift_cart(shift_cart)
    if self._ncs_object:
      self._ncs_object.set_shift_cart(shift_cart)

  def shift_origin(self, desired_origin = (0, 0, 0), log = sys.stdout):
    # shift the origin of all maps/models to desired_origin (usually (0, 0, 0))
    if not self._map_manager:
      print ("No information about origin available", file = log)
      return
    if self._map_manager.map_data().origin() == desired_origin:
      print("Origin is already at %s, no shifts will be applied" %(
       str(desired_origin)), file = log)

    # Figure out shift of model if incoming map and model already had a shift

    if self._ncs_object or self._model:

      # Figure out shift for model and make sure model and map agree
      shift_info = self._map_manager.get_shift_info(
         desired_origin = desired_origin)

      current_origin_shift_cart = self._map_manager.grid_units_to_cart(
       shift_info.current_origin_shift_grid_units)
      expected_model_shift_cart = tuple([-x for x in current_origin_shift_cart])

      shift_to_apply_cart = self._map_manager.grid_units_to_cart(
        shift_info.shift_to_apply)
      new_origin_shift_cart = self._map_manager.grid_units_to_cart(
        shift_info.new_origin_shift_grid_units)
      new_shift_cart = tuple([-x for x in new_origin_shift_cart])
      # shift_to_apply_cart is coordinate shift we are going to apply
      #  new_origin_shift_cart is how to get back to original from new location
      #   current_origin_shift_cart is how to get orig from current location
      assert approx_equal(shift_to_apply_cart, [-(a-b) for a, b in zip(
        new_origin_shift_cart, current_origin_shift_cart)])

      # Get shifts already applied to  model and ncs_object
      #    and check that they match map

      if self._model:
        existing_shift_cart = self._model.shift_cart()
        if existing_shift_cart is not None:
          assert approx_equal(existing_shift_cart, expected_model_shift_cart)

      if self._ncs_object:
        ncs_shift_cart = self._ncs_object.shift_cart()
        assert approx_equal(ncs_shift_cart, expected_model_shift_cart)

      if self._map_manager.origin_is_zero() and \
         expected_model_shift_cart == (0, 0, 0):
        pass # Need to set model shift_cart below

    # Apply shift to model, map and ncs object

    # Shift origin of map_manager
    self._map_manager.shift_origin(desired_origin = desired_origin)

    # Shift origin of model  Note this sets model shift_cart
    if self._model:
      self._model = self.shift_model_to_match_working_map(
        coordinate_shift = shift_to_apply_cart,
        new_shift_cart = new_shift_cart,
        model = self._model, log = log)
    if self._ncs_object:
      self._ncs_object = self.shift_ncs_to_match_working_map(
        coordinate_shift = shift_to_apply_cart,
        new_shift_cart = new_shift_cart,
        ncs_object = self._ncs_object, log = log)

  def shift_ncs_to_match_working_map(self, ncs_object = None, reverse = False,
    coordinate_shift = None,
    new_shift_cart = None,
    log = sys.stdout):
    # Shift an ncs object to match the working map (based
    #    on self._map_manager.origin_shift_grid_units)
    if coordinate_shift is None:
      coordinate_shift = self.get_coordinate_shift(reverse = reverse)

    ncs_object = ncs_object.coordinate_offset(coordinate_shift)
    return ncs_object

  def shift_ncs_to_match_original_map(self, ncs_object = None, log = sys.stdout):
    return self.shift_ncs_to_match_working_map(ncs_object = ncs_object,
      reverse = True, log = log)

  def get_coordinate_shift(self, reverse = False):
    if reverse: # Get origin shift in grid units  ==  position of original origin
                #  on the current grid
      origin_shift = self._map_manager.origin_shift_grid_units
    else:  # go backwards
      a = self._map_manager.origin_shift_grid_units
      origin_shift = [-a[0], -a[1], -a[2]]

    coordinate_shift = []
    for shift_grid_units, spacing in zip(
       origin_shift, self._map_manager.pixel_sizes()):
      coordinate_shift.append(shift_grid_units*spacing)
    return coordinate_shift

  def shift_model_to_match_working_map(self, model = None, reverse = False,
     coordinate_shift = None,
     new_shift_cart = None,
     log = sys.stdout):

    '''
    Shift a model based on the coordinate shift for the working map.

    Optionally specify the shift to apply (coordinate shift) and the
    new value of the shift recorded in the model (new_shift_cart)
    '''

    if coordinate_shift is None:
      coordinate_shift = self.get_coordinate_shift(
       reverse = reverse)
    if new_shift_cart is None:
      new_shift_cart = coordinate_shift

    model.shift_model_and_set_crystal_symmetry(shift_cart = coordinate_shift,
      crystal_symmetry = model.crystal_symmetry())  # keep crystal_symmetry

    # Allow specifying the final shift_cart:
    if tuple(new_shift_cart) !=  tuple(coordinate_shift):
      model.set_shift_cart(new_shift_cart)

    return model

  def shift_model_to_match_original_map(self, model = None, log = sys.stdout):
    # Shift a model object to match the original map (based
    #    on -self._map_manager.origin_shift_grid_units)
    return self.shift_model_to_match_working_map(model = model, reverse = True,
      log = log)

  def generate_map(self,
      output_map_file_name = None,
      map_coeffs = None,
      high_resolution = 3,
      gridding = None,
      origin_shift_grid_units = None,
      low_resolution_fourier_noise_fraction = 0,
      high_resolution_fourier_noise_fraction = 0,
      low_resolution_real_space_noise_fraction = 0,
      high_resolution_real_space_noise_fraction = 0,
      low_resolution_noise_cutoff = None,
      model = None,
      output_map_coeffs_file_name = None,
      scattering_table = 'electron',
      file_name = None,
      n_residues = 10,
      start_res = None,
      b_iso = 30,
      box_buffer = 5,
      space_group_number = 1,
      output_model_file_name = None,
      shake = None,
      random_seed = None,
      log = sys.stdout):

    '''
      Generate a map using generate_model and generate_map_coefficients

      Summary:
      --------

      Calculate a map and optionally add noise to it.  Supply map
      coefficients (miller_array object) and types of noise to add,
      along with optional gridding (nx, ny, nz), and origin_shift_grid_units.
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

      Used in generate_map:
      -----------------------

      output_map_file_name (string, None):  Output map file (MRC/CCP4 format)
      map_coeffs (miller.array object, None) : map coefficients
      high_resolution (float, 3):      high_resolution limit (A)
      gridding (tuple (nx, ny, nz), None):  Gridding of map (optional)
      origin_shift_grid_units (tuple (ix, iy, iz), None):  Move location of
          origin of resulting map to (ix, iy, iz) before writing out
      low_resolution_fourier_noise_fraction (float, 0): Low-res Fourier noise
      high_resolution_fourier_noise_fraction (float, 0): High-res Fourier noise
      low_resolution_real_space_noise_fraction(float, 0): Low-res
          real-space noise
      high_resolution_real_space_noise_fraction (float, 0): High-res
          real-space noise
      low_resolution_noise_cutoff (float, None):  Low resolution where noise
          starts to be added


      Pass-through to generate_map_coefficients (if map_coeffs is None):
      -----------------------
      model (model.manager object, None):    model to use
      output_map_coeffs_file_name (string, None): output model file name
      high_resolution (float, 3):   High-resolution limit for map coeffs (A)
      scattering_table (choice, 'electron'): choice of scattering table
           All choices: wk1995 it1992 n_gaussian neutron electron

      Pass-through to generate_model (used if map_coeffs and model are None):
      -------------------------------

      file_name (string, None):  File containing model (PDB, CIF format)
      n_residues (int, 10):      Number of residues to include
      start_res (int, None):     Starting residue number
      b_iso (float, 30):         B-value (ADP) to use for all atoms
      box_buffer (float, 5):     Buffer (A) around model
      space_group_number (int, 1):  Space group to use
      output_model_file_name (string, None):  File for output model
      shake (float, None):       RMS variation to add (A) in shake
      random_seed (int, None):    Random seed for shake

    '''


    print("\nGenerating new map data\n", file = log)
    if self._map_manager:
      print("NOTE: replacing existing map data\n", file = log)
    if self._model and  file_name:
      print("NOTE: using existing model to generate map data\n", file = log)
      model = self._model
    else:
      model = None

    from iotbx.create_models_or_maps import generate_model, \
       generate_map_coefficients
    from iotbx.create_models_or_maps import generate_map as generate_map_data

    if not model and not map_coeffs:
      model = generate_model(
        file_name = file_name,
        n_residues = n_residues,
        start_res = start_res,
        b_iso = b_iso,
        box_buffer = box_buffer,
        space_group_number = space_group_number,
        output_model_file_name = output_model_file_name,
        shake = shake,
        random_seed = random_seed,
        log = log)

    if not map_coeffs:
      map_coeffs = generate_map_coefficients(model = model,
        high_resolution = high_resolution,
        output_map_coeffs_file_name = output_map_coeffs_file_name,
        scattering_table = scattering_table,
        log = log)

    mm = generate_map_data(
      output_map_file_name = output_map_file_name,
      map_coeffs = map_coeffs,
      high_resolution = high_resolution,
      gridding = gridding,
      origin_shift_grid_units = origin_shift_grid_units,
      low_resolution_fourier_noise_fraction = \
        low_resolution_fourier_noise_fraction,
      high_resolution_fourier_noise_fraction = \
        high_resolution_fourier_noise_fraction,
      low_resolution_real_space_noise_fraction = \
        low_resolution_real_space_noise_fraction,
      high_resolution_real_space_noise_fraction = \
        high_resolution_real_space_noise_fraction,
      low_resolution_noise_cutoff = low_resolution_noise_cutoff,
      log = log)

    mm.show_summary()
    self._map_manager = mm
    self._model = model

  def as_map_and_model(self):

    '''
      Return map_and_model object with contents of this class (not a deepcopy)

    '''
    import iotbx.map_and_model
    mam = iotbx.map_and_model.input(
        map_manager = self.map_manager(),
        model = self.model(),
        ncs_object = self.ncs_object(),
        )
    # Keep track of the gridding and solvent_content (if used) in this boxing.
    return mam
