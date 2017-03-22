#!/usr/bin/env python
#
# nexus.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.


from __future__ import absolute_import, division
try:
  from dxtbx_format_nexus_ext import *
except ImportError:
  pass

class check_dtype(object):
  '''
  A class to check whether the dataset data type matches the expected

  '''

  def __init__(self, dtype):
    self.dtype = dtype

  def __call__(self, dset):
    dtype = dset.dtype
    if not dtype in self.dtype:
      return False, "%s is type %s, expected %s" % (
        dset.name, dtype, ', '.join(self.dtype))
    return True, ""

class check_dims(object):
  '''
  A class to check whether the dataset dimensions matches the expected

  '''

  def __init__(self, dims):
    self.dims = dims

  def __call__(self, dset):
    dims = len(dset.shape)
    if not dims == self.dims:
      return False, '%s has dims %s, expected %s' % (
        dset.name, str(dims), str(self.dims))
    return True, ''

class check_shape(object):
  '''
  A class to check whether the dataset shape matches the expected

  '''

  def __init__(self, shape):
    self.shape = shape

  def __call__(self, dset):
    shape = dset.shape
    if not shape == self.shape:
      return False, '%s has shape %s, expected %s' % (
        dset.name, str(shape), str(self.shape))
    return True, ''

class check_is_scalar(object):
  '''
  A class to check whether the dataset is scalar or not

  '''

  def __init__(self, is_scalar):
    self.is_scalar = is_scalar

  def __call__(self, dset):
    try:
      data = dset[()]
      s = True
    except Exception:
      s = False
    if s != self.is_scalar:
      return False, '%s == scalar is %s, expected %s' % (
        dset.name, s, self.is_scalar)
    return True, ''

class check_dset(object):
  '''
  Check properties of a dataset

  '''

  def __init__(self,
               dtype=None,
               dims=None,
               shape=None,
               is_scalar=None):
    '''
    Set stuff to check
    :param dtype:         The datatype
    :param dims:          The number of dimensions
    :param shape:         The shape of the dataset

    '''
    self.checks = []
    if dtype is not None:
      if not isinstance(dtype, list) and not isinstance(dtype, tuple):
        dtype = [dtype]
      self.checks.append(check_dtype(dtype))
    if dims is not None:
      self.checks.append(check_dims(dims))
    if shape is not None:
      self.checks.append(check_shape(shape))
    if is_scalar is not None:
      self.checks.append(check_is_scalar(is_scalar))

  def __call__(self, dset):
    for check in self.checks:
      passed, errors = check(dset)
      if passed == False:
        raise RuntimeError(errors)


class check_attr(object):
  '''
  Check some properties of an attribute

  '''

  def __init__(self, name, value=None, dtype=None):
    '''
    Set stuff to check
    :param name:  The name of the attribute
    :param value: The value of the attribute
    :param tests: A list of tests to run

    '''
    self.name = name
    self.value = value
    self.dtype = dtype

  def __call__(self, dset):
    if not self.name in dset.attrs.keys():
      raise RuntimeError("'%s' does not have an attribute '%s'" % (
        dset.name, self.name))
    elif self.value is not None and dset.attrs[self.name] != self.value:
      raise RuntimeError("attribute '%s' of %s has value %s, expected %s" % (
        self.name, dset.name, dset.attrs[self.name], self.value))
    elif self.dtype is not None:
      dtype = type(dset.attrs[self.name])
      if not isinstance(dset.attrs[self.name],self.dtype):
        raise RuntimeError("attribute '%s' has type %s, expected %s" % (
          self.name, dset.name, dtype, self.dtype))


def find_entries(nx_file, entry):
  '''
  Find NXmx entries

  '''
  hits = []
  def visitor(name, obj):
    if "NX_class" in obj.attrs.keys():
      if obj.attrs["NX_class"] in ["NXentry", "NXsubentry"]:
        if "definition" in obj.keys():
          if obj["definition"][()] == "NXmx":
            hits.append(obj)
  visitor(entry, nx_file[entry])
  nx_file[entry].visititems(visitor)
  return hits


def find_class(nx_file, nx_class):
  '''
  Find a given NXclass

  '''
  hits = []
  def visitor(name, obj):
    if "NX_class" in obj.attrs.keys():
      if obj.attrs["NX_class"] in [nx_class]:
        hits.append(obj)
  nx_file.visititems(visitor)
  return hits


def convert_units(value, input_units, output_units):
  '''
  Hacky utility function to convert units

  '''
  converters = {
    'm' : {
      'mm'        : lambda x: x * 1e3,
      'microns'   : lambda x: x * 1e6,
      'nm'        : lambda x: x * 1e9
    },
    'mm' : {
      'm'         : lambda x: x * 1e-3,
      'microns'   : lambda x: x * 1e3,
      'nm'        : lambda x: x * 1e6
    },
    'microns' : {
      'm'         : lambda x: x * 1e-6,
      'mm'        : lambda x: x * 1e-3,
      'nm'        : lambda x: x * 1e3
    },
    'nm' : {
      'm'         : lambda x: x * 1e-9,
      'mm'        : lambda x: x * 1e-6,
      'microns'   : lambda x: x * 1e-3,
      'angstroms' : lambda x: x * 10
    }
  }
  if input_units == output_units:
    return value
  try:
    return converters[input_units][output_units](value)
  except Exception, e:
    pass
  raise RuntimeError('Can\'t convert units "%s" to "%s"' % (input_units, output_units))


def visit_dependancies(nx_file, item, visitor = None):
  '''
  Walk the dependency chain and call a visitor function

  '''
  import os.path
  dependency_chain = []
  if os.path.basename(item) == 'depends_on':
    depends_on = nx_file[item][()]
  else:
    depends_on = nx_file[item].attrs['depends_on']
  while not depends_on == ".":
    if visitor is not None:
      visitor(nx_file, depends_on)
    if depends_on in dependency_chain:
      raise RuntimeError("'%s' is a circular dependency" % depends_on)
    try:
      item = nx_file[depends_on]
    except Exception, e:
      raise RuntimeError("'%s' is missing from nx_file" % depends_on)
    dependency_chain.append(depends_on)
    try:
      depends_on = nx_file[depends_on].attrs["depends_on"]
    except Exception:
      raise RuntimeError("'%s' contains no depends_on attribute" % depends_on)


def construct_vector(nx_file, item, vector=None):
  '''
  Walk the dependency chain and create the absolute vector

  '''
  from scitbx import matrix

  class TransformVisitor(object):
    def __init__(self, vector):
      self.vector = matrix.col(vector)

    def __call__(self, nx_file, depends_on):
      from scitbx import matrix
      item = nx_file[depends_on]
      value = item[()]
      units = item.attrs['units']
      ttype = item.attrs['transformation_type']
      vector = matrix.col(item.attrs['vector'])
      if ttype == 'translation':
        value = convert_units(value, units, 'mm')
        self.vector = vector * value + self.vector
      elif ttype == 'rotation':
        if units == 'rad':
          deg = False
        elif units == 'deg':
          deg = True
        else:
          raise RuntimeError('Invalid units: %s' % units)
        self.vector.rotate(axis=vector, angle=value, deg=deg)
      else:
        raise RuntimeError('Unknown transformation_type: %s' % ttype)

    def result(self):
      return self.vector

  if vector is None:
    value = nx_file[item][()]
    units = nx_file[item].attrs['units']
    ttype = nx_file[item].attrs['transformation_type']
    vector = nx_file[item].attrs['vector']
    if ttype == 'translation':
      value = convert_units(value, units, "mm")
      vector = vector * value
  else:
    pass
  visitor = TransformVisitor(vector)

  visit_dependancies(nx_file, item, visitor)

  return visitor.result()


def run_checks(handle, items):
  '''
  Run checks for datasets

  '''
  for item, detail in items.iteritems():
    min_occurs = detail["minOccurs"]
    checks = detail['checks']
    assert(min_occurs in [0, 1])
    try:
      dset = handle[item]
    except Exception:
      dset = None
      if min_occurs != 0:
        raise RuntimeError('Could not find %s in %s' % (item, handle.name))
      else:
        continue
    if dset is not None:
      for check in checks:
        check(dset)


class NXdetector_module(object):
  '''
  A class to hold a handle to NXdetector_module

  '''

  def __init__(self, handle, errors=None):

    self.handle = handle

    items = {
      "data_origin" : {
        "minOccurs" : 1,
        "checks" : [
          check_dset(dtype=["uint32", "uint64", "int32", "int64"], shape=(2,))
        ]
      },
      "data_size" : {
        "minOccurs" : 1,
        "checks" : [
          check_dset(dtype=["int32", "int64", "uint32", "uint64"], shape=(2,))
        ]
      },
      "module_offset" : {
        "minOccurs" : 1,
        "checks" : [
          check_dset(dtype=["float64", "float32", "int64", "int32"], is_scalar=True),
          check_attr("transformation_type"),
          check_attr("vector"),
          check_attr("offset"),
          check_attr("units", dtype=str),
          check_attr("depends_on")
        ]
      },
      "fast_pixel_direction" : {
        "minOccurs" : 1,
        "checks" : [
          check_dset(dtype=["float32", "float64"], is_scalar=True),
          check_attr("transformation_type"),
          check_attr("vector"),
          check_attr("offset"),
          check_attr("units", dtype=str),
          check_attr("depends_on")
        ]
      },
      "slow_pixel_direction" : {
        "minOccurs" : 1,
        "checks" : [
          check_dset(dtype=["float32", "float64"], is_scalar=True),
          check_attr("transformation_type"),
          check_attr("vector"),
          check_attr("offset"),
          check_attr("units", dtype=str),
          check_attr("depends_on"),
        ]
      },
    }

    run_checks(self.handle, items)



class NXdetector(object):
  '''
  A class to handle a handle to NXdetector

  '''

  def __init__(self, handle, errors=None):

    self.handle = handle

    # The items to validate
    items = {
      "depends_on" : {
        "minOccurs" : 1,
        "checks" : []
      },
      "data" : {
        "minOccurs" : 0,
        "checks" : [
          check_dset(dims=3)
        ]
      },
      "description" : {
        "minOccurs" : 1,
        "checks" : []
      },
      "time_per_channel" : {
        "minOccurs" : 0,
        "checks" : []
      },
      "distance" : {
        "minOccurs" : 0,
        "checks" : [
          check_dset(dtype=["float32", "float64"], is_scalar=True)
        ]
      },
      "dead_time" : {
        "minOccurs" : 0,
        "checks" : [
          check_dset(dtype=["float32", "float64"], is_scalar=True)
        ]
      },
      "count_time" : {
        "minOccurs" : 1,
        "checks" : [
          check_dset(dtype=["float32", "float64"], is_scalar=True)
        ]
      },
      "beam_centre_x" : {
        "minOccurs" : 0,
        "checks" : [
          check_dset(dtype=["float32", "float64"], is_scalar=True)
        ]
      },
      "beam_centre_y" : {
        "minOccurs" : 0,
        "checks" : [
          check_dset(dtype=["float32", "float64"], is_scalar=True)
        ]
      },
      "angular_calibration_applied" : {
        "minOccurs" : 0,
        "checks" : [
          check_dset(dtype=['int32', 'int64', 'uint32', 'uint64'], is_scalar=True)
        ]
      },
      "angular_calibration" : {
        "minOccurs" : 0,
        "checks" : [
          check_dset(dtype=["float32", "float64"])
        ]
      },
      "flatfield_applied" : {
        "minOccurs" : 0,
        "checks" : [
          check_dset(dtype=['int32', 'int64', 'uint32', 'uint64'], is_scalar=True)
        ]
      },
      "flatfield" : {
        "minOccurs" : 0,
        "checks" : [
          check_dset(dtype=["float32", "float64"])
        ]
      },
      "flatfield_error" : {
        "minOccurs" : 0,
        "checks" : [
          check_dset(dtype=["float32", "float64"])
        ]
      },
      "pixel_mask_applied" : {
        "minOccurs" : 0,
        "checks" : [
          check_dset(dtype=['int32', 'int64', 'uint32', 'uint64'], is_scalar=True)
        ]
      },
      "pixel_mask" : {
        "minOccurs" : 0,
        "checks" : [
          check_dset(dtype="int32")
        ]
      },
      "countrate_correction_applied" : {
        "minOccurs" : 0,
        "checks" : [
          check_dset(dtype=['int32', 'int64', 'uint32', 'uint64'], is_scalar=True)
        ]
      },
      "bit_depth_readout" : {
        "minOccurs" : 0,
        "checks" : [
          check_dset(dtype=['int32', "int64"], is_scalar=True)
        ]
      },
      "detector_readout_time" : {
        "minOccurs" : 0,
        "checks" : [
          check_dset(dtype=['float32', "float64"], is_scalar=True)
        ]
      },
      "frame_time" : {
        "minOccurs" : 0,
        "checks" : [
          check_dset(dtype=['float32', "float64"], is_scalar=True)
        ]
      },
      "gain_setting" : {
        "minOccurs" : 0,
        "checks" : []
      },
      "saturation_value" : {
        "minOccurs" : 1,
        "checks" : [
          check_dset(dtype=["int32", "int64"], is_scalar=True)
        ]
      },
      "sensor_material" : {
        "minOccurs" : 1,
        "checks" : []
      },
      "sensor_thickness" : {
        "minOccurs" : 1,
        "checks" : [
          check_dset(dtype=["float32", "float64"], is_scalar=True),
          check_attr("units", dtype=str)
        ]
      },
      "threshold_energy" : {
        "minOccurs" : 0,
        "checks" : [
          check_dset(dtype=['float32', "float64"], is_scalar=True)
        ]
      },
      "type" : {
        "minOccurs" : 1,
        "checks" : []
      },
    }

    run_checks(self.handle, items)

    # Find the NXdetector_modules
    self.modules = []
    for entry in find_class(self.handle, "NXdetector_module"):
      try:
        self.modules.append(NXdetector_module(entry, errors=errors))
      except Exception, e:
        if errors is not None:
          errors.append(str(e))

    # Check we've got some stuff
    if len(self.modules) == 0:
      raise RuntimeError('No NXdetector_module in %s' % self.handle.name)


class NXinstrument(object):
  '''
  A class to hold a handle to NXinstrument

  '''

  def __init__(self, handle, errors=None):

    self.handle = handle

    # Find the NXdetector
    self.detectors = []
    for entry in find_class(self.handle, "NXdetector"):
      try:
        self.detectors.append(NXdetector(entry, errors=errors))
      except Exception, e:
        if errors is not None:
          errors.append(str(e))

    # Check we've got stuff
    if len(self.detectors) == 0:
      raise RuntimeError('No NXdetector in %s' % self.handle.name)


class NXbeam(object):
  '''
  A class to hold a handle to NXbeam

  '''

  def __init__(self, handle, errors=None):

    self.handle = handle

    items = {
      "incident_wavelength" : {
        "minOccurs" : 1,
        "checks" : [
          check_dset(dtype=['float32', "float64"], is_scalar=True)
        ]
      },
      "incident_wavelength_spectrum" : {
        "minOccurs" : 0,
        "checks" : []
      },
      "incident_polarization_stokes" : {
        "minOccurs" : 0,
        "checks" : [
          check_dset(dtype=["float32", "float64"], shape=(4,))
        ]
      },
      "flux" : {
        "minOccurs" : 0,
        "checks" : [
          check_dset(dtype=["float32", "float64"], is_scalar=True)
        ]
      },
    }

    run_checks(self.handle, items)


class NXsample(object):
  '''
  A class to hold a handle to NXsample

  '''

  def __init__(self, handle, errors=None):

    self.handle = handle

    items = {
      "name" : {
        "minOccurs" : 0,
        "checks" : []
      },
      "depends_on" : {
        "minOccurs" : 1,
        "checks" : []
      },
      "chemical_formula" : {
        "minOccurs" : 0,
        "checks" : []
      },
      "unit_cell" : {
        "minOccurs" : 0,
        "checks" : [
          check_dset(dtype="float64", dims=2)
        ]
      },
      "unit_cell_class" : {
        "minOccurs" : 0,
        "checks" : []
      },
      "unit_cell_group" : {
        "minOccurs" : 0,
        "checks" : []
      },
      "sample_orientation" : {
        "minOccurs" : 0,
        "checks" : [
          check_dset(dtype="float64", shape=(3,))
        ]
      },
      "orientation_matrix" : {
        "minOccurs" : 0,
        "checks" : [
          check_dset(dtype="float64", dims=3)
        ]
      },
      "temperature" : {
        "minOccurs" : 0,
        "checks" : []
      },
    }

    run_checks(self.handle, items)

    # Find the NXsource
    self.beams = []
    for entry in find_class(self.handle, "NXbeam"):
      try:
        self.beams.append(NXbeam(entry, errors=errors))
      except Exception, e:
        if errors is not None:
          errors.append(str(e))

    # Check we've got stuff
    if len(self.beams) == 0:
      raise RuntimeError('No NXbeam in %s' % self.handle.name)


class NXdata(object):
  '''
  A class to hold a handle to NXdata

  '''

  def __init__(self, handle, errors=None):

    self.handle = handle


class NXmxEntry(object):
  '''
  A class to hold a handle to NXmx entries

  '''

  def __init__(self, handle, errors=None):

    self.handle = handle

    items = {
      'title' : {
        "minOccurs" : 0,
        "checks" : []
      },
      "start_time" : {
        "minOccurs" : 0,
        "checks" : []
      },
      "end_time" : {
        "minOccurs" : 0,
        "checks" : []
      },
    }

    run_checks(self.handle, items)

    # Find the NXinstrument
    self.instruments = []
    for entry in find_class(self.handle, "NXinstrument"):
      try:
        self.instruments.append(NXinstrument(entry, errors=errors))
      except Exception, e:
        if errors is not None:
          errors.append(str(e))

    # Find the NXsample
    self.samples = []
    for entry in find_class(self.handle, "NXsample"):
      try:
        self.samples.append(NXsample(entry, errors=errors))
      except Exception, e:
        if errors is not None:
          errors.append(str(e))

    # Find the NXidata
    self.data = []
    for entry in find_class(self.handle, "NXdata"):
      try:
        self.data.append(NXdata(entry, errors=errors))
      except Exception, e:
        if errors is not None:
          errors.append(str(e))

    # Check we've got some stuff
    if len(self.instruments) == 0:
      raise RuntimeError('No NXinstrument in %s' % self.handle.name)
    if len(self.samples) == 0:
      raise RuntimeError('No NXsample in %s' % self.handle.name)
    if len(self.data) == 0:
      raise RuntimeError('No NXdata in %s' % self.handle.name)


class NXmxReader(object):
  '''
  A hacky class to read an NXmx file

  '''
  def __init__(self, filename=None, handle=None):
    import h5py

    # Get the file handle
    if filename is not None:
      handle = h5py.File(filename, 'r')

    # A list of errors
    self.errors = []

    # Find the NXmx entries
    self.entries = []
    for entry in find_entries(handle, "/"):
      try:
        self.entries.append(NXmxEntry(entry, errors=self.errors))
      except Exception, e:
        self.errors.append(str(e))

    # Check we've got some stuff
    if len(self.entries) == 0:
      raise RuntimeError('''
        Error reading NXmxfile: %s
          No NXmx entries in file

        The following errors occurred:

        %s
      ''' % (filename, "\n".join(self.errors)))

  def print_errors(self):
    '''
    Print any errors that occurred

    '''
    if len(self.errors) > 0:
      print ""
      print "*" * 80
      print "The following errors occurred:\n"
      print "\n".join(self.errors)
      print "*" * 80
      print ""

  def print_description(self):
    '''
    Print a description of the NXmx file

    '''
    print " > Found %d NXmx entries" % len(self.entries)
    for entry in self.entries:
      handle = entry.handle
      instruments = entry.instruments
      samples = entry.samples
      print "  > %s" % handle.name
      for instrument in instruments:
        handle = instrument.handle
        detectors = instrument.detectors
        print "   > %s" % handle.name
        for detector in detectors:
          handle = detector.handle
          modules = detector.modules
          print "    > %s" % handle.name
          for module in modules:
            handle = module.handle
            print "     > %s" % handle.name
      for sample in samples:
        handle = sample.handle
        beams = sample.beams
        print "   > %s" % handle.name
        for beam in beams:
          handle = beam.handle
          print "    > %s" % handle.name


def is_nexus_file(filename):
  '''
  A hacky function to check if this is a nexus file

  '''
  import h5py

  # Get the file handle
  handle = h5py.File(filename, 'r')

  # Find the NXmx entries
  count = 0
  for entry in find_entries(handle, "/"):
    count += 1
  return count > 0


class BeamFactory(object):
  '''
  A class to create a beam model from NXmx stuff

  '''

  def __init__(self, obj):
    from dxtbx.model import Beam

    # Get the items from the NXbeam class
    wavelength = obj.handle['incident_wavelength']
    wavelength_value = wavelength[()]
    wavelength_units = wavelength.attrs['units']

    # Convert wavelength to Angstroms
    wavelength_value = float(convert_units(
      wavelength_value,
      wavelength_units,
      "angstrom"))

    # Construct the beam model
    self.model = Beam(
      direction=(0, 0, -1),
      wavelength=wavelength_value)


class DetectorFactory(object):
  '''
  A class to create a detector model from NXmx stuff

  '''

  def __init__(self, obj, beam):
    from dxtbx.model import Detector, Panel
    from cctbx.eltbx import attenuation_coefficient
    from dxtbx.model import ParallaxCorrectedPxMmStrategy
    from scitbx import matrix

    # Get the handles
    nx_file = obj.handle.file
    nx_detector = obj.handle
    nx_module = obj.modules[0].handle

    # Get the detector name and type
    detector_type = str(nx_detector['type'][()])
    detector_name = str(nx_detector.name)

    # Get the trusted range of pixel values
    trusted_range = (-1, float(nx_detector['saturation_value'][()]))

    # Get the detector thickness
    thickness = nx_detector['sensor_thickness']
    thickness_value = float(thickness[()])
    thickness_units = thickness.attrs['units']
    thickness_value = float(convert_units(
      thickness_value,
      thickness_units,
      "mm"))

    # Get the detector material
    material = str(nx_detector['sensor_material'][()])

    # Get the fast pixel size and vector
    fast_pixel_direction = nx_module['fast_pixel_direction']
    fast_pixel_direction_value = float(fast_pixel_direction[()])
    fast_pixel_direction_units = fast_pixel_direction.attrs['units']
    fast_pixel_direction_vector = fast_pixel_direction.attrs['vector']
    fast_pixel_direction_value = convert_units(
      fast_pixel_direction_value,
      fast_pixel_direction_units,
      "mm")
    fast_axis = matrix.col(fast_pixel_direction_vector).normalize()

    # Get the slow pixel size and vector
    slow_pixel_direction = nx_module['slow_pixel_direction']
    slow_pixel_direction_value = float(slow_pixel_direction[()])
    slow_pixel_direction_units = slow_pixel_direction.attrs['units']
    slow_pixel_direction_vector = slow_pixel_direction.attrs['vector']
    slow_pixel_direction_value = convert_units(
      slow_pixel_direction_value,
      slow_pixel_direction_units,
      "mm")
    slow_axis = matrix.col(slow_pixel_direction_vector).normalize()

    # Get the origin vector
    module_offset = nx_module['module_offset']
    origin = construct_vector(
      nx_file,
      module_offset.name)

    origin = -matrix.col(origin)
    fast_axis = - fast_axis
    slow_axis = - slow_axis

    # Ensure that fast and slow axis are orthogonal
    normal = fast_axis.cross(slow_axis)
    slow_axis = -fast_axis.cross(normal)

    # Compute the attenuation coefficient.
    # This will fail for undefined composite materials
    # mu_at_angstrom returns cm^-1, but need mu in mm^-1
    if material == 'Si':
      pass
    elif material == 'Silicon':
      material = 'Si'
    elif material == 'Sillicon':
      material = 'Si'
    elif material == 'CdTe':
      pass
    elif material == 'GaAs':
      pass
    else:
      raise RuntimeError('Unknown material: %s' % material)
    table = attenuation_coefficient.get_table(material)
    wavelength = beam.get_wavelength()
    mu = table.mu_at_angstrom(wavelength) / 10.0

    # Construct the detector model
    pixel_size = (fast_pixel_direction_value, slow_pixel_direction_value)
    image_size = tuple(map(int, nx_module['data_size']))

    self.model = Detector()
    self.model.add_panel(
      Panel(
        detector_type,
        detector_name,
        tuple(fast_axis),
        tuple(slow_axis),
        tuple(origin),
        pixel_size,
        image_size,
        trusted_range,
        thickness_value,
        material,
        mu))

    # Set the parallax correction
    for panel in self.model:
      panel.set_px_mm_strategy(ParallaxCorrectedPxMmStrategy(mu, thickness_value))
      panel.set_type('SENSOR_PAD')

class GoniometerFactory(object):
  '''
  A class to create a goniometer model from NXmx stuff

  '''
  def __init__(self, obj):
    from dxtbx.model import Goniometer

    # Get the rotation axis
    rotation_axis = construct_vector(
      obj.handle.file,
      obj.handle.file[obj.handle['depends_on'][()]].name)

    # Construct the model
    self.model = Goniometer(
      tuple(rotation_axis))


class ScanFactory(object):
  '''
  A class to create a scan model from NXmx stuff

  '''
  def __init__(self, obj, detector_obj):
    from dxtbx.model import Scan
    from scitbx.array_family import flex

    # Get the image and oscillation range
    phi = obj.handle.file[obj.handle['depends_on'][()]]
    image_range = (1, len(phi))
    if len(phi) > 1:
      oscillation = (float(phi[0]), float(phi[1]-phi[0]))
      is_sweep = True
    else:
      oscillation = (float(phi[0]), 0.0)
      is_sweep = False

    # Get the exposure time
    num_images = len(phi)
    if "frame_time" in detector_obj.handle:
      frame_time = float(detector_obj.handle['frame_time'][()])
      exposure_time = flex.double(num_images, frame_time)
      epochs = flex.double(num_images)
      for i in range(1, len(epochs)):
        epochs[i] = epochs[i-1] + exposure_time[i-1]
    else:
      exposure_time = flex.double(num_images, 0)
      epochs = flex.double(num_images, 0)

    if is_sweep is True:

      # Construct the model
      self.model = Scan(
        image_range,
        oscillation,
        exposure_time,
        epochs)

    else:

      self.model = []
      for i, image in enumerate(range(image_range[0], image_range[1]+1)):
        self.model.append(Scan(
          (image, image),
          oscillation,
          exposure_time[i:i+1],
          epochs[i:i+1]))


class CrystalFactory(object):
  '''
  A class to create a crystal model from NXmx stuff

  '''

  def __init__(self, obj):
    from dxtbx.model import Crystal
    import cctbx.uctbx
    from scitbx import matrix

    # Get the crystal parameters
    unit_cell_parameters = list(obj.handle['unit_cell'][0])
    unit_cell = cctbx.uctbx.unit_cell(unit_cell_parameters)
    U = list(obj.handle['orientation_matrix'][0].flatten())
    U = matrix.sqr(U)
    B = matrix.sqr(unit_cell.fractionalization_matrix()).transpose()
    A = U * B
    Ai = A.inverse()
    real_space_a = Ai[0:3]
    real_space_b = Ai[3:6]
    real_space_c = Ai[6:9]

    # Get the space group symbol
    space_group_symbol = obj.handle['unit_cell_group'][()]

    # Create the model
    self.model = Crystal(
      real_space_a,
      real_space_b,
      real_space_c,
      space_group_symbol)


class DataList(object):
  '''
  A class to make it easier to access the data from multiple datasets.
  FIXME The file should be fixed and this should be removed

  '''

  def __init__(self, obj):
    self.datasets = obj
    self.num_images = 0
    self.lookup = []
    self.offset = [0]
    for i, dataset in enumerate(self.datasets):
      self.num_images += dataset.shape[0]
      self.lookup.extend([i] * dataset.shape[0])
      self.offset.append(self.num_images)
    shape = self.datasets[0].shape
    self.height = shape[1]
    self.width = shape[2]

  def __len__(self):
    return self.num_images

  def __getitem__(self, index):
    from scitbx.array_family import flex
    import numpy as np
    d = self.lookup[index]
    i = index - self.offset[d]

    # Keeping this in for the moment to allow evaluation of speed etc
    # aiming to resolve dials#148
    # mode_148 = True

    # if mode_148:
    #   # allocate empty array, copy data in
    #   data = np.empty((self.height, self.width), dtype='uint32')
    #   self.datasets[d].read_direct(data, np.s_[i,:,:], np.s_[:,:])
    # else:
    #   data = self.datasets[d][i,:,:]
    #   if data.dtype == np.uint16:
    #     data = data.astype(np.uint32)
    # data_as_flex = flex.int(data)
    N, height, width = self.datasets[d].shape
    data_as_flex = dataset_as_flex_int(
      self.datasets[d].id.id,
      (slice(i, i+1, 1),
       slice(0, height, 1),
       slice(0, width, 1)))
    data_as_flex.reshape(flex.grid(data_as_flex.all()[1:]))
    return data_as_flex


class DataFactory(object):
  def __init__(self, obj):
    import h5py
    datasets = []
    for key in sorted(list(obj.handle.iterkeys())):
      if key.startswith("_filename_"):
        continue
      try:
        datasets.append(obj.handle[key])
      except KeyError: # If we cannot follow links due to lack of a write permission
        datasets.append(h5py.File(obj.handle["_filename_" + key].value, "r")["/entry/data/data"])

    self.model = DataList(datasets)


class MaskFactory(object):
  '''
  A class to create an object to hold the pixel mask data

  '''

  def __init__(self, obj):
    handle = obj.handle
    def make_mask(dset):
      from dials.array_family import flex
      height, width = dset.shape
      mask_data = flex.int(dset[:,:].flatten()) == 0
      mask_data.reshape(flex.grid(height, width))
      return mask_data
    self.mask = None
    if "pixel_mask_applied" in handle and handle['pixel_mask_applied']:
      if "pixel_mask" in handle:
        self.mask = make_mask(handle['pixel_mask'])
      elif "detectorSpecific" in handle:
        if "pixel_mask" in handle["detectorSpecific"]:
          self.mask = make_mask(handle["detectorSpecific"]["pixel_mask"])
