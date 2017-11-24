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
  # Workaround for psana build, which doesn't link HDF5 properly
  import os
  if 'SIT_ROOT' not in os.environ:
    raise

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
    if not shape in self.shape:
      return False, '%s has shape %s, expected one of %s' % (
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
    },
    'angstroms' : {
      'angstrom'  : lambda x: x
    }
  }
  if input_units == output_units:
    return value
  try:
    return converters[input_units][output_units](value)
  except Exception:
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
    except Exception:
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
        if hasattr(value, '__iter__') and len(value) == 1:
          value = value[0]
        self.vector = vector * value + self.vector
      elif ttype == 'rotation':
        if hasattr(value, '__iter__') and len(value):
          value = value[0]
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
    if 'offset' in nx_file[item].attrs:
      offset = nx_file[item].attrs['offset']
    else:
      offset = vector * 0.0
    if ttype == 'translation':
      value = convert_units(value, units, "mm")
      try:
        vector = vector * value
      except ValueError:
        vector = vector * value[0]
      vector += offset

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
          check_dset(dtype=["uint32", "uint64", "int32", "int64"], shape=[(2,), (3,)])
        ]
      },
      "data_size" : {
        "minOccurs" : 1,
        "checks" : [
          check_dset(dtype=["int32", "int64", "uint32", "uint64"], shape=[(2,), (3,)])
        ]
      },
      "module_offset" : {
        "minOccurs" : 0,
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
        "minOccurs" : 0,
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
        "minOccurs" : 0,
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

class NXdetector_group(object):
  '''
  A class to hold a handle to NXdetector_group

  '''

  def __init__(self, handle, errors=None):

    self.handle = handle

    items = {
      "group_names" : {
        "minOccurs" : 0,
        "checks" : [
        ]
      },
      "group_index" : {
        "minOccurs" : 0,
        "checks" : [
          check_dset(dtype=['int32', 'int64', 'uint32', 'uint64'], is_scalar=True)
        ]
      },
      "group_parent" : {
        "minOccurs" : 0,
        "checks" : [
          check_dset(dtype=['int32', 'int64', 'uint32', 'uint64'], is_scalar=True)
        ]
      },
      "group_type" : {
        "minOccurs" : 0,
        "checks" : [
          check_dset(dtype=['int32', 'int64', 'uint32', 'uint64'], is_scalar=True)
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
        "minOccurs" : 0,
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
        "minOccurs" : 0,
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
          check_dset(dtype=['uint32', 'int32'])
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
        "minOccurs" : 0,
        "checks" : [
          check_dset(dtype=["int32", "int64"], is_scalar=True)
        ]
      },
      "sensor_material" : {
        "minOccurs" : 0,
        "checks" : []
      },
      "sensor_thickness" : {
        "minOccurs" : 0,
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
        "minOccurs" : 0,
        "checks" : []
      },
    }

    run_checks(self.handle, items)

    # Find the NXdetector_modules
    self.modules = []
    for entry in find_class(self.handle, "NXdetector_module"):
      try:
        self.modules.append(NXdetector_module(entry, errors=errors))
      except Exception as e:
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
      except Exception as e:
        if errors is not None:
          errors.append(str(e))

    # Check we've got stuff
    if len(self.detectors) == 0:
      raise RuntimeError('No NXdetector in %s' % self.handle.name)

    # Find any detector groups
    self.detector_groups = []
    for entry in find_class(self.handle, "NXdetector_group"):
      try:
        self.detector_groups.append(NXdetector_group(entry, errors=errors))
      except Exception as e:
        if errors is not None:
          errors.append(str(e))

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
          check_dset(dtype=["float32", "float64"], shape=[(4,)])
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
        "minOccurs" : 0,
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
          check_dset(dtype="float64", shape=[(3,)])
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
      except Exception as e:
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
      except Exception as e:
        if errors is not None:
          errors.append(str(e))

    # Find the NXsample
    self.samples = []
    for entry in find_class(self.handle, "NXsample"):
      try:
        self.samples.append(NXsample(entry, errors=errors))
      except Exception as e:
        if errors is not None:
          errors.append(str(e))

    # Find the NXidata
    self.data = []
    for entry in find_class(self.handle, "NXdata"):
      try:
        self.data.append(NXdata(entry, errors=errors))
      except Exception as e:
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
      except Exception as e:
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

  def __init__(self, obj, index = None):
    from dxtbx.model import Beam

    # Get the items from the NXbeam class
    wavelength = obj.handle['incident_wavelength']
    wavelength_value = wavelength[()]
    wavelength_units = wavelength.attrs['units']

    if index is not None and hasattr(wavelength_value, '__iter__') and \
          len(wavelength_value) > 1:
      wavelength_value = wavelength_value[index]

    # Convert wavelength to Angstroms
    wavelength_value = float(convert_units(
      wavelength_value,
      wavelength_units,
      "angstrom"))

    # Construct the beam model
    self.model = Beam(
      direction=(0, 0, 1),
      wavelength=wavelength_value)

def get_change_of_basis(transformation):
  '''
  Get the 4x4 homogenous coordinate matrix for a given NXtransformation.

  '''
  from scitbx.matrix import col, sqr
  # Change of basis to convert from NeXus to IUCr/ImageCIF convention
  n2i_cob = sqr((-1,  0,  0,
                  0,  1,  0,
                  0,  0, -1))

  axis_type = transformation.attrs['transformation_type']

  vector = n2i_cob * col(transformation.attrs['vector']).normalize()
  setting = transformation[0]
  units = transformation.attrs['units']

  if 'offset' in transformation.attrs:
    offset = n2i_cob * col(transformation.attrs['offset'])
    offset = convert_units(offset, transformation.attrs['offset_units'], 'mm')
  else:
    offset = col((0,0,0))

  # 4x4 change of basis matrix (homogeneous coordinates)
  cob = None

  if axis_type == "rotation":
    if units == 'rad':
      deg = False
    elif units in ['deg', 'degrees']:
      deg = True
    else:
      raise RuntimeError('Invalid units: %s' % units)
    r3 = vector.axis_and_angle_as_r3_rotation_matrix(setting, deg = deg)
    cob = sqr((r3[0], r3[1], r3[2], offset[0],
               r3[3], r3[4], r3[5], offset[1],
               r3[6], r3[7], r3[8], offset[2],
               0,     0,     0,     1))
  elif axis_type == "translation":
    setting = convert_units(setting, units, 'mm')
    translation = offset + (vector * setting)
    cob = sqr((1,0,0,translation[0],
               0,1,0,translation[1],
               0,0,1,translation[2],
               0,0,0,1))
  else:
    raise Sorry("Unrecognized tranformation type: %d"%axis_type)

  return cob

def get_depends_on_chain_using_equipment_components(transformation):
  '''
  Given an NXtransformation, find the list of dependencies that transformation
  has.  If there are multiple dependencies in a single 'level', as indicated
  by grouping them using equipment_component, then the dependency chain will
  skip the intermediate dependencies, listing only the first at each level.

  '''
  import os
  chain = []
  current = transformation

  while True:
    parent_id = current.attrs['depends_on']

    if parent_id == ".":
      return chain
    parent = current.parent[parent_id]

    if 'equipment_component' in current.attrs:
      eq_comp = current.attrs['equipment_component']
      parent_eq_comp = parent.attrs['equipment_component']
      if eq_comp == parent_eq_comp:
        current = parent
        continue
    chain.append(parent)
    current = parent

def get_cummulative_change_of_basis(transformation):
  ''' Get the 4x4 homogenous coordinate matrix for a given NXtransformation,
  combining it with the change of basis matrices of parent transformations
  with the same equipment component as the given transformation.
  Returns (parent, change of basis matrix), where parent is None if the
  transformation's depends_on is ".".  Parent is the tranformation that the
  top level transformation in this chain of transformations depends on.

  '''

  cob = get_change_of_basis(transformation)

  parent_id = transformation.attrs['depends_on']

  if parent_id == ".":
    return None, cob
  parent = transformation.parent[parent_id]

  if 'equipment_component' in transformation.attrs:
    eq_comp = transformation.attrs['equipment_component']
    parent_eq_comp = parent.attrs['equipment_component']
    if eq_comp == parent_eq_comp:
      non_matching_parent, parent_cob = get_cummulative_change_of_basis(parent)
      return non_matching_parent, parent_cob * cob

  return parent, cob

class DetectorFactoryFromGroup(object):
  '''
  A class to create a detector model from a NXdetector_group

  '''

  def __init__(self, instrument, beam, idx = None):
    from dxtbx.model import Detector, Panel
    from cctbx.eltbx import attenuation_coefficient
    from dxtbx.model import ParallaxCorrectedPxMmStrategy
    from scitbx import matrix
    import os

    if idx is None:
      idx = 0

    assert len(instrument.detector_groups) == 1, "Multiple detectors not supported"

    nx_group = instrument.detector_groups[0].handle
    group_names = nx_group['group_names']
    group_indices = nx_group['group_index']
    group_parent = nx_group['group_parent']

    # Verify the NXdetector objects specified by the detector group are present
    expected_detectors = []
    root_name = None
    for i, parent_id in enumerate(group_parent):
      assert parent_id in [-1, 1], "Hierarchy of detectors not supported. Hierarchy of module components within detector elements is supported"

      if parent_id == -1:
        assert root_name is None, "Multiple roots not supported"
        root_name = group_names[i]
      else:
        expected_detectors.append(group_names[i])

    assert root_name is not None, "Detector root not found"
    assert sorted([os.path.basename(d.handle.name) for d in instrument.detectors]) == sorted(expected_detectors), \
      "Mismatch between detector group names and detectors available"

    root = None

    def set_frame(pg, transformation):
      ''' Function to set the local frame of a panel/panel group given an
      NXtransformation '''
      parent, cob = get_cummulative_change_of_basis(transformation)
      # set up the dxtbx d matrix.  Note use of homogenous coordinates.
      origin = matrix.col((cob * matrix.col((0,0,0,1)))[0:3])
      fast   = matrix.col((cob * matrix.col((1,0,0,1)))[0:3]) - origin
      slow   = matrix.col((cob * matrix.col((0,1,0,1)))[0:3]) - origin

      pg.set_local_frame(
        fast.elems,
        slow.elems,
        origin.elems)

      name = str(os.path.basename(transformation.name))
      pg.set_name(name)

    self.model = Detector()
    for nx_detector in instrument.detectors:
      # Get the detector type
      if 'type' in nx_detector.handle:
        detector_type = str(nx_detector.handle['type'][()])
      else:
        detector_type = "unknown"

      # get the set of leaf modules (handles nested modules)
      modules = []
      for nx_detector_module in nx_detector.modules:
        if len(find_class(nx_detector_module.handle, "NXdetector_module")) == 0:
          modules.append(nx_detector_module)

      n_modules = len(modules)
      # depends_on field for a detector will have NxM entries in it, where
      # N = number of images and M = number of detector modules
      assert len(nx_detector.handle['depends_on']) % n_modules == 0, "Couldn't find NXdetector_modules matching the list of dependencies for this detector"

      for module_number, nx_detector_module in enumerate(modules):
        # Get the depends_on starting point for this module and for this image
        depends_on = nx_detector.handle[nx_detector.handle['depends_on'][(idx//n_modules)+module_number]]
        panel_name = str(os.path.basename(nx_detector_module.handle.name))
        px_fast = int(nx_detector_module.handle['data_size'][0])
        px_slow = int(nx_detector_module.handle['data_size'][1])
        image_size = px_fast, px_slow
        fast_pixel_size = nx_detector_module.handle['fast_pixel_size'][()]
        slow_pixel_size = nx_detector_module.handle['slow_pixel_size'][()]
        pixel_size = fast_pixel_size, slow_pixel_size
        # XXX no units for pixel size in nx example file, plus fast_pixel_size isn't a standardized NXdetector_module field
        #fast_pixel_direction_units = fast_pixel_direction.attrs['units']
        #fast_pixel_size = convert_units(
        #  fast_pixel_direction_value,
        #  fast_pixel_direction_units,
        #  "mm")
        # Get the trusted range of pixel values
        underload = float(nx_detector['undefined_value'][()])  if 'undefined_value'  in nx_detector.handle else -400
        overload  = float(nx_detector['saturation_value'][()]) if 'saturation_value' in nx_detector.handle else 90000
        trusted_range = underload, overload

        # Set up the hierarchical detector by iteraing through the dependencies,
        # starting at the root
        chain = get_depends_on_chain_using_equipment_components(depends_on)
        pg = None
        for transform in reversed(chain):
          name = str(os.path.basename(transform.name))
          if pg is None: # The first transform will be the root of the hiearchy
            if root is None:
              root = self.model.hierarchy()
              set_frame(root, transform)
            else:
              assert root.get_name() == name, "Found multiple roots"
            pg = root
            continue

          # subsequent times through the loop, pg will be the parent of the
          # current transform
          pg_names = [child.get_name() for child in pg]
          if name in pg_names:
            pg = pg[pg_names.index(name)]
          else:
            pg = pg.add_group()
            set_frame(pg, transform)

        # pg is now this panel's parent
        p = pg.add_panel()
        fast = depends_on.attrs['vector']
        fast = matrix.col([-fast[0], fast[1], -fast[2]])
        slow = nx_detector_module.handle[depends_on.attrs['depends_on']].attrs['vector']
        slow = matrix.col([-slow[0], slow[1], -slow[2]])
        parent, cob = get_cummulative_change_of_basis(depends_on)
        origin = matrix.col((cob * matrix.col((0,0,0,1)))[0:3])

        p.set_local_frame(
          fast.elems,
          slow.elems,
          origin.elems)

        p.set_name(panel_name)
        p.set_pixel_size(pixel_size)
        p.set_image_size(image_size)
        p.set_type(detector_type)
        p.set_trusted_range(trusted_range)

        if 'sensor_thickness' in nx_detector.handle:
          # Get the detector thickness
          thickness = nx_detector.handle['sensor_thickness']
          thickness_value = float(thickness[()])
          thickness_units = thickness.attrs['units']
          thickness_value = float(convert_units(
            thickness_value,
            thickness_units,
            "mm"))
          p.set_thickness(thickness_value)

        # Get the detector material
        if 'sensor_material' in nx_detector.handle:
          material = str(nx_detector.handle['sensor_material'][()])
          p.set_material(material)

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
          p.set_mu(table.mu_at_angstrom(wavelength) / 10.0)

        if 'sensor_thickness' in nx_detector.handle and 'sensor_material' in nx_detector.handle:
          p.set_px_mm_strategy(ParallaxCorrectedPxMmStrategy(mu, thickness_value))

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
    if 'type' in nx_detector:
      detector_type = str(nx_detector['type'][()])
    else:
      detector_type = "unknown"
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

    # Change of basis to convert from NeXus to IUCr/ImageCIF convention
    cob = matrix.sqr((-1,  0,  0,
                       0,  1,  0,
                       0,  0, -1))
    origin = cob * matrix.col(origin)
    fast_axis = cob * fast_axis
    slow_axis = cob * slow_axis

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

    def unroll(iterable):
      result = []
      for i in iterable:
        if hasattr(i, '__iter__'):
          result.extend(i)
        else:
          result.append(float(i))
      return result

    # Construct the model - if nonsense present cope by failing over...
    try:
      axis = tuple(unroll(rotation_axis))
      self.model = Goniometer(axis)
    except Exception:
      self.model = Goniometer((1, 0, 0))

def find_goniometer_rotation(obj):
  thing = obj.handle.file[obj.handle['depends_on'][()]]
  tree = get_depends_on_chain_using_equipment_components(thing)
  for t in tree:
    o = obj.handle.file[t.name]
    if o.attrs['transformation_type'] == 'rotation':
      # if this is changing, assume is scan axis
      v = o[()]
      if min(v) < max(v): return o
  raise ValueError('no rotation found')

class ScanFactory(object):
  '''
  A class to create a scan model from NXmx stuff

  '''
  def __init__(self, obj, detector_obj):
    from dxtbx.model import Scan
    from scitbx.array_family import flex

    # Get the image and oscillation range - need to search for rotations
    # in dependency tree
    try:
      phi = find_goniometer_rotation(obj)
    except ValueError:
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

      # if this is not a sweep, then this is stills, in which case... should
      # we have scans? and, even if we do, we have a small problem in that
      # this returns a list of scans which even less works...

      self.model = []
      for i, image in enumerate(range(image_range[0], image_range[1]+1)):
        self.model.append(Scan(
          (image, image),
          oscillation,
          exposure_time[i:i+1],
          epochs[i:i+1]))

      self.model = None


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

  def __init__(self, obj, max_size=0):
    self.datasets = obj
    self.num_images = 0
    self.lookup = []
    self.offset = [0]

    if len(self.datasets) == 1 and max_size:
      self.num_images = max_size
      self.lookup.extend([0] * max_size)
      self.offset.append(max_size)

    else:
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

class DetectorGroupDataList(object):
  '''
  A class to make it easier to access the data from multiple datasets.
  This version brings in all the panels from a detector group with several detectors.

  '''
  def __init__(self, datalists):
    self.datalists = datalists
    lengths = [len(datalist) for datalist in datalists]
    self.num_images = lengths[0]
    assert all([l == self.num_images for l in lengths]), "Not all datasets are the same length"

  def __len__(self):
    return self.num_images

  def __getitem__(self, index):
    data = []
    for datalist in self.datalists:
      data.extend(datalist[index])
    return tuple(data)

def get_detector_module_slices(detector):
  '''
  Helper function to read data_origin and data_size from the NXdetector_modules in a
  NXdetector.  Returns a list of lists, where each sublist is a list of slices in
  slow to fast order.
  Assumes slices are stored in NeXus in fast to slow order.

  '''
  # get the set of leaf modules (handles nested modules)
  modules = []
  for nx_detector_module in detector.modules:
    if len(find_class(nx_detector_module.handle, "NXdetector_module")) == 0:
      modules.append(nx_detector_module)

  all_slices = []
  for module in modules:
    data_origin = module.handle['data_origin']
    data_size = module.handle['data_size']
    all_slices.append(list(reversed([slice(int(start), int(start+step), 1) for start, step in zip(data_origin, data_size)])))
  return all_slices

class MultiPanelDataList(object):
  '''
  A class to make it easier to access the data from multiple datasets.
  Also handles multi-panel data as described in a series of NXdetector_modules

  '''

  def __init__(self, datasets, detector):
    self.datasets = datasets
    self.num_images = 0
    self.lookup = []
    self.offset = [0]
    for i, dataset in enumerate(self.datasets):
      self.num_images += dataset.shape[0]
      self.lookup.extend([i] * dataset.shape[0])
      self.offset.append(self.num_images)

    self.all_slices = get_detector_module_slices(detector)

  def __len__(self):
    return self.num_images

  def __getitem__(self, index):
    from scitbx.array_family import flex
    import numpy as np
    d = self.lookup[index]
    i = index - self.offset[d]

    all_data = []

    for module_slices in self.all_slices:
      slices = [slice(i, i+1, 1)]
      slices.extend(module_slices)
      data_as_flex = dataset_as_flex_int(
        self.datasets[d].id.id,
        tuple(slices))
      data_as_flex.reshape(flex.grid(data_as_flex.all()[-2:])) # handle 3 or 4 dimension arrays
      all_data.append(data_as_flex)
    return tuple(all_data)

class DataFactory(object):
  def __init__(self, obj, max_size=0):
    import h5py
    datasets = []
    for key in sorted(list(obj.handle.iterkeys())):
      if key.startswith("_filename_"):
        continue

      # datasets in this context mean ones which contain diffraction images
      # so must have ndim > 1 - for example omega can also be nexus data set
      # but with 1 dimension...
      if obj.handle[key].ndim == 1:
        continue

      try:
        datasets.append(obj.handle[key])
      except KeyError: # If we cannot follow links due to lack of a write permission
        datasets.append(h5py.File(obj.handle["_filename_" + key].value, "r")["/entry/data/data"])

    self._datasets = datasets

    self.model = DataList(datasets, max_size=max_size)

class DetectorGroupDataFactory(DataFactory):
  """ Class to handle reading data from a detector with a NXdetector_group """
  def __init__(self, obj, instrument):
    import os
    DataFactory.__init__(self, obj)

    # Map NXdetector names to list of datasets
    mapping = {}
    for dataset in self._datasets:
      dataset_name = os.path.basename(dataset.name)
      found_it = False
      for detector in instrument.detectors:
        if dataset_name in detector.handle:
          found_it = True
          detector_name = os.path.basename(detector.handle.name)
          if detector_name in mapping:
            assert dataset_name not in mapping[detector_name]['dataset_names'], "Dataset %s found in > 1 NXdetectors"%dataset_name
            mapping[detector_name]['dataset_names'].append(dataset_name)
            mapping[detector_name]['datasets'].append(dataset)
          else:
            mapping[detector_name] = { 'dataset_names': [dataset_name],
                                       'datasets': [dataset],
                                       'detector': detector }
      assert found_it, "Couldn't match dataset %s to a NXdetector"%dataset_name

    # Create a list of multipanel datalist objects
    datalists = []
    for detector_name in mapping:
      datalists.append(MultiPanelDataList(mapping[detector_name]['datasets'], mapping[detector_name]['detector']))
    self.model = DetectorGroupDataList(datalists)

class MaskFactory(object):
  '''
  A class to create an object to hold the pixel mask data

  '''

  def __init__(self, objects):
    def make_mask(dset):
      from dials.array_family import flex
      mask_data = flex.int(dset[()]) == 0
      mask = []
      for slices in all_slices:
        mask.append(mask_data[slices])
      return tuple(mask)
    self.mask = None
    for obj in objects:
      handle = obj.handle
      if "pixel_mask_applied" in handle and handle['pixel_mask_applied']:
        if self.mask is None:
          self.mask = []
        if "pixel_mask" in handle:
          all_slices = get_detector_module_slices(obj)
          self.mask.extend(list(make_mask(handle['pixel_mask'])))
        elif "detectorSpecific" in handle:
          if "pixel_mask" in handle["detectorSpecific"]:
            all_slices = get_detector_module_slices(obj)
            self.mask.extend(list(make_mask(handle["detectorSpecific"]["pixel_mask"])))
    if self.mask is not None:
      self.mask = tuple(self.mask)
