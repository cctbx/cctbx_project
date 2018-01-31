#!/usr/bin/env python
#
#  experiment_list.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import absolute_import, division
import pkg_resources
from dxtbx.model import Experiment, ExperimentList

class InvalidExperimentListError(RuntimeError):
  """
  Indicates an error whilst validating the experiment list.

  This means that there is some structural problem that prevents the given data
  from representing a well-formed experiment list. This doesn't indicate e.g.
  some problem with the data or model consistency.
  """

class ExperimentListDict(object):
  ''' A helper class for serializing the experiment list to dictionary (needed
  to save the experiment list to JSON format. '''

  def __init__(self, obj, check_format=True):
    ''' Initialise. Copy the dictionary. '''
    from copy import deepcopy
    # Basic check: This is a dict-like object. This can happen if e.g. we
    # were passed a DataBlock list instead of an ExperimentList dictionary
    if isinstance(obj, list) or not hasattr(obj, "get"):
      raise InvalidExperimentListError("Expected dictionary, not {}".format(type(obj)))

    self._obj = deepcopy(obj)
    self._check_format = check_format

  def decode(self):
    ''' Decode the dictionary into a list of experiments. '''

    # If this doesn't claim to be an ExperimentList, don't even try
    if not self._obj.get('__id__', None) == 'ExperimentList':
      raise InvalidExperimentListError(
        "Expected __id__ 'ExperimentList', but found {}".format(repr(self._obj.get("__id__"))))

    # Extract lists of models referenced by experiments
    self._blist = self._extract_models('beam')
    self._dlist = self._extract_models('detector')
    self._glist = self._extract_models('goniometer')
    self._slist = self._extract_models('scan')
    self._clist = self._extract_models('crystal')
    self._plist = self._extract_models('profile')
    self._scalelist = self._extract_models('scaling_model')

    # Go through all the imagesets and make sure the dictionary
    # references by an index rather than a file path. Experiments
    # referencing the same imageset will get different objects
    # due to the fact that we can have different models
    self._ilist = self._extract_imagesets()

    # Extract all the experiments
    return self._extract_experiments()

  def _extract_models(self, name):
    ''' Helper function. Extract the models. '''

    # The from dict function
    from_dict = getattr(self, '_%s_from_dict' % name)

    # Extract all the model list
    mlist = self._obj.get(name, [])

    # Convert the model from dictionary to concreate
    # python class for the model.
    mlist = [from_dict(d) for d in mlist]

    # Dictionaries for file mappings
    mmap = {}

    # For each experiment, check the model is not specified by
    # a path, if it is then get the dictionary of the model
    # and insert it into the list. Replace the path reference
    # with an index
    for eobj in self._obj['experiment']:
      value = eobj.get(name)
      if value is None:
        continue
      elif isinstance(value, str):
        if value not in mmap:
          mmap[value] = len(mlist)
          mlist.append(from_dict(ExperimentListDict._from_file(value)))
        eobj[name] = mmap[value]
      elif not isinstance(value, int):
        raise TypeError('expected int or str, got %s' % type(value))

    # Return the model list
    return mlist

  def _extract_imagesets(self):
    ''' Helper function, extract the imagesets. '''

    # Extract all the model list
    mlist = self._obj.get('imageset', [])

    # Dictionaries for file mappings
    mmap = {}

    # For each experiment, check the imageset is not specified by
    # a path, if it is then get the dictionary of the imageset
    # and insert it into the list. Replace the path reference
    # with an index
    for eobj in self._obj['experiment']:
      value = eobj.get('imageset')
      if value is None:
        continue
      elif isinstance(value, str):
        if value not in mmap:
          mmap[value] = len(mlist)
          mlist.append(ExperimentListDict._from_file(value))
        eobj['imageset'] = mmap[value]
      elif not isinstance(value, int):
        raise TypeError('expected int or str, got %s' % type(value))

    # Return the model list
    return mlist

  def _extract_experiments(self):
    ''' Helper function. Extract the experiments. '''
    from dxtbx.imageset import ImageSweep, ImageSet, ImageGrid
    from dxtbx.serialize.filename import load_path
    import cPickle as pickle
    from dxtbx.format.image import ImageBool, ImageDouble

    # Map of imageset/scan pairs
    imagesets = {}

    # For every experiment, use the given input to create
    # a sensible experiment.
    el = ExperimentList()
    for eobj in self._obj['experiment']:

      # Get the models
      beam = ExperimentListDict.model_or_none(self._blist, eobj, 'beam')
      detector = ExperimentListDict.model_or_none(self._dlist, eobj, 'detector')
      goniometer = ExperimentListDict.model_or_none(self._glist, eobj, 'goniometer')
      scan = ExperimentListDict.model_or_none(self._slist, eobj, 'scan')
      crystal = ExperimentListDict.model_or_none(self._clist, eobj, 'crystal')
      profile = ExperimentListDict.model_or_none(self._plist, eobj, 'profile')
      scaling_model = ExperimentListDict.model_or_none(self._scalelist, eobj, 'scaling_model')
      key = (eobj.get('imageset'), eobj.get('scan'))
      try:
        imageset = imagesets[key]
      except Exception:
        imageset = ExperimentListDict.model_or_none(self._ilist, eobj, 'imageset')

        # Create the imageset from the input data
        if imageset is not None:
          if 'params' in imageset:
            format_kwargs = imageset['params']
          else:
            format_kwargs = {}
          if 'mask' in imageset and imageset['mask'] is not None:
            mask_filename = load_path(imageset['mask'])
            if self._check_format and mask_filename is not "":
              mask = pickle.load(open(mask_filename))
            else:
              mask = None
          else:
            mask_filename = None
            mask = None
          if 'gain' in imageset and imageset['gain'] is not None:
            gain_filename = load_path(imageset['gain'])
            if self._check_format and gain_filename is not "":
              gain = pickle.load(open(gain_filename))
            else:
              gain = None
          else:
            gain_filename = None
            gain = None
          if 'pedestal' in imageset and imageset['pedestal'] is not None:
            pedestal_filename = load_path(imageset['pedestal'])
            if self._check_format and pedestal_filename is not "":
              pedestal = pickle.load(open(pedestal_filename))
            else:
              pedestal = None
          else:
            pedestal_filename = None
            pedestal = None
          if 'dx' in imageset and imageset['dx'] is not None:
            dx_filename = load_path(imageset['dx'])
            if dx_filename is not "":
              dx = pickle.load(open(dx_filename))
            else:
              dx = None
          else:
            dx_filename = None
            dx = None
          if 'dy' in imageset and imageset['dy'] is not None:
            dy_filename = load_path(imageset['dy'])
            if dy_filename is not "":
              dy = pickle.load(open(dy_filename))
            else:
              dy = None
          else:
            dy_filename = None
            dy = None
          if imageset['__id__'] == 'ImageSet':
            imageset = self._make_stills(imageset, format_kwargs=format_kwargs)
          elif imageset['__id__'] == 'ImageGrid':
            imageset = self._make_grid(imageset, format_kwargs=format_kwargs)
          elif imageset['__id__'] == 'ImageSweep':
            imageset = self._make_sweep(
              imageset,
              beam=beam,
              detector=detector,
              goniometer=goniometer,
              scan=scan,
              format_kwargs=format_kwargs)
          elif imageset['__id__'] == 'MemImageSet':
            imageset = self._make_mem_imageset(imageset)
          else:
            raise RuntimeError('Unknown imageset type')

          # Set the external lookup
          if imageset is not None:
            if mask_filename is None:
              mask_filename = ""
            if gain_filename is None:
              gain_filename = ""
            if pedestal_filename is None:
              pedestal_filename = ""
            if dx_filename is None:
              dx_filename = ""
            if dy_filename is None:
              dy_filename = ""
            if mask is None:
              mask = ImageBool()
            else:
              mask = ImageBool(mask)
            if gain is None:
              gain = ImageDouble()
            else:
              gain = ImageDouble(gain)
            if pedestal is None:
              pedestal = ImageDouble()
            else:
              pedestal = ImageDouble(pedestal)
            if dx is None:
              dx = ImageDouble()
            else:
              dx = ImageDouble(dx)
            if dy is None:
              dy = ImageDouble()
            else:
              dy = ImageDouble(dy)
            imageset.external_lookup.mask.data = mask
            imageset.external_lookup.mask.filename = mask_filename
            imageset.external_lookup.gain.data = gain
            imageset.external_lookup.gain.filename = gain_filename
            imageset.external_lookup.pedestal.data = pedestal
            imageset.external_lookup.pedestal.filename = pedestal_filename
            imageset.external_lookup.dx.data = dx
            imageset.external_lookup.dx.filename = dx_filename
            imageset.external_lookup.dy.data = dy
            imageset.external_lookup.dy.filename = dy_filename

          # Update the imageset models
          if isinstance(imageset, ImageSweep):
            imageset.set_beam(beam)
            imageset.set_detector(detector)
            imageset.set_goniometer(goniometer)
            imageset.set_scan(scan)
          elif isinstance(imageset, ImageSet):
            for i in range(len(imageset)):
              imageset.set_beam(beam, i)
              imageset.set_detector(detector, i)
              imageset.set_goniometer(goniometer, i)
              imageset.set_scan(scan, i)
          elif isinstance(imageset, ImageGrid):
            for i in range(len(imageset)):
              imageset.set_beam(beam, i)
              imageset.set_detector(detector, i)
              imageset.set_goniometer(goniometer, i)
              imageset.set_scan(scan, i)
          else:
            pass

          if imageset is not None:
            imageset.update_detector_px_mm_data()

        # Add the imageset to the dict
        imagesets[key] = imageset

      # Append the experiment
      el.append(Experiment(
        imageset=imageset,
        beam=beam,
        detector=detector,
        goniometer=goniometer,
        scan=scan,
        crystal=crystal,
        profile=profile,
        scaling_model=scaling_model))

    # Return the experiment list
    return el

  def _make_mem_imageset(self, imageset):
    ''' Can't make a mem imageset from dict. '''
    return None

  def _make_stills(self, imageset, format_kwargs=None):
    ''' Make a still imageset. '''
    from dxtbx.imageset import ImageSetFactory
    from dxtbx.serialize.filename import load_path
    filenames = [load_path(p) for p in imageset['images']]
    indices = None
    if "single_file_indices" in imageset:
      indices = imageset['single_file_indices']
      assert len(indices) == len(filenames)
    return ImageSetFactory.make_imageset(
      filenames,
      None,
      check_format=self._check_format,
      single_file_indices=indices,
      format_kwargs=format_kwargs)

  def _make_grid(self, imageset, format_kwargs=None):
    ''' Make a still imageset. '''
    from dxtbx.imageset import ImageGrid
    grid_size = imageset['grid_size']
    return ImageGrid.from_imageset(self._make_stills(imageset, format_kwargs=format_kwargs), grid_size)

  def _make_sweep(self,
                  imageset,
                  beam=None,
                  detector=None,
                  goniometer=None,
                  scan=None,
                  format_kwargs=None):
    ''' Make an image sweep. '''
    from dxtbx.sweep_filenames import template_image_range
    from dxtbx.imageset import ImageSetFactory
    from dxtbx.serialize.filename import load_path
    from dxtbx.format.FormatMultiImage import FormatMultiImage

    # Get the template format
    template = load_path(imageset['template'])

    # Get the number of images (if no scan is given we'll try
    # to find all the images matching the template
    if scan is None:
      i0, i1 = template_image_range(template)
    else:
      i0, i1 = scan.get_image_range()

    format_class = None
    if self._check_format is False:
      if "single_file_indices" in imageset:
        format_class = FormatMultiImage

    # Make a sweep from the input data
    return ImageSetFactory.make_sweep(
      template,
      list(range(i0, i1+1)),
      format_class = format_class,
      check_format=self._check_format,
      beam=beam,
      detector=detector,
      goniometer=goniometer,
      scan=scan,
      format_kwargs=format_kwargs)

  @staticmethod
  def model_or_none(mlist, eobj, name):
    ''' Get a model or None. '''
    index = eobj.get(name)
    if index is not None:
      return mlist[index]
    return None

  @staticmethod
  def _beam_from_dict(obj):
    ''' Get a beam from a dictionary. '''
    from dxtbx.model import BeamFactory
    return BeamFactory.from_dict(obj)

  @staticmethod
  def _detector_from_dict(obj):
    ''' Get the detector from a dictionary. '''
    from dxtbx.model import DetectorFactory
    return DetectorFactory.from_dict(obj)

  @staticmethod
  def _goniometer_from_dict(obj):
    ''' Get the goniometer from a dictionary. '''
    from dxtbx.model import GoniometerFactory
    return GoniometerFactory.from_dict(obj)

  @staticmethod
  def _scan_from_dict(obj):
    ''' Get the scan from a dictionary. '''
    from dxtbx.model import ScanFactory
    return ScanFactory.from_dict(obj)

  @staticmethod
  def _crystal_from_dict(obj):
    ''' Get the crystal from a dictionary. '''
    from dxtbx.model import CrystalFactory
    return CrystalFactory.from_dict(obj)

  @staticmethod
  def _profile_from_dict(obj):
    ''' Get the profile from a dictionary. '''
    from dxtbx.model import ProfileModelFactory
    return ProfileModelFactory.from_dict(obj)

  @staticmethod
  def _scaling_model_from_dict(obj):
    ''' Get the scaling model from a dictionary. '''
    for entry_point in pkg_resources.iter_entry_points(
      'dxtbx.scaling_model_ext'):
      if entry_point.name == obj['__id__']:
        return entry_point.load().from_dict(obj)

  @staticmethod
  def _from_file(filename):
    ''' Load a model dictionary from a file. '''
    from dxtbx.serialize.load import _decode_dict
    from dxtbx.serialize.filename import load_path, temp_chdir
    import json
    from os.path import dirname
    filename = load_path(filename)
    try:
      with temp_chdir(dirname(filename)):
        with open(filename, 'r') as infile:
          return json.load(infile, object_hook=_decode_dict)
    except IOError:
      raise IOError('unable to read file, %s' % filename)


class ExperimentListDumper(object):
  ''' A class to help writing JSON files. '''

  def __init__(self, experiment_list):
    ''' Initialise '''
    self._experiment_list = experiment_list

  def as_json(self, filename=None, compact=False, split=False):
    ''' Dump experiment list as json '''
    import json
    from os.path import splitext
    from libtbx.containers import OrderedDict

    # Get the dictionary and get the JSON string
    dictionary = self._experiment_list.to_dict()

    # Split into separate files
    if filename is not None and split:

      # Get lists of models by filename
      basepath = splitext(filename)[0]
      ilist = [('%s_imageset_%d.json' % (basepath, i), d)
                  for i, d in enumerate(dictionary['imageset'])]
      blist = [('%s_beam_%d.json' % (basepath, i), d)
                  for i, d in enumerate(dictionary['beam'])]
      dlist = [('%s_detector_%d.json' % (basepath, i), d)
                  for i, d in enumerate(dictionary['detector'])]
      glist = [('%s_goniometer_%d.json' % (basepath, i), d)
                  for i, d in enumerate(dictionary['goniometer'])]
      slist = [('%s_scan_%d.json' % (basepath, i), d)
                  for i, d in enumerate(dictionary['scan'])]
      clist = [('%s_crystal_%d.json' % (basepath, i), d)
                  for i, d in enumerate(dictionary['crystal'])]
      plist = [('%s_profile_%d.json' % (basepath, i), d)
                  for i, d in enumerate(dictionary['profile'])]
      scalelist = [('%s_scaling_model_%d.json' % (basepath, i), d)
                  for i, d in enumerate(dictionary['scaling_model'])]

      # Get the list of experiments
      edict = OrderedDict([
        ( '__id__', 'ExperimentList' ),
        ( 'experiment', dictionary['experiment'] )
      ])

      # Set paths rather than indices
      for e in edict['experiment']:
        if 'imageset' in e:
          e['imageset'] = ilist[e['imageset']][0]
        if 'beam' in e:
          e['beam'] = blist[e['beam']][0]
        if 'detector' in e:
          e['detector'] = dlist[e['detector']][0]
        if 'goniometer' in e:
          e['goniometer'] = glist[e['goniometer']][0]
        if 'scan' in e:
          e['scan'] = slist[e['scan']][0]
        if 'crystal' in e:
          e['crystal'] = clist[e['crystal']][0]
        if 'profile' in e:
          e['profile'] = plist[e['profile']][0]
        if 'scaling_model' in e:
          e['scaling_model'] = scalelist[e['scaling_model']][0]

      to_write = ilist + blist + dlist + glist + \
                 slist + clist + plist + scalelist + [(filename, edict)]
    else:
      to_write = [(filename, dictionary)]

    for fname, obj  in to_write:
      if compact:
        text = json.dumps(obj, separators=(',',':'), ensure_ascii=True)
      else:
        text = json.dumps(obj, indent=2, ensure_ascii=True)

      # If a filename is set then dump to file otherwise return string
      if fname is not None:
        with open(fname, 'w') as outfile:
          outfile.write(text)
      else:
        return text

  def as_pickle(self, filename=None, **kwargs):
    ''' Dump experiment list as pickle. '''
    import cPickle as pickle

    # Get the pickle string
    text = pickle.dumps(self._experiment_list,
      protocol=pickle.HIGHEST_PROTOCOL)

    # Write the file
    if filename is not None:
      with open(filename, 'wb') as outfile:
        outfile.write(text)
    else:
      return text

  def as_file(self, filename, **kwargs):
    ''' Dump experiment list as file. '''
    from os.path import splitext
    ext = splitext(filename)[1]
    j_ext = ['.json']
    p_ext = ['.p', '.pkl', '.pickle']
    if ext.lower() in j_ext:
      return self.as_json(filename, **kwargs)
    elif ext.lower() in p_ext:
      return self.as_pickle(filename, **kwargs)
    else:
      ext_str = '|'.join(j_ext + p_ext)
      raise RuntimeError('expected extension {%s}, got %s' % (ext_str, ext))


class ExperimentListFactory(object):
  ''' A class to help instantiate experiment lists. '''

  @staticmethod
  def from_args(args, verbose=False, unhandled=None):
    ''' Try to load experiment from any recognised format. '''
    from dxtbx.datablock import DataBlockFactory

    # Create a list for unhandled arguments
    if unhandled is None:
      unhandled = []

    experiments = ExperimentList()
    ## First try as image files
    #experiments = ExperimentListFactory.from_datablock(
      #DataBlockFactory.from_args(args, verbose, unhandled1))

    # Try to load from serialized formats
    for filename in args:
      try:
        experiments.extend(ExperimentListFactory.from_serialized_format(filename))
        if verbose: print 'Loaded experiments from %s' % filename
      except Exception:
        unhandled.append(filename)

    # Return the experiments
    return experiments

  @staticmethod
  def from_imageset_and_crystal(imageset, crystal):
    ''' Load an experiment list from an imageset and crystal. '''
    from dxtbx.imageset import ImageSweep
    if isinstance(imageset, ImageSweep):
      return ExperimentListFactory.from_sweep_and_crystal(imageset, crystal)
    else:
      return ExperimentListFactory.from_stills_and_crystal(imageset, crystal)

  @staticmethod
  def from_sweep_and_crystal(imageset, crystal):
    ''' Create an experiment list from sweep and crystal. '''
    return ExperimentList([
      Experiment(
        imageset=imageset,
        beam=imageset.get_beam(),
        detector=imageset.get_detector(),
        goniometer=imageset.get_goniometer(),
        scan=imageset.get_scan(),
        crystal=crystal)])

  @staticmethod
  def from_stills_and_crystal(imageset, crystal):
    ''' Create an experiment list from stills and crystal. '''
    from itertools import groupby

    # Get a list of models for each image
    beam, detector, gonio, scan = ([], [], [], [])
    for i in xrange(len(imageset)):
      try:
        beam.append(imageset.get_beam(i))
      except Exception:
        beam.append(None)
      try:
        detector.append(imageset.get_detector(i))
      except Exception:
        detector.append(None)
      try:
        gonio.append(imageset.get_goniometer(i))
      except Exception:
        gonio.append(None)
      try:
        scan.append(imageset.get_scan(i))
      except Exception:
        scan.append(None)
    models = zip(beam, detector, gonio, scan)

    # Find subsets where all the models are the same
    experiments = ExperimentList()
    for m, indices in groupby(xrange(len(models)), lambda i: models[i]):
      indices = list(indices)
      experiments.append(Experiment(
        imageset=imageset[indices[0]: indices[-1]+1],
        beam=m[0], detector=m[1],
        goniometer=m[2], scan=m[3], crystal=crystal))

    # Return experiments
    return experiments

  @staticmethod
  def from_datablock_and_crystal(datablock, crystal):
    ''' Load an experiment list from a datablock. '''

    # Initialise the experiment list
    experiments = ExperimentList()

    # If we have a list, loop through
    if isinstance(datablock, list):
      for db in datablock:
        experiments.extend(ExperimentListFactory.from_datablock_and_crystal(
          db, crystal))
      return experiments

    # Add all the imagesets
    for imageset in datablock.extract_imagesets():
      experiments.extend(ExperimentListFactory.from_imageset_and_crystal(
        imageset, crystal))

    # Check the list is consistent
    assert(experiments.is_consistent())

    # Return the experiments
    return experiments

  @staticmethod
  def from_dict(obj, check_format=True):
    ''' Load an experiment list from a dictionary. '''

    # Decode the experiments from the dictionary
    experiments = ExperimentListDict(obj, check_format).decode()

    # Check the list is consistent
    assert(experiments.is_consistent())

    # Return the experiments
    return experiments

  @staticmethod
  def from_json(text, check_format=True):
    ''' Load an experiment list from JSON. '''
    from dxtbx.serialize.load import _decode_dict
    import json
    return ExperimentListFactory.from_dict(
      json.loads(text, object_hook=_decode_dict),
      check_format)

  @staticmethod
  def from_json_file(filename, check_format=True):
    ''' Load an experiment list from a json file. '''
    from dxtbx.serialize.filename import temp_chdir
    from os.path import dirname, abspath
    filename = abspath(filename)
    with temp_chdir(dirname(filename)):
      with open(filename, 'r') as infile:
        return ExperimentListFactory.from_json(infile.read(), check_format)

  @staticmethod
  def from_pickle_file(filename):
    ''' Decode an experiment list from a pickle file. '''
    import cPickle as pickle
    with open(filename, 'rb') as infile:
      obj = pickle.load(infile)
      assert(isinstance(obj, ExperimentList))
      return obj

  @staticmethod
  def from_xds(xds_inp, xds_other):
    ''' Generate an experiment list from XDS files. '''
    from dxtbx.serialize import xds
    from dxtbx.datablock import DataBlockFactory

    # Get the sweep from the XDS files
    sweep = xds.to_imageset(xds_inp, xds_other)

    # Get the crystal from the XDS files
    crystal = xds.to_crystal(xds_other)

    # Create the experiment list
    experiments = ExperimentListFactory.from_imageset_and_crystal(
      sweep, crystal)

    # Set the crystal in the experiment list
    assert(len(experiments) == 1)

    # Return the experiment list
    return experiments

  @staticmethod
  def from_serialized_format(filename, check_format=True):
    ''' Try to load the experiment list from a serialized format. '''

    # First try as a JSON file
    try:
      return ExperimentListFactory.from_json_file(filename, check_format)
    except Exception:
      pass

    # Now try as a pickle file
    return ExperimentListFactory.from_pickle_file(filename)
