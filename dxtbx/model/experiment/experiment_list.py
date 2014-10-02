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
from __future__ import division

class Experiment(object):
  ''' A class to represent what's in an experiment.

  Contains:
    - imageset Access to the image data
    - beam The beam model
    - detector The detector model
    - goniometer The goniometer model
    - scan The scan model
    - crystal The crystal model

  Some of these may be set to "None"

  '''

  def __init__(self, imageset=None, beam=None, detector=None,
               goniometer=None, scan=None, crystal=None):
    ''' Initialise the experiment with the given models. '''
    self.imageset = imageset
    self.beam = beam
    self.detector = detector
    self.goniometer = goniometer
    self.scan = scan
    self.crystal = crystal

  def __contains__(self, item):
    ''' Check if the experiment contains the model. '''
    return (item is self.imageset or
            item is self.beam or
            item is self.detector or
            item is self.goniometer or
            item is self.scan or
            item is self.crystal)

  def __eq__(self, other):
    ''' Check if an experiment is the same as another. '''
    if not isinstance(other, Experiment):
      return False
    return (self.imageset is other.imageset and
            self.beam is other.beam and
            self.detector is other.detector and
            self.goniometer is other.goniometer and
            self.scan is other.scan and
            self.crystal is other.crystal)

  def __ne__(self, other):
    ''' Check if an experiment not equal to another. '''
    return not self.__eq__(other)

  def assert_is_consistent(self):
    ''' If a scan is present, check that it makes sense with the imageset. '''
    return True # XXX FIXME JAMES! 2014_02_04
    #from dxtbx.imageset import ImageSweep
    #if self.scan:
      #if isinstance(self.imageset, ImageSweep):
        #assert(len(self.imageset) == self.scan.get_num_images())
        #assert(self.imageset.get_array_range() == self.scan.get_array_range())
      #elif self.imageset is not None:
        #assert((self.scan.get_num_images() == 1 and
                #self.scan.get_oscillation()[1] == 0.0))
        #assert(len(self.imageset.indices()) == 1)

  def is_consistent(self):
    ''' If a scan is present, check that it makes sense with the imageset. '''
    try:
      self.assert_is_consistent()
      return True
    except Exception:
      return False

class ExperimentList(object):
  ''' The experiment list class. This class is used to manage all the
  experiments and contains methods to get groups of models etc. '''

  def __init__(self, item=None):
    ''' Initialise the list. '''
    if item is not None:
      self._data = list(item)
    else:
      self._data = list()

  def __setitem__(self, index, item):
    ''' Set an experiment. '''
    if isinstance(item, Experiment):
      self._data[index] = item
    else:
      raise TypeError('expected type Experiment, got %s' % type(item))

  def __getitem__(self, index):
    ''' Get an experiment. '''
    if isinstance(index, slice):
      return ExperimentList(self._data[index])
    return self._data[index]

  def __delitem__(self, index):
    ''' Delete an experiment. '''
    del self._data[index]

  def __len__(self):
    ''' Get the number of experiments. '''
    return len(self._data)

  def __iter__(self):
    ''' Iterate through the experiments. '''
    for e in self._data:
      yield e

  def __contains__(self, item):
    ''' Check if an item is contained in the list of experiments.
    Also checks to see if a model is contained in an experiment. '''
    return item in self._data or any(item in e for e in self._data)

  def index(self, item):
    ''' Get the index of an experiment. '''
    return self._data.index(item)

  def append(self, item):
    ''' Add a new experiment to the list. '''
    if isinstance(item, Experiment):
      self._data.append(item)
    else:
      raise TypeError('expected type Experiment, got %s' % type(item))

  def extend(self, other):
    ''' Add another experiment list to this one. '''
    if isinstance(other, ExperimentList):
      self._data.extend(other._data)
    else:
      raise TypeError('expected type ExperimentList, got %s' % type(item))

  def replace(self, a, b):
    ''' Replace all occurances of a with b. '''
    for i in self.indices(a):
      exp = self._data[i]
      if   exp.imageset is a:   exp.imageset = b
      elif exp.beam is a:       exp.beam = b
      elif exp.detector is a:   exp.detector = b
      elif exp.goniometer is a: exp.goniometer = b
      elif exp.scan is a:       exp.scan = b
      elif exp.crystal is a:    exp.crystal = b
      else: raise ValueError('unidentified model %s' % a)

  def remove(self, model):
    ''' Remove all occurances of the model. '''
    self.replace(model, None)

  def indices(self, model):
    ''' Get the indices of the experiments which contains the model. '''
    if isinstance(model, list) or isinstance(model, tuple):
      return list(set.intersection(*[set(self.indices(m)) for m in model]))
    else:
      return [i for i, e in enumerate(self) if model in e]

  def beams(self):
    ''' Get a list of the unique beams (includes None). '''
    from libtbx.containers import OrderedDict
    return OrderedDict([(e.beam, None) for e in self]).keys()

  def detectors(self):
    ''' Get a list of the unique detectors (includes None). '''
    from libtbx.containers import OrderedDict
    return OrderedDict([(e.detector, None) for e in self]).keys()

  def goniometers(self):
    ''' Get a list of the unique goniometers (includes None). '''
    from libtbx.containers import OrderedDict
    return OrderedDict([(e.goniometer, None) for e in self]).keys()

  def scans(self):
    ''' Get a list of the unique scans (includes None). '''
    from libtbx.containers import OrderedDict
    return OrderedDict([(e.scan, None) for e in self]).keys()

  def crystals(self):
    ''' Get a list of the unique crystals (includes None). '''
    from libtbx.containers import OrderedDict
    return OrderedDict([(e.crystal, None) for e in self]).keys()

  def imagesets(self):
    ''' Get a list of the unique imagesets (includes None).

    This returns unique complete sets rather than partial.
    '''
    from libtbx.containers import OrderedDict
    return OrderedDict([(e.imageset, None) for e in self
                        if e.imageset is not None]).keys()
    # temp = OrderedDict([(e.imageset.reader(), i)
    #   for i, e in enumerate(self) if e.imageset is not None])
    # return OrderedDict([(self[i].imageset.complete_set(), None)
    #   for i in temp.itervalues()]).keys()

  def is_consistent(self):
    ''' Check all the models are consistent. '''
    return all([e.is_consistent() for e in self])

  def all_stills(self):
    ''' Check if all the experiments are stills '''
    assert(len(self) > 0)
    result = True
    for e in self:
      if e.goniometer is not None and e.scan is not None:
        result = False
        break
    return result

  def all_sweeps(self):
    ''' Check if all the experiments are from sweeps '''
    assert(len(self) > 0)
    result = True
    for e in self:
      if e.goniometer is None or e.scan is None:
        result = False
        break
    return result

  def to_dict(self):
    ''' Serialize the experiment list to dictionary. '''
    from libtbx.containers import OrderedDict
    from dxtbx.imageset import ImageSet, ImageSweep, MemImageSet
    from dxtbx.serialize.crystal import to_dict as crystal_to_dict
    from dxtbx.format.FormatMultiImage import FormatMultiImage

    # Check the experiment list is consistent
    assert(self.is_consistent())

    # Get the list of unique models
    blist = self.beams()
    dlist = self.detectors()
    glist = self.goniometers()
    slist = self.scans()
    clist = self.crystals()
    ilist = self.imagesets()

    # Create the output dictionary
    result = OrderedDict()
    result['__id__'] = 'ExperimentList'
    result['experiment'] = []

    # Function to find in list by reference
    def find_index(l, m):
      for i, mm in enumerate(l):
        if mm is m:
          return i
      return -1

    # Add the experiments to the dictionary
    for e in self:
      obj = OrderedDict()
      obj['__id__'] = 'Experiment'
      if e.beam is not None:
        obj['beam'] = find_index(blist, e.beam)
      if e.detector is not None:
        obj['detector'] = find_index(dlist, e.detector)
      if e.goniometer is not None:
        obj['goniometer'] = find_index(glist, e.goniometer)
      if e.scan is not None:
        obj['scan'] = find_index(slist, e.scan)
      if e.crystal is not None:
        obj['crystal'] = find_index(clist, e.crystal)
      if e.imageset is not None:
        obj['imageset'] = ilist.index(e.imageset)
        if e.scan is None and not isinstance(e.imageset, ImageSweep):
          if len(e.imageset) != len(e.imageset.complete_set()):
            obj['imageset'] = (obj['imageset'], e.imageset.indices())
      result['experiment'].append(obj)

    def get_template(imset):
      if imset.reader().get_format_class() is None:
        try:
          return imset.get_template()
        except Exception:
          pass
        try:
          return imset.reader().get_path()
        except Exception:
          pass
        return 'none'
      elif issubclass(imset.reader().get_format_class(), FormatMultiImage):
        return imset.reader().get_path()
      else:
        return imset.get_template()

    # Serialize all the imagesets
    result['imageset'] = []
    for imset in ilist:
      if isinstance(imset, ImageSweep):
        # FIXME_HACK
        template = get_template(imset)
        result['imageset'].append(OrderedDict([
          ('__id__', 'ImageSweep'),
          ('template', template)]))
      elif isinstance(imset, MemImageSet):
        result['imageset'].append(OrderedDict([
          ('__id__', 'MemImageSet')]))
      elif isinstance(imset, ImageSet):
        result['imageset'].append(OrderedDict([
          ('__id__', 'ImageSet'),
          ('images', imset.paths())]))
      else:
        raise TypeError('expected ImageSet or ImageSweep, got %s' % type(imset))

    # Extract all the model dictionaries
    result['beam']       = [b.to_dict() for b in blist if b is not None]
    result['detector']   = [d.to_dict() for d in dlist if d is not None]
    result['goniometer'] = [g.to_dict() for g in glist if g is not None]
    result['scan']       = [s.to_dict() for s in slist if s is not None]
    result['crystal']    = [crystal_to_dict(c) for c in clist if c is not None]

    # Return the dictionary
    return result

  def to_datablocks(self):
    ''' Return the experiment list as a datablock list.
    This assumes that the experiment contains 1 datablock.'''
    from dxtbx.datablock import DataBlockFactory

    # Convert the experiment list to dict
    obj = self.to_dict()

    # Convert the dictionary to a datablock dictionary
    obj['__id__'] = 'DataBlock'
    for e in obj['experiment']:
      iid = e['imageset']
      imageset = obj['imageset'][iid]
      if 'beam' in e:
        imageset['beam'] = e['beam']
      if 'detector' in e:
        imageset['detector'] = e['detector']
      if 'goniometer' in e:
        imageset['goniometer'] = e['goniometer']
      if 'scan' in e:
        imageset['scan'] = e['scan']

    # Remove the experiments
    del obj['experiment']

    # Create the datablock
    return DataBlockFactory.from_dict([obj])


class ExperimentListDict(object):
  ''' A helper class for serializing the experiment list to dictionary (needed
  to save the experiment list to JSON format. '''

  def __init__(self, obj, check_format=True):
    ''' Initialise. Copy the dictionary. '''
    from copy import deepcopy
    self._obj = deepcopy(obj)
    self._check_format = check_format

  def decode(self):
    ''' Decode the dictionary into a list of experiments. '''

    # Extract lists of models referenced by experiments
    self._blist = self._extract_models('beam')
    self._dlist = self._extract_models('detector')
    self._glist = self._extract_models('goniometer')
    self._slist = self._extract_models('scan')
    self._clist = self._extract_models('crystal')

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
      value = eobj.get(name, None)
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
      value = eobj.get('imageset', None)
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
    from dxtbx.imageset import ImageSweep

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
      key = (eobj.get('imageset', None), eobj.get('scan',None))
      try:
        imageset = imagesets[key]
      except Exception:
        imageset = ExperimentListDict.model_or_none(self._ilist, eobj, 'imageset')

        # Create the imageset from the input data
        if imageset is not None:
          if imageset['__id__'] == 'ImageSet':
            imageset = self._make_stills(imageset)
          elif imageset['__id__'] == 'ImageSweep':
            imageset = self._make_sweep(imageset, scan)
          elif imageset['__id__'] == 'MemImageSet':
            imageset = self._make_mem_imageset(imageset)
          else:
            raise RuntimeError('Unknown imageset type')

          # Fill in any models if they aren't already there
          if beam is None:
            beam = imageset.get_beam()
          if detector is None:
            detector = imageset.get_detector()
          if goniometer is None:
            goniometer = imageset.get_goniometer()
          if scan is None:
            scan = imageset.get_scan()

          # Update the imageset models
          if isinstance(imageset, ImageSweep):
            imageset.set_beam(beam)
            imageset.set_detector(detector)
            imageset.set_goniometer(goniometer)
            imageset.set_scan(scan)
          else:
            for i in range(len(imageset)):
              imageset.set_beam(beam, i)
              imageset.set_detector(detector, i)
              imageset.set_goniometer(goniometer, i)
              imageset.set_scan(scan, i)

        # Add the imageset to the dict
        imagesets[key] = imageset

      # Append the experiment
      el.append(Experiment(
        imageset=imageset,
        beam=beam,
        detector=detector,
        goniometer=goniometer,
        scan=scan,
        crystal=crystal))

    # Return the experiment list
    return el

  def _make_mem_imageset(self, imageset):
    ''' Can't make a mem imageset from dict. '''
    return None

  def _make_stills(self, imageset):
    ''' Make a still imageset. '''
    from dxtbx.imageset import ImageSetFactory
    from dxtbx.serialize.filename import load_path
    filenames = [load_path(p) for p in imageset['images']]
    return ImageSetFactory.make_imageset(
      filenames, None, check_format=self._check_format)

  def _make_sweep(self, imageset, scan):
    ''' Make an image sweep. '''
    from dxtbx.sweep_filenames import template_image_range
    from dxtbx.imageset import ImageSetFactory
    from dxtbx.serialize.filename import load_path

    # Get the template format
    template = load_path(imageset['template'])

    # Get the number of images (if no scan is given we'll try
    # to find all the images matching the template
    if scan is None:
      i0, i1 = template_image_range(template)
    else:
      i0, i1 = scan.get_image_range()

    # Make a sweep from the input data
    return ImageSetFactory.make_sweep(template,
      list(range(i0, i1+1)), None, check_format=self._check_format)

  @staticmethod
  def model_or_none(mlist, eobj, name):
    ''' Get a model or None. '''
    index = eobj.get(name, None)
    if index is not None:
      return mlist[index]
    return None

  @staticmethod
  def _beam_from_dict(obj):
    ''' Get a beam from a dictionary. '''
    from dxtbx.model import Beam
    return Beam.from_dict(obj)

  @staticmethod
  def _detector_from_dict(obj):
    ''' Get the detector from a dictionary. '''
    from dxtbx.model import Detector, HierarchicalDetector
    if 'hierarchy' in obj:
      return HierarchicalDetector.from_dict(obj)
    else:
      return Detector.from_dict(obj)

  @staticmethod
  def _goniometer_from_dict(obj):
    ''' Get the goniometer from a dictionary. '''
    from dxtbx.model import Goniometer
    return Goniometer.from_dict(obj)

  @staticmethod
  def _scan_from_dict(obj):
    ''' Get the scan from a dictionary. '''
    from dxtbx.model import Scan
    return Scan.from_dict(obj)

  @staticmethod
  def _crystal_from_dict(obj):
    ''' Get the crystal from a dictionary. '''
    from dxtbx.serialize.crystal import from_dict
    return from_dict(obj)

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
          return json.loads(infile.read(), object_hook=_decode_dict)
    except IOError, e:
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

      to_write = ilist + blist + dlist + glist + \
                 slist + clist + [(filename, edict)]
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
    for i in range(len(imageset)):
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
    for m, indices in groupby(range(len(models)), lambda i: models[i]):
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
    except Exception, e:
      pass

    # Now try as a pickle file
    return ExperimentListFactory.from_pickle_file(filename)
