from __future__ import absolute_import, division



class Reader(object):

  _format_class_ = None

  def __init__(self, filenames, **kwargs):
    self.kwargs = kwargs
    self.format_class = Reader._format_class_
    assert len(filenames) == 1
    self._filename = filenames[0]

  def nullify_format_instance(self):
    self.format_class._current_instance_ = None
    self.format_class._current_filename_ = None

  def read(self, index):
    format_instance = self.format_class.get_instance(self._filename, **self.kwargs)
    return format_instance.get_raw_data(index)

  def paths(self):
    return [self._filename]

  def num_images(self):
    format_instance = self.format_class.get_instance(self._filename, **self.kwargs)
    return format_instance.get_num_images()

  def __len__(self):
    return self.num_images()

  def copy(self, filenames, indices=None):
    return Reader(filenames, indices)

  def identifiers(self):
    return ["%s-%d" % (self._filename, index) for index in range(len(self))]

  def is_single_file_reader(self):
    return True

  def master_path(self):
    return self._filename


class Masker(object):

  _format_class_ = None

  def __init__(self, filenames, **kwargs):
    self.kwargs = kwargs
    self.format_class = Masker._format_class_
    assert len(filenames) == 1
    self._filename = filenames[0]

  def get(self, index, goniometer=None):
    format_instance = self.format_class.get_instance(self._filename, **self.kwargs)
    return format_instance.get_mask(index, goniometer)

  def paths(self):
    return [self._filename]

  def num_images(self):
    format_instance = self.format_class.get_instance(self._filename, **self.kwargs)
    return format_instance.get_num_images()

  def __len__(self):
    return self.num_images()

  def copy(self, filenames):
    return Masker(filenames)


class FormatMultiImage(object):

  def __init__(self, **kwargs):
    pass

  def get_num_images(self):
    raise RuntimeError('Overload!')

  def get_goniometer(self, index=None):
    return self._goniometer_instance

  def get_detector(self, index=None):
    return self._detector_instance

  def get_beam(self, index=None):
    return self._beam_instance

  def get_scan(self, index=None):
    return self._scan_instance

  def get_raw_data(self, index=None):
    raise RuntimeError('Overload!')

  def get_mask(self, index=None, goniometer=None):
    return None

  def get_detectorbase(self, index=None):
    raise RuntimeError('Overload!')

  def get_image_file(self, index=None):
    raise RuntimeError('Overload!')

  @classmethod
  def get_reader(Class):
    '''
    Return a reader class

    '''
    obj = Reader
    obj._format_class_ = Class
    return obj

  @classmethod
  def get_masker(Class):
    '''
    Return a reader class

    '''
    obj = Masker
    obj._format_class_ = Class
    return obj

  @classmethod
  def get_imageset(Class,
                   filenames,
                   beam=None,
                   detector=None,
                   goniometer=None,
                   scan=None,
                   as_sweep=False,
                   as_imageset=False,
                   single_file_indices=None,
                   format_kwargs=None):
    '''
    Factory method to create an imageset

    '''
    from dxtbx.imageset import ImageSetData
    from dxtbx.imageset import ImageSet
    from dxtbx.imageset import ImageSweep
    from os.path import abspath

    if isinstance(filenames, str):
      filenames = [filenames]

    # Make filenames absolute
    filenames = map(abspath, filenames)

    # Make it a dictionary
    if format_kwargs is None:
      format_kwargs = {}

    # Get some information from the format class
    reader = Class.get_reader()(filenames, **format_kwargs)
    masker = Class.get_masker()(filenames, **format_kwargs)

    # Get the format instance
    assert len(filenames) == 1
    format_instance = Class(filenames[0], **format_kwargs)

    # Read the vendor type
    vendor = format_instance.get_vendortype()
    params = format_kwargs

    # Check if we have a sweep

    # Make sure only 1 or none is set
    assert [as_imageset, as_sweep].count(True) < 2
    if as_imageset:
      is_sweep = False
    elif as_sweep:
      is_sweep = True
    else:
      if scan is None and format_instance is None:
        raise RuntimeError('''
          One of the following needs to be set
            - as_imageset=True
            - as_sweep=True
            - scan
            - check_format=True
      ''')
      if scan is None:
        test_scan = format_instance.get_scan()
      else:
        test_scan = scan
      if test_scan is not None and test_scan.get_oscillation()[1] != 0:
        is_sweep = True
      else:
        is_sweep = False

    # Create an imageset or sweep
    if not is_sweep:

      # Create the imageset
      iset = ImageSet(
        ImageSetData(
          reader = reader,
          masker = masker,
          properties = {
            "vendor" : vendor,
            "params" : params,
            "format" : Class
          }),
        indices=single_file_indices)

      # If any are None then read from format
      if [beam, detector, goniometer, scan].count(None) != 0:

        # Get list of models
        beam = []
        detector = []
        goniometer = []
        scan = []
        for i in range(format_instance.get_num_images()):
          beam.append(format_instance.get_beam(i))
          detector.append(format_instance.get_detector(i))
          goniometer.append(format_instance.get_goniometer(i))
          scan.append(format_instance.get_scan(i))

      if single_file_indices is None:
        single_file_indices = list(range(format_instance.get_num_images()))

      # Set the list of models
      for i in range(len(single_file_indices)):
        iset.set_beam(beam[single_file_indices[i]], i)
        iset.set_detector(detector[single_file_indices[i]], i)
        iset.set_goniometer(goniometer[single_file_indices[i]], i)
        iset.set_scan(scan[single_file_indices[i]], i)

    else:

      # Get the template
      template = filenames[0]

      # Check indices are sequential
      if single_file_indices is not None:
        assert all(i == j for i, j in zip(
          single_file_indices[:-1],
          single_file_indices[1:]))
        num_images = len(single_file_indices)
      else:
        num_images = format_instance.get_num_images()

      # Check the scan makes sense
      if scan is not None:
        assert scan.get_num_images() == num_images

      # If any are None then read from format
      if beam is None:
        beam = format_instance.get_beam()
      if detector is None:
        detector = format_instance.get_detector()
      if goniometer is None:
        goniometer = format_instance.get_goniometer()
      if scan is None:
        scan = format_instance.get_scan()
        if scan is not None:
          for f in filenames[1:]:
            format_instance = Class(f, **format_kwargs)
            scan += format_instance.get_scan()

      # Create the sweep
      iset = ImageSweep(
        ImageSetData(
          reader     = reader,
          masker     = masker,
          properties = {
            "vendor"   : vendor,
            "params"   : params,
            "format"   : Class,
            "template" : template,
          }),
        beam       = beam,
        detector   = detector,
        goniometer = goniometer,
        scan       = scan,
        indices=single_file_indices)

    # Return the imageset
    return iset
