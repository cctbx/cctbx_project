from __future__ import absolute_import, division
from dxtbx.format.FormatMultiImage import FormatMultiImage
from dxtbx.imageset import ImageSet

class ImageSetLazy(ImageSet):

  def get_detector(self, index=None):
    if index is None: index=0
    detector = super(ImageSetLazy,self).get_detector(index)
    if detector is None:
      format_instance = self.get_format_class()._current_instance_
      #format_instance = self.get_format_class()(self.paths()[self.indices()[index]])
      detector = format_instance.get_detector(self.indices()[index])
      self.set_detector(detector,index)
    return detector


  def get_beam(self, index=None):
    if index is None: index=0
    beam = super(ImageSetLazy,self).get_beam(index)
    if beam is None:
      format_instance = self.get_format_class()._current_instance_
      #format_instance = self.get_format_class()(self.paths()[self.indices()[index]])
      beam = format_instance.get_beam(self.indices()[index])
      self.set_beam(beam,index)
    return beam

  def get_goniometer(self, index=None):
    if index is None: index=0
    goniometer = super(ImageSetLazy,self).get_goniometer(index)
    if goniometer is None:
      format_instance = self.get_format_class()._current_instance_
      #format_instance = self.get_format_class()(self.paths()[self.indices()[index]])
      goniometer = format_instance.get_goniometer(self.indices()[index])
      self.set_goniometer(goniometer,index)
    return goniometer

  def get_scan(self, index=None):
    if index is None: index=0
    scan = super(ImageSetLazy,self).get_scan(index)
    if scan is None:
      format_instance = self.get_format_class()._current_instance_
      #format_instance = self.get_format_class()(self.paths()[self.indices()[index]])
      scan = format_instance.get_scan(self.indices()[index])
      self.set_scan(scan,index)
    return scan

  def __getitem__(self, item):
    if isinstance(item, slice):
      return ImageSetLazy(self.data(), indices = self.indices()[item])
    else:
      # Sets the list for detector, beam etc before being accessed by functions in imageset.h
      self.get_detector(item)
      self.get_beam(item)
      self.get_goniometer(item)
      self.get_scan(item)
    return super(ImageSetLazy,self).__getitem__(item)

class FormatMultiImageLazy(FormatMultiImage):

  def get_goniometer(self, index=None):
    return self._goniometer_instance

  def get_detector(self, index=None):
    return self._detector_instance

  def get_beam(self, index=None):
    return self._beam_instance

  def get_scan(self, index=None):
    return self._scan_instance

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
                   format_kwargs=None,
                   template=None,
                   check_format=True):
    '''
    Factory method to create an imageset

    '''
    from dxtbx.imageset import ImageSetData
    from dxtbx.imageset import ImageSweep
    from os.path import abspath
    from dials.array_family import flex
    if isinstance(filenames, str):
      filenames = [filenames]
    elif len(filenames) > 1:
      assert len(set(filenames)) == 1
      filenames = filenames[0:1]

    # Make filenames absolute
    filenames = map(abspath, filenames)

    # Make it a dictionary
    if format_kwargs is None:
      format_kwargs = {}

    # If we have no specific format class, we need indices for number of images
    if Class == FormatMultiImageLazy:
      assert single_file_indices is not None
      assert min(single_file_indices) >= 0
      num_images = max(single_file_indices) + 1
    else:
      num_images = None

    # Get some information from the format class
    reader = Class.get_reader()(filenames, num_images=num_images, **format_kwargs)
    masker = Class.get_masker()(filenames, num_images=num_images, **format_kwargs)

    # Get the format instance
    assert len(filenames) == 1
    if check_format is True:
      format_instance = Class(filenames[0], **format_kwargs)
    else:
      format_instance = None

    # Read the vendor type
    if check_format is True:
      vendor = format_instance.get_vendortype()
    else:
      vendor = ""

    # Get the format kwargs
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

    if single_file_indices is not None:
      single_file_indices = flex.size_t(single_file_indices)

    # Create an imageset or sweep
    if not is_sweep:

      # Create the imageset
      iset = ImageSetLazy(
        ImageSetData(
          reader = reader,
          masker = masker,
          vendor = vendor,
          params = params,
          format = Class),
        indices=single_file_indices)

      # If any are None then read from format
#      if [beam, detector, goniometer, scan].count(None) != 0:
#
#        # Get list of models
#        beam = []
#        detector = []
#        goniometer = []
#        scan = []
#        for i in range(format_instance.get_num_images()):
#          beam.append(format_instance.get_beam(i))
#          detector.append(format_instance.get_detector(i))
#          goniometer.append(format_instance.get_goniometer(i))
#          scan.append(format_instance.get_scan(i))
#
#      if single_file_indices is None:
#        single_file_indices = list(range(format_instance.get_num_images()))
#
#      # Set the list of models
#      for i in range(len(single_file_indices)):
#        iset.set_beam(beam[single_file_indices[i]], i)
#        iset.set_detector(detector[single_file_indices[i]], i)
#        iset.set_goniometer(goniometer[single_file_indices[i]], i)
#        iset.set_scan(scan[single_file_indices[i]], i)

    else:

      # Get the template
      template = filenames[0]

      # Check indices are sequential
      if single_file_indices is not None:
        assert all(i + 1== j for i, j in zip(
          single_file_indices[:-1],
          single_file_indices[1:]))
        num_images = len(single_file_indices)
      else:
        num_images = format_instance.get_num_images()

      # Check the scan makes sense - we must want to use <= total images
      if scan is not None:
        assert scan.get_num_images() <= num_images

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
          # TODO change this to lazy read of models
          for f in filenames[1:]:
            format_instance = Class(f, **format_kwargs)
            scan += format_instance.get_scan()

      isetdata = ImageSetData(
          reader     = reader,
          masker     = masker,
          vendor     = vendor,
          params     = params,
          format     = Class,
          template   = template)

      # Create the sweep
      iset = ImageSweep(
        isetdata,
        beam       = beam,
        detector   = detector,
        goniometer = goniometer,
        scan       = scan,
        indices    = single_file_indices)

    # Return the imageset
    return iset
