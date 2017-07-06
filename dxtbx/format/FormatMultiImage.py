from __future__ import absolute_import, division

class FormatMultiImage(object):

  def __init__(self, **kwargs):
    pass

  def get_num_images(self):
    raise RuntimeError('Overload!')

  def get_goniometer(self, index=None):
    raise RuntimeError('Overload!')

  def get_detector(self, index=None):
    raise RuntimeError('Overload!')

  def get_beam(self, index=None):
    raise RuntimeError('Overload!')

  def get_scan(self, index=None):
    raise RuntimeError('Overload!')

  def get_raw_data(self, index=None):
    raise RuntimeError('Overload!')

  def get_mask(self, index=None, goniometer=None):
    return None

  def get_detectorbase(self, index=None):
    raise RuntimeError('Overload!')

  def get_image_file(self, index=None):
    raise RuntimeError('Overload!')

  @classmethod
  def get_detectorbase_factory(Class):
    '''
    Return a factory object to create detector base instances

    '''

    class DetectorBaseFactory(object):

      def __init__(self, filenames):
        assert len(filenames) == 1
        self._filename = filenames[0]

      def get(self, index):
        format_instance = Class.get_instance(self._filename)
        return format_instance.get_detectorbase(index)

      def copy(self, filenames):
        return MaskerDetectorBaseFactory(filenames)

      def __len__(self):
        format_instance = Class.get_instance(self._filename)
        return format_instance.get_num_images()

    return DetectorBaseFactory

  @classmethod
  def get_reader(Class):
    '''
    Return a reader class

    '''

    class Reader(object):

      def __init__(self, filenames):
        assert len(filenames) == 1
        self._filename = filenames[0]

      def read(self, index):
        format_instance = Class.get_instance(self._filename)
        return format_instance.get_raw_data(index)

      def paths(self):
        return [self._filename]

      def __len__(self):
        format_instance = Class.get_instance(self._filename)
        return format_instance.get_num_images()

      def copy(self, filenames):
        return Reader(filenames)

      def identifiers(self):
        return ["%s-%d" % (self._filename, index) for index in range(len(self))]

    return Reader

  @classmethod
  def get_masker(Class):
    '''
    Return a reader class

    '''

    class Masker(object):

      def __init__(self, filenames):
        assert len(filenames) == 1
        self._filename = filenames[0]

      def get(self, index):
        format_instance = Class.get_instance(self._filename)
        return format_instance.get_mask(index)

      def paths(self):
        return [self._filename]

      def __len__(self):
        format_instance = Class.get_instance(self._filename)
        return format_instance.get_num_images()

      def copy(self, filenames):
        return Masker(filenames)

    return Masker
