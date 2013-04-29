from __future__ import division

class FormatMultiImage(object):

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

    def get_detectorbase(self, index=None):
        raise RuntimeError('Overload!')

    def get_image_file(self, index=None):
        raise RuntimeError('Overload!')
