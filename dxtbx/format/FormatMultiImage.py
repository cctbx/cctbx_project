from __future__ import division

class FormatMultiImage(Format):

    def __init__(self, image_file):
        Format.__init__(self, image_file)

    def get_num_images(self):
        pass

    def get_goniometer(self, index=None):
        pass

    def get_detector(self, index=None):
        pass

    def get_beam(self, index=None):
        pass

    def get_scan(self, index=None):
        pass

    def get_raw_data(self, index=None):
        pass

    def get_detectorbase(self, index=None):
        pass

    def get_image_file(self, index=None):
        pass
