from __future__ import division
#!/usr/bin/env python
# dxtbx/sweep.py
#
#   Copyright (C) 2013 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Top level interface to the dxtbx: treatment of a sweep of images as a single
# unit, incorporation of the registry and so on... N.B. factory functions will
# be defined which can look up information from the file system etc.

class sweep:
    '''Definition for a sweep of images, defined to be a set of diffraction
    images in files with matching templates, or a volume formatted file.'''

    def __init__(self, list_of_images,
                 beam_model = None,
                 goniometer_model = None,
                 detector_model = None,
                 scan_models = { },
                 check_headers = False):
        '''Construct a sweep object; by default verify that this is accurately
        modelled as one sweep; check consistency in the image headers. N.B.
        here a sweep is defined as contiguous images with uniform oscillation
        widths... Also N.B. generally speaking as a user you are better off with
        one of the factory functions. Finally also please note that the full
        check_headers could be computationally expensive.'''

        assert(len(list_of_images) > 0)

        first_image = list_of_images[0]

        self._beam_model = beam_model
        self._goniometer_model = goniometer_model
        self._detector_model = detector_model
        self._scan_model = scan_models.get(first_image, None)

        from dxtbx.format.Registry import Registry

        format_class = Registry.find(first_image)
        understanding = format_class.understand(first_image)
        format_instance = format_class(first_image)
        
        if not self._beam_model:
            self._beam_model = format_instance.get_beam()

        if not self._goniometer_model:
            self._goniometer_model = format_instance.get_goniometer()

        if not self._detector_model:
            self._detector_model = format_instance.get_detector()

        if not self._scan_model:
            self._scan_model = format_instance.get_scan()

        if check_headers:
            for image in list_of_images[1:]:
                assert(format_class.understand(image) == understanding)
                other_instance = format_class(image)
                assert(other_instance.get_beam() == self._beam_model)
                assert(other_instance.get_goniometer() == self._goniometer_model)
                assert(other_instance.get_detector() == self._detector_model)

                # cache as I will want this later on...
                
                if not image in scan_models:
                    scan_models[image] = other_instance.get_scan()

        # now construct the full scan model from the image headers

        for image in list_of_images[1:]:
            scan = scan_models.get(image, None)
            if not scan:
                scan = format_class(image).get_scan()
            self._scan_model += scan

        self._format_class = format_class

        return
        
    
