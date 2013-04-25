#!/usr/bin/env python
#
# image_set_factory.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

class FilenameAnalyser(object):
    '''Group images by filename into image sets.'''

    def __init__(self):
        '''Initialise the class.'''
        pass

    def __call__(self, filenames):
        '''Group the filenames by imageset.

        Params:
            filenames The list of filenames

        Returns:
            A list of (template, [indices], is_sweep)

        '''
        from sweep_filenames import group_files_by_imageset

        # Analyse filenames to figure out how many imagesets we have
        filelist_per_imageset = group_files_by_imageset(filenames)

        # Label each group as either an imageset or a sweep.
        file_groups = []
        for template, indices in filelist_per_imageset.iteritems():

            # Check if this imageset is a sweep
            is_sweep = self._is_imageset_a_sweep(template, indices)

            # Append the items to the group list
            file_groups.append((template, indices, is_sweep))

        # Return the groups of files
        return file_groups

    def _is_imageset_a_sweep(self, template, indices):
        ''' Return True/False if the imageset is a sweep or not.

        Where more than 1 image that follow sequential numbers are given
        the images are catagorised as belonging to a sweep, otherwise they
        belong to an image set.

        '''
        if len(indices) <= 1:
            return False
        else:
            indices = sorted(indices)
            if self._indices_sequential(indices):
                return True
            else:
                return False

    def _indices_sequential(self, indices):
        ''' Determine if indices are sequential.'''
        prev = indices[0]
        for curr in indices[1:]:
            if curr != prev + 1:
                return False
            prev = curr

        return True


class ImageSetFactory(object):
    ''' Factory to create imagesets and sweeps. '''

    @staticmethod
    def new(filenames, check_headers=False):
        ''' Create an imageset or sweep

        Params:
            filenames A list of filenames
            check_headers Check the headers to ensure all images are valid

        Returns:
            A list of imagesets

        '''
        from sweep_filenames import find_matching_images

        # Ensure we have enough images
        if isinstance(filenames, list):
            assert(len(filenames) > 0)
        elif isinstance(filenames, str):
            filenames = [filenames]
        else:
            raise RuntimeError, 'unknown argument passed to ImageSetFactory'

        # Analyse the filenames and group the images into imagesets.
        analyse_files = FilenameAnalyser()
        filelist_per_imageset = analyse_files(filenames)

        # For each file list denoting an image set, create the imageset
        # and return as a list of imagesets. N.B sweeps and image sets are
        # returned in the same list.
        return [ImageSetFactory._create_imageset_or_sweep(
            filelist, check_headers) for filelist in filelist_per_imageset]

    @staticmethod
    def _create_imageset_or_sweep(filelist, check_headers):
        '''Create either an imageset of sweep.'''
        if filelist[2] == True:
            return ImageSetFactory._create_sweep(filelist, check_headers)
        else:
            return ImageSetFactory._create_imageset(filelist, check_headers)

    @staticmethod
    def _create_imageset(filelist, check_headers):
        '''Create an image set'''
        from dxtbx.format.Registry import Registry
        from dxtbx.image_set import ImageSet, BufferedImageSetReader

        # Extract info from filelist
        template, indices, is_sweep = filelist

        # Get the template format
        count = template.count('#')
        if count > 0:
            pfx = template.split('#')[0]
            sfx = template.split('#')[-1]
            template_format = '%s%%0%dd%s' % (pfx, template.count('#'), sfx)
            filenames = [template_format % index for index in indices]
        else:
            filenames = [template]

        # Sort the filenames
        filenames = sorted(filenames)

        # Get the format object
        format_class = Registry.find(filenames[0])

        # Create the image set object
        image_set = ImageSet(BufferedImageSetReader(format_class, filenames))

        # Check the image set is valid
        if check_headers and not image_set.is_valid():
            raise RuntimeError('Invalid ImageSet')

        # Return the image set
        return image_set

    @staticmethod
    def _create_sweep(filelist, check_headers):
        '''Create a sweep'''
        import os
        from dxtbx.format.Registry import Registry
        from dxtbx.sweep import Sweep, BufferedSweepReader

        # Extract info from filelist
        template, indices, is_sweep = filelist

        # Get the template format
        count = template.count('#')
        if count > 0:
            pfx = template.split('#')[0]
            sfx = template.split('#')[-1]
            template_format = '%s%%0%dd%s' % (pfx, template.count('#'), sfx)
            filenames = [template_format % index for index in indices]
        else:
            filenames = [template]

        # Sort the filenames
        filenames = sorted(filenames)

        # Get the format object
        format_class = Registry.find(filenames[0])

        # Get the first image and our understanding
        first_image = filenames[0]

        # Get the directory and first filename and set the template format
        directory, first_image_name = os.path.split(first_image)
        first_image_number = indices[0]

        # Get the template format
        pfx = template.split('#')[0]
        sfx = template.split('#')[-1]
        template_format = '%s%%0%dd%s' % (pfx, template.count('#'), sfx)
        template = os.path.join(directory, template)

        # Set the image range
        image_range = (min(indices), max(indices))

        # Create the sweep object
        sweep = Sweep(BufferedSweepReader(
            format_class, template, image_range))

        # Check the sweep is valid
        if check_headers and not sweep.is_valid():
            raise RuntimeError('Invalid sweep of images')

        # Return the sweep
        return sweep












if __name__ == '__main__':
    pass
#    filenames = ['shot-s.pickle',
#                 'shot-ss.pickle',
#                 'shot-s00_20130311224756207.pickle',
#                 'shot-s0020130311224758207.pickle',
#                 'shot-s00-201-Y-30311224758207.pickle',
#
#                 'shot-s00-20130311224753249.pickle',
#                 'shot-s00-20130311224756207.pickle',
#                 'shot-s00-20130311224758207.pickle',
#                 'shot-s00-20130311224837983.pickle',
#                 'shot-s00-20130311224842441.pickle',
#                 'shot-s01-20130311224501171.pickle',
#                 'shot-s01-20130311224501254.pickle',
#                 'shot-s01-20130311224502379.pickle',
#                 'shot-s01-20130311224505379.pickle',
#                 'shot-s02-20130311224602371.pickle',
#                 'shot-s02-20130311224603163.pickle',
#                 'shot-s02-20130311224603413.pickle',
#                 'shot-s02-20130311224609164.pickle',
#                 'shot-s02-20130311224609789.pickle',
#                 'shot-s02-20130311224610622.pickle',
#                 'shot-s02-20130311224632835.pickle',
#                 'shot-s02-20130311224636002.pickle',
#                 'shot-s02-20130311224639336.pickle',
#                 'shot-s04-20130311224524986.pickle',
#                 'shot-s04-20130311224526569.pickle',
#                 'shot-s04-20130311224526944.pickle',
#                 'shot-s04-20130311224532611.pickle',
#                 'shot-s04-20130311224533319.pickle',
#                 'shot-s04-20130311224535069.pickle',
#                 'shot-s04-20130311224535402.pickle',
#                 'shot-s04-20130311224537903.pickle',
#                 'shot-s04-20130311224538903.pickle']
#
#    filenames = [
#            'SIM_MX_mod_001.cbf',
#            'thaumatin_die_M1S5_1_asc_0041.cbf',
#            'q315r_lyso_1_001.img',
#            'als501_q4_1_001.img',
#            '200mMNaCl5pcGlyc_400.edf',
#            'q210_lyso_1_101.img',
#            'q315r_lyso_001.img',
#            'q315_1_001.img',
#            'q210_1_001.img',
#            'APS1_0.0001',
#            'q315_unbinned_a.0001.img',
#            'mar300.0001',
#            'mar300_1_E1.0001',
#            'pilatus_1_0001.cbf',
#            'q315_1_001.img',
#            'mar225_2_E0_0001.img',
#            'mar345_01_001.mar2300',
#            'q210_2_001.img',
#            's01f0001.osc',
#            'mar165_001.mccd',
#            'mar225_1_001.mccd',
#            'q315r_7_001.img',
#            'lys_001.osc',
#            '000.pickle',
#            'shot-s00-2011-12-02T21_07Z29.723_00569.pickle',
#            'shot-s04-20111204004533388.pickle',
#            'reallysurprise_001.ipf',
#            'fake_00001.img',
#            'test1_lysozyme_0111060001.osc',
#            'lyso_00001.img',
#            'mar225_2_001.img',
#            'pilatus6m_1_00001.cbf',
#            'Xtal17-2phi_3_015.cbf',
#            'lys_00001.img',
#            'mar225_001.img',
#            'q4_1_001.img',
#            'mar325_1_001.mccd',
#            'q315_1_001.img',
#            'lyziph6p5_01_0001.sfrm',
#            'XPARM.XDS']
#
#    imagesets = ImageSetFactory.new(filenames)

#    for s in imagesets:
#        print s
