from __future__ import division
#!/usr/bin/env python
# scan_helpers.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Helpers for the scan class, which are things for handling e.g. filenames,
# templates and so on.

import os
import re
import string
import math
import copy

def image2template(filename):
    '''Return a template to match this filename.'''

    assert(not '#' in filename)

    # the patterns in the order I want to test them - these should be kept
    # in an extensible singleton somewhere

    pattern_keys = [r'([^\.]*)\.([0-9]+)\Z',
                    r'([^\.]*)\.([0-9]+)(.*)',
                    r'(.*)_([0-9]*)\.(.*)',
                    r'(.*?)([0-9]*)\.(.*)']

    # patterns is a dictionary of possible regular expressions with
    # the format strings to put the file name back together

    patterns = {r'([^\.]*)\.([0-9]+)\Z':'%s.%s%s',
                r'([^\.]*)\.([0-9]+)(.*)':'%s.%s%s',
                r'(.*)_([0-9]*)\.(.*)':'%s_%s.%s',
                r'(.*?)([0-9]*)\.(.*)':'%s%s.%s'}

    for pattern in pattern_keys:
        match = re.compile(pattern).match(filename)

        if match:
            prefix = match.group(1)
            number = match.group(2)
            try:
                exten = match.group(3)
            except: # intentional
                exten = ''

            for digit in string.digits:
                number = number.replace(digit, '#')

            return patterns[pattern] % (prefix, number, exten)

    raise RuntimeError, 'filename %s not understood as a template' % \
          filename

def image2image(filename):
    '''Return an integer for the template to match this filename.'''

    # check that the file name doesn't contain anything mysterious
    if filename.count('#'):
        raise RuntimeError, '# characters in filename'

    # the patterns in the order I want to test them

    pattern_keys = [r'([^\.]*)\.([0-9]+)\Z',
                    r'([^\.]*)\.([0-9]+)(.*)',
                    r'(.*)_([0-9]*)\.(.*)',
                    r'(.*?)([0-9]*)\.(.*)']

    for pattern in pattern_keys:
        match = re.compile(pattern).match(filename)

        if match:
            prefix = match.group(1)
            number = match.group(2)
            try:
                exten = match.group(3)
            except: # intentional
                exten = ''

            return int(number)

    raise RuntimeError, 'filename %s not understood as a template' % \
          filename

def image2template_directory(filename):
    '''Separate out the template and directory from an image name.'''

    directory = os.path.dirname(filename)

    if not directory:

        # then it should be the current working directory
        directory = os.getcwd()

    image = os.path.split(filename)[-1]
    template = image2template(image)

    return template, directory

def find_matching_images(template, directory):
    '''Find images which match the input template in the directory
    provided.'''

    files = os.listdir(directory)

    # to turn the template to a regular expression want to replace
    # however many #'s with EXACTLY the same number of [0-9] tokens,
    # e.g. ### -> ([0-9]{3})

    # change 30/may/2008 - now escape the template in this search to cope with
    # file templates with special characters in them, such as "+" -
    # fix to a problem reported by Joel B.

    length = template.count('#')
    regexp_text = re.escape(template).replace(
        '\\#' * length, '([0-9]{%d})' % length)
    regexp = re.compile(regexp_text)

    images = []

    for f in files:
        match = regexp.match(f)

        if match:
            images.append(int(match.group(1)))

    images.sort()

    return images

def template_directory_number2image(template, directory, number):
    '''Construct the full path to an image from the template, directory
    and image number.'''

    # FIXME why does this duplicate code shown below??

    length = template.count('#')

    # check that the number will fit in the template

    if (math.pow(10, length) - 1) < number:
        raise RuntimeError, 'number too big for template'

    # construct a format statement to give the number part of the
    # template
    format = '%%0%dd' % length

    # construct the full image name
    image = os.path.join(directory,
                         template.replace('#' * length,
                                          format % number))

    return image

def template_number2image(template, number):
    '''Construct the an image from the template and image number.'''

    length = template.count('#')

    # check that the number will fit in the template

    if (math.pow(10, length) - 1) < number:
        raise RuntimeError, 'number too big for template'

    format = '%%0%dd' % length

    image = template.replace('#' * length, format % number)

    return image

class scan_helper_image_files:
    '''A helper class which handles things like image names, making templates,
    finding matching images and so on. Currently this just provides aliases
    to existing functions elsewhere, but ultimately it would be good if they
    were all encapsulated herein.'''

    @staticmethod
    def image_to_template(filename):
        '''From an image name, return a file template which should match.'''
        return image2template(filename)

    @staticmethod
    def image_to_index(filename):
        '''From an image name, determine the index within the scan for this
        image, complementary to the image_to_template method above.'''
        return image2image(filename)

    @staticmethod
    def image_to_template_directory(filename):
        '''From a full path to an image, return the filename template and
        directory.'''
        return image2template_directory(filename)

    @staticmethod
    def template_directory_to_indices(template, directory):
        '''For a given template and directory, return a list of image indices
        which match. Also complementary with image_to_template_directory.'''
        return find_matching_images(template, directory)

    @staticmethod
    def template_directory_index_to_image(template, directory, index):
        '''Construct the full image name from the template, directory and
        file index.'''
        return template_directory_number2image(template, directory, index)

    @staticmethod
    def template_index_to_image(template, index):
        '''Construct the image file name from the template and file index.'''
        return template_number2image(template, index)

class scan_helper_image_formats:
    '''A helper class which enxapsulates the allowed and supported image
    formats namely CBF, TIFF, SMV, RAXIS, MAR. N.B. there will be some
    crosstalk between this class and the _image_format classes.'''

    FORMAT_CBF = 'FORMAT_CBF'
    FORMAT_TIFF = 'FORMAT_TIFF'
    FORMAT_SMV = 'FORMAT_SMV'
    FORMAT_RAXIS = 'FORMAT_RAXIS'
    FORMAT_MAR = 'FORMAT_MAR'

    @staticmethod
    def check_format(format):
        if format in [scan_helper_image_formats.FORMAT_CBF,
                      scan_helper_image_formats.FORMAT_TIFF,
                      scan_helper_image_formats.FORMAT_SMV,
                      scan_helper_image_formats.FORMAT_RAXIS,
                      scan_helper_image_formats.FORMAT_MAR]:
            return True
        print "Format %s may not be supported--contact Nick Sauter"%format
        return True
