from __future__ import division

def template_regex(filename):
  '''Try a bunch of templates to work out the most sensible. N.B. assumes
  that the image index will be the last digits found in the file name.'''

  import re

  # filename template code stolen from xia2...

  if not hasattr(template_regex, "pattern"):
    # Compile regular expressions once and store compiled version

    # The filename is reversed for pattern evaluation.
    # These are reversed patterns.

    patterns = ['()([0-9]+)(\..*)',
                   # filename ends with numbers
                   #  img.0815
                '([a-zA-Z-]+\.)([0-9]+)(\..*)',
                   # filename ends with numbers followed by simple extension
                   #  img.0815.cbf
                '(.*?\.)([0-9]+)(_.*)',
                   # last number in the filename standing between _ and .
                   #  NO2_0100.sweep.18keV
                '(.*?\.)([0-9]+)(.*)']
                   # last number in the filename before a .
                   #  NO2.00100.sweep.18keV or image-00001.cbf

    # patterns are evaluated in order, the first to match will be used
    # each pattern needs to describe the entire filename with three braces:
    #  (suffix)(imagenumber)(prefix)
    # if more braces are used in a pattern they MUST start with (?:...) so that
    # they do not interfere with the matching
    template_regex.pattern = [re.compile(pattern) for pattern in patterns]

  rfilename = filename[::-1]

  template = None
  digits = 0
  for j, cp in enumerate(template_regex.pattern):
    match = cp.match(rfilename)
    if not match:
      continue
    groups = match.groups()

    exten = groups[0][::-1]
    digits = groups[1][::-1]
    prefix = groups[2][::-1]

    template = prefix + ('#' * len(digits)) + exten
    break

# print "File name template:\n    %s\n -> %s (%d)" % (filename, template, int(digits))
  return template, int(digits)


def group_files_by_imageset(filenames):
  '''Group filenames by supposed imageset.

  Get the template for each file in the list. Then add to a dictionary
  indexed by template containing a list of indices within that template.
  For files that do not match a template, these are added by filename
  instead.

  '''
  from collections import defaultdict

  # Calculate the template for each image. If the template is None
  # (i.e. there are no numbers to identify the filename, add the
  # filename itself.
  template = []
  for f in filenames:
    t = template_regex(f)
    if t[0] == None:
      template.append((f, None))
    else:
      template.append(t)

  # Loop through all the templates and add the new item to a dictionary
  # with a list of files per template.
  matched = defaultdict(list)
  for t in template:
    matched[t[0]].append(t[1])

  # Return the matched filenames
  return matched


def find_matching_images(image_name):
  '''Search in the directory in which this image is for images which share
  the same template: return this list.'''

  import os

  directory, filename = os.path.split(image_name)
  if directory is None or directory == '':
    directory = '.'

  template, digits = template_regex(filename)

  if template:
    len_digits = template.count('#')
    pfx = template.split('#')[0]
    sfx = template.split('#')[-1]
    template_str = pfx + '%%0%dd' % len_digits + sfx

    files_in_directory = os.listdir(directory)

    matching_images = []

    for j in range(0, 10 ** len_digits):
      if template_str % j in files_in_directory:
        matching_images.append(os.path.join(directory, template_str % j))

  else:
    matching_images = [image_name]

  return matching_images

def replace_template_format_with_hash(match):
  '''Replace the format match with hashes'''
  return '#' * len(match.group(0) % 0)

def template_format_to_string(template):
  '''Convert the template format to the template string.'''
  import re
  return re.sub(r'%0[0-9]+d', replace_template_format_with_hash, template)

def template_string_to_glob_expr(template):
  '''Convert the template to a glob expression.'''
  pfx = template.split('#')[0]
  sfx = template.split('#')[-1]
  return '%s%s%s' % (pfx, '[0-9]'*template.count('#'), sfx)

def template_string_number_index(template):
  '''Get the number idex of the template.'''
  pfx = template.split('#')[0]
  sfx = template.split('#')[-1]
  return len(pfx), len(pfx) + template.count('#')

def locate_files_matching_template_string(template):
  '''Return all files matching template.'''
  from glob import glob
  return glob(template_string_to_glob_expr(template))

def template_image_range(template):
  '''Return the image range of files with this template.'''

  # Find the files matching the template
  filenames = locate_files_matching_template_string(template)
  filenames = sorted(filenames)

  # Check that the template matches some files
  if len(filenames) == 0:
    raise ValueError('Template {0} doesn\'t match any files.'.format(
        template))

  # Get the templete format
  index = slice(*template_string_number_index(template))

  # Get the first and last indices
  first = int(filenames[0][index])
  last  = int(filenames[-1][index])

  # Reutrn the image range
  return (first, last)
