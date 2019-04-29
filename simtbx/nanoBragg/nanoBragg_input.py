from __future__ import division, print_function

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 11/03/2017
Description : SIMTBX (nanoBragg) I/O module. Reads PHIL input.
'''

import iotbx.phil as ip
import inspect

from simtbx.nanoBragg import nanoBragg

def isprop(v):
  ''' Test if attribute is a property '''
  return isinstance(v, property)

def flatten(L):
  for item in L:
    try:
      for i in flatten(item):
        yield i
    except TypeError:
      yield item

def generate_simtbx_phil():
  ''' Generates a PHIL object from simTBX params

  @return: sim_phil: a PHIL object
  '''

  ui_phil_string = """
description = None
  .type = str
  .help = Run description (optional).
  .optional = True
output = None
  .type = path
  .help = Base output directory, current directory in command-line, can be set in GUI
  .optional = False
image_prefix = synthetic_image
  .type = str
  .help = Prefix for synthetic image filenames (numbers will be added)
image_format = *img cbf mccd pickle
  .type = choice
  .help = Will be format in the future, just extension for now
reference_coordinates = None
  .type = path
  .help = Coordinates (PDB) to generate FCalc
reference_FCalc = None
  .type = path
  .help = FCalc (MTZ) that can be used to generate the diffraction pattern
reference_image = None
  .type = path
  .help = Diffraction image that can be used to match orientation
radial_average_background = None
  .type = path
  .help = Radial average background (sin(theta) / lambda, STOL)
dataset
  .help = Synthetic dataset options
{
  start_phi = 0
    .type = float
    .help = phi angle for the first image in the dataset
  finish_phi = 90
    .type = float
    .help = phi angle for the last image in the dataset
  oscillation = 1
    .type = float
    .help = phi oscillation per dataset image
}
"""

  sim = nanoBragg()
  params = [name for (name, value) in inspect.getmembers(nanoBragg, isprop)]

  param_lines = ['simtbx',
                 '{']
  for i in params:
    # The below properties result in boost errors (Fbg_vs_stol), TypeError
    # exceptions (xray_beams), or yield objects rather than values (the rest)
    if i in ('Fbg_vs_stol', 'Fhkl_tuple', 'amplitudes', 'indices', 'raw_pixels',
             'xray_source_XYZ', 'xray_source_intensity_fraction', 'xray_beams',
             'xray_source_wavelengths_A', 'unit_cell_tuple', 'progress_pixel'):
      continue

    try:
      param_value = getattr(sim, i)
      if type(param_value).__name__ == 'tuple':
        element_type = type(list(flatten(param_value))[0]).__name__
        param_type = "{}s".format(element_type)
        param_value = ' '.join(str(i) for i in list(flatten(param_value)))
        param_doc = getattr(nanoBragg, i).__doc__
      else:
        param_type = type(param_value).__name__
        param_doc = getattr(nanoBragg, i).__doc__
        if param_type == 'unit_cell':
          param_value = ' '.join(str(i) for i in param_value.parameters())

      # Correct a number of weird vartypes
      if param_type in ('convention', 'pivot', 'shapetype'):
        param_type = 'str'


      # PHIL interprets ';' as a line break
      param_doc = param_doc.replace(';', '.')

      param_line = '  {} = {} \n    .type = {}\n    .help = {}' \
                   ''.format(i, param_value, param_type, param_doc)
      param_lines.append(param_line)

    except TypeError as e:
      print(i, e)
      pass

  param_lines.append('}')
  sim_phil_string = '\n'.join(param_lines)
  combo_string = ui_phil_string + sim_phil_string

  sim_phil = ip.parse(combo_string)
  return sim_phil
