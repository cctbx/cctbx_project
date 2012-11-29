from __future__ import division
from __future__ import print_function
from cctbx.eltbx import wavelengths
from cctbx.iotbx.fullprof.write_pcr import write_pcr

""" The purpose of this module is to interface cctbx with the FullProf progamm
suite.

FullProf can be obtained from http://www.ill.eu/sites/fullprof/
"""

def cctbx_xray_structure_from(cls, file=None, filename=None):
  # XXX: a pcr reader has still to be implemented
  from iotbx import builders
  #builder = builders.crystal_structure_builder(
  #  set_grad_flags=set_grad_flags,
  #  min_distance_sym_equiv=min_distance_sym_equiv)
  #stream = command_stream(file=file, filename=filename)
  #stream = crystal_symmetry_parser(stream, builder)
  #stream = atom_parser(stream.filtered_commands(), builder, strictly_shelxl)
  #stream.parse()
  return None #builder.structure

def rietfeld_refine_structure(crystalstructure, Iobs):
  # XXX: todo implement
  pass

def simulate_powder_pattern(crystalstructure,
                     wavelength=wavelengths.characteristic("CU").as_angstrom()):
    return None

def run_fullprof(pcrfile, verbose=0):
  from libtbx import easy_run
  import sys, os
  """Run fullprof on a prepared .pcr file

  :param verbose: be verbose
  :type verbose: integer
  """
  sys.stdout.flush()
  sys.stderr.flush()
  try: os.unlink(pcrfile + ".out")
  except KeyboardInterrupt: raise
  except Exception: pass
  fullprof_out = easy_run.fully_buffered(command="fp2k " + pcrfile) \
    .raise_if_errors() \
    .stdout_lines
  if (0 or verbose):
    for l in fullprof_out: print l
  f = open(pcrfile + ".out", "r")
  fullprof_out = f.readlines()
  f.close()
  sys.stderr.flush()
  if (0 or verbose):
    for l in fullprof_out: print l[:-1]
  sys.stdout.flush()
