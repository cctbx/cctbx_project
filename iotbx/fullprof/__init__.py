from __future__ import division
from __future__ import print_function
from cctbx.eltbx import wavelengths

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
                     wavelength=wavelengths.characteristic("CU").as_angstrom(),
                     filename="",
                     keep_results=False):
  """
  Get integrated intensities and a a simulated XRD profile calculated with
  FullProf (has to be installed and callable via "fp2k").

  :param crystalstructure: a crystal structure to calculate the intensities for
  :type crystalstructure: cctbx.xray.structure
  :param wavelength: x-ray wavelength in angstroms
  :type wavelength: float
  :param filename: a filepath to save the in- and output of FullProf to
  :type filename: string
  :param keep_results: keep the (temporary) files from FullProf for a later \
  manual inspection
  :type keep_results: boolean

  :returns: calculated integral intensities, calculated profile
  :rtype: cctbx.miller, list(tuple(float,int))

  XXX Todo: implement extraction of calculated profile
  """
  from write_pcr import write_pcr
  from iotbx.reflection_file_reader import any_reflection_file
  import tempfile
  import os
  # write pcr file and execute FullProf
  try:
    if filename == "":
      f = tempfile.NamedTemporaryFile(suffix=".pcr", delete=False)
    else:
      f = open(filename, "w")
  except IOError: raise
  pcrfile = f.name
  basepath = os.path.splitext(pcrfile)[0]
  write_pcr(f, crystalstructure, jobtype=2, wavelength=wavelength)
  f.close()
  run_fullprof(pcrfile, verbose=0)

  # fix hkl file for hkl reader
  hklfile = basepath + ".fou"
  f = open(hklfile, "r")
  lines = f.readlines()[1:]
  f.close()
  hklfile = basepath + ".hkl"
  if keep_results == True:
    os.rename(hklfile, basepath + ".hkldata")   # backup old hkl file
  f = open(hklfile, "w")
  f.writelines(lines)
  f.close()

  # extract intensities
  f_calc = None
  try:
    rf = any_reflection_file(file_name=hklfile)
    f_calc = rf.as_miller_arrays(
                          crystal_symmetry=crystalstructure.crystal_symmetry(),
                          assume_shelx_observation_type_is="intensities")[0]
  except KeyboardInterrupt: raise
  except Exception: raise
  profile = None

  # clean up
  if keep_results == False:
    for ext in [".pcr", ".sym", ".fou", ".ins", ".prf", ".sim", ".sum", ".out",
                ".hkl", ".hkldata", "1.fst", "1.sub"]:
      try:
        os.unlink(basepath + ext)
      except KeyboardInterrupt: raise
      except Exception: pass

  return f_calc, profile

def run_fullprof(pcrfile, verbose=0):
  from libtbx import easy_run
  import sys, os
  """Run fullprof on a prepared .pcr file

  :param verbose: be verbose
  :type verbose: integer
  """
  sys.stdout.flush()
  sys.stderr.flush()
  if pcrfile.lower().endswith(".pcr"):
    pcrfile = os.path.splitext(pcrfile)[0]
  try: os.unlink(pcrfile + ".out")
  except KeyboardInterrupt: raise
  except Exception: pass
  fullprof_out = easy_run.fully_buffered(command="fp2k " + pcrfile) \
    .raise_if_errors() \
    .stdout_lines
  if (0 or verbose):
    for l in fullprof_out: print(l)
  f = open(pcrfile + ".out", "r")
  fullprof_out = f.readlines()
  f.close()
  sys.stderr.flush()
  if (0 or verbose):
    for l in fullprof_out: print(l[:-1])
  sys.stdout.flush()

if __name__ == '__main__':
  # just a little test for debugging
  from cctbx import sgtbx
  from cctbx.development import random_structure
  xrs = random_structure.xray_structure(
        space_group_info=sgtbx.space_group_info(number=1),
        elements=["C"]*10,u_iso=0.005)
  print(list(simulate_powder_pattern(xrs)[0]))#, filename="/tmp/test.pcr", keep_results=True))
