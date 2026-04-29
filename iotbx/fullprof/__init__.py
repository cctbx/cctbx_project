""" Interface cctbx with the FullProf program suite.

FullProf can be obtained from http://www.ill.eu/sites/fullprof/
"""
from __future__ import absolute_import, division, print_function
from cctbx.eltbx import wavelengths

def rietveld_refine_structure(crystalstructure,
                      wavelength=wavelengths.characteristic("CU").as_angstrom(),
                      I_obs=None, Profile=None, ProfileFile=None):
  """This function tries to rietveld refine a structure using FullProf.

  If I_obs is given neither a Profile or a ProfileFile may be specified.
  'Profile' and 'ProfileFile' are also exclusive.

  :param crystalstructure: the starting model for the refinement
  :type crystalstructure: cctbx.xray.structure
  :param wavelength: the x-ray wavelength for the given intensity data
  :type wavelength: float
  :param I_obs: observed Intensities
  :type I_obs: cctbx.miller
  :param Profile: a XRD profile given as a list/tuple of (2-theta, intensity)-tuples
  :type Profile: list/tuple(tuple(2theta,I))
  :param ProfileFile: path to a XRD profile as FullProf .prf file
  :type ProfileFile: string
  """
  # Check preconditions
  if [I_obs, Profile, ProfileFile].count(None) != 2:
    raise ValueError("You may only pass one of I_obs, Profile and ProfileFile")
  # start work
  from iotbx.write_pcr import write_pcr
  import tempfile
  import shutil
  import os
  # write pcr file and execute FullProf
  try:
    f = tempfile.NamedTemporaryFile(suffix=".pcr", delete=False)
  except IOError: raise
  pcrfile = f.name
  basepath = os.path.splitext(pcrfile)[0]
  write_pcr(f, crystalstructure, jobtype=0,
            wavelength=wavelength, I_obs=I_obs)
  f.close()
  try:
    if ProfileFile is not None:
      shutil.copyfile(ProfileFile, basepath+".dat")
    elif Profile is not None:
      # write out profile file for FullProf
      # XXX: todo implement
      pass
  except IOError: raise
  run_fullprof(pcrfile, verbose=0)
  # XXX Todo: extract refined structure from resulting .pcr/.new
  # XXX Todo: extract Rwp of refined structure from resulting .pcr/.new
  return None, None



def simulate_powder_pattern(crystalstructure,
                     wavelength=wavelengths.characteristic("CU").as_angstrom(),
                     filename="",
                     keep_results=False,
                     scale_down=1.0):
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
  :param scale_down: factor to divide intensities by (to avoid overflows)
  :type scale_down: float

  :returns: calculated integral intensities, calculated profile
  :rtype: cctbx.miller, list(tuple(float,int))

  XXX Todo: implement extraction of calculated profile
  """
  from iotbx.write_pcr import write_pcr
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
  write_pcr(f, crystalstructure, jobtype=2,
            wavelength=wavelength, scale_down=scale_down)
  f.close()
  if scale_down > 10000: raise
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
  except Exception:
    # an overflow occured
    return simulate_powder_pattern(crystalstructure, wavelength=wavelength,
                                   filename=filename, keep_results=keep_results,
                                   scale_down=scale_down*10.0)
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
  if not os.path.exists(pcrfile + ".pcr"):
    raise IOError(pcrfile + ".pcr not found!")
  pcrfile = os.path.abspath(pcrfile)
  old_cwd = os.getcwd()
  workdir = os.path.split(pcrfile)[0]
  try:
    os.chdir(workdir)
    os.unlink(pcrfile + ".out")
  except KeyboardInterrupt: raise
  except Exception: pass
  fullprof_out = easy_run.fully_buffered(command="fp2k "+os.path.split(pcrfile)[1]) \
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
  os.chdir(old_cwd)

if __name__ == '__main__':
  # just a little test for debugging
  from cctbx import sgtbx
  from cctbx.development import random_structure
  xrs = random_structure.xray_structure(
        space_group_info=sgtbx.space_group_info(number=1),
        elements=["C"]*10,u_iso=0.005)
  print(list(simulate_powder_pattern(xrs)[0]))#, filename="/tmp/test.pcr", keep_results=True))
