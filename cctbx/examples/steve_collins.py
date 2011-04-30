# > What we are interested in is a standard crystallographic library that we
# > can call from Jython as well as data analysis packages such as
# > Matlab/Mathematica. The kinds of things we would like would include:
#
# General note: the cctbx is a Python-based library with C++ extensions.
# The extensions are essential for all crystallographic algorithms.
# I've never tested how this plays with Jython.
#
# > Low-level:
# > Generate space groups (in matrix/vector form) based on spacegroup number
# > (names are a pain)
# > Generate conditions for allowed reflections etc
# > Generate point group of space group in matrix form
# > Generate point-group matrices for given coordinate
# > Access databases for form factors, anomalous scattering factors,
#
# See the attached example script.
# This script should work with one of the bundles from here:
#
#   http://cci.lbl.gov/cctbx_build/
#
# If you can, start with a binary bundle. These install in seconds.
# If you start from sources it may take 30 minutes or more.
#
# > Higher level:
# > Read/(write?) CIF file
#
# PyCifRW is included in the cctbx bundles. See the PyCifRW documentation:
#
#   http://www.ansto.gov.au/natfac/ANBF/CIF/
#
# See also newsletter No. 5:
#
#   http://cci.lbl.gov/publications/download/iucrcompcomm_jan2005.pdf
#
# > Calculate reflection lists with intensities
#
# See newsletters No. 1 (quartz example) and No. 5 (PyCifRW example):
#
#   http://cci.lbl.gov/publications/download/iucrcompcomm_jan2003.pdf
#   http://cci.lbl.gov/publications/download/iucrcompcomm_jan2005.pdf

def examples():
  # Generate space groups (in matrix/vector form) based on spacegroup number
  # (names are *not* a pain)
  # See also: http://cctbx.sourceforge.net/current/c_plus_plus/classcctbx_1_1sgtbx_1_1space__group__symbols.html#_details
  from cctbx import sgtbx
  for s in sgtbx.space_group_info(symbol="I41/amd").group():
    print s # in "xyz" notation
    print s.r().as_rational().mathematica_form(), \
          s.t().as_rational().transpose().mathematica_form()
  print

  # now with a space group number
  space_group_info = sgtbx.space_group_info(number=123)
  space_group_info.show_summary()
  print

  # Generate conditions for allowed reflections etc
  from cctbx import crystal
  from cctbx import miller
  crystal_symmetry = crystal.symmetry(
    unit_cell=(10,10,13,90,90,90),
    space_group_info=space_group_info)
  miller_set = miller.build_set(
    crystal_symmetry=crystal_symmetry,
    anomalous_flag=False,
    d_min=4)
  # change the space group in order to get a few systematic absences
  miller_set = miller_set.customized_copy(
    space_group_info=sgtbx.space_group_info(symbol="I41/amd"))
  sys_absent_flags = miller_set.sys_absent_flags()
  for h,f in zip(sys_absent_flags.indices(), sys_absent_flags.data()):
    print h, f
  print
  # try also (from the command line): libtbx.help cctbx.miller

  # Generate point group of space group in matrix form
  point_group = miller_set.space_group().build_derived_point_group()
  point_group_info = sgtbx.space_group_info(group=point_group)
  point_group_info.show_summary()
  for s in point_group:
    print s
  print

  # Generate point-group matrices for given coordinate
  # first we have to define what we consider as special position
  special_position_settings = crystal.special_position_settings(
    crystal_symmetry=miller_set,
    min_distance_sym_equiv=0.5) # <<<<< here
  site_symmetry = special_position_settings.site_symmetry(
    site=(0,0.48,0))
  print "special position operator:", site_symmetry.special_op_simplified()
  print "distance to original site:", site_symmetry.distance_moved()
  print "point group of the special position:"
  for s in site_symmetry.matrices():
    print s
  print
  # See also: http://cci.lbl.gov/~rwgk/my_papers/iucr/au0265_reprint.pdf

  # Access database for form factors
  from cctbx.eltbx import xray_scattering
  si_form_factor = xray_scattering.it1992("Si")
  gaussians = si_form_factor.fetch()
  for stol in [0, 0.01, 0.02, 0.5]:
    print stol, gaussians.at_stol(stol)
  print

  # anomalous scattering factors: Sasaki tables
  from cctbx.eltbx import sasaki
  si_table = sasaki.table("Si")
  for wavelength in [0.5, 0.8, 0.9]:
    data = si_table.at_angstrom(wavelength)
    print wavelength, data.fp(), data.fdp()
  print

  # anomalous scattering factors: Henke tables
  from cctbx.eltbx import henke
  si_table = henke.table("Si")
  for wavelength in [0.5, 0.8, 0.9]:
    data = si_table.at_angstrom(wavelength)
    print wavelength, data.fp(), data.fdp()
  print
  print "OK"

if (__name__ == "__main__"):
  examples()
