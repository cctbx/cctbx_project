from __future__ import division
from scitbx.array_family import flex
from scitbx.matrix import sqr
from simtbx.nanoBragg import testuple
from simtbx.nanoBragg import pivot
from simtbx.nanoBragg import shapetype
from simtbx.nanoBragg import convention
from simtbx.nanoBragg import nanoBragg
import libtbx.load_env # possibly implicit
from cctbx import crystal
from cctbx import miller

pdb_lines = """HEADER TEST
CRYST1   50.000  100.000  150.000  90.00  90.00  90.00 P 1
ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00 20.00           O
ATOM      1 SED  MSE A   1       1.000   2.000   3.000  1.00 20.00          SE
END
"""

def fcalc_from_pdb(resolution,algorithm=None,wavelength=0.9):
  from iotbx import pdb
  pdb_inp = pdb.input(source_info=None,lines = pdb_lines)
  xray_structure = pdb_inp.xray_structure_simple()
  #
  # take a detour to insist on calculating anomalous contribution of every atom
  scatterers = xray_structure.scatterers()
  for sc in scatterers:
    from cctbx.eltbx import sasaki, henke
    #expected_sasaki = sasaki.table(sc.element_symbol()).at_angstrom(wavelength)
    expected_henke = henke.table(sc.element_symbol()).at_angstrom(wavelength)
    sc.fp = expected_henke.fp()
    sc.fdp = expected_henke.fdp()
  # how do we do bulk solvent?
  primitive_xray_structure = xray_structure.primitive_setting()
  P1_primitive_xray_structure = primitive_xray_structure.expand_to_p1()
  fcalc = P1_primitive_xray_structure.structure_factors(
    d_min=resolution, anomalous_flag=True, algorithm=algorithm).f_calc()
  return fcalc.amplitudes()


def run_sim2smv(fileout):
  SIM = nanoBragg(detpixels_slowfast=(1000,1000),Ncells_abc=(5,5,5),verbose=99)
  SIM.pixel_size_mm=0.1
  SIM.distance_mm=100
  SIM.oversample=1
  SIM.wavelength_A=1
  SIM.unit_cell_tuple=(50,50,50,90,90,90)
  SIM.default_F=100
  #SIM.missets_deg = (10,20,30)
  #sfall = fcalc_from_pdb(resolution=2.,algorithm="direct",wavelength=SIM.wavelength_A)
  #
  #SIM.Fhkl=sfall
  SIM.xtal_shape=shapetype.Tophat
  SIM.progress_meter=False
  #SIM.printout_pixel_fastslow=(0,0)
  #SIM.printout=True
  SIM.show_params()
  SIM.add_nanoBragg_spots()
  SIM.to_smv_format(fileout="intimage_001.img")
  SIM.verbose=99
  # rough approximation to water: interpolation points for sin(theta/lambda) vs structure factor
  bg = flex.vec2_double([(0,2.57),(0.0365,2.58),(0.07,2.8),(0.12,5),(0.162,8),(0.2,6.75),(0.18,7.32),(0.216,6.75),(0.236,6.5),(0.28,4.5),(0.3,4.3),(0.345,4.36),(0.436,3.77),(0.5,3.17)])
  SIM.Fbg_vs_stol = bg
  SIM.amorphous_sample_thick_mm = 0.1
  SIM.amorphous_density_gcm3 = 1
  SIM.amorphous_molecular_weight_Da = 18
  SIM.flux=1e12
  SIM.beamsize_mm=0.1
  SIM.exposure_s=0.1
  SIM.add_background()
  SIM.to_smv_format(fileout="intimage_002.img")
  # rough approximation to air
  bg = flex.vec2_double([(0,14.1),(0.045,13.5),(0.174,8.35),(0.35,4.78),(0.5,4.22)])
  SIM.Fbg_vs_stol = bg
  SIM.amorphous_sample_thick_mm = 35 # between beamstop and collimator
  SIM.amorphous_density_gcm3 = 1.2e-3
  SIM.amorphous_sample_molecular_weight_Da = 28 # nitrogen = N2
  SIM.add_background()
  SIM.detector_psf_kernel_radius_pixels=5;
  SIM.detector_psf_fwhm_mm=0.2;
  SIM.detector_psf_type=shapetype.Fiber
  SIM.apply_psf()
  print SIM.raw[500000]
  SIM.to_smv_format(fileout="intimage_003.img")
  SIM.detector_psf_fwhm_mm=0
  SIM.add_noise()

  #fileout = "intimage_001.img"
  SIM.to_smv_format(fileout="noiseimage_001.img",intfile_scale=1)


#run_sim2smv("intimage_001.img")


def tst_all():
  F = testuple()
  assert F == (1,2,3,4)
  #
  fileout = "intimage_001.img"
  run_sim2smv(fileout)
  import os
  assert os.path.isfile(fileout)

  exit()
  #simulation is complete, now we'll autoindex the image fragment and verify
  # that the indexed cell is similar to the input cell.

  if (not libtbx.env.has_module("annlib")):
    print "Skipping some tests: annlib not available."
    return
  # 1. Analysis of the image to identify the Bragg peak centers.
  #    This step uses an inefficient algorithm and implementation and
  #    is most time consuming; but the code is only for testing, not production
  from rstbx.diffraction.fastbragg.tst_utils_clustering import specific_libann_cluster
  M=specific_libann_cluster(scale_factor*data,intensity_cutoff = 25,distance_cutoff=17)
  # M is a dictionary of peak intensities indexed by pixel coordinates

  # 2. Now autoindex the pattern
  from rstbx.diffraction.fastbragg.tst_utils_clustering import index_wrapper
  SIM.C = C
  SIM.D = D
  ai,ref_uc = index_wrapper(M.keys(), SIM, PSII)
  tst_uc = ai.getOrientation().unit_cell()
  #print ref_uc  # (127.692, 225.403, 306.106, 90, 90, 90)
  #print tst_uc  # (106.432, 223.983, 303.102, 90.3185, 91.5998, 90.5231)

  # 3. Final assertion.  In the given orientation,
  #  the unit cell A vector is into the beam and is not well sampled,
  #  so tolerances have to be fairly relaxed, 5%.
  #  Labelit does better with the target_cell restraint, but this improved
  #  algorithm is not used here for the test
  assert ref_uc.is_similar_to(tst_uc, relative_length_tolerance=0.20,
                                     absolute_angle_tolerance= 2.0)


if __name__=="__main__":
  tst_all()
  print "OK"
