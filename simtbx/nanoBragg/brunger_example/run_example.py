from __future__ import division
from __future__ import print_function
from scitbx.array_family import flex # import dependency
from simtbx.nanoBragg import shapetype
from simtbx.nanoBragg import convention
from simtbx.nanoBragg import nanoBragg
import libtbx.load_env # possibly implicit
from cctbx import crystal


# get the structure factor of spots
from iotbx.reflection_file_reader import any_reflection_file
mtz_file = any_reflection_file("./model_nophase.mtz")
#mtz_file = any_reflection_file("./shit.hkl=amplitudes")
#mtz_file = any_reflection_file("./justone.mtz")
Fhkl = mtz_file.as_miller_arrays()[0]

# get the structure factors of the background
stolfile = "./bg.stol"
Fbg_vs_stol = []
with open(stolfile, "rb") as fp:
  for i in fp.readlines():
    tmp = i.split(" ")
    try:
      Fbg_vs_stol.append((float(tmp[0]), float(tmp[1])))
    except Exception:pass
# now Fbg_vs_stol is a list of stol,Fbg tuples

# open the existing diffraction image: we need it for the background profile
import dxtbx
img = dxtbx.load("./F4_0_00008.mccd.gz")
panel = img.get_detector().to_dict()['panels'][0]
pixel_size_mm = panel['pixel_size'][0]
distance_mm = -panel['origin'][2]

# create the simulation
SIM = nanoBragg(detpixels_slowfast=(4096,4096),pixel_size_mm=0.079346,verbose=9)
SIM.Fhkl = Fhkl
SIM.Fbg_vs_stol = Fbg_vs_stol
SIM.close_distance_mm=299.83
SIM.wavelength_A=1.304735
SIM.polarization=0.99
SIM.beamsize_mm=0.03
#SIM.fluence=4.28889e+18
# fluence scaled to make crystal look bigger
SIM.fluence=1.03e+27
SIM.beamcenter_convention=convention.Custom
SIM.beam_center_mm=( 160.53, 182.31 )
SIM.dispersion_pct = 0.5
SIM.dispsteps=6
print("dispsteps=",SIM.dispsteps)
SIM.divergence_hv_mrad = ( 0.02, 0.02 )
SIM.divsteps_hv = ( 2 , 2 )
print(SIM.divsteps_hv)
SIM.round_div=True
print(SIM.divsteps_hv)
#SIM.detector_thick_mm = 0.037
SIM.detector_thick_mm = 0.
SIM.detector_thicksteps = 1
# override mtz unit cell
SIM.unit_cell_tuple = ( 68.78, 169.26, 287.42, 90, 90, 90 )
#SIM.Ncells_abc = ( 1, 1, 1 )
SIM.Ncells_abc = ( 14, 6, 4 )
#SIM.Ncells_abc = ( 35, 15, 10 )
print("Ncells_abc=",SIM.Ncells_abc)
SIM.xtal_shape=shapetype.Tophat
print("xtal_size_mm=",SIM.xtal_size_mm)
SIM.interpolate=0
SIM.progress_meter=True
SIM.mosaic_spread_deg = 0.2
SIM.mosaic_domains = 30
SIM.oversample = 1
SIM.detector_psf_type=shapetype.Fiber
SIM.adc_offset_adu = 10
SIM.readout_noise_adu = 1.5

# speedups, comment out for realism
#SIM.divergence_hv_mrad = ( 0,0 )
#SIM.dispersion_pct = 0
#SIM.mosaic_spread_deg = 0
#SIM.region_of_interest=((1450,1850),(1550,1950))
SIM.printout_pixel_fastslow=(1782,1832)
# set this to 0 or -1 to trigger automatic radius.  could be very slow with bright images
#SIM.detector_psf_kernel_radius_pixels=5;


SIM.amorphous_sample_thick_mm = 0.1
SIM.amorphous_density_gcm3 = 7e-7
SIM.amorphous_sample_molecular_weight_Da = 18 # default

print("dispsteps=",SIM.dispsteps)
print("divsteps=",SIM.divsteps_hv)
print("oversample=",SIM.oversample)
SIM.add_background(oversample=1)
print("mid_sample=",SIM.raw[1782,1832])
print("dispsteps=",SIM.dispsteps)
print("divsteps=",SIM.divsteps_hv)
print("oversample=",SIM.oversample)
SIM.to_smv_format(fileout="intimage_001.img",intfile_scale=1)

# three clusters of mosaic domains
SIM.fluence /= 3
SIM.missets_deg = ( 96.9473, -52.0932, -32.518 )
#SIM.missets_deg = ( 96.544, -51.9673, -32.4243 )
SIM.add_nanoBragg_spots()
SIM.to_smv_format(fileout="intimage_002.img",intfile_scale=1)
SIM.missets_deg = ( 97.5182, -52.3404, -32.7289 )
SIM.add_nanoBragg_spots()
SIM.to_smv_format(fileout="intimage_003.img",intfile_scale=1)
SIM.missets_deg = ( 97.1251, -52.2242, -32.751 )
SIM.add_nanoBragg_spots()
SIM.to_smv_format(fileout="intimage_004.img",intfile_scale=1)


SIM.detector_psf_fwhm_mm=0.08;
SIM.detector_psf_type=shapetype.Fiber

# get same noise each time this test is run
SIM.seed = 1
print("seed=",SIM.seed)
print("calib_seed=",SIM.calib_seed)
print("quantum_gain=",SIM.quantum_gain)
print("adc_offset_adu=",SIM.adc_offset_adu)
print("detector_calibration_noise_pct=",SIM.detector_calibration_noise_pct)
print("flicker_noise_pct=",SIM.flicker_noise_pct)
print("readout_noise_adu=",SIM.readout_noise_adu)
print("detector_psf_type=",SIM.detector_psf_type)
print("detector_psf_fwhm_mm=",SIM.detector_psf_fwhm_mm)
print("detector_psf_kernel_radius_pixels=",SIM.detector_psf_kernel_radius_pixels)
SIM.show_params()
SIM.add_noise()

print("raw=",SIM.raw)
SIM.to_smv_format(fileout="noiseimage_001.img",intfile_scale=1)

print("mosaic_domains=",SIM.mosaic_domains)
print("mosaic_spread_deg=",SIM.mosaic_spread_deg)
print("dispersion_pct=",SIM.dispersion_pct)
print("dispsteps=",SIM.dispsteps)
print("divergence_hv_mrad=",SIM.divergence_hv_mrad)
print("divergence_hv=",SIM.divsteps_hv)
