from __future__ import division
from __future__ import print_function
from scitbx.array_family import flex
from simtbx.nanoBragg import nanoBragg, shapetype, convention, pivot, testuple
#from simtbx.nanoBragg import shapetype
#from simtbx.nanoBragg import pivot
#from simtbx.nanoBragg import convention
#from simtbx.nanoBragg import nanoBragg
import dxtbx
from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus
from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus
import libtbx.load_env # possibly implicit
from cctbx import crystal
from cctbx import miller
assert miller

pdb_lines = """HEADER TEST
CRYST1   50.000   60.000   70.000  90.00  90.00  90.00 P 21 21 21
ATOM      1  O   HOH A   1      56.829   2.920  55.702  1.00 20.00           O
ATOM      2  O   HOH A   2      49.515  35.149  37.665  1.00 20.00           O
ATOM      3  O   HOH A   3      52.667  17.794  69.925  1.00 20.00           O
ATOM      4  O   HOH A   4      40.986  20.409  18.309  1.00 20.00           O
ATOM      5  O   HOH A   5      46.896  37.790  41.629  1.00 20.00           O
ATOM      6 SED  MSE A   6       1.000   2.000   3.000  1.00 20.00          SE
END
"""

# traces of F vs sin(theta)/lambda for various materials
# rough approximation to water: interpolation points for sin(theta/lambda) vs structure factor
Fbg_water = flex.vec2_double([(0,2.57),(0.0365,2.58),(0.07,2.8),(0.12,5),(0.162,8),(0.2,6.75),(0.18,7.32),(0.216,6.75),(0.236,6.5),(0.28,4.5),(0.3,4.3),(0.345,4.36),(0.436,3.77),(0.5,3.17)])
# rough approximation to air
Fbg_air = flex.vec2_double([(0,14.1),(0.045,13.5),(0.174,8.35),(0.35,4.78),(0.5,4.22)])
# structure factors of nanocrystalline ice Ic
Fbg_nanoice = flex.vec2_double([(0.0,0.677252),(0.0707968,1.43379),(0.0981619,2.3679),(0.114124,3.75888),(0.120588,5.12927),(0.127821,9.97093),(0.129198,11.9564),(0.132627,21.7772),(0.134056,27.7952),(0.13532,31.3918),(0.137142,30.3891),(0.138649,24.5024),(0.142137,13.6271),(0.144115,10.6367),(0.146059,8.67322),(0.147749,7.72956),(0.156921,5.57316),(0.166321,4.98756),(0.189445,4.45263),(0.210479,5.46593),(0.212578,5.78441),(0.220728,12.9117),(0.221969,15.9529),(0.224001,13.2016),(0.225279,10.8369),(0.226852,9.07622),(0.232851,7.00939),(0.242413,6.05003),(0.250842,5.50808),(0.257907,7.96463),(0.259429,10.0142),(0.262963,6.37791),(0.269814,3.89966),(0.271936,3.76553),(0.326095,3.71447),(0.328613,3.77097),(0.340562,4.56185),(0.343282,5.10066),(0.347891,4.34943),(0.349662,4.39573),(0.450378,2.30382),(0.5,1.31048)])

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


def run_sim(seed=1,wavelength=0.9,distance=500,random_orientation=False,phi=0,osc=0):
  SIM = nanoBragg(detpixels_slowfast=(2527,2463),pixel_size_mm=0.172,Ncells_abc=(15,15,15),verbose=0)
  # get same noise each time this test is run
  SIM.seed = seed
  # user may select random orientation
  if random_orientation :
    SIM.randomize_orientation()
  else:
    # fixed orientation, for doing rotation data sets
    SIM.missets_deg= (10,20,30)
  SIM.distance_mm = distance
  SIM.wavelength_A = wavelength
  # option to render a phi swell
  SIM.phi_deg = phi
  SIM.osc_deg = osc
  if osc >= 0 :
    # need quite a few steps/deg for decent Rmeas
    SIM.phistep_deg = 0.001
  SIM.oversample = 1
  SIM.polarization = 1
  # this will become F000, marking the beam center
  SIM.F000 = 100
  print("mosaic_seed=",SIM.mosaic_seed)
  print("seed=",SIM.seed)
  print("calib_seed=",SIM.calib_seed)
  print("missets_deg =", SIM.missets_deg)
  # re-set the detector to be in
  SIM.beamcenter_convention=convention.DIALS
  SIM.beam_center_mm=(211,214)
  sfall = fcalc_from_pdb(resolution=1.6,algorithm="direct",wavelength=SIM.wavelength_A)
  # use crystal structure to initialize Fhkl array
  SIM.Fhkl=sfall
  # fastest option, least realistic
  SIM.xtal_shape=shapetype.Tophat
  # only really useful for long runs
  SIM.progress_meter=False
  # prints out value of one pixel only.  will not render full image!
  #SIM.printout_pixel_fastslow=(500,500)
  #SIM.printout=True
  # flux is always in photons/s
  SIM.flux=1e12
  # assumes round beam
  SIM.beamsize_mm=0.1
  SIM.exposure_s=0.1
  # Pilatus detector is thick
  SIM.detector_thick_mm=0.450
  # get attenuation depth in mm
  from cctbx.eltbx import attenuation_coefficient
  table = attenuation_coefficient.get_table("Si")
  SIM.detector_attenuation_length_mm = 10.0/(table.mu_at_angstrom(SIM.wavelength_A))
  print("detector_attenuation_length =",SIM.detector_attenuation_length_mm)
  #SIM.detector_thick_mm=0
  SIM.detector_thicksteps=3
  #
  # will simulate module mis-alignment by adjusting beam center of each module
  beam_center0 = SIM.beam_center_mm
  # simulated crystal is only 3375 unit cells (42 nm wide)
  # amplify spot signal to simulate physical crystal of 100 um
  # this is equivalent to adding up 28e9 mosaic domains in exactly the same orientation, but a lot faster
  aggregate_xtal_volume_mm3 = 0.1*0.1*0.1
  mosaic_domain_volume_mm3 = SIM.xtal_size_mm[0]*SIM.xtal_size_mm[1]*SIM.xtal_size_mm[2]
  SIM.spot_scale = aggregate_xtal_volume_mm3/mosaic_domain_volume_mm3
  #
  # option for detector rotations to be about beam center or sample position
  #SIM.detector_pivot=pivot.Sample
  SIM.detector_pivot=pivot.Beam
  print("pivoting detector about",SIM.detector_pivot)
  # make the module roi list
  # and also make up per-module mis-alignments shifts and rotations
  import random
  # make sure shifts are the same from image to image
  random.seed(12345)
  module_size = (487, 195)
  gap_size = (7, 17)
  modules = []
  for mx in range(0,5):
   for my in range(0,12):
    # define region of interest for this module
    fmin = mx*(module_size[0]+gap_size[0])
    fmax = fmin+module_size[0]
    smin = my*(module_size[1]+gap_size[1])
    smax = smin+module_size[1]
    # assume misalignments are uniform
    # up to ~0.5 pixel and 0.0 deg in each direction
    dx = 0.07*(random.random()-0.5)*2
    dy = 0.07*(random.random()-0.5)*2
    rotx = 0.0*(random.random()-0.5)*2;
    roty = 0.0*(random.random()-0.5)*2;
    rotz = 0.0*(random.random()-0.5)*2;
    modules.append((((fmin,fmax),(smin,smax)),(dx,dy),(rotx,roty,rotz)))
  #
  # defeat module-by-module rendering by uncommenting next 3 lines
  #dim=SIM.detpixels_fastslow
  #modules = [(( ((0,dim[0]),(0,dim[1])) ),(0,0),(0,0,0))]
  #print modules
  i=-1
  for roi,(dx,dy),drot in modules:
    i=i+1
    print("rendering module:",i,"roi=",roi)
    SIM.region_of_interest = roi
    # shift the beam center
    x = beam_center0[0]+dx
    y = beam_center0[1]+dy
    SIM.beam_center_mm = (x,y)
    # also tilt the module by up to 0.05 deg in all directions
    SIM.detector_rotation_deg = drot
    print("beam center: ",SIM.beam_center_mm,"rotations:",drot)
    # now actually burn up some CPU
    SIM.add_nanoBragg_spots()
    # add water contribution
    SIM.Fbg_vs_stol = Fbg_water
    SIM.amorphous_sample_thick_mm = 0.1
    SIM.amorphous_density_gcm3 = 1
    SIM.amorphous_molecular_weight_Da = 18
    SIM.add_background()
    # add air contribution
    SIM.Fbg_vs_stol = Fbg_air
    SIM.amorphous_sample_thick_mm = 35 # between beamstop and collimator
    SIM.amorphous_density_gcm3 = 1.2e-3
    SIM.amorphous_sample_molecular_weight_Da = 28 # nitrogen = N2
    SIM.add_background()
    # add nanocrystalline cubic ice contribution (unscaled)
    SIM.Fbg_vs_stol = Fbg_nanoice
    SIM.amorphous_sample_thick_mm = 0.05 # between beamstop and collimator
    SIM.amorphous_density_gcm3 = 0.95
    SIM.amorphous_sample_molecular_weight_Da = 18 # H2O
    SIM.add_background()
  #
  # set this to 0 or -1 to trigger automatic radius.  could be very slow with bright images
  SIM.detector_psf_kernel_radius_pixels=5;
  # turn off the point spread function for Pilatus
  SIM.detector_psf_fwhm_mm=0.0;
  SIM.detector_psf_type=shapetype.Gauss
  SIM.adc_offset_adu=0
  SIM.readout_noise_adu=0
  dim=SIM.detpixels_fastslow
  SIM.region_of_interest = ( (0,dim[0]),(0,dim[1]) )
  SIM.add_noise()
  return SIM


# function for setting bad pixels
def mask_pixels(raw_pixels):
  # randomly select bad pixels
  import random
  # make sure bad pixel mask is always the same
  random.seed(12345)
  maxy,maxx = raw_pixels.all()
  for i in range(1,600):
    # assume any pixel is equally likely to be bad
    x,y = maxx,maxy
    x *= random.random()
    y *= random.random()
    # use relative position on pixel below
    ix,rx = int(x),x%1
    iy,ry = int(y),y%1
    # avoid out-of-range errors
    ix = ( ix if ix<=maxx-2 else maxx-2 )
    iy = ( iy if iy<=maxy-2 else maxy-2 )
    # do not put -2 in the missing-pixels areas
    if raw_pixels[iy,ix] == -1 :
      continue
    # flag selected pixel as bad
    raw_pixels[iy,ix]=-2
    # allow for sets of up to four baddies in a bunch
    if rx>0.5 :
      raw_pixels[iy,ix+1]=-2
    if ry>0.5 :
      raw_pixels[iy+1,ix]=-2
    if ( rx*ry > 0.25 ):
      raw_pixels[iy+1,ix+1]=-2
  #
  # define the mask (f0, f1, s0, s1)
  from dxtbx.format.FormatPilatusHelpers import pilatus_6M_mask
  for pane in pilatus_6M_mask():
    print("negating pane:",pane)
    fs=pane[0]-1
    fe=pane[1]
    ss=pane[2]-1
    se=pane[3]
    # default code for missing pixel is -1, will overwrite anything else
    for y in range(ss,se):
      for x in range(fs,fe):
        raw_pixels[y,x]=-1
  # return the whole sim object
  return raw_pixels

def write_cbf(SIM,fileout="noiseimage_00001.cbf"):
  # first and foremost, convert the pixel data to integers
  dxData = SIM.raw_pixels.iround()
  # now get the image writing goodies
  import dxtbx
  from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus
  # apparently, cant initialize a detector witout doing a panel first
  dxPanel = dxtbx.model.Panel()
  dxPanel.set_name('Panel')
  dxPanel.set_type('SENSOR_PAD')
  dxPanel.set_image_size(SIM.detpixels_fastslow)
  dxPanel.set_pixel_size((SIM.pixel_size_mm,SIM.pixel_size_mm))
  dxPanel.set_gain(SIM.quantum_gain)
  dxPanel.set_trusted_range((-1,1.3e6))
  dxPanel.set_frame(SIM.fdet_vector,SIM.sdet_vector,SIM.dials_origin_mm)
  dxPanel.set_material('Si')
  dxPanel.set_thickness(SIM.detector_thick_mm)
  dxPanel.set_mu(SIM.detector_attenuation_length_mm*10)
  # use this panel to create a detector
  dxDet = dxtbx.model.Detector(dxPanel)
  # now create a beam from nothing
  dxBeam = dxtbx.model.Beam()
  #dxBeam.set_direction(tuple(flex.double(SIM.beam_vector)*-1))
  dxBeam.set_s0(SIM.beam_vector)
  dxBeam.set_wavelength(SIM.wavelength_A)
  dxBeam.set_flux(SIM.flux)
  dxBeam.set_divergence(SIM.divergence_hv_mrad[0])
  dxBeam.set_polarization_normal(tuple(flex.double(SIM.polar_Bvector)*-1))
  dxBeam.set_polarization_fraction(SIM.polarization)
  # or cheat and choose first internal sub-beam
  #dxBeam = SIM.xray_beams[0]
  # create a goniomater from nothing
  dxGoni = dxtbx.model.Goniometer()
  dxGoni.set_rotation_axis(SIM.spindle_axis)
  # define the scan from nothing
  dxScan = dxtbx.model.Scan()
  dxScan.set_image_range((1,1));
  dxScan.set_oscillation([SIM.phi_deg,SIM.osc_deg]);
  dxScan.set_exposure_times([SIM.exposure_s]);
  # now, finally, do the image write
  FormatCBFMiniPilatus.as_file(
    detector=dxDet,beam=dxBeam,gonio=dxGoni,scan=dxScan,
    data=dxData,path=fileout)


def tst_all():
  # make sure NAN stuff is working
  F = testuple()
  assert F == (1,2,3,4)
  #
  # image number
  img_num = 1
  # check for any command-line options
  import sys
  if len(sys.argv)>1:
    img_num = int(sys.argv[1])
  num = format(img_num,"05d")
  fileout = "noiseimage_"+num+".cbf"
  SIM=run_sim(seed=img_num,phi=(img_num-1)/10,osc=0.1)
  # mark the window panes
  mask_pixels(SIM.raw_pixels)
  # output the file
  write_cbf(SIM,fileout)

  # how do you do a destructor in Python?
  SIM.free_all()
  import os
  assert os.path.isfile(fileout)

  exit()
  #simulation is complete, now we'll autoindex the image fragment and verify
  # that the indexed cell is similar to the input cell.


if __name__=="__main__":
  tst_all()
  print("OK")
