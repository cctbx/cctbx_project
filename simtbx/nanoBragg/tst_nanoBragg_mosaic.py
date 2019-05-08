from __future__ import division, print_function
# test that mosaic spread is working
#from scitbx.array_family import flex
from simtbx.nanoBragg import shapetype
#from simtbx.nanoBragg import convention
from simtbx.nanoBragg import nanoBragg

# create the simulation object, all parameters have sensible defaults
SIM = nanoBragg(detpixels_slowfast=(150,150),pixel_size_mm=0.1,verbose=0)

# fairly big cell
SIM.unit_cell_tuple = (100,100,100,90,90,90)

# nice, sharp spots
SIM.Ncells_abc = (50,50,50)
SIM.xtal_shape=shapetype.Tophat

# dont bother with importing a structure, we just want spots
SIM.default_F = 1
SIM.F000 = 0.03

# the usual wavelength
SIM.wavelength_A = 1.0

# detector far away
SIM.distance_mm = 1000

# beam center near one edge
SIM.beam_center_mm = (1,7)
#print SIM.XDS_ORGXY
#print SIM.dials_origin_mm
print(SIM.beamcenter_convention)

# default orientation is with a axis down the beam, lets rotate it a bit
# tilt to Bragg angle
SIM.phi_deg = 89.7

# specify a broad enough mosaicity so we can measure it
SIM.mosaic_spread_deg = 10
# display randomly-picked missetting angles
SIM.mosaic_domains = 1000

# no need to over-do realism
SIM.interpolate = 0
SIM.oversample = 1

# show all parameters
#SIM.show_params()

# now actually run the simulation
SIM.add_nanoBragg_spots()

# check test points
test1 = (           SIM.raw_pixels[9,70] > 0 )
test1 = ( test1 or  SIM.raw_pixels[10,70] > 0 )
test1 = ( test1 or  SIM.raw_pixels[11,70] > 0 )
test1 = ( test1 or  SIM.raw_pixels[10,69] > 0 )
test1 = ( test1 or  SIM.raw_pixels[10,71] > 0 )
test1 = ( test1 and SIM.raw_pixels[8,70] == 0 )
test1 = ( test1 and SIM.raw_pixels[12,70] == 0 )
test1 = ( test1 and SIM.raw_pixels[10,68] == 0 )
test1 = ( test1 and SIM.raw_pixels[10,72] == 0 )

# left edge
test2 = (           SIM.raw_pixels[108,52] > 0 )
test2 = ( test2 or  SIM.raw_pixels[109,52] > 0 )
test2 = ( test2 and SIM.raw_pixels[108,50] == 0 )
test2 = ( test2 and SIM.raw_pixels[109,50] == 0 )

# right edge
test3 = (          SIM.raw_pixels[108,88] > 0 )
test3 = ( test3 or SIM.raw_pixels[109,88] > 0 )
test3 = ( test3 and SIM.raw_pixels[108,90] == 0 )
test3 = ( test3 and SIM.raw_pixels[109,90] == 0 )

# all together now
test = ( test1 and test2 and test3 )

# write out a file on arbitrary scale, header contains beam center in various conventions
#SIM.to_smv_format(fileout="10deg_mosaic_spread.img")

if not test:
  print("test failed")

if test:
  print("OK")

