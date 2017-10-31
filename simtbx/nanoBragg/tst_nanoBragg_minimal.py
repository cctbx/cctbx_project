# absolute bare-minimum diffraction image simulation

from simtbx.nanoBragg import nanoBragg

# create the simulation object, all parameters have sensible defaults
SIM = nanoBragg()

# dont bother with importing a structure, we just want spots
SIM.default_F = 1

# default is one unit cell, lets to 125
SIM.Ncells_abc = (5,5,5)

# default orientation is with a axis down the beam, lets pick a random one
SIM.randomize_orientation()
# display randomly-picked missetting angles
print SIM.missets_deg
# or an Arndt-Wonacott A matrix (U*B), same as used by mosflm
print SIM.Amatrix

# show all parameters
SIM.show_params()

# now actually run the simulation
SIM.add_nanoBragg_spots()

# write out a file on arbitrary scale, header contains beam center in various conventions
SIM.to_smv_format(fileout="intimage_001.img")

# now apply a scale to get more than 1 photon/pixel and add noise
SIM.raw_pixels*=2000
SIM.add_noise()

SIM.to_smv_format(fileout="noiseimage_001.img")


if __name__=="__main__":
  print "OK"
