from __future__ import absolute_import, division, print_function
from rstbx.diffraction.fastbragg import detector,standard_camera,crystal_structure
from rstbx.simulation.constants import electron_radius
from rstbx.simulation import xfel1
from rstbx.simulation.sim_pdf import PDF
from rstbx.simulation.sim_utils import generate_random_rotation

# Instantiation of the Detector, array size & pixel edge in meters
Det = detector(slowpixels=1516,fastpixels=1516,pixel_size=0.00022)
Det.show()

Cam = standard_camera(detector = Det,
                      mean_xtd_distance_m = 0.150,
                      mean_xray_wavelength_m = 1.24E-10)
Cam.show()

Pdb = crystal_structure(pdb_code='2alp',standard_camera=Cam)
Pdb.show()

square_crystal_edge = 5.E-6 # in meters
square_focus_edge = 5.E-6 # in meters
vol_crystal = square_crystal_edge**3 # volume of the crystal in meters**3
I_beam = 2E12 #number of photons per LCLS pulse
I_beam_flux = I_beam/(square_focus_edge**2)
darwin_numerator = I_beam_flux * (electron_radius**2) * vol_crystal

class Params:pass
Sim = Params()
Sim.mosaicity = 0.10 #degrees
Sim.bandpass = 0.002
Sim.tracing_impacts = 1000
Sim.darwin_factor = darwin_numerator/(Pdb.p1_cell().volume()/1E30)#vol:meters**3

if __name__=="__main__":

    simulation1 = xfel1(detector=Det,camera=Cam,structure=Pdb,simulation=Sim)
    simulation1.show()
    plotter = PDF("./xsim1_%s_m%03d_xtd%03d_res%02d.pdf"%(Pdb.pdb_code,
      int(1000*Sim.mosaicity),
      int(Cam.distance),
      int(10*Pdb.limiting_resolution)
      ))

    #all indices for the asymmetric unit, and all intensities
    idx,intensity = (simulation1.indices_all, simulation1.intensities_all)

    for shot in generate_random_rotation(1):
      #shot = ((0.1551649112858198, 0.987395705748476, -0.031202092480033007),
      #        2.4190024383592474)

      simulation1.compute(axis = shot[0], angle = shot[1])

      ##### State information that is available for each image ####
      #miller indices that are near or at reflection condition on this shot
      simulation1.selected
      #raw counts for these selected miller indices (the partial recorded intensity)
      simulation1.signals
      #the fraction of spot recorded for these selected miller indices
      simulation1.partialities
      #the recorded image as flex<int> with flex.grid(2)
      simulation1.image
      #pixel coordinates of each spots (floating point x,y values; not sure if offset by 1/2 pxl)
      #not actually tested yet
      simulation1.spots

      plotter.make_image_plots_detail(simulation1)

"""
where to go from here 6/30/11:
apply point spread function to make things easier to see
test whether the formula for proximity is really inclusive enough
test mosaicity, bandpass
figure out what to do about contrast level; apply these things at the C++ level
make things consistent in terms of what things are called & exactly
where they are stored.  Use boost python injector?
"""
