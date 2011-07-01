import boost.python
ext = boost.python.import_ext("rstbx_simulation_ext")
from rstbx_simulation_ext import *
from scitbx import matrix
from cctbx import crystal_orientation
from rstbx.simulation.constants import eV_per_inv_meter
import math

class xfel1(ext.xfel1):

  def show(self):
    print "In XFEL1 simulation"

  def __init__(self,detector=None,camera=None,structure=None,simulation=None):
    ext.xfel1.__init__(self)

    self.detector=detector
    self.camera=camera
    self.structure=structure
    self.sim=simulation

    self.uc = self.structure.p1_cell()

    self.bmat = matrix.sqr(self.uc.orthogonalization_matrix()).inverse().transpose()
    self.Ori = crystal_orientation.crystal_orientation( self.bmat,
      crystal_orientation.basis_type.reciprocal )

    energy_in_eV = eV_per_inv_meter / (self.camera.lambda0) # lambda in meters

    # top hat function for dispersion
    self.full_pass_eV = energy_in_eV * matrix.col([1.-(self.sim.bandpass/2.),
                                              1.+(self.sim.bandpass/2.)])

    self.full_pass_lambda = eV_per_inv_meter * matrix.col((1./self.full_pass_eV[0],
                                                           1./self.full_pass_eV[1]))

    intensities = self.structure.p1_intensities()
    self.set_indices(intensities.indices())
    self.set_intensities(intensities.data())

  def compute(self,axis,angle):
    shotOri = self.Ori.rotate_thru(axis,angle)
    beam_vector_B = matrix.col([0, 0, 1 / (self.camera.lambda0)])
    self.selected = self.select_proximal_indices(half_edge=int(self.detector.raw.focus()[0]/2),
                 detector_distance_m=self.camera.distance,
                 pixel_size_m=self.detector.pixel_sz,
                 orientation = shotOri,
                 mosaicity_full_width  = self.sim.mosaicity*(math.pi/180.),
                 bandpass_full_width = self.sim.bandpass,
                 wavelength_m = self.camera.lambda0,
                 limiting_resolution_Ang = self.structure.limiting_resolution)
    print "Selected %d reflections close to the Ewald sphere in this orientation"%len(self.selected)

    self.detector.raw.fill(0.0)
    self.image = self.raw_diffraction(selection = self.selected,
                            pixels = self.detector.raw,
                            mosaic_domains = self.sim.tracing_impacts,
                            detector_distance_m=self.camera.distance,
                            pixel_size_m=self.detector.pixel_sz,
                            darwin_factor = self.sim.darwin_factor)
