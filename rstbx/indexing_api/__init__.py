from __future__ import division
import boost.python
boost.python.import_ext("rstbx_indexing_api_ext")
from rstbx_indexing_api_ext import *
import rstbx_indexing_api_ext as ext

class _(boost.python.injector, ext.dps_extended):

  def set_beam_vector(self,beam):
    self.beam = beam
    self.beam_vector = beam # will be deprecated soon XXX
    self.inv_wave = self.beam.length() # will be deprecated soon XXX
    self.wavelength_set = 1./self.inv_wave

  def set_rotation_axis(self,axis):
    self.rotation_vector = axis
    self.axis = axis # will be deprecated soon XXX
    assert axis.length() == 1.0

  def set_detector(self,input_detector):
    from scitbx.matrix import col
    self.origin = col(input_detector.get_origin())
    self.d1 = col(input_detector.get_fast_axis())
    self.d2 = col(input_detector.get_slow_axis())
    self.detector = input_detector

  def set_detector_position(self,origin,d1,d2): # optional, alternate form
    from dxtbx.model.detector import detector_factory
    self.detector = detector_factory.make_detector(
      stype = "indexing",
      fast_axis = d1,
      slow_axis = d2,
      origin = origin,
      pixel_size = (1.0,1.0),  #not actually using pixels for indexing
      image_size = (100,100),  #not using pixels
      )
    self.origin = origin
    self.d1 = d1 # detector fast axis
    self.d2 = d2 # detector slow axis

  #def raw_spot_input_to_reciprocal_space( # to be implemented

  def model_likelihood(self,separation_mm):
    TOLERANCE = 0.5
    fraction_properly_predicted = 0.0

    #help(self.detector)
    #print self.detector[0]
    #help(self.detector[0])
    panel = self.detector[0]
    from scitbx import matrix
    Astar = matrix.sqr(self.getOrientation().reciprocal_matrix())

    import math
    xyz = self.getXyzData()

    # step 1.  Deduce fractional HKL values from the XyzData.  start with x = A* h
    #          solve for h:  h = (A*^-1) x
    Astarinv = Astar.inverse()
    from cctbx.array_family import flex
    Hint = flex.vec3_double()
    for x in xyz:
      H = Astarinv * x
      #print "%7.3f %7.3f %7.3f"%(H[0],H[1],H[2])
      Hint.append((round(H[0],0), round(H[1],0), round(H[2],0)))
    xyz_miller = flex.vec3_double()
    from rstbx.diffraction import rotation_angles
    ra = rotation_angles(limiting_resolution=1.0,orientation = Astar,
                         wavelength = self.wavelength_set, axial_direction = self.rotation_vector)
    for ij,hkl in enumerate(Hint):
      xyz_miller.append( Astar * hkl ) # figure out how to do this efficiently on vector data
      if ra(hkl):
        omegas = ra.get_intersection_angles()
        rotational_diffs = [ abs((-omegas[omegaidx] * 180./math.pi)-self.raw_spot_input[ij][2])
                             for omegaidx in [0,1] ]
        min_diff = min(rotational_diffs)
        min_index = rotational_diffs.index(min_diff)
        omega = omegas[min_index]
        rot_mat = self.rotation_vector.axis_and_angle_as_r3_rotation_matrix(omega)

        Svec = (rot_mat * Astar) * hkl + self.beam_vector
#        print panel.get_ray_intersection(Svec), self.raw_spot_input[ij]
        calc = matrix.col(panel.get_ray_intersection(Svec))
        pred = matrix.col(self.raw_spot_input[ij][0:2])
#        print (calc-pred).length(), separation_mm * TOLERANCE
        if ((calc-pred).length() < separation_mm * TOLERANCE):
          fraction_properly_predicted += 1./ self.raw_spot_input.size()
    print "fraction properly predicted",fraction_properly_predicted,"with spot sep (mm)",separation_mm
    return fraction_properly_predicted
