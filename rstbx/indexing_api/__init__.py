from __future__ import absolute_import, division, print_function
from six.moves import range
import boost_adaptbx.boost.python as bp
import rstbx.dps_core # import dependency
bp.import_ext("rstbx_indexing_api_ext")
from rstbx_indexing_api_ext import *
import rstbx_indexing_api_ext as ext
from rstbx.array_family import flex
from scitbx.matrix import col
from rstbx_ext import * # gets us SpotClass

@bp.inject_into(ext.dps_extended)
class _():

  def set_beam_vector(self,beam):
    # input vector "beam" points from crystal to source.
    # S0 = -beam
    self.S0_vector = -beam
    self.inv_wave = self.S0_vector.length() # will be deprecated soon XXX

  def set_rotation_axis(self,axis):
    self.axis = axis
    assert axis.length() == 1.0

  def set_detector(self,input_detector):
    self.detector = input_detector

  @staticmethod
  def raw_spot_positions_mm_to_S1_vector( raw_spot_input, # as vec3_double
      detector, inverse_wave,
      panelID=None
      ):

    if panelID is None:
      panelID = flex.int(len(raw_spot_input),0)

    reciprocal_space_vectors = flex.vec3_double()

    # tile surface to laboratory transformation
    for n in range(len(raw_spot_input)):
      pid = panelID[n]
      lab_direct = col(detector[pid].get_lab_coord(raw_spot_input[n][0:2]))

    # laboratory direct to reciprocal space xyz transformation
      lab_recip = (lab_direct.normalize() * inverse_wave)

      reciprocal_space_vectors.append ( lab_recip )
    return reciprocal_space_vectors

  @staticmethod
  def raw_spot_positions_mm_to_reciprocal_space( raw_spot_input, # as vec3_double
      detector, inverse_wave, beam, axis, # beam, axis as scitbx.matrix.col
      panelID=None
      ):

    if panelID is None:
      panelID = flex.int(len(raw_spot_input),0)

    if axis is None:
      return raw_spot_positions_mm_to_reciprocal_space_xyz (
          raw_spot_input, detector, inverse_wave, beam, panelID )
    else:
      return raw_spot_positions_mm_to_reciprocal_space_xyz (
          raw_spot_input, detector, inverse_wave, beam, axis, panelID )

    """Assumptions:
    1) the raw_spot_input is in the same units of measure as the origin vector (mm).
       they are given in physical length, not pixel units
    2) the raw_spot centers of mass are given with the same corner/center convention
       as the origin vector.  E.g., spotfinder assumes that the mm scale starts in
       the middle of the lower-corner pixel.
    """

    reciprocal_space_vectors = flex.vec3_double()

    # tile surface to laboratory transformation
    for n in range(len(raw_spot_input)):
      pid = panelID[n]
      lab_direct = col(detector[pid].get_lab_coord(raw_spot_input[n][0:2]))

    # laboratory direct to reciprocal space xyz transformation
      lab_recip = (lab_direct.normalize() * inverse_wave) - beam

      reciprocal_space_vectors.append ( lab_recip.rotate_around_origin(
        axis=axis, angle=raw_spot_input[n][2], deg=True)
        )
    return reciprocal_space_vectors

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
    Hint = flex.vec3_double()
    for x in xyz:
      H = Astarinv * x
      #print "%7.3f %7.3f %7.3f"%(H[0],H[1],H[2])
      Hint.append((round(H[0],0), round(H[1],0), round(H[2],0)))
    xyz_miller = flex.vec3_double()
    from rstbx.diffraction import rotation_angles
    # XXX limiting shell of 1.0 angstroms probably needs to be changed/ removed.  How?
    ra = rotation_angles(limiting_resolution=1.0,orientation = Astar,
                         wavelength = 1./self.inv_wave, axial_direction = self.axis)
    for ij,hkl in enumerate(Hint):
      xyz_miller.append( Astar * hkl ) # figure out how to do this efficiently on vector data
      if ra(hkl):
        omegas = ra.get_intersection_angles()
        rotational_diffs = [ abs((-omegas[omegaidx] * 180./math.pi)-self.raw_spot_input[ij][2])
                             for omegaidx in [0,1] ]
        min_diff = min(rotational_diffs)
        min_index = rotational_diffs.index(min_diff)
        omega = omegas[min_index]
        rot_mat = self.axis.axis_and_angle_as_r3_rotation_matrix(omega)

        Svec = (rot_mat * Astar) * hkl + self.S0_vector
#        print panel.get_ray_intersection(Svec), self.raw_spot_input[ij]
        if self.panelID is not None: panel = self.detector[ self.panelID[ij] ]
        calc = matrix.col(panel.get_ray_intersection(Svec))
        pred = matrix.col(self.raw_spot_input[ij][0:2])
#        print (calc-pred).length(), separation_mm * TOLERANCE
        if ((calc-pred).length() < separation_mm * TOLERANCE):
          fraction_properly_predicted += 1./ self.raw_spot_input.size()
    #print "fraction properly predicted",fraction_properly_predicted,"with spot sep (mm)",separation_mm
    return fraction_properly_predicted

  def get_predicted_spot_positions_and_status(self, old_status=None): #similar to above function; above can be refactored
    panel = self.detector[0]
    from scitbx import matrix
    Astar = matrix.sqr(self.getOrientation().reciprocal_matrix())
    # must have set the basis in order to generate Astar matrix.  How to assert this has been done???

    import math,copy
    xyz = self.getXyzData()
    if old_status is None:
      spot_status = [SpotClass.GOOD]*len(xyz)
    else:
      assert len(old_status)==len(xyz)
      spot_status = copy.copy( old_status ) # valid way of copying enums w/o reference
    self.assigned_hkl= [(0,0,0)]*len(xyz)

    # step 1.  Deduce fractional HKL values from the XyzData.  start with x = A* h
    #          solve for h:  h = (A*^-1) x
    Astarinv = Astar.inverse()
    Hint = flex.vec3_double()
    results = flex.vec3_double()
    for x in xyz:
      H = Astarinv * x
      Hint.append((round(H[0],0), round(H[1],0), round(H[2],0)))
    xyz_miller = flex.vec3_double()
    from rstbx.diffraction import rotation_angles
    ra = rotation_angles(limiting_resolution=1.0,orientation = Astar,
                         wavelength = 1./self.inv_wave, axial_direction = self.axis)
    for ij,hkl in enumerate(Hint):
      xyz_miller.append( Astar * hkl ) # figure out how to do this efficiently on vector data
      if ra(hkl):
        omegas = ra.get_intersection_angles()
        rotational_diffs = [ abs((-omegas[omegaidx] * 180./math.pi)-self.raw_spot_input[ij][2])
                             for omegaidx in [0,1] ]
        min_diff = min(rotational_diffs)
        min_index = rotational_diffs.index(min_diff)
        omega = omegas[min_index]
        rot_mat = self.axis.axis_and_angle_as_r3_rotation_matrix(omega)

        Svec = (rot_mat * Astar) * hkl + self.S0_vector
        if self.panelID is not None: panel = self.detector[ self.panelID[ij] ]
        xy = panel.get_ray_intersection(Svec)
        results.append((xy[0],xy[1],0.0))
        self.assigned_hkl[ij]=hkl
      else:
        results.append((0.0,0.0,0.0))
        spot_status[ij]=SpotClass.NONE
    return results,spot_status

  def get_hkl(self,idx):
    return self.assigned_hkl[idx]
