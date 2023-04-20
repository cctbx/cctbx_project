from __future__ import division

def basis_from_geo(geo, use_z = True):
  """ Given a psana GeometryObject, construct a basis object """
  rotx = matrix.col((1,0,0)).axis_and_angle_as_r3_rotation_matrix(
    geo.rot_x + geo.tilt_x, deg=True)
  roty = matrix.col((0,1,0)).axis_and_angle_as_r3_rotation_matrix(
    geo.rot_y + geo.tilt_y, deg=True)
  rotz = matrix.col((0,0,1)).axis_and_angle_as_r3_rotation_matrix(
    geo.rot_z + geo.tilt_z, deg=True)

  rot = (rotx*roty*rotz).r3_rotation_matrix_as_unit_quaternion()

  if use_z:
    trans = matrix.col((geo.x0/1000, geo.y0/1000, geo.z0/1000))
  else:
    trans = matrix.col((geo.x0/1000, geo.y0/1000, 0))

  return basis(orientation = rot, translation = trans)

def evt_wavelength(evt, delta_k=0):
  """The evt_wavelength() function returns the wavelength in Ångström
  of the event pointed to by @p evt.  From Margaritondo & Rebernik
  Ribic (2011): the dimensionless relativistic γ-factor is derived
  from beam energy in MeV and the electron rest mass, K is a
  dimensionless "undulator parameter", and L is the macroscopic
  undulator period in Ångström.  See also
  https://people.eecs.berkeley.edu/~attwood/srms/2007/Lec10.pdf

  @param evt     Event data object, a configure object
  @param delta_k Optional K-value correction
  @return        Wavelength, in Ångström
  """

  if evt is not None:
    ebeam = get_ebeam(evt)

    if hasattr(ebeam, 'fEbeamPhotonEnergy') and ebeam.fEbeamPhotonEnergy > 0:
      # pyana
      return 12398.4187 / ebeam.fEbeamPhotonEnergy
    if hasattr(ebeam, 'ebeamPhotonEnergy') and ebeam.ebeamPhotonEnergy() > 0:
      # psana
      return 12398.4187 / ebeam.ebeamPhotonEnergy()

    if hasattr(ebeam, 'fEbeamL3Energy') and ebeam.fEbeamL3Energy > 0:
      # pyana
      gamma = ebeam.fEbeamL3Energy / 0.510998910
    elif hasattr(ebeam, 'ebeamL3Energy') and ebeam.ebeamL3Energy() > 0:
      # psana
      gamma = ebeam.ebeamL3Energy() / 0.510998910
    else:
      return None
    K = 3.5 + delta_k
    L = 3.0e8
    return L / (2 * gamma**2) * (1 + K**2 / 2)
  return None

def env_distance(address, env, offset):
  """The env_distance() function returns the distance between the
  sample and the detector with the given address string in mm.  The
  distance between the sample and the the detector's zero-point can
  vary by an inch or more between different LCLS runs.  According to
  Sébastien Boutet the offset should be stable to within ±0.5 mm
  during a normal experiment.

  @param address Full data source address of the DAQ device
  @param env     Environment object
  @param offset  Detector-sample offset in mm, corresponding to
                 longest detector-sample distance
  @return        Detector-sample distance, in mm
  """

  detz = env_detz(address, env)
  if detz is not None:
    return detz + offset
  return None
