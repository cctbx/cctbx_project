from __future__ import division

from scitbx import matrix
from serialtbx.detector import basis
from cctbx import factor_ev_angstrom

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

def old_address_to_new_address(address):
  """ Change between old and new style detector addresses.
  I.E. CxiDs1-0|Cspad-0 becomes CxiDs1.0:Cspad.0
  @param address detector address to convert
  """
  return address.replace('-','.').replace('|',':')

def address_split(address, env=None):
  """The address_split() function splits an address into its four
  components.  Address strings are on the form
  detector-detectorID|device-deviceID, where the detectors must be in
  dir(xtc.DetInfo.Detector) and device must be in
  (xtc.DetInfo.Device).
  @param address Full data source address of the DAQ device
  @param env     Optional env to dereference an alias into an address
  @return        Four-tuple of detector name, detector ID, device, and
                 device ID
  """

  import re

  # pyana
  m = re.match(
    r"^(?P<det>\S+)\-(?P<det_id>\d+)\|(?P<dev>\S+)\-(?P<dev_id>\d+)$", address)
  if m is not None:
    return (m.group('det'), m.group('det_id'), m.group('dev'), m.group('dev_id'))

  # psana
  m = re.match(
    r"^(?P<det>\S+)\.(?P<det_id>\d+)\:(?P<dev>\S+)\.(?P<dev_id>\d+)$", address)
  if m is not None:
    return (m.group('det'), m.group('det_id'), m.group('dev'), m.group('dev_id'))

  # psana DetInfo string
  m = re.match(
    r"^DetInfo\((?P<det>\S+)\.(?P<det_id>\d+)\:(?P<dev>\S+)\.(?P<dev_id>\d+)\)$", address)
  if m is not None:
    return (m.group('det'), m.group('det_id'), m.group('dev'), m.group('dev_id'))

  if env is not None:
    # Try to see if this is a detector alias, and if so, dereference it. Code from psana's Detector/PyDetector.py
    amap = env.aliasMap()
    alias_src = amap.src(address) # string --> DAQ-style psana.Src

    # if it is an alias, look up the full name
    if amap.alias(alias_src) != '':         # alias found
      address = str(alias_src)
      return address_split(address)

  return (None, None, None, None)

def get_ebeam(evt):
  try:
    # pyana
    ebeam = evt.getEBeam()
  except AttributeError:
    try:
      from psana import Source, Bld
      src = Source('BldInfo(EBeam)')
      ebeam = evt.get(Bld.BldDataEBeamV6, src)
      if ebeam is None:
        ebeam = evt.get(Bld.BldDataEBeamV5, src)
      if ebeam is None:
        ebeam = evt.get(Bld.BldDataEBeamV4, src)
      if ebeam is None:
        ebeam = evt.get(Bld.BldDataEBeamV3, src)
      if ebeam is None:
        ebeam = evt.get(Bld.BldDataEBeamV2, src)
      if ebeam is None:
        ebeam = evt.get(Bld.BldDataEBeamV1, src)
      if ebeam is None:
        ebeam = evt.get(Bld.BldDataEBeamV0, src)
      if ebeam is None:
        ebeam = evt.get(Bld.BldDataEBeam, src) # recent version of psana will return a V7 event or higher if this type is asked for
    except ImportError:
      # psana2
      try:
        ebeam = evt.run().Detector('ebeamh')
      except KeyError: # UED has no ebeam
        return None
  return ebeam

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
      return factor_ev_angstrom / ebeam.fEbeamPhotonEnergy
    if hasattr(ebeam, 'ebeamPhotonEnergy') and ebeam.ebeamPhotonEnergy() > 0:
      # psana
      return factor_ev_angstrom / ebeam.ebeamPhotonEnergy()

    if hasattr(ebeam, 'fEbeamL3Energy') and ebeam.fEbeamL3Energy > 0:
      # pyana
      gamma = ebeam.fEbeamL3Energy / 0.510998910
    elif hasattr(ebeam, 'ebeamL3Energy') and ebeam.ebeamL3Energy() > 0:
      # psana
      gamma = ebeam.ebeamL3Energy() / 0.510998910
    elif hasattr(ebeam, 'raw'):
      # psana2
      return factor_ev_angstrom / ebeam.raw.ebeamPhotonEnergy(evt)
    else:
      return None
    K = 3.5 + delta_k
    L = 3.0e8
    return L / (2 * gamma**2) * (1 + K**2 / 2)
  return None

def env_detz(address, env):
  """The env_detz() function returns the position of the detector with
  the given address string on the z-axis in mm.  The zero-point is as
  far away as possible from the sample, and values decrease as the
  detector is moved towards the sample.
  @param address Full data source address of the DAQ device
  @param env     Environment object
  @return        Detector z-position, in mm
  """

  if env is not None:
    detector = address_split(address, env)[0]
    if detector is None:
      return None
    elif detector == 'CxiDs1':
      pv = env.epicsStore().value('CXI:DS1:MMS:06.RBV')
      if pv is None:
        # Even though potentially unsafe, fall back on the commanded
        # value if the corresponding read-back value cannot be read.
        # According to Sébastien Boutet, this particular motor has not
        # caused any problem in the past.
        pv = env.epicsStore().value('CXI:DS1:MMS:06')
      if pv is None:
        # Try the other detector. These are sometimes inconsistent
        pv = env.epicsStore().value('CXI:DS2:MMS:06.RBV')
    elif detector == 'CxiDsd' or detector == 'CxiDs2':
      # XXX Note inconsistency in naming: Dsd vs Ds2!
      pv = env.epicsStore().value('CXI:DS2:MMS:06.RBV')
      if pv is None:
        # Try the other detector. These are sometimes inconsistent
        pv = env.epicsStore().value('CXI:DS1:MMS:06.RBV')
    elif detector == 'XppGon':
      # There is no distance recorded for the XPP's CSPAD on the robot
      # arm.  Always return zero to allow the distance to be set using
      # the offset.
      return 0
    elif detector == 'XppEndstation' or \
         detector == 'MfxEndstation':
      # There is no distance recorded for the XPP's or MFX's Rayonix
      # on the robot arm.  Always return zero to allow the distance to
      # be set using the offset.
      return 0
    else:
      return None

    if pv is None:
      return None

    if hasattr(pv, "values"):
      if len(pv.values) == 1:
        return pv.values[0]
      else:
        return None
    return pv

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
