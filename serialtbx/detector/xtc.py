from __future__ import division

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
