from __future__ import absolute_import, division, print_function

from scitbx.matrix import col, sqr
from serialtbx.mono_simulation import max_like
import numpy as np
from scipy import constants
ENERGY_CONV = 1e10*constants.c*constants.h / constants.electron_volt


def extract_mosaic_parameters_using_lambda_spread(experiment, refls, verbose=True, as_eV=True,
                                                  Deff_init=3000, eta_init=0.01):
  """
  Estimates a wavelength spread from a set of indexed Bragg spots
  Recomputes delta_psi, Deffective (mosaic domain size), and eta (mosaic spread)
  using the estimated wavelength spread following Acta Cryst. (2014). D70, 3299–3309

  :param experiment: dxtbx.model.Experiment object
  :param refls: dials.array_family.flex reflection_table object
  :param verbose: bool, verbosity flag
  :param as_eV: bool, if True, retrun all_spot_wavelen in electron volts, else leave as Angstroms
  :param Deff_init: float, initial guess value for Deffective (mosaic domain size in Angstroms)
  :param eta_init: float, initial guess value for eta (mosaic spread in degrees)
  :return: tuple all_spot_wavelen, delpsi_new, Deff_opt, eta_opt
      all_spot_wavelen: list of floats, wavelengths estimated from spot predictor
      delpsi_new: list of floats, new value for delta psi, computed using per-spot wavelength estimates
      Deff_opt: float, mosaic domain size (in Angstroms), optimized using protocol outlined in above reference
      eta_opt: float, mosaic spread (in degrees), optimized along with Deff_opt
  """
  crystal = experiment.crystal
  beam = experiment.beam
  det = experiment.detector

  if crystal is None:
    print("Crystal required in expeirment")
    return
  if beam is None:
    print("beam required in expeirment")
    return
  if det is None:
    print("detector required in expeirment")
    return

  Amat = crystal.get_A()
  wave = beam.get_wavelength()
  s0_nominal = col(beam.get_s0())

  required_cols = 'panel', 'xyzobs.mm.value', 'xyzcal.mm', 'rlp', 'miller_index'
  nmisses = 0
  for column in required_cols:
    if column not in refls:
      print("Column %s missing from refl table" % column)
      nmisses += 1
  if nmisses > 0:
    print("Please format refl table to include above missing columns. Exiting.")
    return

  all_spot_wavelen = []
  all_new_delpsi = []
  d_i = []  # resolution
  twotheta = []

  try:
    delpsi_table_vals = refls["delpsical.rad"]
  except KeyError:
    delpsi_table_vals = [None] * len(refls)

  for i_ref in range(len(refls)):
    delpsi_table = delpsi_table_vals[i_ref]

    pid = refls[i_ref]['panel']
    panel = det[pid]

    # is this necessary ?
    #panel.set_px_mm_strategy(SimplePxMmStrategy())
    #panel.set_mu(0)
    #panel.set_thickness(0)

    # Extract an estimated wavelength for this spot
    res = 1/np.linalg.norm(refls[i_ref]['rlp'])
    d_i.append(res)
    cent = panel.get_beam_centre_lab(beam.get_unit_s0())
    x,y,_ = refls[i_ref]['xyzobs.mm.value']
    xcal, ycal, _ = refls[i_ref]['xyzcal.mm']
    xyz = panel.get_lab_coord((x,y))
    rad_dist = np.sqrt(np.sum((np.array(xyz) - np.array(cent))**2))
    det_dist = panel.get_distance()
    theta = 0.5*np.arctan( rad_dist / det_dist)
    twotheta.append( 2*theta)
    spot_wavelen = 2*res*np.sin(theta)  # extracted per-spot wavelength
    all_spot_wavelen.append(spot_wavelen)

    # get the rlp
    hkl = col(refls[i_ref]['miller_index'])
    qhkl = sqr(Amat)*hkl
    qnorm2 = qhkl.length()**2
    q0 = qhkl.normalize()

    # get the s1 vector with nominal shot wavelength
    xyzcal = col(panel.get_lab_coord((xcal,ycal)))
    s1_nominal = xyzcal.normalize() / wave

    # r-vector (Q) with shot wavelength
    qcal = s1_nominal - s0_nominal

    # s0 vector with new wavelength
    s0_recalc = s0_nominal.normalize() / spot_wavelen

    # this is just a shortcut to repredict the spot position with new wavelength, and then compute the new delta psi
    e1 = (q0.cross(s0_recalc)).normalize()
    c0 = (s0_recalc.cross(e1)).normalize()
    q1 = q0.cross(e1)
    a = qnorm2*spot_wavelen / 2.
    b = np.sqrt(qnorm2 - a*a)
    qcal_recalc = -a*(s0_nominal.normalize()) + b*c0

    delpsi_nominal = -np.arctan2(qcal.dot(q1), qcal.dot(q0))
    delpsi_new = -np.arctan2(qcal_recalc.dot(q1), qcal_recalc.dot(q0))
    all_new_delpsi.append( delpsi_new)
    if verbose:
      print("\n<><><><><><><><><><><><><><><><><><><>")
      print("Delta psi spot %d/%d" % (i_ref+1, len(refls)))
      print("<><><><><><><><><><><><><><><><><><><>")
      print("nominal wavelen: %.6f" % wave)
      print("per-spot wavelen: %.6f" % spot_wavelen)
      if delpsi_table is not None:
        print("Value in table: %.6f " % delpsi_table)
      print("recompute table value (sanity test): %.6f" % delpsi_nominal)
      print("compute updated value with per-spot wavelen: %.6f" % delpsi_new)

  delpsi_vals = flex.double(all_new_delpsi)
  out = max_like.minimizer(flex.double(d_i), delpsi_vals, eta_init*np.pi/180., Deff_init)
  Deff_opt = 2/out.x[0]
  eta_opt = out.x[1] * 180 / np.pi

  all_spot_wavelen = flex.double(all_spot_wavelen)
  all_spot_energies = flex.double([ENERGY_CONV/wave for wave in all_spot_wavelen])

  mean_en = np.mean(all_spot_energies)
  sig_en = np.std(all_spot_energies)

  if verbose:
    print("\n<><><><>")
    print("Results:")
    print("<><><><>")
    print("Deff with per-spot wavelen: %.6f Angstrom." % Deff_opt)
    print("eta with per-spot wavelen: %.6f deg." % eta_opt)
    print("Per-spot wavelengths: %.4f +- %.4f eV. \n" % (mean_en, sig_en))

  poly_data = all_spot_wavelen
  if as_eV:
    poly_data = all_spot_energies

  return poly_data, delpsi_vals,  Deff_opt, eta_opt


if __name__ == "__main__":
  import sys
  from dxtbx.model.experiment_list import ExperimentListFactory
  from dials.array_family import flex
  exp_fname = sys.argv[1]
  refl_fname = sys.argv[2]
  exper = ExperimentListFactory.from_json_file(exp_fname, check_format=False)[0]
  refl = flex.reflection_table.from_file(refl_fname)
  _ = extract_mosaic_parameters_using_lambda_spread(exper, refl)
