from __future__ import division

import dxtbx
from xfel.cxi.cspad_ana import cspad_tbx
from matplotlib import pyplot as plt
import numpy as np
from cctbx import factor_kev_angstrom

def gather_ebeam_fee_data(locfiles, runs=None, max_events=None, verbose=True):
  """Read every event of each loc file (up to max_events) that has both a FEE
  spectrum and an ebeam record. Returns one dict per loc file with the run
  label and the per-event energies and wavelengths from both sources; the
  lists are empty for a run with no usable events. This is the slow part of
  the ebeam calibration and does not touch any figure, so it can run on a
  worker thread."""
  if runs is None:
    runs = list(range(len(locfiles)))
  rundata = []
  for i, locfile in enumerate(locfiles):
    data = {'run': runs[i], 'ebeams_eV': [], 'ebeams_wav': [], 'fee_coms_eV': [], 'fee_coms_wav': []}
    rundata.append(data)
    if locfile is None:
      continue

    img = dxtbx.load(locfile)
    n_img = img.get_num_images()
    if max_events is not None:
      n_img = min(n_img, max_events)

    for j in range(n_img):
      if not img.get_spectrum(j):
        continue # no FEE
      ewav = cspad_tbx.evt_wavelength(img._get_event(j), run=img.get_run_from_index(j))
      if not ewav:
        continue # no ebeam
      fwav = img.get_beam(j).get_wavelength()
      feV = factor_kev_angstrom*1000./fwav
      eeV = factor_kev_angstrom*1000./ewav
      data['fee_coms_wav'].append(fwav)
      data['fee_coms_eV'].append(feV)
      data['ebeams_wav'].append(ewav)
      data['ebeams_eV'].append(eeV)
      if verbose:
        print(f'{i}: {int(feV)} eV FEE / {int(eeV)} eV Ebeam')

    if len(data['ebeams_eV']) == 0:
      print(f"No events found with both FEE and eBeam for run {runs[i]}")
  return rundata

def ebeam_fee_offsets(rundata):
  """Average over runs of the per-run median FEE minus ebeam difference, in eV
  and in Angstroms. Returns (None, None) when no run had usable events."""
  ebeam_eV_offsets = []
  ebeam_wav_offsets = []
  for data in rundata:
    if len(data['ebeams_eV']) == 0:
      continue
    diffs_eV = np.array(data['fee_coms_eV']) - np.array(data['ebeams_eV'])
    ebeam_eV_offsets.append(np.median(diffs_eV))
    diffs_wav = np.array(data['fee_coms_wav']) - np.array(data['ebeams_wav'])
    ebeam_wav_offsets.append(np.median(diffs_wav))
  if len(ebeam_eV_offsets) == 0:
    return None, None
  return (sum(ebeam_eV_offsets)/len(ebeam_eV_offsets),
          sum(ebeam_wav_offsets)/len(ebeam_wav_offsets))

def plot_ebeam_fee_comparison(rundata, fig):
  """Histogram the ebeam and calibrated FEE energies of every run in rundata
  on the given figure."""
  ax = fig.subplots()
  for data in rundata:
    if len(data['ebeams_eV']) == 0:
      ax.plot([], [], label=f"No data for run {data['run']}")
      continue
    median_fee_eV = np.median(data['fee_coms_eV'])
    median_ebeam_eV = np.median(data['ebeams_eV'])
    ax.hist(data['ebeams_eV'], alpha=0.5, bins=40, label=f"run {data['run']} ebeams ({round(median_ebeam_eV)} eV)")
    ax.hist(data['fee_coms_eV'], alpha=0.5, bins=40, label=f"run {data['run']} FEE COMs ({round(median_fee_eV)} eV)")
  ax.legend()
  ax.set_xlabel('Energy (eV)')
  ax.set_ylabel('Counts')
  ax.set_title('Ebeam vs Calibrated FEE')

def compare_ebeams_with_fees(locfiles, runs=None, plot=True, use_figure=None, max_events=None):
  rundata = gather_ebeam_fee_data(locfiles, runs=runs, max_events=max_events)
  if plot:
    fig = use_figure or plt.figure()
    plot_ebeam_fee_comparison(rundata, fig)
    if not use_figure:
      plt.show()
  return ebeam_fee_offsets(rundata)
