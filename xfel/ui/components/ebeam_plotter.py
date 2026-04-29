from __future__ import division

import dxtbx
from xfel.cxi.cspad_ana import cspad_tbx
from matplotlib import pyplot as plt
import numpy as np
from cctbx import factor_kev_angstrom

def compare_ebeams_with_fees(locfiles, runs=None, plot=True, use_figure=None, max_events=None):
  if plot:
    fig = use_figure or plt.figure()
    ax = fig.subplots()

  ebeam_eV_offsets = []
  ebeam_wav_offsets = []

  for i in range(len(locfiles)):
    if locfiles[i] is None:
      if plot:
        ax.plot([],[], label='No data for run {runs[i]}')
      continue

    img = dxtbx.load(locfiles[i])
    n_img = img.get_num_images()
    if max_events is not None:
      n_img = min(n_img, max_events)

    ebeams_eV = []
    ebeams_wav = []
    fee_coms_eV = []
    fee_coms_wav = []

    for j in range(n_img):
      if not img.get_spectrum(j):
        continue # no FEE
      ewav = cspad_tbx.evt_wavelength(img._get_event(j))
      if not ewav:
        continue # no ebeam
      fee_coms_wav.append(fwav:=img.get_beam(j).get_wavelength())
      fee_coms_eV.append(feV:=factor_kev_angstrom*1000./fwav)
      ebeams_wav.append(ewav)
      ebeams_eV.append(eeV:=factor_kev_angstrom*1000./ewav)
      print(f'{i}: {int(feV)} eV FEE / {int(eeV)} eV Ebeam')

    if len(ebeams_eV) == 0:
      print("No events found with both FEE and eBeam")
      return None, None
    fee_coms_eV = np.array(fee_coms_eV)
    ebeam_eV = np.array(ebeams_eV)
    fee_coms_wav = np.array(fee_coms_wav)
    ebeams_wav = np.array(ebeams_wav)

    diffs_eV = fee_coms_eV - ebeam_eV
    ebeam_eV_offsets.append(np.median(diffs_eV))
    diffs_wav = fee_coms_wav - ebeams_wav
    ebeam_wav_offsets.append(np.median(diffs_wav))

    median_fee_eV = np.median(fee_coms_eV)
    median_ebeam_eV = np.median(ebeam_eV)

    if plot:
      ax.hist(ebeams_eV, alpha=0.5, bins=40, label=f'run {runs[i]} ebeams ({round(median_ebeam_eV)} eV)')
      ax.hist(fee_coms_eV, alpha=0.5, bins=40, label=f'run {runs[i]} FEE COMs ({round(median_fee_eV)} eV)')

  if plot:
    ax.legend()
    ax.set_xlabel('Energy (eV)')
    ax.set_ylabel('Counts')
    ax.set_title('Ebeam vs Calibrated FEE')
    if not use_figure:
      plt.show()

  return (sum(ebeam_eV_offsets)/len(ebeam_eV_offsets),
          sum(ebeam_wav_offsets)/len(ebeam_wav_offsets))

