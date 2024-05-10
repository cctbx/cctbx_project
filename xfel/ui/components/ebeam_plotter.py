from __future__ import division

import dxtbx
from xfel.cxi.cspad_ana import cspad_tbx
from matplotlib import pyplot as plt
import numpy as np
from simtbx.nanoBragg.utils import ENERGY_CONV

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
      fee_coms_eV.append(feV:=ENERGY_CONV/fwav)
      ebeams_wav.append(ewav)
      ebeams_eV.append(eeV:=ENERGY_CONV/ewav)
      print(f'{i}: {int(feV)} eV FEE / {int(eeV)} eV Ebeam')

    if plot:
      ax.hist(ebeams_eV, alpha=0.5, bins=40, label=f'run {runs[i]} ebeams ({int(eeV)} eV)')
      ax.hist(fee_coms_eV, alpha=0.5, bins=40, label=f'run {runs[i]} FEE COMs ({int(feV)} eV)')

    diffs_eV = np.array(fee_coms_eV) - np.array(ebeams_eV)
    ebeam_eV_offsets.append(sum(diffs_eV)/len(diffs_eV))
    diffs_wav = np.array(fee_coms_wav) - np.array(ebeams_wav)
    ebeam_wav_offsets.append(sum(diffs_wav)/len(diffs_wav))

  if plot:
    ax.legend()
    ax.set_xlabel('Energy (eV)')
    ax.set_ylabel('Counts')
    ax.set_title('Ebeam vs Calibrated FEE')
    if not use_figure:
      plt.show()

  return (sum(ebeam_eV_offsets)/len(ebeam_eV_offsets),
          sum(ebeam_wav_offsets)/len(ebeam_wav_offsets))

