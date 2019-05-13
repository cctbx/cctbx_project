from __future__ import division
from __future__ import print_function
import psana
import numpy as np
from libtbx import easy_pickle
from xfel.cxi.cspad_ana import cspad_tbx
from libtbx.phil import parse
from libtbx.utils import Sorry

phil_scope = parse("""
  spectra_filter {
    name = None
      .type = str
      .help = Name of this set of filters
    detector_address = None
      .type = str
      .help = Address for fee spectrometer, eg FeeHxSpectrometer.0:Opal1000.1
    roi = None
      .type = ints(size=4)
      .help = Spectra region of interest in pixels (xmin, xmax, ymin, ymax)
    background_roi = None
      .type = ints(size=4)
      .help = Background region of interest in pixels (xmin, xmax, ymin, ymax)
    peak_range = None
      .type = ints(size=2)
      .help = Range of pixels for the peak region (xmin, xmax). ymin and ymax \
              come from the roi.
    background_path = None
      .type = str
      .help = Path to background pickle file
    filter
      .multiple = True
      .help = Set of available filters
    {
      name = None
        .type = str
        .help = Name of this filter
      flux_min = None
        .type = float
        .help = Minimum flux to accept an event, where flux is defined as \
                sum(roi pixels-background pickle pixels) - mean(background roi pixels - background pickle pixels) \
                None: do not filter with this parameter
      flux_max = None
        .type = float
        .help = Maximum flux to accept an event. None: do not filter with this parameter
      peak_ratio_min = None
        .type = float
        .help = Minimum peak ratio to accept an event, where the peak ratio is the peak_sum/flux, and peak_sum is \
                sum(peak pixels-background pickle pixels) - mean(background roi pixels - background pickle pixels) \
                None: do not filter with this parameter
      peak_ratio_max = None
        .type = float
        .help = Maximum peak ratio to accept an event. None: do not filter with this parameter
    }
  }
""")

class spectra_filter(object):
  def __init__(self, params):
    self.src = psana.Source(params.spectra_filter.detector_address)
    self.roi = params.spectra_filter.roi
    self.bg_roi = params.spectra_filter.background_roi
    if params.spectra_filter.background_path is None:
      self.dark_pickle = None
    else:
      self.dark_pickle = easy_pickle.load(params.spectra_filter.background_path)['DATA'].as_numpy_array()
    self.peak_range = params.spectra_filter.peak_range
    self.params = params

  def filter_event(self, evt, filter_name):
    filter_phils = [f for f in self.params.spectra_filter.filter if f.name == filter_name]
    if len(filter_phils) == 0:
      raise Sorry("Filter %s not found in phil file"%filter_name)
    elif len(filter_phils) > 1:
      raise Sorry("Ambiguous filter name: %s"%filter_name)

    filter_phil = filter_phils[0]
    return self._filter_event(evt, filter_phil.flux_min, filter_phil.flux_max, filter_phil.peak_ratio_min, filter_phil.peak_ratio_max)

  def _filter_event(self, evt, flux_min, flux_max, peak_ratio_min, peak_ratio_max):
    all_data = evt.get(psana.Camera.FrameV1, self.src)
    if all_data is None:
      print("No data")
      return False, None, None, None, None, None
    print(cspad_tbx.evt_timestamp(cspad_tbx.evt_time(evt)), end=' ')
    data = np.array(all_data.data16().astype(np.int32))
    if self.dark_pickle is None:
      self.dark_pickle = np.zeros(data.shape)
    if self.bg_roi is None:
      dc_offset = None
      dc_offset_mean = 0.0
    else:
      xmin, xmax, ymin, ymax = self.bg_roi
      dc_offset = data[ymin:ymax,xmin:xmax] - self.dark_pickle[ymin:ymax,xmin:xmax] #dc: direct current
      dc_offset_mean = np.mean(dc_offset)
    if self.roi is None:
      xmin = 0; ymin = 0
      ymax, xmax = data.shape
    else:
      xmin, xmax, ymin, ymax = self.roi
    data_signal = data[ymin:ymax,xmin:xmax] - self.dark_pickle[ymin:ymax,xmin:xmax]
    data = data_signal - dc_offset_mean

    # make a 1D trace
    spectrum = np.sum(data, 0)
    flux = np.sum(spectrum)
    if self.peak_range is None:
      print("Run", evt.run(), "flux:", flux)
    else:
      peak = spectrum[self.peak_range[0]-xmin:self.peak_range[1]-xmin]
      peak_max = np.sum(peak)
      print("Run", evt.run(), "peak max:", peak_max, "at", np.argmax(peak)+xmin, "px", "flux:", flux, "m/f:", peak_max/flux)

      if flux_min is not None and flux < flux_min: return False, None, None, None, None, None
      if flux_max is not None and flux >= flux_max: return False, None, None, None, None, None
      if peak_ratio_min is not None and peak_max/flux < peak_ratio_min: return False, None, None, None, None, None
      if peak_ratio_max is not None and peak_max/flux >= peak_ratio_max: return False, None, None, None, None, None

    all_data = all_data.data16().astype(np.int32)
    return True, data, spectrum, dc_offset, all_data - self.dark_pickle, all_data
