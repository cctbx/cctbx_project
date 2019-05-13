# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#

''' Filters shots from FEE spectrometer that are two color and flags xtc
    stream events as being a two color event or not.
'''
from __future__ import division
from __future__ import print_function
from six.moves import range
from xfel.cxi.cspad_ana import cspad_tbx
from xfel.cxi.cspad_ana import skip_event_flag
import numpy as np
from libtbx import easy_pickle
from psana import *

class mod_spectrum_filter_remote_only(object):
  def __init__(self,
               address,
               peak_one_range_min,
               peak_one_range_max,
               peak_two_range_min,
               peak_two_range_max,
               peak_ratio = 0.4,
               normalized_peak_to_noise_ratio = 0.4,
               spectrometer_dark_path = None,
               metal_edge_px_position = None,
               metal_edge_eV = None,
               forbidden_range_eV=20):
    """The mod_spectrum_filter class constructor stores the parameters passed
    from the psana configuration file in instance variables.

    @param address Full data source address of the FEE device
    @param peak_one_range_min the minimum x coordinate in pixel units
     of the first peak range on the detector
    @param peak_one_range_max the maximum x coordinate in pixel units
     of the first peak range on the detector
    @param peak_two_range_min the minimum x coordinate in pixel units
     of the second peak range on the detector
    @param peak_two_range_max the maximum x coordinate in pixel units
     of the second peak range on the detector
    @param peak_ratio the ratio of the two peak heights
    @param normalized_peak_to_noise_ratio ratio of the normalized integrated
      peak to the normalized integrated regions between the peaks
    @param spectrometer_dark_path path to pickle file of FEE dark
      if None then no dark is subtracted from the spectrum
    @param metal_edge_position the position in pixels of the metal edge
      if used to absorb photons in experiment if set to None it is not used as a filtering parameter
    @param percent_eV_range range in percent of metal edge region
      used as +/- a percent of the given metal edge if set to None is is not used as a filtering parameter
    """
    #self.logger = logging.getLogger(self.__class__.__name__)
    #self.logger.setLevel(logging.INFO)
    #self.logging = logging

    self.src = Source('%s'%address)
    if spectrometer_dark_path is not None:
      self.dark = easy_pickle.load(spectrometer_dark_path)
    else:
      self.dark = None
    self.peak_one_range_min = int(peak_one_range_min)
    self.peak_one_range_max = int(peak_one_range_max)
    self.peak_two_range_min = int(peak_two_range_min)
    self.peak_two_range_max = int(peak_two_range_max)
    self.normalized_peak_to_noise_ratio = float(normalized_peak_to_noise_ratio)
    self.peak_ratio = float(peak_ratio)
    self.forbidden_range_eV = int(forbidden_range_eV)
    if metal_edge_px_position is not None:
      self.metal_edge_px_position = int(metal_edge_px_position)
    else:
      self.metal_edge_px_position = None

    if metal_edge_eV is not None:
      self.metal_edge_eV = int(metal_edge_eV)
    else:
      self.metal_edge_eV = None
    self.ntwo_color = 0
    self.nnodata = 0
    self.nshots = 0
    self.naccepted= 0

  def beginjob(self, evt, env):
    pass

  def event(self,evt,evn):
    """The event() function puts a "skip_event" object with value @c
    True into the event if the shot is to be skipped.

    @param evt Event data object, a configure object
    @param env Environment object
    """
    #import pdb; pdb.set_trace()
    if (evt.get("skip_event")):
      return
    # check if FEE data is one or two dimensional
    data = evt.get(Camera.FrameV1, self.src)
    if data is None:
      one_D = True
      data = evt.get(Bld.BldDataSpectrometerV1, self.src)
    else:
      one_D = False
    # get event timestamp
    timestamp = cspad_tbx.evt_timestamp(cspad_tbx.evt_time(evt)) # human readable format

    if data is None:
      self.nnodata +=1
      #self.logger.warning("event(): No spectrum data")
      evt.put(skip_event_flag(),"skip_event")

    if timestamp is None:
      evt.put(skip_event_flag(),"skip_event")
      #self.logger.warning("event(): No TIMESTAMP, skipping shot")

    elif data is not None:
      self.nshots +=1
      # get data as array and split into two half to find each peak
      if one_D:
        # filtering out outlier spikes in FEE data
        data = np.array(data.hproj().astype(np.float64))
        for i in range(len(data)):
          if data[i]>1000000000:
            data[i]=data[i]-(2**32)
        if self.dark is not None:
          data = data - self.dark
        spectrum = data
        spectrum1 = data[:data.shape[0]//2]
        spectrum2 = data[data.shape[0]//2:]
      else:
        data = np.array(data.data16().astype(np.int32))
        if self.dark is not None:
          data = data - self.dark
        data = np.double(data)
        data_split1 = data[:,:data.shape[1]//2]
        data_split2 = data[:,data.shape[1]//2:]
        # make a 1D trace of entire spectrum and each half to find peaks
        spectrum  = np.sum(data,0)/data.shape[0]
        spectrum1 = np.sum(data_split1,0)/data_split1.shape[0]
        spectrum2 = np.sum(data_split2,0)/data_split2.shape[0]
      if not one_D:
        # the x-coordinate of the weighted center of peak region
        weighted_peak_one_positions = []
        for i in range(self.peak_one_range_min,self.peak_one_range_max):
          weighted_peak_one_positions.append(spectrum[i]*i)
        weighted_sum_peak_one = np.sum(weighted_peak_one_positions)
        weighted_peak_one_center_position = weighted_sum_peak_one/np.sum(spectrum[self.peak_one_range_min:self.peak_one_range_max])

        weighted_peak_two_positions = []
        for i in range(self.peak_two_range_min,self.peak_two_range_max):
          weighted_peak_two_positions.append(spectrum[i]*i)
        weighted_sum_peak_two = np.sum(weighted_peak_two_positions)
        weighted_peak_two_center_position = weighted_sum_peak_two/np.sum(spectrum[self.peak_two_range_min:self.peak_two_range_max])

        # normalized integrated regions between the peaks
        #int_left_region = np.sum(spectrum[weighted_peak_one_center_position+len(spectrum[self.peak_one_range_min:self.peak_one_range_max])/2:(weighted_peak_two_center_position-len(spectrum[self.peak_two_range_min:self.peak_two_range_max])/2)])
        int_left_region = np.sum(spectrum[:weighted_peak_two_center_position/2])

        #int_left_region_norm = np.sum(spectrum[weighted_peak_one_center_position+len(spectrum[self.peak_one_range_min:self.peak_one_range_max])/2:(weighted_peak_two_center_position-len(spectrum[self.peak_two_range_min:self.peak_two_range_max])/2)])/len(spectrum[weighted_peak_one_center_position+len(spectrum[self.peak_one_range_min:self.peak_one_range_max])/2:(weighted_peak_two_center_position-len(spectrum[self.peak_two_range_min:self.peak_two_range_max])/2)])
        int_left_region_norm = np.sum(spectrum[:weighted_peak_two_center_position/2])/len(spectrum[:weighted_peak_two_center_position/2])

        int_right_region = np.sum(spectrum[self.peak_two_range_max:])

        int_right_region_norm = np.sum(spectrum[self.peak_two_range_max:])/len(spectrum[self.peak_two_range_max:])

        # normalized integrated peaks
        int_peak_one = np.sum(spectrum[(weighted_peak_one_center_position-len(spectrum[self.peak_one_range_min:self.peak_one_range_max])/2):(weighted_peak_one_center_position+len(spectrum[self.peak_one_range_min:self.peak_one_range_max])/2)])

        int_peak_one_norm = np.sum(spectrum[(weighted_peak_one_center_position-len(spectrum[self.peak_one_range_min:self.peak_one_range_max])/2):(weighted_peak_one_center_position+len(spectrum[self.peak_one_range_min:self.peak_one_range_max])/2)])/len(spectrum[(weighted_peak_one_center_position-len(spectrum[self.peak_one_range_min:self.peak_one_range_max])/2):(weighted_peak_one_center_position+len(spectrum[self.peak_one_range_min:self.peak_one_range_max])/2)])

        int_peak_two = np.sum(spectrum[(weighted_peak_two_center_position-len(spectrum[self.peak_two_range_min:self.peak_two_range_max])/2):(weighted_peak_two_center_position+len(spectrum[self.peak_two_range_min:self.peak_two_range_max])/2)])

        int_peak_two_norm = np.sum(spectrum[(weighted_peak_two_center_position-len(spectrum[self.peak_two_range_min:self.peak_two_range_max])/2):(weighted_peak_two_center_position+len(spectrum[self.peak_two_range_min:self.peak_two_range_max])/2)])/len(spectrum[(weighted_peak_two_center_position-len(spectrum[self.peak_two_range_min:self.peak_two_range_max])/2):(weighted_peak_two_center_position+len(spectrum[self.peak_two_range_min:self.peak_two_range_max])/2)])

      if not one_D:
        if int_peak_one_norm/int_peak_two_norm > self.peak_ratio:
          print("event(): inflection peak too high")
          evt.put(skip_event_flag(), "skip_event")
          return
        if int_left_region_norm > self.normalized_peak_to_noise_ratio*int_peak_two_norm:
          print("event(): noisy left of low energy peak")
          evt.put(skip_event_flag(), "skip_event")
          return
        if int_right_region_norm > self.normalized_peak_to_noise_ratio*int_peak_two_norm:
          print("event(): noisy right of high energy peak")
          evt.put(skip_event_flag(), "skip_event")
          return
      #self.logger.info("TIMESTAMP %s accepted" %timestamp)
      self.naccepted += 1
      self.ntwo_color += 1
      print("%d Remote shot"  %self.ntwo_color)
      print("%s Remote timestamp" %timestamp)
  def endjob(self, obj1, obj2=None):
    """
    @param evt Event object (psana only)
    @param env Environment object
    """

    if obj2 is None:
      env = obj1
    else:
      evt = obj1
      env = obj2

    #self.logger.info("Saw %d shots, two_color %d, nodata %d " % (self.nshots, self.naccepted, self.nnodata))

  #def __del__(self):
    #logging.shutdown()
