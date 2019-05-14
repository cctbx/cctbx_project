# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#

''' Filters shots from FEE spectrometer that are two color and flags xtc
    stream events as being a two color event or not.
'''
from __future__ import absolute_import, division, print_function
from six.moves import range
from xfel.cxi.cspad_ana import cspad_tbx
from xfel.cxi.cspad_ana import skip_event_flag
import numpy as np
from libtbx import easy_pickle
from psana import *

class mod_spectrum_filter(object):
  def __init__(self,address,
               peak_one_position_min,
               peak_one_position_max,
               peak_two_position_min,
               peak_two_position_max,
               peak_one_width,
               peak_two_width,
               peak_ratio = 0.4,
               normalized_peak_to_noise_ratio = 0.4,
               spectrometer_dark_path = None,
               iron_edge_position = None):
    """The mod_spectrum_filter class constructor stores the parameters passed
    from the psana configuration file in instance variables.

    @param address Full data source address of the FEE device
    @param peak_one_position_min the minimum x coordinate in pixel units
     of the first peak position on the detector
    @param peak_one_position_max the maximum x coordinate in pixel units
     of the first peak position on the detector
    @param peak_two_position_min the minimum x coordinate in pixel units
     of the second peak position on the detector
    @param peak_two_position_max the maximum x coordinate in pixel units
     of the second peak position on the detector
    @param peak_one_width the width in pixels of the first peak
    @param peak_two_width the width in pixels of the second peak
    @param peak_ratio the ratio of the two peak heights
    @param normalized_peak_to_noise_ratio ratio of the normalized integrated
      peak to the normalized integrated regions between the peaks
    @param spectrometer_dark_path path to pickle file of FEE dark
      if None then no dark is subtracted from the spectrum
    @param iron_edge_position the position in pixels of the iron edge if absorbing
      iron foil is used in the experiment if None this is not used as a filtering parameter
    """
    #self.logger = logging.getLogger(self.__class__.__name__)
    #self.logger.setLevel(logging.INFO)
    #self.logging = logging

    self.src = Source('%s'%address)
    if spectrometer_dark_path is not None:
      self.dark = easy_pickle.load(spectrometer_dark_path)
    else:
      self.dark = None
    self.peak_one_position_min = int(peak_one_position_min)
    self.peak_one_position_max = int(peak_one_position_max)
    self.peak_two_position_min = int(peak_two_position_min)
    self.peak_two_position_max = int(peak_two_position_max)
    self.peak_one_width = int(peak_one_width)
    self.peak_two_width = int(peak_two_width)
    self.normalized_peak_to_noise_ratio = float(normalized_peak_to_noise_ratio)
    self.peak_ratio = float(peak_ratio)
    if iron_edge_position is not None:
      self.iron_edge_position = int(iron_edge_position)
    else:
      self.iron_edge_position = None

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
      self.logger.warning("event(): No TIMESTAMP, skipping shot")

    elif data is not None:
      self.nshots +=1
      # get data as array and split into two half to find each peak
      if one_D:
        data = np.array(data.hproj().astype(np.float64))
        if 'dark' in locals():
          data = data - self.dark
        spectrum = data
        spectrum1 = data[:data.shape[0]//2]
        spectrum2 = data[data.shape[0]//2:]
      else:
        data = np.array(data.data16().astype(np.float64))
        if 'dark' in locals():
          data = data - self.dark
        data_split1 = data[:,:data.shape[1]//2]
        data_split2 = data[:,data.shape[1]//2:]
        # make a 1D trace of entire spectrum and each half to find peaks
        spectrum  = np.sum(data,0)/data.shape[0]
        spectrum1 = np.sum(data_split1,0)/data_split1.shape[0]
        spectrum2 = np.sum(data_split2,0)/data_split2.shape[0]

      peak_one = np.max(spectrum1)
      peak_two = np.max(spectrum2)
      peak_one_position = np.argmax(spectrum1)
      peak_two_position = np.argmax(spectrum2) + len(spectrum2)
    # define the limits of the regions between the two peaks
      peak_one_lower_limit = self.peak_one_position_min - self.peak_one_width
      peak_one_upper_limit = self.peak_one_position_max + self.peak_one_width
      peak_two_lower_limit = self.peak_two_position_min - self.peak_two_width
      peak_two_upper_limit = self.peak_two_position_max + self.peak_two_width

      # the x-coordinate of the weighted center of peak region
      weighted_peak_one_positions = []
      for i in range(peak_one_lower_limit,peak_one_upper_limit):
        weighted_peak_one_positions.append(spectrum[i]*i)

      weighted_sum_peak_one = sum(weighted_peak_one_positions)
      weighted_peak_one_center_position = weighted_sum_peak_one//sum(spectrum[peak_one_lower_limit:peak_one_upper_limit])

      weighted_peak_two_positions = []
      for i in range(peak_two_lower_limit,peak_two_upper_limit):
        weighted_peak_two_positions.append(spectrum[i]*i)

      weighted_sum_peak_two = sum(weighted_peak_two_positions)
      weighted_peak_two_center_position = weighted_sum_peak_two//sum(spectrum[peak_two_lower_limit:peak_two_upper_limit])
    # normalized integrated peaks
      int_peak_one = np.sum(spectrum[peak_one_lower_limit:self.peak_one_position_max])/len(spectrum[peak_one_lower_limit:self.peak_one_position_max])
      int_peak_two = np.sum(spectrum[peak_two_lower_limit:self.peak_two_position_max])/len(spectrum[peak_two_lower_limit:self.peak_two_position_max])
    # normalized integrated regions between the peaks
      int_left_region = np.sum(spectrum[0:peak_one_lower_limit])/len(spectrum[0:peak_one_lower_limit])
      int_middle_region = np.sum(spectrum[peak_one_upper_limit:peak_two_lower_limit])/len(spectrum[peak_one_upper_limit:peak_two_lower_limit])
      int_right_region = np.sum(spectrum[peak_two_upper_limit:])/len(spectrum[:peak_two_upper_limit])
      # now to do the filtering
      if peak_one/peak_two < self.peak_ratio or peak_one/peak_two > 1/self.peak_ratio:
        print("event(): too low")
        evt.put(skip_event_flag(), "skip_event")
        return
      if weighted_peak_two_center_position < self.peak_two_position_min or weighted_peak_two_center_position > self.peak_two_position_max:
        print("event(): out of range high energy peak")
        evt.put(skip_event_flag(), "skip_event")
        return
      if weighted_peak_one_center_position < self.peak_one_position_min or weighted_peak_one_center_position > self.peak_one_position_max:
        print("event(): out of range low energy peak")
        evt.put(skip_event_flag(), "skip_event")
        return
      if not one_D and (int_left_region/int_peak_one > self.normalized_peak_to_noise_ratio):
        print("event(): noisy left of low energy peak")
        evt.put(skip_event_flag(), "skip_event")
        return
      if not one_D and (int_middle_region/int_peak_one > self.normalized_peak_to_noise_ratio):
        print("event(): noisy middle")
        evt.put(skip_event_flag(), "skip_event")
        return
      if not one_D and (int_middle_region/int_peak_two > self.normalized_peak_to_noise_ratio):
        print("event(): noisy middle")
        evt.put(skip_event_flag(), "skip_event")
        return
      if not one_D and (int_right_region/int_peak_two > self.normalized_peak_to_noise_ratio):
        print("event(): noisy right of high energy peak")
        evt.put(skip_event_flag(), "skip_event")
        return
      #iron edge at 738 pixels on FFE detetor
      if one_D and (spectrum[self.iron_edge_position]>=spectrum[weighted_peak_one_center_position]):
        print("event(): peak at iron edge")
        evt.put(skip_event_flag(), "skip_event")
        return

      if one_D and (spectrum[self.iron_edge_position]>=spectrum[weighted_peak_two_center_position]):
        print("event(): peak at iron edge")
        evt.put(skip_event_flag(), "skip_event")
        return

      #self.logger.info("TIMESTAMP %s accepted" %timestamp)
      self.naccepted += 1
      self.ntwo_color += 1
      print("%d Two Color shots"  %self.ntwo_color)

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
