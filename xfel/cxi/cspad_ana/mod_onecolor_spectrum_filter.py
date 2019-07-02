# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#

''' Filters shots from FEE spectrometer '''
from __future__ import absolute_import, division, print_function
from libtbx.phil import parse
from xfel.cxi.cspad_ana import cspad_tbx
from xfel.cxi.cspad_ana import skip_event_flag
from xfel.cxi.spectra_filter import spectra_filter, phil_scope

class mod_onecolor_spectrum_filter(object):
  def __init__(self,
               phil_file = None,
               selected_filter = None):
    """The mod_onecolor_spectrum_filter class constructor stores the parameters passed
    from the psana configuration file in instance variables.

    @param phil_file Filter parameters. See cxi/spectra_filter.py
    @param selected_filter Which of the filters in phil_file to apply
    """
    self.params = phil_scope.fetch(parse(file_name = phil_file)).extract()
    self.filter = spectra_filter(self.params)
    self.selected_filter = selected_filter
    self.n_accepted = 0
    self.n_total = 0

  def beginjob(self, evt, env):
    pass

  def event(self,evt,evn):
    """The event() function puts a "skip_event" object with value @c
    True into the event if the shot is to be skipped.

    @param evt Event data object, a configure object
    @param env Environment object
    """
    if (evt.get("skip_event")):
      return
    self.n_total += 1
    ts = cspad_tbx.evt_timestamp(cspad_tbx.evt_time(evt))

    accept = self.filter.filter_event(evt, self.selected_filter)[0]
    if not accept:
      print("Skipping event", ts, ": didn't pass filter", self.selected_filter)
      evt.put(skip_event_flag(), "skip_event")
      return

    print("Accepting event", ts, ": passed filter", self.selected_filter)
    self.n_accepted += 1

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

    print("Accepted %d of %d shots"%(self.n_accepted, self.n_total))
