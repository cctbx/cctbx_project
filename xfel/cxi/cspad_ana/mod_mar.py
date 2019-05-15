# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# $Id$

from __future__ import absolute_import, division, print_function

import logging, os

from xfel.cxi.cspad_ana import cspad_tbx
from xfel.cxi.cspad_ana import skip_event_flag
from scitbx.array_family import flex
from six.moves import range

class mod_mar(object):
  def __init__(self, address, directory, beam_x = None, beam_y = None, template = None):
    """
    @param address     Address string XXX Que?!
    @param directory   Directory portion of the MAR image paths
    @param beam_x      Beam center x in pixels. Uses dxtbx beam center if not specified.
    @param beam_y      Beam center y in pixels. Uses dxtbx beam center if not specified.
    @param template    Simple template for filtering files processed. If @param template
                       is in the file name, the file is accepted
    """
    self._logger = logging.getLogger(self.__class__.__name__)
    self._logger.setLevel(logging.INFO)

    # This is needed as a key for the image data.
    self._address     = cspad_tbx.getOptString(address)
    self._directory   = cspad_tbx.getOptString(directory)

    # Save the beam center and image template
    self._beam_x      = cspad_tbx.getOptInteger(beam_x)
    self._beam_y      = cspad_tbx.getOptInteger(beam_y)
    self._template    = cspad_tbx.getOptString(template)

    # At the moment, we cannot use multiprocessing for this anyway.
    # This has to be an absolute path to the image.
    self._path = None
    self._mccd_name = None


  def begincalibcycle(self, evt, env):
    # This is all very specific to the XPP experiments, so don't worry
    # so much about hardcoded stuff.

    from os.path import join
    from pypdsdata import xtc

    ctrl_config = env.getConfig(xtc.TypeId.Type.Id_ControlConfig)
    if ctrl_config is None:
      return

    # This is for the L650 experiment and beyond
    path_d = {
      'filename':"",
      'directory_part1':"",
      'directory_part2':""
    }

    for i in range(ctrl_config.npvLabels()): # interface change for psana: used to add 1 to the return value of npvLabels()
      try:
        pv = ctrl_config.pvLabel(i)
      except AttributeError:
        # interface change for psana
        pv = ctrl_config.pvLabels()[i]

      if pv.name() in path_d:
        path_d[pv.name()] = pv.value()

    if path_d['filename'] != "":
      # directory_part1 contains two parts, including where it came from and a user specific part.
      # Example: /data/blctl/G78/test/A1. We need only 'test/A1' from that string, but 'G78' varies.
      path_d['directory_part1'] = os.sep.join(path_d['directory_part1'].split(os.sep)[4:])
      self._path = join(self._directory,
                        path_d['directory_part1'],
                        path_d['directory_part2'],
                        path_d['filename'] + '.mccd')
      self._mccd_name  = path_d['filename']
      return

    # The value of this control must be of integer type.  If that is
    # not the case, the error will be caught during access.  This is
    # for the L748 experiment.
    for i in range(ctrl_config.npvControls() + 1):
      if not hasattr(ctrl_config, 'pvControl'): continue
      pv = ctrl_config.pvControl(i)
      if pv is None or pv.name() != 'marccd_filenumber':
        continue
      self._path = join(
        self._directory, 'xpp74813_image_%06d.mccd' % int(round(pv.value())))
      return


  def beginjob(self, evt, env):
    self._fmt = None


  def event(self, evt, env):
    from dxtbx.format.Registry import Registry
    from os.path import exists
    from time import sleep

    # Nop if there is no image.  For experiments configured to have
    # exactly one event per calibration cycle, this should never
    # happen.
    if self._path is None:
      evt.put(skip_event_flag(), "skip_event")
      return

    # Skip this event if the template isn't in the path
    if self._template is not None and not True in [t in self._path for t in self._template.split(',')]:
      evt.put(skip_event_flag(), "skip_event")
      return

    if "phi" in self._path:
      evt.put(skip_event_flag(), "skip_event")
      return

    # Wait for the image to appear in the file system, probing for it
    # at exponentially increasing delays.
    t = 1
    t_tot = 0

    if not exists(self._path):
      self._logger.info("Waiting for path %s"%self._path)

    while not exists(self._path):
      if t_tot > 1:
        self._logger.info("Timeout waiting for path %s"%self._path)
        evt.put(skip_event_flag(), "skip_event")
        self._logger.info("Image not found:  %s"%self._path)
        return
      sleep(t)
      t_tot += t
      t *= 2

    # Find a matching Format object and instantiate it using the
    # given path.  If the Format object does not understand the image,
    # try determining a new format.  XXX Emits "Couldn't create a
    # detector model for this image".
    if self._fmt is None:
      self._fmt = Registry.find(self._path)
      if self._fmt is None:
        evt.put(skip_event_flag(), "skip_event")
        return

    img = self._fmt(self._path)
    if img is None:
      self._fmt = Registry.find(self._path)
      if self._fmt is None:
        evt.put(skip_event_flag(), "skip_event")
        return
      img = self._fmt(self._path)
      if img is None:
        evt.put(skip_event_flag(), "skip_event")
        return

    self._logger.info(
      "Reading %s using %s" % (self._path, self._fmt.__name__))

    # Get the raw image data and convert to double precision floating
    # point array.  XXX Why will img.get_raw_data() not work, like it
    # does in print_header?
    db = img.get_detectorbase()
    db.readHeader()
    db.read()
    data = db.get_raw_data().as_double()

    # Get the pixel size and store it for common_mode.py
    detector = img.get_detector()[0]
    ps = detector.get_pixel_size()
    assert ps[0] == ps[1]
    pixel_size = ps[0]
    evt.put(ps[0],"marccd_pixel_size")
    evt.put(detector.get_trusted_range()[1],"marccd_saturated_value")
    evt.put(detector.get_distance(),"marccd_distance")

    # If the beam center isn't provided in the config file, get it from the
    # image.  It will probably be wrong.
    if self._beam_x is None or self._beam_y is None:
      self._beam_x, self._beam_y = detector.get_beam_centre_px(img.get_beam().get_s0())
      self._beam_x = int(round(self._beam_x))
      self._beam_y = int(round(self._beam_y))

    # Crop the data so that the beam center is in the center of the image
    maxy, maxx = data.focus()

    minsize = min([self._beam_x,self._beam_y,maxx-self._beam_x,maxy-self._beam_y])

    data = data[self._beam_y-minsize:self._beam_y+minsize,self._beam_x-minsize:self._beam_x+minsize]
    evt.put((minsize,minsize),"marccd_beam_center")

    evt.put(flex.int([0,0,minsize*2,minsize*2]),"marccd_active_areas")

    # Store the image in the event.
    evt.put(data, self._address)
    # Store the .mmcd file name in the event
    evt.put(self._mccd_name, "mccd_name")

    evt_time = cspad_tbx.evt_time(evt) # tuple of seconds, milliseconds
    timestamp = cspad_tbx.evt_timestamp(evt_time) # human readable format
    self._logger.info("converted  %s to pickle with timestamp %s" %(self._path, timestamp))

    # This should not be necessary as the machine is configured with
    # one event per calibration cycle.
    self._path = None

  #signature for pyana:
  #def endjob(self, env):

  #signature for psana:
  #def endjob(self, evt, env):

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

    pass
