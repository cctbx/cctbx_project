# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# $Id$

# XXX Better dark generation: are the bright pixels really due to
# cosmic radiation?

# Does it even make sense to run this on more than one core?  NO IT DOES NOT!

from __future__ import absolute_import, division, print_function

import math
import numpy

from pypdsdata import xtc


from scitbx.array_family import flex

from xfel.cxi.cspad_ana import common_mode
from xfel.cxi.cspad_ana import cspad_tbx

from collections import deque
from six.moves import range


def _summed_area_table(f, m, n):
  """The _summed_area_table() function calculates the summed-area tables
  required for the denominator in the expansion of the normalised
  correlation coefficient.

  @param f Two-dimensional intensity image
  @param m Height of the template, in pixels
  @param n Width of the template, in pixels
  @return  The energy of @f as a function of the template position
  """

  from scitbx.array_family import flex

  s = flex.double(flex.grid(f.focus()[0] + 2 * m,
                            f.focus()[1] + 2 * n))
  s.matrix_paste_block_in_place(block=f, i_row=m, i_column=n)

  # Treat the first row and the first column as special cases, then
  # fill in the rest of matrix.
  for u in range(1, s.focus()[0]):
    s[u, 0] += s[u - 1, 0]
  for v in range(1, s.focus()[1]):
    s[0, v] += s[0, v - 1]
  for u in range(1, s.focus()[0]):
    for v in range(1, s.focus()[1]):
      s[u, v] += s[u - 1, v] + s[u, v - 1] - s[u - 1, v - 1]

  # Note the off-by-one discrepancy with respect to Lewis (1995), due
  # to the padding.
  e = flex.double(flex.grid(f.focus()[0] + m - 1, f.focus()[1] + n - 1))
  e[0, 0] = s[m, n]
  for u in range(1, e.focus()[0]):
    e[u, 0] = s[u + m, n] - s[u, n]
  for v in range(1, e.focus()[1]):
    e[0, v] = s[m, v + n] - s[m, v]
  e[1:e.focus()[0], 1:e.focus()[1]] = \
    s[(m + 1):(e.focus()[0] + m ), (n + 1):(e.focus()[1] + n)] - \
    s[    (1):(e.focus()[0]),      (n + 1):(e.focus()[1] + n)] - \
    s[(m + 1):(e.focus()[0] + m),      (1):(e.focus()[1])]     + \
    s[    (1):(e.focus()[0]),          (1):(e.focus()[1])]

  return e


def lewis(template, img):
  """The lewis() function computes the normalised cross-correlation
  (NCC) of an image, @p img, and a template, @p template.  Both image
  and template must be two-dimensional, real, and finite.  The image
  must be larger than the template, which in turn should contain more
  than one pixel, and the template must have positive variance.  The
  function returns the correlation coefficients in the range [-1, +1].
  See Lewis, J. P. (1995) "Fast Template Matching", Vision Interface,
  120-123.

  @note This function should be equivalent to MATLAB's normxcorr2()
        function.

  @param img      Two-dimensional intensity image
  @param template Two-dimensional intensity template
  @return         Correlation coefficients
  """

  import math
  import numpy
  from sys import float_info

  from scitbx import fftpack
  from scitbx.array_family import flex

  # Assert that image and template are two-dimensional arrays, and
  # that the template is no larger than image.  Assert that template
  # is not flat.  XXX Check for real and finite, too?
  assert len(img.focus()) == 2 \
    and  len(template.focus()) == 2
  assert img.focus()[0] >= template.focus()[0] \
    and  img.focus()[1] >= template.focus()[1]
  assert template.sample_standard_deviation() > 0

  # For conformance with MATLAB's normxcorr2() and geck 342320: for
  # numerical robustness, ensure that both image and template are
  # always non-negative.
  img_nn = img - min(0, flex.min(img))
  template_nn = template - min(0, flex.min(template))

  # Calculate the terms of the denominator of gamma.  Must guard
  # against negative variance of the image due to inaccuracies in the
  # one-pass formula.
  img_sum = _summed_area_table(
    img_nn, template_nn.focus()[0], template_nn.focus()[1])
  img_ssq = _summed_area_table(
    flex.pow2(img_nn), template_nn.focus()[0], template_nn.focus()[1])

  f_sigma = (img_ssq - img_sum * img_sum /
             (template_nn.focus()[0] * template_nn.focus()[1]))
  f_sigma.set_selected(f_sigma < 0, 0)
  f_sigma = flex.sqrt(f_sigma)
  t_sigma = (template_nn - flex.mean(template_nn)).norm()

  gamma_denominator =  f_sigma * t_sigma

  # Zero-pad the image to permit partial overlap of template and
  # image, and embed the time-reversed template in a zero-padded array
  # of the same size.  Zero-padding means the entire template is
  # always overlapping the image, and terms involving the template
  # mean cancel in the expansion of the numerator of gamma.
  #
  # Note: the NCC demands the template to be time-reversed, which can
  # be accomplished by conjugation in the frequency domain.  An
  # implementation following that approach would however require
  # special care to be taken for the first rows and columns:
  #
  #   from numpy import roll
  #   t_embed.matrix_paste_block_in_place(
  #     block=template_nn,
  #     i_row=full[0] - template_nn.focus()[0],
  #     i_column=full[1] - template_nn.focus()[1])
  #   t_embed = flex.double(roll(
  #     roll(t_embed.as_numpy_array(), 1, axis=0), 1, axis=1))
  #
  # Calculate correlation in frequency domain.  XXX Could use spatial
  # domain calculation in cases where it's faster (see MATLAB's
  # implementation).
  full = (img_nn.focus()[0] + template_nn.focus()[0] - 1,
          img_nn.focus()[1] + template_nn.focus()[1] - 1)

  f_embed = flex.double(flex.grid(full))
  f_embed.matrix_paste_block_in_place(
    block=img_nn, i_row=0, i_column=0)
  f_prime = flex.complex_double(
    reals=f_embed, imags=flex.double(flex.grid(full)))

  t_embed = flex.double(flex.grid(full))
  t_embed.matrix_paste_block_in_place(
    block=template_nn.matrix_rot90(2), i_row=0, i_column=0)
  t_prime = flex.complex_double(
    reals=t_embed, imags=flex.double(flex.grid(full)))

  fft = fftpack.complex_to_complex_2d(full)
  fft.forward(f_prime)
  fft.forward(t_prime)
  gamma_numerator = f_prime * t_prime
  fft.backward(gamma_numerator)
  gamma_numerator = flex.real(gamma_numerator) / (fft.n()[0] * fft.n()[1]) \
                    - img_sum * flex.mean(template_nn)

  # For conformance with MATLAB: set the NCC to zero in regions where
  # the image has zero variance over the extent of the template.  If,
  # due to small variances in the image or the template, a correlation
  # coefficient falls outside the range [-1, 1], set it to zero to
  # reflect the undefined 0/0 condition.
  tol = math.sqrt(math.ldexp(
      float_info.epsilon,
      math.frexp(flex.max(flex.abs(gamma_denominator)))[1] - 1))
  sel = gamma_denominator <= tol
  gamma = gamma_numerator.set_selected(sel, 0) / \
          gamma_denominator.set_selected(sel, 1)
  gamma.set_selected(flex.abs(gamma) > 1 + math.sqrt(float_info.epsilon), 0)

  return gamma


class ring_buffer(deque):
  """Better named history_buffer?
  """

  def __init__(self, maxlen=1024):
    """Must ensure to call the constructor with maxlen other than @c
    None.  The maxlen parameter requires Python 2.6 or later.
    """
    assert maxlen is not None
    super(ring_buffer, self).__init__(maxlen=maxlen)


  def frequency(self):
    """If time is in seconds, returns the frequency, in Hz, of the
    value.
    """
    dt = self[-1][0] - self[0][0]
    if dt <= 0:
      return 0
    return (len(self) - 1) / dt


  def push(self, time, value):
    if len(self) == 0 or self[-1][1] != value:
      self.append((time, value))


class mod_ledge(common_mode.common_mode_correction):

  def __init__(
      self, address, display, mat_path, table_path, template_path=None, **kwds):
    """The mod_average class constructor stores the parameters passed
    from the pyana configuration file in instance variables.  All
    parameters, except @p address are optional, and hence need not be
    defined in pyana.cfg.

    @param address Address string XXX Que?!
    """

    super(mod_ledge, self).__init__(address=address, **kwds)

    # XXX Should be forced to false if graphics unavailable (e.g. when
    # running on the cluster).
    self._display = cspad_tbx.getOptBool(display)

    # Use line buffering to allow tailing the output in real time.
    # XXX Make name configurable.
    self._mat_path = cspad_tbx.getOptString(mat_path)
    self._table_path = cspad_tbx.getOptString(table_path)

    # Ring buffers for history-keeping.
    self._history_injector_xyz = ring_buffer()
    self._history_spectrometer_xyz = ring_buffer()
    self._history_energy = ring_buffer()

    # Get the template for image registration.
    if template_path is not None:
      from libtbx import easy_pickle
      self._template = easy_pickle.load(template_path)
    else:
      self._template = None

    # Optionally, initialise the viewer with a custom callback, see
    # mod_view.  The regions of interest are communicated to the
    # viewer through a shared multiprocessing array.  XXX This needs a
    # bit more thought to come up with a sensible interface.
    if self._display:
      from .mod_view import _xray_frame_process
      from multiprocessing import Array, Process, Manager
      from rstbx.viewer import display

      manager = Manager()
      self._queue = manager.Queue()
      self._proc = Process(
        target=_xray_frame_process, args=(self._queue, True, 0))
      self._roi = Array('i', 15 * 4, lock=False) # XXX Make variable!
      display.user_callback = self._callback
      self._proc.start()


  def _callback(self, dc, wxpanel, wx):
    # Draws a rectangle with the given top left corner, and with the
    # given size.  The current pen is used for the outline and the
    # current brush for filling the shape.
    #
    # These are like in wxPython, (x, y, width, height).  Origin is at
    # top left corner, width increases to the left, height downwards.
    # Log it!

    dc.SetBrush(wx.TRANSPARENT_BRUSH)
    dc.SetPen(wx.Pen('red'))
    for i in range(0, len(self._roi), 4):
      roi = self._roi[i:i + 4]

      # XXX Figure out what's going on here--subtract one because it's
      # inclusive?!
      tl = wxpanel._img.image_coords_as_screen_coords(
        roi[0] - 0.5, roi[1] - 0.5)
      br = wxpanel._img.image_coords_as_screen_coords(
        roi[0] + roi[2], roi[1] + roi[3])
      dc.DrawRectangle(tl[0], tl[1], br[0] - tl[0], br[1] - tl[1])

      """
      tl = wxpanel._img.image_coords_as_screen_coords(
        roi[0] + roi[2] - 0.5, roi[1] - 0.5)
      br = wxpanel._img.image_coords_as_screen_coords(
        roi[0] - 0.5 + roi[2] + 10, roi[1] - 0.5 + 10)
      dc.DrawRectangle(tl[0], tl[1], br[0] - tl[0], br[1] - tl[1])

      tl = wxpanel._img.image_coords_as_screen_coords(
        roi[0] - 0.5, roi[1] - 0.5 + roi[3])
      br = wxpanel._img.image_coords_as_screen_coords(
        roi[0] - 0.5 + 10, roi[1] - 0.5 + roi[3] + 10)
      dc.DrawRectangle(tl[0], tl[1], br[0] - tl[0], br[1] - tl[1])
      """


  def _reset_counters(self):
    self._I0 = flex.double()
    self._fee_before = flex.double()
    self._fee_after = flex.double()
    self._hit = flex.bool()

    self._injector_plus_current = flex.double()
    self._injector_plus_voltage = flex.double()
    self._injector_minus_current = flex.double()
    self._injector_minus_voltage = flex.double()

    self._injector_micos_xyz = flex.vec3_double()
    self._injector_rough_xyz = flex.vec3_double()
    self._repetition_rate = flex.double()
    self._timestamp = flex.double()
    self._spectrometer_xyz = flex.vec3_double()
    self._energy = flex.double()

    self._acq_apd_integral = flex.double()
    self._acq_opto_diode_integral = flex.double()

#    self._injector_x = flex.double()
#    self._injector_y = flex.double()
#    self._injector_z = flex.double()

#    self._spectrometer_x = flex.double()
#    self._spectrometer_y = flex.double()
#    self._spectrometer_z = flex.double()


  @staticmethod
  def _filtered_stats(function, iterable):
    """The _filtered_stats() computes first- and second-order statistics
    on an array, after first applying a filter.

    @param function Function which returns @c True for those elements
                    of @p iterable which should be excluded in the
                    returned statistics, or @c None to include all
                    data
    @param iterable An iterable sequence of data items
    @return         A four-tuple of mean, standard deviation,
                    effective sample size, and number of rejected
                    samples
    """

    filtered_data = filter(function, iterable)
    if len(filtered_data) == 0:
      return (0, 0, 0, 0)

    stats = flex.mean_and_variance(flex.double(filtered_data))
    mean = stats.mean()
    if len(filtered_data) > 1:
      stddev = stats.unweighted_sample_standard_deviation()
    else:
      stddev = 0

    return (mean,
            stddev,
            len(filtered_data),
            len(iterable) - len(filtered_data))


  def beginjob(self, evt, env):
    """The beginjob() function does one-time initialisation from
    event- or environment data.  It is called at an XTC configure
    transition.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    from os import makedirs, path

    super(mod_ledge, self).beginjob(evt, env)
    self._reset_counters()
    self._nframes = 0

    # XXX Don't really want to store the stream, but how to properly
    # close it on exit?
    self._table_path = cspad_tbx.pathsubst(self._table_path, evt, env)
    if not path.isdir(path.dirname(self._table_path)):
      makedirs(path.dirname(self._table_path))

    self._stream_table = open(self._table_path, mode='wb', buffering=1)

    # Initial camera settings requested by Rolf.  The Andor's region
    # of interest is defined by width, height, orgX, and orgY,
    # independent of the binning.  The camera has model number
    # DO936N-M0W-#BN.  Note that the configuration may change on an
    # event-basis.  XXX Does the width really correspond to X, and the
    # height to Y?
    config = cspad_tbx.getConfig(self.address, env)
    if config is not None:
      self.logger.info("Initial effective frame size: %dx%d" %
                       (config.width()  // config.binX(),
                        config.height() // config.binY()))
      self.logger.info("Initial region of interest: (%d, %d) - (%d, %d)" % \
                         (config.orgX(),
                          config.orgY(),
                          config.orgX() + config.width(),
                          config.orgY() + config.height()))

      if config.highCapacity() >= 0 and config.highCapacity() <= 1:
        self.logger.info(
          "Initial output mode: %s" % \
            ("high sensitivity", "high capacity")[config.highCapacity()])
      if config.readoutSpeedIndex() >= 0 and config.readoutSpeedIndex() <= 3:
        self.logger.info(
          "Initial pixel readout rate: %s" % \
            ("5 MHz", "3 MHz", "1 MHz", "50 kHz")[config.readoutSpeedIndex()])
      if config.gainIndex() >= 0 and config.gainIndex() <= 2:
        self.logger.info(
          "Initial pre-amplifier gain: %s" % \
            ("1x", "2x", "4x")[config.gainIndex()])


  def event(self, evt, env):
    """The event() function is called for every L1Accept transition.

    For now, log error and set bogus value to allow stuff to continue
    -- must check for the bogosity later

    XXX The dead time of the detector complicates checking how often
    things are updated!  Move this to the ring buffer?

    @param evt Event data object, a configure object
    @param env Environment object
    """

    from pyana.event import Event
    from acqiris_ext import acqiris_integrate, apd_hitfind

    super(mod_ledge, self).event(evt, env)
    if evt.status() != Event.Normal:
      pass # XXX return -- Never skip because arrays will end up
           # different length, so ignore this?

    # Get the time of the event, in fractional seconds since the
    # epoch.  This is needed for all subsequent history-keeping, and
    # is hence determined first.  XXX Is history-keeping even
    # justified?
    time = cspad_tbx.evt_time(evt)
    if time is None:
      time = float('nan')
    else:
      time = time[0] + time[1] / 1e3
    self._timestamp.append(time)

    # The repetition rate is currently just used for sanity checking.
    repetition_rate = cspad_tbx.evt_repetition_rate(evt)
    if repetition_rate is None:
      repetition_rate = float('nan')
    self._repetition_rate.append(repetition_rate)

    # Get the I0.  No need to warn about it here, it will be done once
    # the image is written out.
    I0 = cspad_tbx.evt_pulse_energy(evt)
    if I0 is None:
      I0 = float('nan')
    self._I0.append(I0)

    # Get the FEE energy.  Average the two readings before and after
    # attenuation separately.  XXX What are the units?  It look like
    # it could be mJ?
    fee_before = 0.5 * sum(evt.getFeeGasDet()[0:2])
    if fee_before is None:
      fee_before = float('nan')
    self._fee_before.append(fee_before)

    fee_after = 0.5 * sum(evt.getFeeGasDet()[2:4])
    if fee_after is None:
      fee_after = float('nan')
    self._fee_after.append(fee_after)

    # XXX Just a check: this is what xtcexplorer does:
    fee_energy = evt.get(xtc.TypeId.Type.Id_FEEGasDetEnergy)
    if fee_energy is not None:
      assert evt.getFeeGasDet()[0] == fee_energy.f_11_ENRC \
        and  evt.getFeeGasDet()[1] == fee_energy.f_12_ENRC \
        and  evt.getFeeGasDet()[2] == fee_energy.f_21_ENRC \
        and  evt.getFeeGasDet()[3] == fee_energy.f_22_ENRC

    """
    # For Bill: expect 84240 data points for r0054
    #
    # grep "^BILL_POINT" | cut -d' ' -f2,3,4,5,6 > t.dat
    # gnuplot> m=0.1 ; k=-0.01e-8; f(x) = k * x + m
    # gnuplot> fit f(x) "t.dat" using ($3):($5) via k,m
    if not hasattr(self, '_gmd_seqno'):
      self._gmd_seqno = 0
    gmd = evt.get(key=xtc.TypeId.Type.Id_GMD)
    if gmd is None:
      return
    acq_apd = evt.getAcqValue('SxrEndstation-0|Acqiris-1', 0, env)
    if acq_apd is not None and acq_apd.waveform() is not None:
      w = acq_apd.waveform()
      baseline = numpy.mean(w[0:(w.shape[0] / 5)])
      peak = numpy.min(w[(w.shape[0] / 5):w.shape[0]])
      self._gmd_seqno += 1
      print "BILL_POINT %d %s %s %s %s" % (self._gmd_seqno,
                                           repr(gmd.fBgValuePerSample),
                                           repr(gmd.fCorrectedSumPerPulse),
                                           repr(gmd.fRelativeEnergyPerPulse),
                                           repr(peak - baseline))
    return
    """

    """
    # XXX Record injector motion--note that they cannot be added--see
    # Ray's email.
    injector_micos_xyz = cspad_tbx.env_pv3_get(
      env,
      ['SXR:EXP:MZM:%02d:ENCPOSITIONGET' % i for i in [1, 2, 3]])
    if injector_micos_xyz is None:
      self.logger.error("No micos injector motor positions")
      injector_micos_xyz = (float('nan'), float('nan'), float('nan'))
    self._injector_micos_xyz.append(injector_micos_xyz)

    injector_rough_xyz = cspad_tbx.env_pv3_get(
      env,
      ['SXR:EXP:MMS:%02d.RBV' % i for i in [1, 2, 3]])
    if injector_rough_xyz is None:
      self.logger.error("No rough injector motor positions")
      injector_rough_xyz = (float('nan'), float('nan'), float('nan'))
    self._injector_rough_xyz.append(injector_rough_xyz)

    # Injector power supplies XXX There is a third PSU, no?
    #
    # The -5kV supply
    # SXR:EXP:SHV:VHS6:CH0:VoltageMeasure
    # SXR:EXP:SHV:VHS6:CH0:CurrentMeasure
    #
    # The plus 5kV supply
    # SXR:EXP:SHV:VHS2:CH0:VoltageMeasure
    # SXR:EXP:SHV:VHS2:CH0:CurrentMeasure
    injector_plus_current = cspad_tbx.env_pv1_get(
      env, 'SXR:EXP:SHV:VHS6:CH0:CurrentMeasure')
    if injector_plus_current is None:
      self.logger.error("No plus-motor current")
      injector_plus_current = -1
    self._injector_plus_current.append(injector_plus_current)

    injector_plus_voltage = cspad_tbx.env_pv1_get(
      env, 'SXR:EXP:SHV:VHS6:CH0:VoltageMeasure')
    if injector_plus_voltage is None:
      self.logger.error("No plus-motor voltage")
      injector_plus_voltage = -1
    self._injector_plus_voltage.append(injector_plus_voltage)

    injector_minus_current = cspad_tbx.env_pv1_get(
      env, 'SXR:EXP:SHV:VHS2:CH0:CurrentMeasure')
    if injector_minus_current is None:
      self.logger.error("No minus-motor current")
      injector_minus_current = -1
    self._injector_minus_current.append(injector_minus_current)

    injector_minus_voltage = cspad_tbx.env_pv1_get(
      env, 'SXR:EXP:SHV:VHS2:CH0:VoltageMeasure')
    if injector_minus_voltage is None:
      self.logger.error("No minus-motor voltage")
      injector_minus_voltage = -1
    self._injector_minus_voltage.append(injector_minus_voltage)
    """

    """
    # The spectrometer motor positions are just used for sanity
    # checking.
    spectrometer_xyz = cspad_tbx.env_spectrometer_xyz_sxr(env)
    if spectrometer_xyz is None:
      self.logger.error("No spectrometer motor positions")
      spectrometer_xyz = (float('nan'), float('nan'), float('nan'))
    self._spectrometer_xyz.append(spectrometer_xyz)
    """

    # Get the pulse energy after monochromator, and fall back on the
    # pre-monochromator energy if the former is absent.  Record in
    # list for mean and stddev.  XXX Verify that the wavelength after
    # the monochromator is updated at around 1 Hz.
    #
    # For the publication an offset and scale were calibrated.
    wavelength = cspad_tbx.env_wavelength_sxr(evt, env)
    if wavelength is None:
      wavelength = cspad_tbx.evt_wavelength(evt)
    if wavelength is None:
      energy = float('nan')
    else:
      energy = 12398.4187 / wavelength
    self._energy.append(energy)
    self._history_energy.push(time, energy) # XXX Not necessary?!

    """
    # Laser shutters XXX need to sort out laser numbering XXX Laser
    # power stuff? XXX Position of polarizer/analyser
    shutters = cspad_tbx.env_laser_shutters(env)
    #print "Got shutters", shutters
    """

    # Read out the diode traces from the via the Acqiris.  XXX In any
    # case, the APD and the more sensitive Opto Diode in the monitor
    # tank (i.e. the transmission diode) should be anti-correlated, so
    # check it!  The entire trace always covers 10 us.  XXX Could this
    # be figured out from xtc.TypeId.Type.Id_AcqConfig?
    #
    # XXX This appears to be suboptimal: look at the
    # skewness-transform for the APD to sort this out.
    acq_apd = evt.getAcqValue('SxrEndstation-0|Acqiris-1', 0, env)
    acq_apd_integral = float('nan')
    if acq_apd is not None:
      waveform = acq_apd.waveform()
      if waveform is not None:
        # With a 40k-point trace, one should integrate from 18200 to
        # 18400.
        waveform = waveform.flatten()
        nmemb = len(waveform) // 200
        if nmemb > 0:
          acq_apd_integral = acqiris_integrate(
            flex.double(waveform), 91 * nmemb, 100 * nmemb, nmemb)
    self._acq_apd_integral.append(acq_apd_integral)

    if evt.expNum() == 208:
      # Opto diode address for L632.
      acq_opto_diode = evt.getAcqValue('SxrEndstation-0|Acqiris-1', 1, env)
    elif evt.expNum() == 363:
      # Opto diode address for LB68.
      acq_opto_diode = evt.getAcqValue('SxrEndstation-0|Acqiris-2', 2, env)
    acq_opto_diode_integral = float('nan')
    if acq_opto_diode is not None:
      waveform = acq_opto_diode.waveform()
      if waveform is not None:
        # With a 40k-point trace, one should integrate from 16000 to
        # 24000.  With a 20k-point trace, a suitable integration
        # region is bounded by 8000 and 12000.  There is no need for
        # thresholding, because the integral of the Opto Diode will
        # not be used for hit finding.  XXX What are the "misses" we
        # record on the Opto Diode?  XXX The direct beam is completely
        # gone after it hits the sample, because soft X-rays.
        waveform = waveform.flatten()
        nmemb = len(waveform) // 5
        if nmemb > 0:
          acq_opto_diode_integral = acqiris_integrate(
            flex.double(waveform), 2 * nmemb, 4 * nmemb, nmemb)
    self._acq_opto_diode_integral.append(acq_opto_diode_integral)

    # Sanity check: verify that the timestamps for the two Acqiris
    # traces are similar enough.
    if acq_apd is not None and acq_opto_diode is not None:
      assert \
          len(acq_apd.timestamps()) == len(acq_opto_diode.timestamps()) and \
          numpy.any(numpy.abs(acq_apd.timestamps() -
                              acq_opto_diode.timestamps())) < 1e-6

    #self.logger.info("DIODE INTEGRALS: %f %f %f" % (I0, acq_apd_integral, acq_opto_diode_integral))

    """
    import matplotlib.pyplot as plt

    hit_array_apd = apd_hitfind(
      flex.double(acq_apd.waveform()),
      len(acq_apd.waveform()) // 5)
    hit_array_opto_diode = apd_hitfind(
      flex.double(acq_opto_diode.waveform()),
      len(acq_opto_diode.waveform()) // 5)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    #ax.plot(
    #  range(len(acq_apd.timestamps())), acq_apd.waveform())
    ax.plot(
      range(len(acq_opto_diode.timestamps())), acq_opto_diode.waveform()[0, :])
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    #ax.plot(
    #  acq_apd.timestamps()[0:len(hit_array_apd)], hit_array)
    ax.plot(
      acq_opto_diode.timestamps()[0:len(hit_array_opto_diode)], hit_array)
    plt.show()
    """

    # Determine whether the beam hit the sample, and register the
    # outcome.  If not using any diodes for hit-finding, every shot is
    # assumed to be a hit.  XXX Unfortunately, this crucial piece is
    # very unreliable.  The threshold for the APD needs to be
    # verified--inspect all the histograms.  XXX hitfind_flags is
    # probable better as a module parameter.
#    hitfind_flags = 0x3
    hitfind_flags = 0
    hit = False
    if not hitfind_flags:
      hit = True
    elif hitfind_flags & 0x1 and acq_apd_integral > 0.2:
      hit = True
    self._hit.append(hit)

    # Always proceed all the way through (even if some shots have
    # invalid values of e.g. I0) because images are precious.  XXX
    # Must reset counters before returning!  XXX What about skipping
    # all of the above if display is True?
    if self.cspad_img is not None:
      self._nframes += 1

      """
      # The spectrometer should not move!
      t = (self._spectrometer_xyz -
           self._spectrometer_xyz.mean()).rms_length()
      print "Spectrometer displacement", t

      # Fine/rough motor position deviations from the mean.  See Ray's
      # email.
      t = (self._injector_micos_xyz -
           self._injector_micos_xyz.mean()).rms_length()
      print "Injector micos displacement", t

      t = (self._injector_rough_xyz -
           self._injector_rough_xyz.mean()).rms_length()
      print "Injector rough displacement", t

      # Injector motor position means and deviations
      if self._injector_plus_current.size() > 1:
        t = flex.mean_and_variance(self._injector_plus_current)
        print "Injector plus current mean %10e stddev %10e" % \
            (t.mean(), t.unweighted_sample_standard_deviation())
      if self._injector_plus_voltage.size() > 1:
        t = flex.mean_and_variance(self._injector_plus_voltage)
        print "Injector plus voltage mean %10e stddev %10e" % \
            (t.mean(), t.unweighted_sample_standard_deviation())

      if self._injector_minus_current.size() > 1:
        t = flex.mean_and_variance(self._injector_minus_current)
        print "Injector minus current mean %10e stddev %10e" % \
            (t.mean(), t.unweighted_sample_standard_deviation())
      if self._injector_minus_voltage.size() > 1:
        t = flex.mean_and_variance(self._injector_minus_voltage)
        print "Injector minus voltage mean %10e stddev %10e" % \
            (t.mean(), t.unweighted_sample_standard_deviation())

      """

      # Energy statistics are collected from all shots, regardless of
      # whether they are hits or not.  Since this statistic mentions
      # the frame number, it should be reported first.  XXX The energy
      # should have a really small standard deviation.  Check
      # self._energy.size() and self._history_energy.frequency() XXX
      # verify that it works for one data point.
      (energy_mean, energy_stddev, energy_nmemb, n) = self._filtered_stats(
        lambda x: not math.isnan(x) and x > 0, self._energy)
      if n > 0:
        self.logger.warning("%d shots have undefined energy" % n)

      (I0_mean, I0_stddev, I0_nmemb, n) = self._filtered_stats(
        lambda x: not math.isnan(x), self._I0)
      if n > 0:
        self.logger.warning("%d shots have undefined I0" % n)

      self.logger.info(
        "Frame %d: E=%.3f+/-%.3f (N=%d) I0=%.0f+/-%.0f (N=%d)" %
        (self._nframes,
         energy_mean, energy_stddev, energy_nmemb,
         I0_mean, I0_stddev, I0_nmemb))

      # Sanity check: unless changed while integrating the frame, the
      # repetition rate should have a standard deviation of zero.
      dt = self._timestamp[-1] - self._timestamp[0]
      rr_mean = rr_observed = rr_stddev = 0
      if dt > 0:
        rr_observed = (len(self._timestamp) - 1) / dt
        rr = filter(
          lambda x: not math.isnan(x) and x > 0, self._repetition_rate)
        if len(rr) > 1:
          rr_stats = flex.mean_and_variance(flex.double(rr))
          rr_mean = rr_stats.mean()
          rr_stddev = rr_stats.unweighted_sample_standard_deviation()
      self.logger.info(
        "Repetition rate: %.3f Hz (observed), %.3f+/-%.3f Hz (expected)" %
        (rr_observed, rr_mean, rr_stddev))

      # Compare observed and configured exposure time.
      config = cspad_tbx.getConfig(self.address, env)
      exposure_time = 0
      if config is not None and dt > 0 and len(self._timestamp) > 0:
        exposure_time = dt * (len(self._timestamp) + 1) / len(self._timestamp)
      self.logger.info(
        "Exposure time: %.3f s (observed), %.3f s (configured)" %
        (exposure_time, config.exposureTime()))

      # Compute the leading dead time, the time between starting the
      # readout of the previous frame and the arrival of the shot
      # immediately following it.  This is an interesting statistic,
      # no matter what.  XXX Maybe look at its distribution?
      dead_time = 0
      if rr_observed > 0 and hasattr(self, '_previous_readout_time'):
        dead_time = \
            self._timestamp[0] - self._previous_readout_time - 1 / rr_observed
        if math.isnan(dead_time):
          dead_time = 0
      self.logger.info("Dead time: %.3f s" % dead_time)
      self._previous_readout_time = self._timestamp[-1]

      assert time == self._timestamp[-1] # XXX ZAP once one run survives it!

      # Flag blank images (i.e. images that had no hits), because
      # these may interesting for background subtraction.
      hits = self._hit.count(True)
      self.logger.info("Hit rate: %d/%d (%.2f%%)" %
                       (hits, self._hit.size(), 100 * hits / self._hit.size()))
      if hits == 0:
        self.logger.info("Frame %d is blank" % self._nframes)

      # Get the normalisation factor by summing up I0 for all hits.
      # Invalid and non-positive values of I0 are treated as zeroes.
      # XXX Make this kind of summing a function of its own.
      I0 = sum(filter(lambda x: not math.isnan(x) and x > 0,
                      self._I0.select(self._hit)))
      I0_all = sum(filter(lambda x: not math.isnan(x) and x > 0,
                          self._I0))

      fee_before_all = sum(filter(lambda x: not math.isnan(x) and x > 0,
                                  self._fee_before))
      fee_after_all = sum(filter(lambda x: not math.isnan(x) and x > 0,
                                 self._fee_after))

      # Register the template to the image and locate the regions of
      # interest based on the registration parameters.  XXX Should
      # also give contrast: fit 2D-Gaussian to peak and report its
      # standard deviations and fit?
      if self._template is not None:
        gamma = lewis(self._template, self.cspad_img)
        p = flex.max_index(gamma)
        peak = (p // gamma.focus()[1] - self._template.focus()[0] + 1,
                p % gamma.focus()[1] - self._template.focus()[1] + 1)

        #"""
        ### REFERENCE CHECK ###
        from os.path import dirname, isdir, join
        from scipy import io

        mat_dirname = dirname(cspad_tbx.pathsubst(
          self._mat_path, evt, env, frame_number=self._nframes))
        if not isdir(mat_dirname):
          makedirs(mat_dirname)

        io.savemat(
          file_name=join(mat_dirname, 'cross-check-%05d.mat' % self._nframes),
          mdict=dict(
            image=self.cspad_img.as_numpy_array(),
            template=self._template.as_numpy_array(),
            gamma=gamma.as_numpy_array(),
            peak=numpy.array(peak)),
          appendmat=False,
          do_compression=True,
          oned_as='column')

        return
        ### REFERENCE CHECK ###
        #"""
      else:
        # Alternative: position everything with respect to the frame
        # origin.
        peak = (0, 0)

      # XXX Come up with a better way to handle the offsets!  They
      # really do depend on the template, and should therefore be
      # packaged with it.
      self.logger.info("Template registration anchor point (%d, %d)" %
                       (peak[0], peak[1]))

      roi = []
      if evt.expNum() == 208:
        # Regions of interest for L632 (experiment number 208).  XXX
        # Could perhaps migrate the template matching here instead?

        # The left, middle, and right manganese signals.  XXX Extend the
        # rightmost ROI three pixels in upward direction (see runs 145
        # and onwards, also note narrower slit)?
        roi.append((peak[0] + 59, peak[1] - 24, 12, 5))
        roi.append((peak[0] + 61, peak[1] + 28, 12, 4))
        roi.append((peak[0] + 61, peak[1] + 79, 12, 5))

        # Two background regions between the manganese spots, with the
        # same total area as the signal.
        roi.append((peak[0] + 62, peak[1] +  1, 8, 8))
        roi.append((peak[0] + 63, peak[1] + 51, 8, 8))

        # The left and right direct reflections from the Si substrate
        # (i.e. the areas between the zone plates).  These were the
        # features used for template registration.
        roi.append((peak[0], peak[1],      40, 10))
        roi.append((peak[0], peak[1] + 50, 40,  9))

        # Spot between the direct reflections.  XXX What is this?
        roi.append((peak[0] + 1, peak[1] + 23, 22, 13))

        # The horizontal slit, where the direct reflection occurs.  This
        # is fixed.  XXX Verify this!
        roi.append((22, 0, 41, 128))

        # Background stripe, below the manganese spots.  This is fixed
        # to the bottom of the detector.
        roi.append((104, 0, 20, 128))

      elif evt.expNum() == 363:
        # Regions of interest for LB68 (experiment number 363).
        # 0-pixel are active, 255-pixel are inactive
        from scipy.misc import imread

        # Dec 5, 2013 (09:00 - 21:00): initial estimates from r0010
        """
        roi.append((peak[0] +  14, peak[1] + 138 + 23, 25, 50 - 25))
        roi.append((peak[0] +  45, peak[1] + 138 + 23, 25, 50 - 25))
        roi.append((peak[0] +  78, peak[1] + 137 + 23, 25, 50 - 25))
        roi.append((peak[0] + 111, peak[1] + 137 + 23, 25, 50 - 25))
        roi.append((peak[0] + 144, peak[1] + 137 + 23, 25, 50 - 25))
        roi.append((peak[0] + 177, peak[1] + 136 + 23, 25, 50 - 25))
        roi.append((peak[0] + 210, peak[1] + 136 + 23, 25, 50 - 25))
        roi.append((peak[0] + 243, peak[1] + 136 + 23, 25, 50 - 25))
        roi.append((peak[0] + 278, peak[1] + 135 + 23, 25, 50 - 25))
        roi.append((peak[0] + 312, peak[1] + 135 + 23, 25, 50 - 25))
        roi.append((peak[0] + 344, peak[1] + 135 + 23, 25, 50 - 25))
        roi.append((peak[0] + 376, peak[1] + 135 + 23, 25, 50 - 25))
        roi.append((peak[0] + 408, peak[1] + 135 + 23, 25, 50 - 25))
        roi.append((peak[0] + 442, peak[1] + 135 + 23, 25, 50 - 25))
        roi.append((peak[0] + 475, peak[1] + 135 + 23, 25, 50 - 25))
        """

        # Dec 6, 2013 (09:00 - 21:00): rough estimates
        """
        roi.append((peak[0] + 0, peak[1] +  25, 512,  25)) # bkg
        roi.append((peak[0] + 0, peak[1] + 135, 512,  25)) # oxygen
        roi.append((peak[0] + 0, peak[1] + 160, 512,  25)) # signal
        roi.append((peak[0] + 0, peak[1] + 300, 512, 130)) # zeroth order
        """

        # Dec 7, 2013 (09:00 - 21:00): overlap between oxygen and
        # signal.  Will loose some signal.
        """
        roi.append((peak[0] + 0, peak[1] +  25, 512,  25)) # bkg
        roi.append((peak[0] + 0, peak[1] + 135, 512,  50)) # oxygen
        roi.append((peak[0] + 0, peak[1] + 185, 512,  40)) # signal
        roi.append((peak[0] + 0, peak[1] + 270, 512, 170)) # zeroth order
        """

        """
        # Dec 7 2013 (09:00 - 21:00): binary masks stored in PNG
        # images.

        roi.append((peak[0] + 0, peak[1] +  25, 512,  25)) # bkg
        roi.append((peak[0] + 0, peak[1] + 135, 512,  25)) # oxygen

        #roi_image = flex.float(
        #  imread('/reg/neh/home1/hattne/myrelease/LB68-r0039-max-mask.png',
        #         flatten=True))
        #roi_image = flex.float(
        #  imread('/reg/neh/home1/hattne/myrelease/LB68-r0039-std-mask.png',
        #         flatten=True))
        roi_image = flex.float(
          imread('/reg/neh/home1/hattne/myrelease/LB68-r0052-avg-mask.png',
                 flatten=True))
        roi_image = (255 - roi_image)

        #roi.append((0, 0, self.cspad_img.focus()[0], self.cspad_img.focus()[1]))
        roi.append(roi_image)

        roi.append((peak[0] + 0, peak[1] + 270, 512, 170)) # zeroth order
        """

        # Dec 9, 2013 (09:00 - 21:00)
        #"""
        roi.append((peak[0] + 0, peak[1] +  25, 512,  25)) # bkg
        roi.append((peak[0] + 0, peak[1] + 135, 512,  25)) # oxygen
        #roi.append((peak[0] + 0, peak[1] + 160, 512,  25)) # signal
        roi_image = flex.float(
          imread('/reg/neh/home1/hattne/myrelease/LB68-r0067-max-mask.png',
                 flatten=True))
        roi.append(roi_image)

        roi.append((peak[0] + 0, peak[1] + 240, 512, 180)) # zeroth order
        #"""

      else:
        self.logger.error(
          "No regions of interest for %s (experiment number %d)" % (
            env.experiment(), evt.expNum()))

      # Clip the regions of interest to the actual image.  If the ROI
      # does not overlap with the image at all, set its width and
      # height to zero.  XXX Do the integration here as well?
      for i in range(len(roi)):
        if not isinstance(roi[i], tuple):
          continue

        r = roi[i]
        if    r[0] + r[2] < 0 or r[0] >= self.cspad_img.focus()[0] or \
              r[1] + r[3] < 0 or r[1] >= self.cspad_img.focus()[1]:
          roi[i] = (r[0], r[1], 0, 0)
          continue

        r = roi[i]
        if r[0] < 0:
          roi[i] = (0, r[1], r[2] + r[0], r[3])

        r = roi[i]
        if r[1] < 0:
          roi[i] = (r[0], 0, r[2], r[3] + r[1])

        r = roi[i]
        if r[0] + r[2] > self.cspad_img.focus()[0]:
          roi[i] = (r[0], r[1], self.cspad_img.focus()[0] - r[0], r[3])

        r = roi[i]
        if r[1] + r[3] > self.cspad_img.focus()[1]:
          roi[i] = (r[0], r[1], r[2], self.cspad_img.focus()[1] - r[1])

      # Sum up intensities in all regions of interest, and keep track
      # of the actual number of pixels summed.  The common_mode module
      # takes care of dark-subtraction.  XXX Would like to estimate
      # sigma for spot, like in spotfinder/LABELIT.
      I = flex.double(len(roi))
      I_nmemb = flex.int(len(roi))
      for i in range(len(roi)):
        if isinstance(roi[i], flex.float):
          sel = roi[i].as_1d() < 128
          I[i] = flex.sum(self.cspad_img.as_1d().select(sel))
          I_nmemb[i] = sel.count(True)
          continue

        if roi[i][2] <= 0 or roi[i][3] <= 0:
          I[i] = 0
          I_nmemb[i] = 0
        else:
          I[i] = flex.sum(self.cspad_img.matrix_copy_block(
              i_row=roi[i][0],
              i_column=roi[i][1],
              n_rows=roi[i][2],
              n_columns=roi[i][3]))
          I_nmemb[i] = roi[i][2] * roi[i][3]
          """
          # Sanity check: white out the region of interest.
          self.cspad_img.matrix_paste_block_in_place(
            block=flex.double(flex.grid(roi[i][2], roi[i][3])),
            i_row=roi[i][0],
            i_column=roi[i][1])
          """

      acq_apd_sum = sum(
        filter(lambda x: not math.isnan(x) and x > 0,
               self._acq_apd_integral.select(self._hit)))
      acq_opto_diode_sum = sum(
        filter(lambda x: not math.isnan(x) and x > 0,
               self._acq_opto_diode_integral.select(self._hit)))

      acq_apd_sum_all = sum(
        filter(lambda x: not math.isnan(x) and x > 0,
               self._acq_apd_integral))
      acq_opto_diode_sum_all = sum(
        filter(lambda x: not math.isnan(x) and x > 0,
               self._acq_opto_diode_integral))

      # Append the data point to the stream: shots, hits, energy, and
      # I.  XXX OrderedDict requires Python 2.7, could fall back on
      # regular Dict at the price of non-deterministic column order.
      from collections import OrderedDict
      csv_dict = OrderedDict([
        ('n_frames', self._hit.size()),
        ('n_hits', hits),
        ('I0', I0),
        ('I0_all', I0_all),
        ('fee_before_all', fee_before_all),
        ('fee_after_all', fee_after_all),
        ('energy_mean', energy_mean),
        ('acq_apd_sum', acq_apd_sum),
        ('acq_apd_sum_all', acq_apd_sum_all),
        ('acq_opto_diode_sum', acq_opto_diode_sum),
        ('acq_opto_diode_sum_all', acq_opto_diode_sum_all)])
      for (i, item) in enumerate(zip(roi, I, I_nmemb)):
        key = 'roi_' + ('bkg', 'oxygen', 'manganese', 'zeroth_order')[i]
        csv_dict['%s_nmemb' % key] = item[2]

        if isinstance(item[0], tuple):
          csv_dict['%s_ss_start' % key] = item[0][0]
          csv_dict['%s_fs_start' % key] = item[0][1]
          csv_dict['%s_ss_size' % key] = item[0][2]
          csv_dict['%s_fs_size' % key] = item[0][3]
        else:
          csv_dict['%s_ss_start' % key] = 0
          csv_dict['%s_fs_start' % key] = 0
          csv_dict['%s_ss_size' % key] = item[0].focus()[0]
          csv_dict['%s_fs_size' % key] = item[0].focus()[1]

        csv_dict['%s_I' % key] = item[1]

      # XXX assert that keys match up with what's in the file already?
      # Or exploit the error-reporting mechanism already implemented?
      # Write the header.  XXX How to control the order of the
      # columns?
      if not hasattr(self, '_csv'):
        from csv import DictWriter
        self._csv = DictWriter(self._stream_table, list(csv_dict.keys()))
        self._csv.writerow({key: key for key in csv_dict.keys()})
      self._csv.writerow(csv_dict)

      # Output the non-normalised image and all other relevant data to
      # a binary MATLAB file.  XXX What if scipy is not available?
      from os import makedirs, path
      from scipy import io

      mat_path = cspad_tbx.pathsubst(
        self._mat_path, evt, env, frame_number=self._nframes)
      if not path.isdir(path.dirname(mat_path)):
        makedirs(path.dirname(mat_path))

      io.savemat(
        file_name=mat_path,
        mdict=dict(
          DATA=self.cspad_img.as_numpy_array(),
          DIODES=numpy.array((acq_apd_sum, acq_apd_sum_all,
                              acq_opto_diode_sum, acq_opto_diode_sum_all)),
          ENERGY=energy_mean,
          HITS=numpy.array((hits, self._hit.size())),
          I0=numpy.array((I0, I0_all)),
          INTENSITIES=numpy.array(I),
          ROIS=numpy.array([r for r in roi if isinstance(r, tuple)])),
        appendmat=False,
        do_compression=True,
        oned_as='column')

      # Optionally update the image in the viewer.  See mod_view.
      if self._display:
        from time import localtime, strftime

        # Copy over regions of interest to shared multiprocessing
        # array.  XXX Flip to honour wxPython convention.
        for i in range(len(roi)):
          if not isinstance(roi[i], tuple):
            continue
          self._roi[4 * i + 0] = roi[i][1]
          self._roi[4 * i + 1] = roi[i][0]
          self._roi[4 * i + 2] = roi[i][3]
          self._roi[4 * i + 3] = roi[i][2]

        time_str = strftime("%H:%M:%S", localtime(evt.getTime().seconds()))
        title = "r%04d@%s: frame %d on %s" \
            % (evt.run(), time_str, self._nframes, self.address)

        # XXX No distance in the Andor experiment.  So don't bother
        # with the fictional beam center, distance, and saturation
        # value?  See also mod_average.endjob()
        img_obj = (dict(BEAM_CENTER=(0, 0),
                        DATA=self.cspad_img,
                        DETECTOR_ADDRESS=self.address,
                        DISTANCE=10, # XXX Evil kludge to keep dxtbx happy!
                        PIXEL_SIZE=13.5e-3, # XXX Hard-coded, again!
                        SATURATED_VALUE=10000,
                        TIME_TUPLE=cspad_tbx.evt_time(evt),
                        WAVELENGTH=12398.4187 / energy),
                   title)

        while not self._queue.empty():
          if not self._proc.is_alive():
            evt.setStatus(Event.Stop)
            return
        while True:
          try:
            self._queue.put(img_obj, timeout=1)
            break
          except Exception: #Queue.Full:
            pass

      self._reset_counters()
      return


  #signature for pyana:
  #def endjob(self, env):

  #signature for psana:
  #def endjob(self, evt, env):

  def endjob(self, obj1, obj2=None):
    """The endjob() function terminates the viewer process by sending
    it a @c None object, and waiting for it to finish.

    @param evt Event object (psana only)
    @param env Environment object
    """

    if obj2 is None:
      env = obj1
    else:
      evt = obj1
      env = obj2

    # Close the stream.
    self._stream_table.close()

    self.logger.info("XXX We gotta get out into space")

    # Optionally, wait for the viewer to shut down, see mod_view.
    if self._display:
      try:
        self._queue.put(None)
      except Exception:
        pass
      self._proc.join()
