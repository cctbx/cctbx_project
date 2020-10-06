from __future__ import absolute_import, division, print_function
#from dxtbx.format.FormatXTC import
from dxtbx.format.FormatXTCCspad import FormatXTCCspad, cspad_locator_scope
from dxtbx.format.FormatXTCEpix import FormatXTCEpix, epix_locator_scope
from dxtbx.format.FormatXTCJungfrau import FormatXTCJungfrau, jungfrau_locator_scope
from dxtbx.format.Format import Format
from dxtbx.format.FormatMultiImageLazy import FormatMultiImageLazy
from dxtbx.format.FormatStill import FormatStill
from libtbx.phil import parse
from dials.command_line.stills_process import Processor, phil_scope, sync_geometry, Script as DSPScript
from dxtbx.imageset import ImageSet, ImageSetData, MemReader
from dxtbx.model.experiment_list import ExperimentListFactory
from dials.util import log

"""
API for using cctbx.xfel/DIALS processing of events without using cctbx/DIALS commands directly
"""

class FormatXTCCspadSingleEvent(FormatXTCCspad):
  """ Class to stub out and override FormatXTCCspad methods to allow in-memory processing """
  def __init__(self, params, run, detector, **kwargs):
    FormatMultiImageLazy.__init__(self, **kwargs)
    FormatStill.__init__(self, None, **kwargs)
    Format.__init__(self, None, **kwargs)

    self.event = None
    self._psana_runs = {run.run():run}
    self._psana_det = {run.run():detector}
    self.params = params
    self.params.detector_address = [str(detector.name)]
    self._cache_psana_pedestals()
    self._cache_psana_gain()

  def get_raw_data(self, index=None):
    if index is None: index = 0
    return super(FormatXTCCspadSingleEvent, self).get_raw_data(index)

  def _get_event(self, index):
    return self.event

class FormatXTCEpixSingleEvent(FormatXTCEpix):
  """ Class to stub out and override FormatXTCEpix methods to allow in-memory processing """
  def __init__(self, params, run, detector, **kwargs):
    FormatMultiImageLazy.__init__(self, **kwargs)
    FormatStill.__init__(self, None, **kwargs)
    Format.__init__(self, None, **kwargs)

    self.event = None
    self._psana_det = {run.run():detector}
    self.params = params
    self.params.detector_address = [str(detector.name)]
    self._cached_detector = {}
    self._cached_psana_detectors = {run.run():detector}
    self._run = run

  def get_raw_data(self, index=None):
    if index is None: index = 0
    return super(FormatXTCEpixSingleEvent, self).get_raw_data(index)

  def _get_event(self, index):
    return self.event

  def get_run_from_index(self, index):
    return self._run

class FormatXTCJungfrauSingleEvent(FormatXTCJungfrau):
  """ Class to stub out and override FormatXTCJungfrau methods to allow in-memory processing """
  def __init__(self, params, run, detector, **kwargs):
    FormatMultiImageLazy.__init__(self, **kwargs)
    FormatStill.__init__(self, None, **kwargs)
    Format.__init__(self, None, **kwargs)

    self.event = None
    self.params = params
    self.params.detector_address = [str(detector.name)]
    self._cached_detector = {}
    self._cached_psana_detectors = {run.run():detector}
    self._beam_index = None
    self._beam_cache = None
    self._run = run

  def get_raw_data(self, index=None):
    if index is None: index = 0
    return super(FormatXTCJungfrauSingleEvent, self).get_raw_data(index)

  def _get_event(self, index):
    return self.event

  def get_run_from_index(self, index):
    return self._run

class SimpleScript(DSPScript):
  """ Override dials.stills_process Script class to use its load_reference_geometry function """
  def __init__(self, params):
    self.params = params

class CctbxPsanaEventProcessor(Processor):
  """ Processor class for psana events """
  def __init__(self, params_filename, output_tag, logfile = None):
    """
    @param params_filename cctbx.xfel/DIALS parameter file for processing
    @output_tag String that will prefix output files
    @logfile File name for logging
    """
    self.parsed_params = parse(file_name='params.phil')
    dials_params = phil_scope.fetch(self.parsed_params).extract()
    super(CctbxPsanaEventProcessor, self).__init__(dials_params, output_tag)
    simple_script = SimpleScript(dials_params)
    simple_script.load_reference_geometry()
    self.reference_detector = simple_script.reference_detector
    self.output_tag = output_tag
    self.detector_params = None

    if logfile is not None:
      log.config(logfile=logfile)

  def setup_run(self, run, psana_detector):
    """ Initialize processing for a given run
    @param run psana Run object
    @param psana_detector psana Detector object
    """
    if psana_detector.is_cspad():
      format_class = FormatXTCCspadSingleEvent
      detector_scope = cspad_locator_scope
    elif psana_detector.is_epix10ka2m():
      format_class = FormatXTCEpixSingleEvent
      detector_scope = epix_locator_scope
    elif psana_detector.is_jungfrau():
      format_class = FormatXTCJungfrauSingleEvent
      detector_scope = jungfrau_locator_scope
    else:
      raise RuntimeError('Unrecognized detector %s'%psana_detector.name)

    detector_params = detector_scope.fetch(self.parsed_params).extract()
    self.dxtbx_img = format_class(detector_params, run, psana_detector)
    self.imageset = ImageSet(ImageSetData(MemReader([self.dxtbx_img]), None))

  def process_event(self, event, event_tag):
    """ Process a single psana event
    @param event psana Event object
    @param event_tag string identifying the event
    """
    experiments = self.experiments_from_event(event)

    self.process_experiments('%s_%s'%(self.output_tag, event_tag), experiments)

  def experiments_from_event(self, event):
    """ Create an ExperimentList from a psana Event
    @param event psana Event object
    """
    self.dxtbx_img.event = event
    self.imageset.set_beam(self.dxtbx_img.get_beam())
    self.imageset.set_detector(self.dxtbx_img.get_detector())
    experiments = ExperimentListFactory.from_imageset_and_crystal(self.imageset, None)

    if self.reference_detector is not None:
      experiment = experiments[0]
      sync_geometry(self.reference_detector.hierarchy(), self.imageset.get_detector().hierarchy())
      experiment.detector = self.imageset.get_detector()
    return experiments

