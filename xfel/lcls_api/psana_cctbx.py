from __future__ import absolute_import, division, print_function
from dxtbx.format.FormatXTCCspad import FormatXTCCspad, cspad_locator_scope
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

class FormatXTCCspadExtern(FormatXTCCspad):
  """ Class to stub out and override FormatXTCCspad methods to allow in-memory processing """
  def __init__(self, run, detector, **kwargs):
    FormatMultiImageLazy.__init__(self, **kwargs)
    FormatStill.__init__(self, None, **kwargs)
    Format.__init__(self, None, **kwargs)

    self.event = None
    self._psana_runs = {run.run():run}
    self._psana_det = {run.run():detector}
    self.params = cspad_locator_scope.extract()
    self.params.detector_address = [str(detector.name)]
    self.params.cspad.detz_offset = 580
    #self.params.wavelength_offset = -0.02386
    self._cache_psana_pedestals()
    self._cache_psana_gain()

  def get_raw_data(self, index=None):
    if index is None: index = 0
    return super(FormatXTCCspadExtern, self).get_raw_data(index)

  def _get_event(self, index):
    return self.event

class SimpleScript(DSPScript):
  """ Override dials.stills_process Script class to use its load_reference_geometry function """
  def __init__(self, params):
    self.params = params

class PsanaEventProcessor(object):
  """ Processor class for psana events """
  def __init__(self, params_filename, output_tag, logfile = None):
    """
    @param params_filename cctbx.xfel/DIALS parameter file for processing
    @output_tag String that will prefix output files
    @logfile File name for logging
    """
    dials_params = phil_scope.fetch(parse(file_name='params.phil')).extract()
    simple_script = SimpleScript(dials_params)
    simple_script.load_reference_geometry()
    self.reference_detector = simple_script.reference_detector
    self.processor = Processor(dials_params, output_tag)
    self.output_tag = output_tag

    if logfile is None:
      logfile = output_tag + ".log"
    log.config(logfile=logfile)

  def setup_run(self, run, psana_detector):
    """ Initialize processing for a given run
    @param run psana Run object
    @param psana_detector psana Detector object
    """
    if psana_detector.is_cspad():
      self.dxtbx_img = FormatXTCCspadExtern(run, psana_detector)
    else:
      raise RuntimeError('Unrecognized detector %s'%psana_detector.name)
    self.imageset = ImageSet(ImageSetData(MemReader([self.dxtbx_img]), None))

  def process_event(self, event, event_tag):
    """ Process a single psana event
    @param event psana Event object
    @param event_tag string identifying the event
    """
    self.dxtbx_img.event = event
    self.imageset.set_beam(self.dxtbx_img.get_beam())
    self.imageset.set_detector(self.dxtbx_img.get_detector())
    experiments = ExperimentListFactory.from_imageset_and_crystal(self.imageset, None)

    if self.reference_detector is not None:
      experiment = experiments[0]
      sync_geometry(self.reference_detector.hierarchy(), self.imageset.get_detector().hierarchy())
      experiment.detector = self.imageset.get_detector()

    self.processor.process_experiments('%s_%s'%(self.output_tag, event_tag), experiments)

  def finalize(self):
    """ Write out final output if composite_mode is True, plus do any other final steps """
    self.processor.finalize()
