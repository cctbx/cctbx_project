from __future__ import absolute_import, division, print_function
import psana
from xfel.lcls_api.psana_cctbx import CctbxPsanaEventProcessor

def simple_example(experiment, run_number, detector_address, params_file, event_num):
  """ Demo using the cctbx/lcls api
  @param experiment LCLS experiment string
  @param run_number Run number
  @param params_file cctbx/DIALS parameter file for processing
  @param event_num Index for specific event to process
  """
  output_tag = '%s_run%d'%(experiment, run_number)
  print("Getting datasource")
  ds = psana.DataSource('exp=%s:run=%d'%(experiment, run_number))

  processor = CctbxPsanaEventProcessor(params_file, output_tag, logfile = output_tag + ".log")

  for run in ds.runs():
    print("Getting detector")
    det = psana.Detector(detector_address)
    processor.setup_run(run, det)
    for event_id, event in enumerate(ds.events()):
      print(event_id)
      if event_num is not None and event_id != event_num: continue
      processor.process_event(event, str(event_id))
      break
    break
  processor.finalize()

def full_api_example(experiment, run_number, detector_address, params_file, event_num):
  """ Demo using the cctbx/lcls api
  @param experiment LCLS experiment string
  @param run_number Run number
  @param params_file cctbx/DIALS parameter file for processing
  @param event_num Index for specific event to process
  """
  output_tag = '%s_run%d'%(experiment, run_number)
  print("Getting datasource")
  ds = psana.DataSource('exp=%s:run=%d'%(experiment, run_number))

  processor = CctbxPsanaEventProcessor(params_file, output_tag) # note, logfile already initialized in this demo, so don't do it twice

  for run in ds.runs():
    print("Getting detector")
    det = psana.Detector(detector_address)
    processor.setup_run(run, det)
    for event_id, event in enumerate(ds.events()):
      print(event_id)
      if event_num is not None and event_id != event_num: continue

      tag = '%s_%s'%(output_tag, str(event_id))
      experiments = processor.experiments_from_event(event)
      processor.tag = tag
      processor.setup_filenames(tag)
      try:
        processor.pre_process(experiments)
        observed = processor.find_spots(experiments)
        experiments, indexed = processor.index(experiments, observed)
        experiments, indexed = processor.refine(experiments, indexed)
        integrated = processor.integrate(experiments, indexed)
        print("Integrated %d spots on %d lattices"%(len(integrated), len(experiments)))
      except Exception as e:
        print("Couldn't process event %d"%event_id, str(e))
      break
    break
  processor.finalize()

if __name__ == '__main__':
  import sys
  experiment, run_number, detector_address, params_file, event_num = sys.argv[1:6]
  simple_example(experiment, int(run_number), detector_address, params_file, int(event_num))
  full_api_example(experiment, int(run_number), detector_address, params_file, int(event_num))
