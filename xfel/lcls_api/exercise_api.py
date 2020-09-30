from __future__ import absolute_import, division, print_function
import psana
from xfel.lcls_api.psana_cctbx import PsanaEventProcessor

def run_psana_cctbx_processor(experiment, run, params_file, event_num = None):
  """ Demo using the cctbx/lcls api
  @param experiment LCLS experiment string
  @param run Run number
  @param params_file cctbx/DIALS parameter file for processing
  @param event_num Optional index for specific event to process
  """
  output_tag = '%s_run%d'%(experiment, run)
  print("Getting datasource")
  ds = psana.DataSource('exp=%s:run=%d'%(experiment, run))
  print("Getting detector")

  processor = PsanaEventProcessor(params_file, output_tag, logfile = output_tag + ".log")

  for run in ds.runs():
    print(run.run())
    det = psana.Detector('DsdCsPad')
    processor.setup_run(run, det)
    for event_id, event in enumerate(ds.events()):
      print(event_id)
      if event_num is not None and event_id != event_num: continue
      processor.process_event(event, str(event_id))
      if event_num is not None: break
    if event_num is not None: break
  processor.finalize()

if __name__ == '__main__':
  import sys
  experiment, run, params_file = sys.argv[1:4]
  if len(sys.argv) == 5:
    event_num = int(sys.argv[4])
  else:
    event_num = None
  run_psana_cctbx_processor(experiment, int(run), params_file, event_num)
