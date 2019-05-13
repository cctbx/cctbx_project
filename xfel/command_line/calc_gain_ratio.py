from __future__ import division
from __future__ import print_function
from six.moves import range
# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cxi.calc_gain_ratio
#
import dxtbx, sys
from libtbx import phil
from libtbx.utils import Sorry
import psana
from scitbx.array_family import flex

help_message = """
Program to compute the gain ratio for the CSPAD detector in mixed gain mode. For each pair of adjacent pixels, if one is in high gain mode and the other in low gain mode, then the ratio between them is computed. These ratios are averaged for each panel and the average over the whole detector is reported.

Inputs: requires an averaged run in CBF format and access to psana (to extract the original gain mask for the run)

Example usage:
calc_gain_ratio.py average=cxid9114_avg-r0097.cbf experiment=cxid9114 run=97 address=CxiDs2.0:Cspad.0
"""

phil_str = """
  experiment = None
    .type = str
    .help = LCLS experiment name
  run=None
    .type = int
    .help = Run number
  average = None
    .type = str
    .help = Averaged image (can be maximum projection instead)
  address = None
    .type = str
    .help = Detector address, eg CxiDs1.0:Cspad.0
"""

phil_scope = phil.parse(phil_str)

def run(args):
  if "-c" in args or "-h" in args or "--help" in args:
    print(help_message)
  user_phil = []
  for arg in args :
    try :
      user_phil.append(phil.parse(arg))
    except RuntimeError as e :
      raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))
  params = phil_scope.fetch(sources=user_phil).extract()

  img = dxtbx.load(params.average)
  dataset_name = "exp=%s:run=%s:idx"%(params.experiment,params.run)
  ds = psana.DataSource(dataset_name)
  run = next(ds.runs())

  psana_det = psana.Detector(params.address, ds.env())
  psana_gain_mask = psana_det.gain_mask()
  psana_gain_mask = flex.bool(psana_gain_mask==1)

  gain_masks = []
  assert psana_gain_mask.focus() == (32, 185, 388)
  for i in range(32):
    gain_masks.append(psana_gain_mask[i:i+1,:,:194])
    gain_masks[-1].reshape(flex.grid(185,194))
    gain_masks.append(psana_gain_mask[i:i+1,:,194:])
    gain_masks[-1].reshape(flex.grid(185,194))

  ratios = flex.double()
  counts = flex.int()
  for panel_id, (data, mask) in enumerate(zip(img.get_raw_data(), gain_masks)):
    if mask.all_eq(True) or mask.all_eq(False):
      continue

    panel_sum = 0
    panel_count = 0
    for s in range(data.focus()[1]):
      for f in range(data.focus()[0]):
        if f+1 == data.focus()[0]:
          continue
        if (not mask[f,s]) and mask[f+1,s] and data[f+1,s] != 0:
          panel_sum += data[f,s]/data[f+1,s]
          panel_count += 1
        elif mask[f,s] and not mask[f+1,s] and data[f,s] != 0:
          panel_sum += data[f+1,s]/data[f,s]
          panel_count += 1
    if panel_count > 0:
      ratio = panel_sum/panel_count
      ratios.append(ratio)
      counts.append(panel_count)
      print("Panel", panel_id, "ratio:", ratio, "N pairs", panel_count)

  if len(ratios) <= 1:
    return
  print("Mean:", flex.mean(ratios))
  print("Standard deviation", flex.mean_and_variance(ratios).unweighted_sample_standard_deviation())

  stats = flex.mean_and_variance(ratios, counts.as_double())
  print("Weighted mean:", stats.mean())
  print("Weighted standard deviation", stats.gsl_stats_wsd())

if __name__ == "__main__":
  run(sys.argv[1:])
