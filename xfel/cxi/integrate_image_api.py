from __future__ import division
from __future__ import print_function
import os

def integrate_one_image(data, **kwargs):
  from xfel.cxi.display_spots import run_one_index_core
  from labelit.dptbx.error import NoAutoIndex
  from libtbx.utils import Sorry
  from spotfinder.exception import SpotfinderError
  from labelit.exception import AutoIndexError
  from iotbx.detectors.cspad_detector_formats import detector_format_version as detector_format_function
  from iotbx.detectors.cspad_detector_formats import reverse_timestamp

  basename = kwargs.get("integration_basename")
  if (basename is None):
    basename = ""

  dirname  = kwargs.get("integration_dirname")
  if (dirname is None):
    dirname = "integration"
  if (not os.path.isdir(dirname)):
    import errno
    try:
      os.makedirs(dirname)
    except OSError as exc:
      if exc.errno==errno.EEXIST: pass
  path = os.path.join(dirname, basename          \
                        +      data['TIMESTAMP'] \
                        +      ("_%05d.pickle" % data['SEQUENCE_NUMBER']))

  args = ["indexing.data=dummy",
          "beam_search_scope=0.5",
          "lepage_max_delta = 3.0",
          "spots_pickle = None",
          "subgroups_pickle = None",
          "refinements_pickle = None",
          "rmsd_tolerance = 5.0",
          "mosflm_rmsd_tolerance = 5.0",
          "indexing.completeness_pickle=%s"%path,
          "difflimit_sigma_cutoff=2.0",
          #"indexing.open_wx_viewer=True"
          ]

  detector_format_version = detector_format_function(
    data['DETECTOR_ADDRESS'], reverse_timestamp(data['TIMESTAMP'])[0])
  args += ["distl.detector_format_version=%s" % detector_format_version]

  from xfel.phil_preferences import load_cxi_phil
  horizons_phil = load_cxi_phil(data["xtal_target"], args)
  horizons_phil.indexing.data = data
  print("XFEL processing: %s"%path)
  try:
    return run_one_index_core(horizons_phil)
  except NoAutoIndex as e:
    print("NoAutoIndex", data['TIMESTAMP'], e)
    info = e.info
  except AutoIndexError as e:
    print("FailedAutoIndex", data['TIMESTAMP'], e)
    info = e.info
  except Sorry as e:
    print("Sorry", data['TIMESTAMP'], e)
    info = e.info
  except ZeroDivisionError as e:
    print("ZeroDivisionError", data['TIMESTAMP'], e)
    info = e.info
  except SpotfinderError as e:
    print("Too few spots from Spotfinder", data['TIMESTAMP'], e)
    info = e.info
  except Exception as e:
    print("ANOTHER exception", data['TIMESTAMP'], e)
    import traceback
    traceback.print_exc()
    info = e.info

  # return number of spotfinder spots
  try:
    return len(info.organizer.S.images[info.organizer.frames[0]]['spots_total'])
  except Exception as e:
    print("Couldn't find spotfinding results", data['TIMESTAMP'])

if __name__=="__main__":
  pass
