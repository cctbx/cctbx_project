import os

def integrate_one_image(data, **kwargs):
  from display_spots import run_one_index_core
  from labelit.dptbx.error import NoAutoIndex
  from libtbx.utils import Sorry
  from spotfinder.exception import SpotfinderError
  from labelit.exception import AutoIndexError
  from cxi_user.xfel_targets import targets

  basename = kwargs.get("integration_basename")
  if (basename is None):
    basename = ""

  dirname  = kwargs.get("integration_dirname")
  if (dirname is None):
    dirname = "integration"
  if (not os.path.isdir(dirname)):
    os.makedirs(dirname)

  path = os.path.join(dirname, basename          \
                        +      data['TIMESTAMP'] \
                        +      ("_%05d.pickle" % data['SEQUENCE_NUMBER']))

  args = ["indexing.data=dummy",
          "distl.detector_format_version=CXI 5.1",
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
          ] + targets[data["xtal_target"]]

  from spotfinder.applications.xfel import cxi_phil
  horizons_phil = cxi_phil.cxi_versioned_extract(args)
  horizons_phil.indexing.data = data

  try:
    return run_one_index_core(horizons_phil)
  except NoAutoIndex,e:
    print "NoAutoIndex"
    print e
  except AutoIndexError,e:
    print "FailedAutoIndex"
    print e
  except Sorry,e:
    print "Sorry"
    print e
  except ZeroDivisionError,e:
    print "ZeroDivisionError"
    print e
  except SpotfinderError,e:
    print "Too few spots from Spotfinder"
    print e
  except Exception,e:
    print "ANOTHER exception"
    print e

if __name__=="__main__":
  pass
