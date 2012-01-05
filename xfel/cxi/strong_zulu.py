import os

def integrate_one_image(data):
  from display_spots import run_one_index_core
  from labelit.dptbx.error import NoAutoIndex
  from libtbx.utils import Sorry

  path = "strong/" + data['TIMESTAMP'] + "_" + "%05d" % data['SEQUENCE_NUMBER'] + ".pickle"
  args = ["indexing.data=dummy",
          "distl.detector_format_version=CXI 3.2",
          "force_method2_resolution_limit=2.2",
          "distl_highres_limit=2.2",
          "distl.res.outer=2.0",
          "beam_search_scope=0.5",
          "target_cell=38,78,78,90,90,90",
          "lepage_max_delta = 3.0",
          "spots_pickle = None",
          "subgroups_pickle = None",
          "refinements_pickle = None",
          "rmsd_tolerance = 5.0",
          "mosflm_rmsd_tolerance = 5.0",
          "known_setting = 9",
          "mosaicity_limit=1.0",
          "indexing.completeness_pickle=%s"%path,
          "difflimit_sigma_cutoff=2.0",
          "indexing.verbose_cv=True",
          #"indexing.open_wx_viewer=True"
          ]

  from spotfinder.applications.xfel import cxi_phil
  horizons_phil = cxi_phil.cxi_versioned_extract(args)
  horizons_phil.indexing.data = data

  try:
    run_one_index_core(horizons_phil)
  except NoAutoIndex,e:
    print "NoAutoIndex"
    print e
  except Sorry,e:
    print e
  except ZeroDivisionError,e:
    print "ZeroDivisionError"
    print e
  except Exception,e:
    print "ANOTHER exception"
    print e

if __name__=="__main__":

  dirname = "/net/sunbird/raid1/sauter/rawdata/L291/lyso/sept02/new_test/"
  files = """2011-02-20T19:15Z39:482_00734.pickle
2011-02-20T19:18Z25:271_00708.pickle
2011-02-20T19:20Z12:239_00917.pickle
2011-02-20T19:14Z42:365_00021.pickle
2011-02-20T19:17Z33:519_00155.pickle
2011-02-20T19:17Z40:287_00358.pickle
2011-02-20T19:18Z34:937_00998.pickle
2011-02-20T19:13Z34:768_00993.pickle
2011-02-20T19:14Z02:108_00813.pickle
2011-02-20T19:20Z30:520_00465.pickle
2011-02-20T19:20Z47:900_00987.pickle
2011-02-20T19:15Z21:798_00204.pickle
2011-02-20T19:15Z16:423_00042.pickle
2011-02-20T19:16Z35:167_00405.pickle
2011-02-20T19:18Z33:220_00946.pickle
2011-02-20T19:20Z31:811_00504.pickle
2011-02-20T19:18Z53:079_00542.pickle
2011-02-20T19:17Z52:237_00717.pickle
2011-02-20T19:14Z58:115_00493.pickle
2011-02-20T19:20Z42:534_00826.pickle
2011-02-20T19:20Z42:534_00826.pickle
2011-02-20T19:14Z24:074_00472.pickle
2011-02-20T19:15Z46:566_00947.pickle
2011-02-20T19:19Z50:902_00277.pickle
2011-02-20T19:14Z37:340_00870.pickle
2011-02-20T19:20Z17:638_00079.pickle
2011-02-20T19:20Z19:013_00120.pickle
2011-02-20T19:13Z15:453_00413.pickle
2011-02-20T19:11Z17:086_00862.pickle
2011-02-20T19:19Z39:021_00920.pickle
2011-02-20T19:18Z13:570_00357.pickle
2011-02-20T19:15Z46:974_00959.pickle
2011-02-20T19:19Z26:906_00557.pickle
2011-02-20T19:15Z22:856_00235.pickle
2011-02-20T19:14Z51:049_00281.pickle
2011-02-20T19:19Z17:505_00275.pickle
2011-02-20T19:18Z06:395_00141.pickle
2011-02-20T19:19Z01:329_00789.pickle
2011-02-20T19:18Z27:104_00763.pickle
2011-02-20T19:18Z01:195_00985.pickle
2011-02-20T19:18Z56:721_00651.pickle
2011-02-20T19:19Z28:822_00614.pickle
2011-02-20T19:14Z40:707_00971.pickle""".split("\n")



  def get_paths_single():
      for file in files:
          path = os.path.join(dirname,"out_shot%s"%file)
          yield path


  get_paths = get_paths_single

  for path in get_paths():
    print path
    from libtbx import easy_pickle
    L = easy_pickle.load(path)
    integrate_one_image(L)
