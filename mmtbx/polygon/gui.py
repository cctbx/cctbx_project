from mmtbx import polygon, model_vs_data
import iotbx.phil

def get_histogram_data (d_min) :
  polygon_params = iotbx.phil.parse("""\
polygon {
  keys_to_show = *r_work_pdb *r_free_pdb *bonds_rmsd *angles_rmsd *adp_mean
  number_of_histogram_slots = 10
  filter {
    key = *d_min
    value_min = %.1f
    value_max = %.1f
  }
}""" % (d_min - 0.1, d_min + 0.1))
  params = polygon.master_params.fetch(sources=[polygon_params])
  return polygon.polygon(params=params.extract(),
                         d_min=d_min,
                         show_histograms=False,
                         extract_gui_data=True)

# XXX: not pickle-able - run this in GUI thread
def convert_histogram_data (polygon_result) :
  histograms = {}
  for (stat, data) in polygon_result :
    histograms[stat] = polygon.convert_to_histogram(data=data,
                                                    n_slots=10)
  return histograms

def get_stats_and_histogram_data (mvd_object) :
  pdb_file = mvd_object.pdb_file
  fmodel = mvd_object.fmodel
  d_min = fmodel.info().d_min
  model = mvd_object.models[0]
  x = model.xray_structure_stat
  g = model.model_statistics_geometry
  stats = { "r_work_pdb" : fmodel.r_work(),
            "r_free_pdb" : fmodel.r_free(),
            "adp_mean" : float(x.b_mean),
            "bonds_rmsd" : g.b_mean,
            "angles_rmsd" : g.a_mean }
  histograms = get_histogram_data(d_min=d_min)
  return stats, histograms
