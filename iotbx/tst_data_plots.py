
from iotbx import data_plots
import libtbx.load_env
import os

def exercise () :
  loggraph1 = """\
$TABLE: Resolution shell statistics
$GRAPHS
:R-free vs. resolution
:A:1,3:
:FOM vs. resolution
:A:1,4:
$$
1/resol^2  Nrefl      R-free     FOM       $$
$$
0.02       2004       0.25       0.89
0.04       2084       0.23       0.88
0.06       2037       0.27       0.83
0.08       1949       0.28       0.75
0.1        1783       0.38       *
$$
"""
  t = data_plots.table_data(None)
  t.import_loggraph(loggraph1)
#  print t.format_loggraph()
#  print "---"
#  print loggraph1
  assert len(t.data) == 4
  assert t.data[0] == [0.02, 0.04, 0.06, 0.08, 0.10]
  assert t.data[3][4] is None
  assert t.format_loggraph() == loggraph1
  f = open("_tst_data_plots.log", "w")
  f.write("\nRandom non-loggraph text\n\n")
  f.write(loggraph1)
  f.write("\n\n")
  f.write(loggraph1)
  f.close()
  tables = data_plots.import_ccp4i_logfile("_tst_data_plots.log")
  assert len(tables) == 2
  assert tables[0].format_loggraph() == loggraph1
  assert tables[0].format_loggraph() == tables[1].format_loggraph()
  t2 = data_plots.table_data(
    title = "Resolution shell statistics",
    column_labels = ["1/resol^2", "Nrefl", "R-free", "FOM"],
    graph_names = ["R-free vs. resolution", "FOM vs. resolution"],
    graph_columns = [[0,2], [0,3]],
    data = [[0.02, 0.04, 0.06, 0.08, 0.10],
            [2004, 2084, 2037, 1949, 1783],
            [0.25, 0.23, 0.27, 0.28, 0.38],
            [0.89, 0.88, 0.83, 0.75, None]]
  )
  #print loggraph1
  #print "---"
  #print t2.format_loggraph()
  assert t2.format_loggraph() == loggraph1
  g1 = t2.get_graph("R-free vs. resolution")
  g2 = t2.get_graph("FOM vs. resolution")
  assert len(g1.data) == 2 and len(g2.data) == 2
  p1 = g1.get_plots(fill_in_missing_y=None)
  p2 = g2.get_plots(fill_in_missing_y=None)
  assert len(p1) == 1 and len(p2) == 1
  (plot_x, plot_y) = p1.pop()
  assert len(plot_x) == 5
  (plot_x, plot_y) = p2.pop()
  assert len(plot_x) == 4
  formatted_table = """\
  -------------------------------------------------
  | Resolution shell statistics                   |
  |-----------------------------------------------|
  | 1/resol^2 | Nrefl     | R-free    | FOM       |
  |-----------------------------------------------|
  | 0.02      | 2004      | 0.25      | 0.89      |
  | 0.04      | 2084      | 0.23      | 0.88      |
  | 0.06      | 2037      | 0.27      | 0.83      |
  | 0.08      | 1949      | 0.28      | 0.75      |
  | 0.1       | 1783      | 0.38      | *         |
  -------------------------------------------------
"""
  #print t2.format()
  assert t2.format(indent=2) == formatted_table
  simple_table = """\
  Resolution shell statistics
  1/resol^2 Nrefl     R-free    FOM
  0.02      2004      0.25      0.89
  0.04      2084      0.23      0.88
  0.06      2037      0.27      0.83
  0.08      1949      0.28      0.75
  0.1       1783      0.38      *
"""
  # as above, but without indentation
  simpler_table = """\
Resolution shell statistics
1/resol^2 Nrefl     R-free    FOM
0.02      2004      0.25      0.89
0.04      2084      0.23      0.88
0.06      2037      0.27      0.83
0.08      1949      0.28      0.75
0.1       1783      0.38      *
"""
  assert t2.format_simple(indent=2) == simple_table
  assert str(t2) == simpler_table
  log_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/tracking/scala.log",
    test=os.path.isfile)
  if log_file is not None :
    tables = data_plots.import_ccp4i_logfile(log_file)
    assert ([ t.title for t in tables ] ==
      ['>>> Scales v rotation range, Unspecified ',
       'Analysis against Batch, Unspecified ',
       'Analysis against resolution , Unspecified ',
       'Analysis against intensity, Unspecified ',
       'Completeness, multiplicity, Rmeas v. resolution, Unspecified ',
       'Correlations within dataset, Unspecified ',
       'Axial reflections, axis h, Unspecified ',
       'Axial reflections, axis k, Unspecified ',
       'Axial reflections, axis l, Unspecified ',
       'Run     1, standard deviation v. Intensity, Unspecified '])
  print "OK"

if __name__ == "__main__" :
  exercise()
#---end
