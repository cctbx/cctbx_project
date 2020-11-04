
from __future__ import absolute_import, division, print_function
from iotbx import data_plots
import libtbx.load_env
import json
import os

def exercise_inline():
  loggraph1 = """\
$TABLE: Resolution shell statistics
$GRAPHS
:R-free vs. resolution:A:1,3:
:FOM vs. resolution:A:1,4:
$$
1/resol^2  Nrefl      R-free     FOM       $$
$$
0.02       2004       0.25       0.89
0.04       2084       0.23       0.88
0.06       *          0.27       nan
0.08       1949       0.28       0.75
0.1        1783       0.38       *
$$
"""
  # Different newline arrangement, otherwise identical
  loggraph_input = """\
$TABLE: Resolution shell statistics: $GRAPHS :R-free vs. resolution:A:1,3: :FOM vs. resolution:A:1,4:
$$
1/resol^2  Nrefl      """+"""
R-free     FOM       $$ $$
0.02       2004       """+"""
0.25       0.89
0.04       2084       0.23
       0.88
0.06       nan        0.27       nan
0.08       1949       0.28       0.75 0.1        1783       0.38       * $$
"""
  t = data_plots.table_data(None)
  t.import_loggraph(loggraph_input)
#  print t.format_loggraph()
#  print "---"
#  print loggraph1
  assert (len(t.data) == 4)
  assert (t.data[0] == [0.02, 0.04, 0.06, 0.08, 0.10])
  assert (t.data[3][4] is None)
  assert (t.format_loggraph() == loggraph1), t.format_loggraph()
  assert (t.export_rows()[-1] == ['0.1', '1783', '0.38', '*'])
  json_t = t.export_json_table()
  json_d = json.loads(json_t)
  assert (json_d['rows'] == [
    ["1/resol^2", "Nrefl", "R-free", "FOM"],
    ["0.02", "2004", "0.25", "0.89"],
    ["0.04", "2084", "0.23", "0.88"],
    ["0.06", "*", "0.27", "nan"],
    ["0.08", "1949", "0.28", "0.75"],
    ["0.1", "1783", "0.38", "*"]]), json_d['rows']
  assert (json_d['title'] == "Resolution shell statistics"), json_d['title']

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
            [2004, 2084, None, 1949, 1783],
            [0.25, 0.23, 0.27, 0.28, 0.38],
            [0.89, 0.88, float('NaN'), 0.75, None]]
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
  | 0.06      | *         | 0.27      | nan       |
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
  0.06      *         0.27      nan
  0.08      1949      0.28      0.75
  0.1       1783      0.38      *
"""
  # as above, but without indentation
  simpler_table = """\
Resolution shell statistics
1/resol^2 Nrefl     R-free    FOM
0.02      2004      0.25      0.89
0.04      2084      0.23      0.88
0.06      *         0.27      nan
0.08      1949      0.28      0.75
0.1       1783      0.38      *
"""
  formatted = t2.format_simple(indent=2)
  assert (formatted == simple_table), formatted
  assert str(t2) == simpler_table
  json_str = t2.export_json()
  json_dict = json.loads(json_str)
  expected_dict = {"graph_types": ["A", "A"], "graph_columns": [[0, 2], [0, 3]], "title": "Resolution shell statistics", "column_labels": ["1/resol^2", "Nrefl", "R-free", "FOM"], "data": [[0.02, 0.04, 0.06, 0.08, 0.1], [2004, 2084, None, 1949, 1783], [0.25, 0.23, 0.27, 0.28, 0.38], [0.89, 0.88, float('nan'), 0.75, None]], "graph_names": ["R-free vs. resolution", "FOM vs. resolution"], "x_is_inverse_d_min": False}
  for key in expected_dict:
    if key != 'data':
      assert (key in json_dict), key
      assert (json_dict[key] == expected_dict[key]), \
              (key, json_dict[key], expected_dict[key])
  assert ('"data": [[0.02, 0.04, 0.06, 0.08, 0.1], [2004, 2084, null, 1949, 1783], [0.25, 0.23, 0.27, 0.28, 0.38], [0.89, 0.88, NaN, 0.75, null]]' in
          json_str), json_str


def exercise_column_formats():
  t = data_plots.table_data(
    title = "Resolution shell statistics",
    column_labels = ["1/resol^2", "Nrefl", "R-free", "FOM"],
    column_formats = ["%.2g", "%i", "%.2f", "%.2f"],
    data = [[0.02, 0.04, 0.06, 0.08, 0.10],
            [2004, 2084, None, 1949, 1783],
            [0.25, 0.23, 0.27, 0.28, 0.38],
            [0.89, 0.88, float('NaN'), 0.75, None]]
  )
  assert t.format() == """\
-------------------------------------------------
| Resolution shell statistics                   |
|-----------------------------------------------|
| 1/resol^2 | Nrefl     | R-free    | FOM       |
|-----------------------------------------------|
| 0.02      | 2004      | 0.25      | 0.89      |
| 0.04      | 2084      | 0.23      | 0.88      |
| 0.06      | *         | 0.27      | nan       |
| 0.08      | 1949      | 0.28      | 0.75      |
| 0.1       | 1783      | 0.38      | *         |
-------------------------------------------------
"""


def exercise_logfile():
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
  loggraph3 = """\
$TABLE: Matthews coefficients:
$GRAPHS: Protein crystal computed at resolution of 2.450 :A:1,2,3,4
$$ Nmol/asym Matthews_Coeff sovlent_frac P(2.450) $$
$$
    1         3.13            0.61         0.99
    2         1.56            0.21         0.01
$$
"""
  t3 = data_plots.table_data(None)
  t3.import_loggraph(loggraph3)
  g3 = t3.get_graph("Protein crystal computed at resolution of 2.450")
  p3 = g3.get_plots()

if __name__ == "__main__" :
  exercise_column_formats()
  exercise_inline()
  exercise_logfile()
  print("OK")
