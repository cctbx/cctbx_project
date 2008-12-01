
from iotbx import data_plots

def exercise () :
  loggraph = """
$TABLE: Resolution shell statistics
$GRAPHS
:R-free vs. resolution
:A:1,2:
$$
1/resol^2   R-free   FOM $$
$$
0.02        0.25     0.89
0.04        0.23     0.88
0.06        0.27     0.83
0.08        0.28     0.75
$$
"""
  t = data_plots.table_data(None)
  t.import_loggraph(loggraph)
  

if __name__ == "__main__" :
  exercise()
