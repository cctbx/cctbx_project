from __future__ import absolute_import, division, print_function
from six.moves import range
# LIBTBX_SET_DISPATCHER_NAME distl.find_active_area
import os, sys
from iotbx.detectors.npy import NpyImage
from spotfinder.core_toolbox import find_active_area

#special import of pickled NumPy array: CXI/CSPad data file
def ImageFactory(filename):
  if os.path.isfile(filename):
    I = NpyImage(filename)
    I.readHeader()
    return I

class graph_tracker:
  def has_one(self,graph):
    for key in graph.keys():
      if len(graph[key])==1:
        self.key = key
        self.item_sink = graph[key][0]
        return True
    return False
  def prune(self,graph):
    for key in graph.keys():
      try:
        graph[key].remove(self.item_sink)
      except ValueError: pass

def run_one(path, display):
  image = ImageFactory(path)
  image.read()
  data = image.linearintdata
  PC = find_active_area(data)

  sources = []; sinks = []
  for x in range(0,len(PC),2):
    if PC[x]>=0:
      sources.append((PC[x],PC[x+1]))
    else:
      sinks.append((-PC[x],-PC[x+1]))
  print(len(sources),len(sinks))
  assert len(sources)==len(sinks)
  graph = {}
  final_graph = {}
  for src in sources:
    item_sinks = [i for i in sinks if i[0]>src[0] and i[1]>src[1]]
    graph[src]=item_sinks

  G = graph_tracker()
  while G.has_one(graph):
    print(G.key, G.item_sink)
    final_graph[G.key]=G.item_sink
    del graph[G.key]
    G.prune(graph)

  assert len(graph)==0

if __name__ == "__main__":
  for arg in sys.argv[1:]:
    run_one(arg, display=True)
