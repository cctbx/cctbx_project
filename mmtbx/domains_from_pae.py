from __future__ import division, print_function
import sys
from scitbx.array_family import flex
from libtbx import group_args
from mmtbx.process_predicted_model import get_indices_as_ranges

##############################################################################
#################  domains_from_pae  by Tristan Croll#########################
##############################################################################
"""
This license applies to the routines: parse_pae_file,
domains_from_pae_matrix_networkx, domains_from_pae_matrix_igraph
in this file

These routines are from:
https://github.com/tristanic/pae_to_domains

MIT License

Copyright (c) 2021 Tristan Croll

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
def parse_pae_file(pae_json_file):
    import json, numpy

    with open(pae_json_file, 'rt') as f:
        data = json.load(f)[0]

    r1, d = data['residue1'],data['distance']

    size = max(r1)

    matrix = numpy.empty((size,size))

    matrix.ravel()[:] = d

    return matrix

def domains_from_pae_matrix_networkx(pae_matrix, pae_power=1,
       pae_cutoff=5, graph_resolution=1):
    '''
    Takes a predicted aligned error (PAE) matrix representing the predicted error in distances between each
    pair of residues in a model, and uses a graph-based community clustering algorithm to partition the model
    into approximately rigid groups.

    Arguments:

        * pae_matrix: a (n_residues x n_residues) numpy array. Diagonal elements should be set to some non-zero
          value to avoid divide-by-zero warnings
        * pae_power (optional, default=1): each edge in the graph will be weighted proportional to (1/pae**pae_power)
        * pae_cutoff (optional, default=5): graph edges will only be created for residue pairs with pae<pae_cutoff
        * graph_resolution (optional, default=1): regulates how aggressively the clustering algorithm is. Smaller values
          lead to larger clusters. Value should be larger than zero, and values larger than 5 are unlikely to be useful.

    Returns: a series of lists, where each list contains the indices of residues belonging to one cluster.
    '''
    try:
        import networkx as nx
    except ImportError:
        print('ERROR: This method requires NetworkX (>=2.6.2) to be installed. Please install it using "pip install networkx" '
            'in a Python >=3.7 environment and try again.')
        import sys
        sys.exit()
    import numpy
    weights = 1/pae_matrix**pae_power

    g = nx.Graph()
    size = weights.shape[0]
    g.add_nodes_from(range(size))
    edges = numpy.argwhere(pae_matrix < pae_cutoff)
    sel_weights = weights[edges.T[0], edges.T[1]]
    wedges = [(i,j,w) for (i,j),w in zip(edges,sel_weights)]
    g.add_weighted_edges_from(wedges)

    from networkx.algorithms import community

    try:
      clusters = community.greedy_modularity_communities(g, weight='weight',
          resolution=graph_resolution)
    except Exception as e: # run without resolution
      clusters = community.greedy_modularity_communities(g, weight='weight')
    return clusters

def domains_from_pae_matrix_igraph(pae_matrix, pae_power=1, pae_cutoff=5, graph_resolution=1):
    '''
    Takes a predicted aligned error (PAE) matrix representing the predicted error in distances between each
    pair of residues in a model, and uses a graph-based community clustering algorithm to partition the model
    into approximately rigid groups.

    Arguments:

        * pae_matrix: a (n_residues x n_residues) numpy array. Diagonal elements should be set to some non-zero
          value to avoid divide-by-zero warnings
        * pae_power (optional, default=1): each edge in the graph will be weighted proportional to (1/pae**pae_power)
        * pae_cutoff (optional, default=5): graph edges will only be created for residue pairs with pae<pae_cutoff
        * graph_resolution (optional, default=1): regulates how aggressively the clustering algorithm is. Smaller values
          lead to larger clusters. Value should be larger than zero, and values larger than 5 are unlikely to be useful.

    Returns: a series of lists, where each list contains the indices of residues belonging to one cluster.
    '''
    try:
        import igraph
    except ImportError:
        print('ERROR: This method requires python-igraph to be installed. Please install it using "pip install python-igraph" '
            'in a Python >=3.6 environment and try again.')
        import sys
        sys.exit()
    import numpy
    weights = 1/pae_matrix**pae_power

    g = igraph.Graph()
    size = weights.shape[0]
    g.add_vertices(range(size))
    edges = numpy.argwhere(pae_matrix < pae_cutoff)
    sel_weights = weights[edges.T[0], edges.T[1]]
    g.add_edges(edges)
    g.es['weight']=sel_weights

    vc = g.community_leiden(weights='weight', resolution_parameter=graph_resolution/100, n_iterations=-1)
    membership = numpy.array(vc.membership)
    from collections import defaultdict
    clusters = defaultdict(list)
    for i, c in enumerate(membership):
        clusters[c].append(i)
    clusters = list(sorted(clusters.values(), key=lambda l:(len(l)), reverse=True))
    return clusters
##############################################################################
#################  end of domains_from_pae  by Tristan Croll##################
##############################################################################


def cluster_as_selection(c, first_resno = None):
  # first_resno is residue number for selections corresponding to index of zero
  c = sorted(c)
  ranges = get_indices_as_ranges(c)
  selection_string = ""
  if first_resno is not None:
    offset = first_resno
  else:
    offset = 0
  for r in ranges:
    r_start = r.start + offset
    r_end = r.end + offset
    if not selection_string:
      selection_string = "(resseq %s:%s)"  %(r_start, r_end)
    else:
      selection_string = "%s or (resseq %s:%s)"  %(
        selection_string,r_start, r_end)
  return selection_string

def get_domain_selections_from_pae_matrix(pae_matrix = None,
      pae_file = None,
      pae_power = None, pae_cutoff = None, resolution = None,
      first_resno = None):

    # first_resno is residue number of residue with index of zero

    if pae_matrix is None:
      pae_matrix = parse_pae_file(pae_file)
    clusters = domains_from_pae_matrix_networkx(pae_matrix,
      pae_power=pae_power, pae_cutoff=pae_cutoff,
        graph_resolution=resolution) # NOTE graph_resolution not used in 2.7

    new_clusters = []
    for c in clusters:
      new_clusters.append(sorted(c))
    clusters = sorted(new_clusters,
        key = lambda c: c[0])
    selections = []
    for c in clusters:
      selections.append(cluster_as_selection(c, first_resno = first_resno))
    return selections

def write_pae_file(pae_matrix, file_name, range_to_keep = None):

  pae_as_list_of_lists = pae_matrix.tolist()
  if range_to_keep:
    n = range_to_keep.end + 1 - range_to_keep.start
    shape = (n, n)
    pae_1d = []
    for l in pae_as_list_of_lists[range_to_keep.start:range_to_keep.end+1]:
      new_l = l[range_to_keep.start:range_to_keep.end+1]
      pae_1d += new_l
  else:
    assert type(pae_matrix.tolist()[0][0])==type(float(1))

    shape=tuple(pae_matrix.shape)

    # Flatten it out
    pae_1d =pae_matrix.flatten().tolist()

  # Read in to flex array
  flex_array=flex.float(pae_1d)

  from scitbx.array_family.flex import grid
  flex_grid=grid(shape)

  n,n = shape

  # Reshape the flex array
  flex_array.reshape(flex_grid)

  # Write out array to text file as json
  residues_1 = []
  residues_2 = []
  distances = []
  for i in range(n):
    ii = i + 1
    for j in range(n):
      jj= j + 1
      residues_1.append(ii)
      residues_2.append(jj)
      distances.append(float("%.2f" %(flex_array[i,j])))

  residue_dict = {"residue1":residues_1,
                   "residue2":residues_2,
                  "distance":distances,
                  "max_predicted_aligned_error":0}
  values = [residue_dict]

  text = str(values).replace(" ","").replace("'",'"')

  f = open(file_name, 'w')
  print(text, file = f)
  f.close()
  print("Wrote json file to %s" %(file_name))


if __name__ == '__main__':
    input_file_name = sys.argv[1]
    args = group_args(
      group_args_type = 'parameters',
      pae_file = input_file_name,
      library = 'networkx',
      pae_power = 1.0,
      pae_cutoff = 5.0,
      resolution = 1.0,
      select_range = False)

    if args.select_range:
      range_to_keep = group_args(
        group_args_type  = 'range to keep',
        start = 283,  # starts with 0
        end = 364,
         )
      pae_matrix = parse_pae_file(args.pae_file)
      write_pae_file(pae_matrix, "pae.json", range_to_keep = range_to_keep)
    else:
      pae_matrix = parse_pae_file(args.pae_file)
      selections = get_domain_selections_from_pae_matrix(
        pae_matrix = pae_matrix,
        pae_power = args.pae_power, pae_cutoff = args.pae_cutoff,
        resolution = args.resolution,)
      print("Selections:")
      for s in selections:
       print(s)
