""" This module is designed to provide tools to deal with groups of serial
crystallography images, treated as a graph problem.

**Author:**   Oliver Zeldin <zeldin@stanford.edu>
"""
from __future__ import absolute_import, division, print_function
from functools import reduce

__author__ = 'zeldin'

import os
import logging
from xfel.graph_proc.components import ImageNode, Edge
from xfel.clustering.cluster import Cluster
import numpy as np
import matplotlib.pyplot as plt
import random

class Graph(Cluster):
  """
  .. note::
    Developmental code. Do not use without contacting zeldin@stanford.edu

  Class for treating a graph of serial still XFEL shots as a graph.
  """
  def __init__(self, vertices, min_common_reflections=10):
    """
    Extends the constructor from cluster.Cluster to describe the cluster as a
    graph.

    :param min_common_reflections: number of reflections two images must have in
    common for an edge to be created.
    """
    Cluster.__init__(self, vertices, "Graph cluster", "made as a graph")
    self.common_miller_threshold = min_common_reflections
    self.edges = self._make_edges()


  def _make_edges(self):
    all_edges = []
    for a, vertex1 in enumerate(self.members):
      for b, vertex2 in enumerate(self.members[:a]):
        miller1 = vertex1.miller_array
        miller2 = vertex2.miller_array
        edge_weight = miller1.common_set(miller2,
                                      assert_is_similar_symmetry=False).size()
        if edge_weight >= self.common_miller_threshold:
          logging.debug("Making edge: ({}, {}) {}, {}"
                        .format(a, b, vertex1.name, vertex2.name))
          this_edge = Edge(vertex1, vertex2, edge_weight)
          logging.debug("Edge weight: {}".format(edge_weight))
          all_edges.append(this_edge)
          vertex1.edges.append(this_edge)
          vertex2.edges.append(this_edge)

    assert len(all_edges) == sum([len(v.edges) for v in self.members])/2
    return all_edges


  @classmethod
  def from_directories(cls, folder, **kwargs):
    """
    Creates a Graph class instance from a directory of files (recursively)

    :param folder: the root directory containing integration pickles.
    :return: A Graph instance.
    """
    def add_frame(all_vertices, dirpath, filename, **kwargs):
      path = os.path.join(dirpath, filename)
      this_frame = ImageNode(path, filename, **kwargs)
      if hasattr(this_frame, 'miller_array'):
        all_vertices.append(this_frame)
      else:
        logging.info('skipping file {}'.format(filename))

    nmax = kwargs.get('nmax', None)

    all_vertices = []
    for arg in folder:
      for (dirpath, dirnames, filenames) in os.walk(arg):
        for filename in filenames:
          if not nmax:
            add_frame(all_vertices, dirpath, filename, **kwargs)
          else:
            if len(all_vertices) <= nmax:
              add_frame(all_vertices, dirpath, filename)
    return cls(all_vertices)

  def make_nx_graph(self, saveas=None, show_singletons=True, force_res=False,
                    **kwargs):
    """
    Create a networkx Graph, and visualise it.

    :param saveas: optional filename to save the graph as a pdf as. Displays graph if not specified.
    :param pos: an optional networkx position dictionairy for plotting the graph
    :param force_res: if True, use the edge residual field to color them. Otherwise use edge labels if available.
    :return: the position dictionairy, so that mutiple plots can be made using the same set of positions.
    """
    import networkx as nx
    import brewer2mpl
    cols = brewer2mpl.get_map('BrBG', 'Diverging', 3).mpl_colors

    nx_graph = nx.Graph()
    # Add vertices to graph
    if show_singletons:
      nx_graph.add_nodes_from(self.members)
    else:
      nx_graph.add_nodes_from([vert for vert in self.members if vert.edges])
    logging.debug("Vertices added to graph.")

    # Add edges to graph
    for edge in self.edges:
      nx_graph.add_edge(edge.vertex_a, edge.vertex_b,
                        {'n_conn': edge.weight,
                         'residual': np.mean(np.abs(edge.residuals())),
                         'intra': edge.intra})
    logging.debug("Edges added to graph.")

    # Position vertices
    if 'pos' not in kwargs :
      kwargs['pos'] = nx.spring_layout(nx_graph, iterations=1000)

    # Organise Edge colors
    if all([e.intra is None for e in self.edges]) or force_res:
      # color by residual
      e_cols = [edge[2]["residual"]
                     for edge in nx_graph.edges_iter(data=True)]
      e_cmap = plt.get_cmap('jet')
      e_vmax = max(e_cols)
      e_vmin = 0
    else:
      # color by edge
      e_cmap = None
      e_vmin = None
      e_vmax = None
      e_cols = []
      for edge in nx_graph.edges_iter(data=True):
        intra = edge[2]["intra"]
        if intra is None:
          e_cols.append('0.5')
        elif intra is True:
          e_cols.append('black')
        elif intra is False:
          e_cols.append('red')
        else:
          raise Exception("something bad happened here")

    # Organise Vertex colors
    for v in nx_graph:
      if v.label is None:
        v.col = 0.7
      elif v.label == 0:
        v.col = cols[0]
      elif v.label == 1:
        v.col = cols[2]

    if all([v.source is not None for v in nx_graph]):
      nxlabels = {v: v.source for v in nx_graph}
    else:
      nxlabels = None

    fig = plt.figure(figsize=(12, 9))
    im_vertices = nx.draw_networkx_nodes(nx_graph, with_labels=False,
                                         node_color=[v.col for v in nx_graph],
                                         **kwargs)

    im_edges = nx.draw_networkx_edges(nx_graph,
                                      vmin=e_vmin,
                                      vmax=e_vmax,
                                      edge_color=e_cols,
                                      edge_cmap=e_cmap,
                                      **kwargs)

    # Organize labels
    if all([v.source is not None for v in nx_graph]):
      nx.draw_networkx_labels(nx_graph,labels=nxlabels, **kwargs)

    if e_cmap:
      cb = plt.colorbar(im_edges)
      cmin, cmax = cb.get_clim()
      ticks = np.linspace(cmin, cmax, 7)
      cb.set_ticks(ticks)
      cb.set_ticklabels(['{:.3f}'.format(t) for t in ticks])
      cb.set_label("Edge Residual")

    plt.axis('off')
    plt.tight_layout()
    plt.title("Network Processing\nThickness shows number of connections")
    if saveas is not None:
      plt.savefig(saveas, bbox_inches='tight')
    else:
      plt.show()
    return kwargs['pos']

  def merge_dict(self, estimate_fullies=True):
    """ Make a dict of Miller indices with  ([list of intensities], resolution)
    value tuples for each miller index. Use the fully-corrected equivalent.

    :param estimate_fullies: Boolean for if the partialities and scales should be used.
    """
    intensity_miller = {}
    for vertex in self.members:
      if vertex.edges:
        miller_array = vertex.miller_array
        if estimate_fullies:
          miller_array = miller_array / (vertex.partialities
                                         * vertex.scales)
        #miller_array = miller_array.map_to_asu()
        d_spacings = list(miller_array.d_spacings().data())
        miller_indeces = list(miller_array.indices())
        miller_intensities = list(miller_array.data())

        for observation in zip(miller_indeces, miller_intensities, d_spacings):
          try:
            intensity_miller[observation[0]][0].append(observation[1])
          except KeyError:
            intensity_miller[observation[0]] = ([observation[1]],
                                                observation[2])
    return intensity_miller

  def cc_half(self, nbins=10, estimate_fullies=True):
    """
    Calculate the correlation coefficients for the vertices of the graph.

    :param nbins: number of bins to use for binned correlations.
    :param estimate_fullies: if True, use current orientation, spot profile and scale to estimate the 'full' intensity of each reflection.

    :return cc_half, p_value: overall CC1/2 for all reflections in the cluster.
    :return pretty_string: string with the CC by bin. Lowest resolution bin is extended if data cannot be split into n_bins evenly.
    """
    from scipy.stats import pearsonr
    import operator

    intensity_miller = self.merge_dict(estimate_fullies=estimate_fullies)

    # Remove miller indices that are only measured once.
    multiple_millers = [values for values in intensity_miller.values()
                        if len(values[0]) > 1]

    # Sort, since we will be binning later
    multiple_millers.sort(key=lambda x: x[1])

    # Avoid corner case where number of millers goes exactly into the bins
    if len(multiple_millers) % nbins == 0:
      nbins += 1

    # Figure out bin sizes:
    bin_size = int(len(multiple_millers) / nbins)
    bin_edges = [multiple_millers[i][1] for i in range(0, len(multiple_millers),
                                                       bin_size)]

    # Extend the last bin does not cover the whole resolution range
    bin_edges[-1] = multiple_millers[-1][1]

    assert bin_edges[0] == multiple_millers[0][1]
    assert bin_edges[-1] == multiple_millers[-1][1]
    assert len(bin_edges) == nbins + 1

    # For bins of size bin_size, split each miller index into two, and add half
    # the observations to first_half and half to second_half
    first_half = [[] for _ in range(nbins)]
    second_half = [[] for _ in range(nbins)]
    for indx, intensities in enumerate(multiple_millers):
      # Extending the last bin if necesarry
      bin_id = indx // bin_size if indx // bin_size < nbins else nbins - 1
      obs = intensities[0]
      random.shuffle(obs)
      middle = len(obs) // 2
      first_half[bin_id].append(np.mean(obs[:middle]))
      second_half[bin_id].append(np.mean(obs[middle:]))

    # Calculate the CC for the two halves of each resolution bin.
    logging.info("Calculating CC1/2 by resolution bin. "
                 "{} miller indices per bin".format(bin_size))
    pretty_string = ''
    for bin_id in reversed(range(nbins)):
      bin_cc, bin_p = pearsonr(first_half[bin_id], second_half[bin_id])
      bin_str = "{:6.3f} - {:6.3f}: {:4.3f}".format(bin_edges[bin_id + 1],
                                                    bin_edges[bin_id],
                                                    bin_cc)
      logging.info(bin_str)
      pretty_string += bin_str + "\n"

    # Combine all the split millers, and take the CC of everything
    first_half_all = reduce(operator.add, first_half)
    second_half_all = reduce(operator.add, second_half)
    cc_half, p_value = pearsonr(first_half_all, second_half_all)
    logging.info("CC 1/2 calculated: {:4.3f}, p-value: {}".format(cc_half,
                                                                  p_value))

    return cc_half, p_value, pretty_string


  def k_means_edges(self):
    """
    Preforms k-means clustering on the edge residuals, updating the Edge.intra attribute depending on if it is in the 1st (edge.intra = True) or 2nd (edge.intra=False) cluster.
    """
    from scipy.cluster.vq import kmeans2

    # 0. Make an array of edge residuals
    edge_residuals = np.array([e.mean_residual() for e in self.edges])

    # 1. split into in-cluster and between-cluster using k-means
    labels = kmeans2(edge_residuals, 2)

    # 2. find which cluster corresponds to intra-group edges (ie. is smallest)
    if labels[0][0] < labels[0][1]:
      intra = 0
    else:
      intra = 1

    # 2. assign these labels to the edges
    for edge, group in zip(self.edges, labels[1]):
      if group == intra:
        edge.intra = True
      else:
        edge.intra = False

    logging.info(('initial edge assignment complete: {} intra edges, '
                  '{} inter edges').format(
                                  sum([e.intra == True for e in self.edges]),
                                  sum([e.intra == False for e in self.edges])))


  def label_vertices(self, max_iter=100):
    """
    Labels the vertices in the graph as members of cluster 0 or 1.
    Currently just does a BFS to get through the graph, and updates each node
    the first time it comes across is.

    Next step would be to take a majority vote for each vertex, update, and
    repeat this process until convergence.

    :param max_iter: Maximum number of rounds of 'majority vote' updates for
    vertex labeling.

    """
    from collections import deque

    if not all([e.intra is not None for e in self.edges]):
       logging.warning("k_means_edges has not been run -- edges are not \
           labled. Running this now.")
       self.k_means_edges()
    assert all([e.intra is not None for e in self.edges])

    def majority_vote(vertex):
      """ Reasign the vertex id based on the nearest neighbours.
      """
      from scipy.stats import mode
      votes = []
      for e in vertex.edges:
        other_class = e.other_vertex(vertex).label
        if e.intra:
          votes.append(other_class)
        else:
          if other_class == 0:
            votes.append(1)
          elif vertex.label == 1:
            votes.append(0)

      new_label, frequency = mode(votes)
      if new_label != vertex.label:
        logging.debug("updating vertex to {}, with {}/{} votes" \
                              .format(new_label[0], frequency[0], len(votes)))
        vertex.label = int(new_label)

    # 0. Get the most connected node.
    most_conn = max(self.members, key=lambda x: len(x.edges))
    most_conn.label = 0
    q = deque([most_conn])

    # 1. BFS, calling the near_edges each time a node is popped.
    # ToDo: when we move to a mixed model for edges, use a priority queue to
    # pick the best edge to move along.
    for v in self.members:
      v.visited = False
    curr_node = most_conn
    curr_node.visited = True
    while q:
      curr_node = q.popleft()  # Take a node of the stack
      curr_node.visited = True
      """ Take a vertex, and updates the labels of all the
      vertices around it using the argument vertex's edge labels. If a vertex
      at the other end of the edge already has a label, skip it.
      """
      assert curr_node.label is not None, "Vertex must already have a label."
      for e in curr_node.edges:
        other_v = e.other_vertex(curr_node)
        # Is the curr_node on the other end already assigned? if yes, skip
        if other_v.label is not None:
          continue
        if e.intra:
          other_v.label = curr_node.label
        elif not e.intra:
          if curr_node.label == 0:
            other_v.label = 1
          elif curr_node.label == 1:
            other_v.label = 0
      # add this node's connections to the stack
      q.extend([e.other_vertex(curr_node) for e in curr_node.edges
                if not e.other_vertex(curr_node).visited])

    logging.info(('initial vertex assignment complete: {} vertices in group 0, '
                  '{} in group 1').format(
                                    sum([v.label == 0 for v in self.members]),
                                    sum([v.label == 1 for v in self.members])))

    if logging.Logger.root.level <= logging.DEBUG:  # Extreme debug!
      self.make_nx_graph()

    # 2. Majority vote
    original_labels = [v.label for v in self.members]
    current_labels = [v.label for v in self.members]
    v_iter = 0
    while v_iter <= max_iter:
      v_iter += 1
      for v in self.members:
        majority_vote(v)
      new_labels = [v.label for v in self.members]
      if current_labels == new_labels:
        logging.info("Majority vote converged after {} iterations. {} "
                     "vertices out of {} have been changed." \
                         .format(v_iter, sum([new_labels[i] != original_labels[i]
                                              for i in range(len(new_labels))]),
                               len(new_labels)))
        break
      if iter == max_iter:
        logging.info("Majority voting did no make vertex labels "
                     "converge. final state: {}/{} different" \
                     .format(sum([new_labels[i] != original_labels[i]
                                  for i in range(len(new_labels))]),
                             len(new_labels)))

      current_labels = new_labels

