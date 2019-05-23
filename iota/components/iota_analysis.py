# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from six.moves import range, zip

'''
Author      : Lyubimov, A.Y.
Created     : 04/07/2015
Last Changed: 02/15/2019
Description : Analyzes integration results and outputs them in an accessible
              format. Includes (optional) unit cell analysis by hierarchical
              clustering (Zeldin, et al., Acta Cryst D, 2013). In case of
              multiple clusters outputs a file with list of integrated pickles
              that comprise each cluster. (The clustering module requires scipy
              and is thus currently suspended.) Populates a PHIL file for PRIME
              with information from integration results (e.g. unit cell,
              resolution, data path, etc.)
'''

import os
import numpy as np
from collections import Counter
import math

from libtbx import easy_pickle as ep
from cctbx import crystal, uctbx, statistics
from cctbx.sgtbx.lattice_symmetry import metric_subgroups

from iota import iota_version, now
import iota.components.iota_utils as util

from prime.postrefine.mod_mx import mx_handler
from prime.postrefine import mod_input


def isprop(v):
  """ Test if attribute is a property """
  return isinstance(v, property)


class AnalysisResult(object):
  pass


class Plotter(object):

  def __init__(self, params, info):

    self.info = info
    self.params = params
    self.final_objects = self.info.get_final_objects()

    self.hm_file = os.path.join(self.info.viz_base, 'heatmap.pdf')
    self.hi_file = os.path.join(self.info.viz_base, 'res_histogram.pdf')
    self.xy_file = os.path.join(self.info.viz_base, 'beamXY.pdf')

    self.font = {'fontfamily': 'sans-serif', 'fontsize': 12}

  def plot_spotfinding_heatmap(self, write_files=False):

    import matplotlib.pyplot as plt

    hlist = [i.final['sph'] for i in self.final_objects]
    alist = [i.final['spa'] for i in self.final_objects]

    ch = max(hlist) - min(hlist) + 1
    ca = max(alist) - min(alist) + 1
    ints = [(i.final['sph'], i.final['spa']) for i in self.final_objects]
    ic = Counter(ints)

    hm_data = np.zeros((ch, ca))
    for i in ic.items():
      hm_data[i[0][0] - min(hlist), i[0][1] - min(alist)] = i[1]

    rows = range(min(hlist), max(hlist) + 1)
    cols = range(min(alist), max(alist) + 1)
    row_labels = [str(i) for i in rows]
    col_labels = [str(j) for j in cols]

    fig, ax = plt.subplots()
    fig.canvas.draw()
    heatmap = plt.pcolor(hm_data, cmap='Reds')

    ax.set_yticks(np.arange(len(rows)) + .5, minor=False)
    ax.set_xticks(np.arange(len(cols)) + .5, minor=False)
    ax.set_yticklabels(row_labels, minor=False)
    ax.set_xticklabels(col_labels, minor=False)
    ax.set_xlabel('Spot area')
    ax.set_ylabel('Spot height')

    plt.gca().set_xlim(0, len(cols))
    plt.gca().set_ylim(0, len(rows))

    # Annotate
    for y in range(hm_data.shape[0]):
      for x in range(hm_data.shape[1]):
        plt.text(x + 0.5, y + 0.5, '%3d' % hm_data[y, x],
                 horizontalalignment='center',
                 verticalalignment='center',
                 )

    if write_files:
      fig.savefig(self.hm_file, format='pdf', bbox_inches=0)
    else:
      plt.show()

  def calculate_beam_xy(self):
    """ calculates beam xy and other parameters """
    info = []

    # Import relevant info
    pixel_size = self.info.pixel_size
    for i in [j.final for j in self.final_objects]:
      try:
        info.append([i, i['beamX'], i['beamY'], i['wavelength'], i['distance'],
                     (i['a'], i['b'], i['c'], i['alpha'], i['beta'],
                      i['gamma'])])
      except IOError as e:
        print('IOTA ANALYSIS ERROR: BEAMXY failed! ', e)
        pass

    # Calculate beam center coordinates and distances
    beamX = [i[1] for i in info]
    beamY = [j[2] for j in info]
    beam_dist = [math.hypot(i[1] - np.median(beamX), i[2] - np.median(beamY))
                 for i in info]
    beam_dist_std = np.std(beam_dist)
    img_list = [[i[0], i[1], i[2], i[3], i[4], i[5], j] for i, j in
                list(zip(info, beam_dist))]

    # Separate out outliers
    outliers = [i for i in img_list if i[3] > 2 * beam_dist_std]
    clean = [i for i in img_list if i[3] <= 2 * beam_dist_std]
    cbeamX = [i[1] for i in clean]
    cbeamY = [j[2] for j in clean]
    obeamX = [i[1] for i in outliers]
    obeamY = [j[2] for j in outliers]

    # Calculate median wavelength, detector distance and unit cell params from
    # non-outliers only
    wavelengths = [i[3] for i in clean]
    distances = [i[4] for i in clean]
    cells = [i[5] for i in clean]

    wavelength = np.median(wavelengths)
    det_distance = np.median(distances)
    a = np.median([i[0] for i in cells])
    b = np.median([i[1] for i in cells])
    c = np.median([i[2] for i in cells])

    # Calculate predicted L +/- 1 misindexing distance for each cell edge
    aD = det_distance * math.tan(2 * math.asin(wavelength / (2 * a)))
    bD = det_distance * math.tan(2 * math.asin(wavelength / (2 * b)))
    cD = det_distance * math.tan(2 * math.asin(wavelength / (2 * c)))

    return beamX, beamY, cbeamX, cbeamY, obeamX, obeamY, beam_dist, \
           [i[4] for i in info], aD, bD, cD, pixel_size

  def plot_beam_xy(self, write_files=False, return_values=False, threeD=False):
    """ Plot beam center coordinates and a histogram of distances from the median
        of beam center coordinates to each set of coordinates. Superpose a
        predicted mis-indexing shift by L +/- 1 (calculated for each axis).
    """

    import matplotlib.pyplot as plt

    # Get values
    beamX, beamY, cbeamX, cbeamY, obeamX, obeamY, \
    beam_dist, distances, aD, bD, cD, pixel_size = self.calculate_beam_xy()

    # Plot figure
    if threeD:
      fig = plt.figure(figsize=(8, 8))
      ax1 = fig.add_subplot(111, projection='3d')
    else:
      fig = plt.figure(figsize=(9, 13))
      gsp = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
      ax1 = fig.add_subplot(gsp[0, :], aspect='equal')

    # Calculate axis limits of beam center scatter plot
    ax1_delta = np.ceil(np.max(beam_dist))
    xmax = round(np.median(beamX) + ax1_delta)
    xmin = round(np.median(beamX) - ax1_delta)
    ymax = round(np.median(beamY) + ax1_delta)
    ymin = round(np.median(beamY) - ax1_delta)
    zmax = round(np.ceil(np.max(distances)))
    zmin = round(np.floor(np.min(distances)))

    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(ymin, ymax)
    if threeD:
      ax1.set_zlim(zmin, zmax)

    # Plot beam center scatter plot
    if threeD:
      ax1.scatter(beamX, beamY, distances, alpha=1, s=20, c='grey', lw=1)
      ax1.plot([np.median(beamX)], [np.median(beamY)], [np.median(distances)],
               markersize=8, marker='o', c='yellow', lw=2)
    else:
      ax1.scatter(cbeamX, cbeamY, alpha=1, s=20, c='grey', lw=1)
      ax1.scatter(obeamX, obeamY, alpha=1, s=20, c='red', lw=1)
      ax1.plot(np.median(beamX), np.median(beamY), markersize=8, marker='o',
               c='yellow', lw=2)

      # Plot projected mis-indexing limits for all three axes
      circle_a = plt.Circle((np.median(beamX), np.median(beamY)), radius=aD,
                            color='r', fill=False, clip_on=True)
      circle_b = plt.Circle((np.median(beamX), np.median(beamY)), radius=bD,
                            color='g', fill=False, clip_on=True)
      circle_c = plt.Circle((np.median(beamX), np.median(beamY)), radius=cD,
                            color='b', fill=False, clip_on=True)
      ax1.add_patch(circle_a)
      ax1.add_patch(circle_b)
      ax1.add_patch(circle_c)

    # Set labels
    ax1.set_xlabel('BeamX (mm)', fontsize=15)
    ax1.set_ylabel('BeamY (mm)', fontsize=15)
    if threeD:
      ax1.set_zlabel('Distance (mm)', fontsize=15)
      ax1.set_title('Beam XYZ Coordinates')
    else:
      ax1.set_title('Beam XY Coordinates')

    if not threeD:
      # Plot histogram of distances to each beam center from median
      ax2 = fig.add_subplot(gsp[1, :])
      ax2_n, ax2_bins, ax2_patches = plt.hist(beam_dist, 20, facecolor='b',
                                              alpha=0.75, histtype='stepfilled')
      ax2_height = (np.max(ax2_n) + 9) // 10 * 10
      ax2.axis([0, np.max(beam_dist), 0, ax2_height])
      ax2.set_xlabel('Distance from median (mm)', fontsize=15)
      ax2.set_ylabel('No. of images', fontsize=15)

    if write_files:
      fig.savefig(self.xy_file, format='pdf', bbox_inches=0)
    else:
      plt.show()

    if return_values:
      return np.median(beamX), np.median(beamY), pixel_size

  def plot_res_histogram(self, write_files=False):

    import matplotlib.pyplot as plt

    # Get resolution values
    hres = [i.final['res'] for i in self.final_objects]
    lres = [i.final['lres'] for i in self.final_objects]

    # Plot figure
    fig = plt.figure(figsize=(9, 13))
    gsp = gridspec.GridSpec(2, 1)
    hr = fig.add_subplot(gsp[0, :])
    hr_n, hr_bins, hr_patches = plt.hist(hres, 20, facecolor='b', alpha=0.75,
                                         histtype='stepfilled')
    hr_height = (np.max(hr_n) + 9) // 10 * 10
    hr.axis([np.min(hres), np.max(hres), 0, hr_height])
    reslim = 'High Resolution Limit ({})'.format(r'$\AA$')
    hr.set_xlabel(reslim, fontsize=15)
    hr.set_ylabel('No. of frames', fontsize=15)

    lr = fig.add_subplot(gsp[1, :])
    lr_n, lr_bins, lr_patches = plt.hist(lres, 20, facecolor='b', alpha=0.75,
                                         histtype='stepfilled')
    lr_height = (np.max(lr_n) + 9) // 10 * 10
    lr.axis([np.min(lres), np.max(lres), 0, lr_height])
    reslim = 'Low Resolution Limit ({})'.format(r'$\AA$')
    lr.set_xlabel(reslim, fontsize=15)
    lr.set_ylabel('No. of frames', fontsize=15)

    if write_files:
      fig.savefig(self.hi_file, format='pdf', bbox_inches=0)
    else:
      plt.show()


class Analyzer(object):
  """ Class to analyze integration results """

  def __init__(self, info=None, params=None, gui_mode=False):

    self.info = info
    self.params = params
    self.gui_mode = gui_mode

    # Attributes for LivePRIME override
    self.best_pg = None
    self.best_uc = None

  def get_results(self, finished_objects=None):
    if not finished_objects:
      finished_objects = self.info.get_finished_objects()
      if not finished_objects:
        return False
    final_objects = []

    if self.gui_mode:
      self.info.unplotted_stats = {}
      for key in self.info.stats:
        self.info.unplotted_stats[key] = dict(lst=[])

    for obj in finished_objects:
      if len(self.info.unprocessed) > 0:
        for item in self.info.unprocessed:
          if item[0] == obj.img_index:
            self.info.unprocessed.remove(item)
            break

      if len(self.info.categories['not_processed'][0]) > 0:
        self.info.categories['not_processed'][0].remove(obj.img_path)

      if obj.fail:
        key = obj.fail.replace(' ', '_')
        if key in self.info.categories:
          self.info.categories[key][0].append(obj.img_path)
      else:
        self.info.categories['integrated'][0].append(obj.final['final'])
        self.info.final_objects.append(obj.obj_file)
        final_objects.append(obj)

      if not obj.fail or 'triage' not in obj.fail:
        self.info.categories['have_diffraction'][0].append(obj.img_path)

    # Calculate processing stats from final objects
    if final_objects:
      self.info.pixel_size = final_objects[0].final['pixel_size']

      # Get observations from file
      try:
        all_obs = ep.load(self.info.idx_file)
      except Exception:
        all_obs = None

      # Collect image processing stats
      for obj in final_objects:
        for key in self.info.stats:
          if key in obj.final:
            stat_tuple = (obj.img_index, obj.img_path, obj.final[key])
            self.info.stats[key]['lst'].append(stat_tuple)

            if self.gui_mode:
              if key not in self.info.unplotted_stats:
                self.info.unplotted_stats[key] = dict(lst=[])
              self.info.unplotted_stats[key]['lst'].append(stat_tuple)

        # Unit cells and space groups (i.e. cluster iterable)
        self.info.cluster_iterable.append(
          [float(obj.final['a']),
           float(obj.final['b']),
           float(obj.final['c']),
           float(obj.final['alpha']),
           float(obj.final['beta']),
           float(obj.final['gamma']),
           str(obj.final['sg'])
           ])

        # Get observations from this image
        obs = None
        if 'observations' in obj.final:
          obs = obj.final['observations'].as_non_anomalous_array()
        else:
          pickle_path = obj.final['final']
          if os.path.isfile(pickle_path):
            try:
              pickle = ep.load(pickle_path)
              obs = pickle['observations'][0].as_non_anomalous_array()
            except Exception as e:
              print('IMAGE_PICKLE_ERROR for {}: {}'.format(pickle_path, e))

        with util.Capturing():
          if obs:
            # Append observations to combined miller array
            obs = obs.expand_to_p1()
            if all_obs:
              all_obs = all_obs.concatenate(obs,
                                            assert_is_similar_symmetry=False)
            else:
              all_obs = obs

            # Get B-factor from this image
            try:
              mxh = mx_handler()
              asu_contents = mxh.get_asu_contents(500)
              observations_as_f = obs.as_amplitude_array()
              observations_as_f.setup_binner(auto_binning=True)
              wp = statistics.wilson_plot(observations_as_f, asu_contents,
                                          e_statistics=True)
              b_factor = wp.wilson_b
            except RuntimeError as e:
              b_factor = 0
              print('B_FACTOR_ERROR: ', e)
            self.info.b_factors.append(b_factor)

      # Save collected observations to file
      if all_obs:
        ep.dump(self.info.idx_file, all_obs)

      # Calculate dataset stats
      for k in self.info.stats:
        stat_list = list(zip(*self.info.stats[k]['lst']))[2]
        stats = dict(lst=self.info.stats[k]['lst'],
                     median=np.median(stat_list),
                     mean=np.mean(stat_list),
                     std=np.std(stat_list),
                     max=np.max(stat_list),
                     min=np.min(stat_list),
                     cons=Counter(stat_list).most_common(1)[0][0])
        self.info.stats[k].update(stats)
      return True
    else:
      return False

  def print_results(self, final_table=None):
    """ Prints diagnostics from the final integration run. """

    assert self.info

    if not final_table:
      final_table = ["\n\n{:-^80}\n".format('ANALYSIS OF RESULTS')]

    if not self.info.categories['integrated']:
      final_table.append('NO IMAGES INTEGRATED!')
    else:
      label_lens = [len(v['label']) for k, v in self.info.stats.items()]
      max_label = int(5 * round(float(np.max(label_lens)) / 5)) + 5
      for k, v in self.info.stats.items():
        if k in ('lres', 'res', 'beamX', 'beamY'):
          continue
        line = '{: <{l}}:  max = {:<6.2f} min = {:<6.2f} ' \
               'avg = {:<6.2f} ({:<6.2f})' \
               ''.format(v['label'], v['max'], v['min'],
                         v['mean'], v['std'], l=max_label)
        final_table.append(line)

      # TODO: Figure out what to do with summary charts
      # # If more than one integrated image, plot various summary graphs
      # if len(self.info.categories['integrated']) > 1:
      #   plot = Plotter(self.params, self.info)
      #   if self.params.analysis.summary_graphs:
      #     if ( self.params.advanced.processing_backend == 'ha14' and
      #             self.params.cctbx_ha14.grid_search.type is not None
      #     ):
      #       plot.plot_spotfinding_heatmap(write_files=True)
      #     plot.plot_res_histogram(write_files=True)
      #     med_beamX, med_beamY, pixel_size = plot.plot_beam_xy(write_files=True,
      #                                                        return_values=True)
      #   else:
      #     with warnings.catch_warnings():
      #       # To catch any 'mean of empty slice' runtime warnings
      #       warnings.simplefilter("ignore", category=RuntimeWarning)
      #       beamXY_info = plot.calculate_beam_xy()
      #       beamX, beamY = beamXY_info[:2]
      #       med_beamX = np.median(beamX)
      #       med_beamY = np.median(beamY)
      #       pixel_size = beamXY_info[-1]

      final_table.append("{: <{l}}:  X = {:<4.2f}, Y = {:<4.2f}"
                         "".format('Median Beam Center',
                                   self.info.stats['beamX']['mean'],
                                   self.info.stats['beamY']['mean'],
                                   l=max_label))

      # Special entry for resolution last
      v = self.info.stats['res']
      final_table.append('{: <{l}}:  low = {:<6.2f} high = {:<6.2f} ' \
                         'avg = {:<6.2f} ({:<6.2f})' \
                         ''.format(v['label'], v['max'], v['min'],
                                   v['mean'], v['std'], l=max_label))

    for item in final_table:
      util.main_log(self.info.logfile, item, False)
    self.info.update(final_table=final_table)

  def unit_cell_analysis(self):
    """ Calls unit cell analysis module, which uses hierarchical clustering
        (Zeldin, et al, Acta D, 2015) to split integration results according to
        detected morphological groupings (if any). Most useful with preliminary
        integration without target unit cell specified. """

    # Will not run clustering if only one integration result found or if turned off
    if not self.info.categories['integrated']:
      util.main_log(self.info.logfile,
                    "\n\n{:-^80}\n".format(' UNIT CELL ANALYSIS '), True)
      util.main_log(self.info.logfile,
                    '\n UNIT CELL CANNOT BE DETERMINED!', True)

    elif len(self.info.categories['integrated']) == 1:
      unit_cell = (self.info.cluster_iterable[0][:5])
      point_group = self.info.cluster_iterable[0][6]
      util.main_log(self.info.logfile,
                    "\n\n{:-^80}\n".format(' UNIT CELL ANALYSIS '), True)
      uc_line = "{:<6} {:^4}:  {:<6.2f}, {:<6.2f}, {:<6.2f}, {:<6.2f}, " \
                "{:<6.2f}, {:<6.2f}".format('(1)', point_group,
                                            unit_cell[0], unit_cell[1],
                                            unit_cell[2],
                                            unit_cell[3], unit_cell[4],
                                            unit_cell[5])
      util.main_log(self.info.logfile, uc_line, True)

      self.info.best_pg = str(point_group)
      self.info.best_uc = unit_cell

    else:
      uc_table = []
      uc_summary = []

      if self.params.analysis.run_clustering:
        # run hierarchical clustering analysis
        from xfel.clustering.cluster import Cluster

        counter = 0
        self.info.clusters = []

        threshold = self.params.analysis.cluster_threshold
        cluster_limit = self.params.analysis.cluster_limit
        final_pickles = self.info.categories['integrated'][0]

        pickles = []
        if self.params.analysis.cluster_n_images > 0:
          import random

          for i in range(len(self.params.analysis.cluster_n_images)):
            random_number = random.randrange(0, len(final_pickles))
            if final_pickles[random_number] in pickles:
              while final_pickles[random_number] in pickles:
                random_number = random.randrange(0, len(final_pickles))
              pickles.append(final_pickles[random_number])
        else:
          pickles = final_pickles

        # Cluster from files (slow, but will keep for now)
        ucs = Cluster.from_files(pickle_list=pickles)

        # Do clustering
        clusters, _ = ucs.ab_cluster(threshold=threshold,
                                     log=False, write_file_lists=False,
                                     schnell=False, doplot=False)
        uc_table.append("\n\n{:-^80}\n" \
                        "".format(' UNIT CELL ANALYSIS '))

        # extract clustering info and add to summary output list
        if cluster_limit is None:
          if len(pickles) / 10 >= 10:
            cluster_limit = 10
          else:
            cluster_limit = len(pickles) / 10

        for cluster in clusters:
          sorted_pg_comp = sorted(cluster.pg_composition.items(),
                                  key=lambda x: -1 * x[1])
          pg_nums = [pg[1] for pg in sorted_pg_comp]
          cons_pg = sorted_pg_comp[np.argmax(pg_nums)]

          if len(cluster.members) > cluster_limit:
            counter += 1

            # Write to file
            cluster_filenames = [j.path for j in cluster.members]
            if self.params.analysis.cluster_write_files:
              output_file = os.path.join(self.info.int_base,
                                         "uc_cluster_{}.lst".format(counter))
              for fn in cluster_filenames:
                with open(output_file, 'a') as scf:
                  scf.write('{}\n'.format(fn))

              mark_output = os.path.basename(output_file)
            else:
              mark_output = '*'
              output_file = None

          else:
            mark_output = ''
            output_file = None

          # Populate clustering info for GUI display
          uc_init = uctbx.unit_cell(cluster.medians)
          symmetry = crystal.symmetry(unit_cell=uc_init,
                                      space_group_symbol='P1')
          groups = metric_subgroups(input_symmetry=symmetry, max_delta=3)
          top_group = groups.result_groups[0]
          best_sg = str(groups.lattice_group_info()).split('(')[0]
          best_uc = top_group['best_subsym'].unit_cell().parameters()
          # best_sg = str(top_group['best_subsym'].space_group_info())

          uc_no_stdev = "{:<6.2f} {:<6.2f} {:<6.2f} " \
                        "{:<6.2f} {:<6.2f} {:<6.2f} " \
                        "".format(best_uc[0], best_uc[1], best_uc[2],
                                  best_uc[3], best_uc[4], best_uc[5])
          cluster_info = {'number': len(cluster.members),
                          'pg': best_sg,
                          'uc': uc_no_stdev,
                          'filename': mark_output}
          self.info.clusters.append(cluster_info)

          # format and record output
          # TODO: How to propagate stdevs after conversion from Niggli?
          # uc_line = "{:<6} {:^4}:  {:<6.2f} ({:>5.2f}), {:<6.2f} ({:>5.2f}), "\
          #           "{:<6.2f} ({:>5.2f}), {:<6.2f} ({:>5.2f}), "\
          #           "{:<6.2f} ({:>5.2f}), {:<6.2f} ({:>5.2f})   "\
          #           "{}".format('({})'.format(len(cluster.members)), cons_pg[0],
          #                                 cluster.medians[0], cluster.stdevs[0],
          #                                 cluster.medians[1], cluster.stdevs[1],
          #                                 cluster.medians[2], cluster.stdevs[2],
          #                                 cluster.medians[3], cluster.stdevs[3],
          #                                 cluster.medians[4], cluster.stdevs[4],
          #                                 cluster.medians[5], cluster.stdevs[5],
          #                                 mark_output)
          # uc_table.append(uc_line)
          uc_table.append("{:<6}:  {} {}".format(len(cluster.members),
                                                 uc_no_stdev, mark_output))
          lattices = ', '.join(['{} ({})'.format(i[0], i[1]) for
                                i in sorted_pg_comp])
          # uc_info = [len(cluster.members), cons_pg[0], cluster.medians,
          #            output_file, uc_line, lattices]
          uc_info = [len(cluster.members), best_sg, best_uc, output_file,
                     uc_no_stdev, lattices]
          uc_summary.append(uc_info)

      else:
        # generate average unit cell
        uc_table.append("\n\n{:-^80}\n" \
                        "".format(' UNIT CELL AVERAGING (no clustering) '))
        uc_a, uc_b, uc_c, uc_alpha, \
        uc_beta, uc_gamma, uc_sg = list(zip(*self.info.cluster_iterable))
        cons_pg = Counter(uc_sg).most_common(1)[0][0]
        all_pgs = Counter(uc_sg).most_common()
        unit_cell = (np.median(uc_a), np.median(uc_b), np.median(uc_c),
                     np.median(uc_alpha), np.median(uc_beta),
                     np.median(uc_gamma))

        # Populate clustering info for GUI display
        uc_init = uctbx.unit_cell(unit_cell)
        symmetry = crystal.symmetry(unit_cell=uc_init,
                                    space_group_symbol='P1')
        groups = metric_subgroups(input_symmetry=symmetry, max_delta=3)
        top_group = groups.result_groups[0]
        best_sg = str(groups.lattice_group_info()).split('(')[0]
        best_uc = top_group['best_subsym'].unit_cell().parameters()
        # best_sg = str(top_group['best_subsym'].space_group_info())

        uc_no_stdev = "{:<6.2f} {:<6.2f} {:<6.2f} " \
                      "{:<6.2f} {:<6.2f} {:<6.2f} " \
                      "".format(best_uc[0], best_uc[1], best_uc[2],
                                best_uc[3], best_uc[4], best_uc[5])
        cluster_info = {'number': len(self.info.cluster_iterable),
                        'pg': best_sg,
                        'uc': uc_no_stdev,
                        'filename': None}
        self.info.clusters.append(cluster_info)

        # uc_line = "{:<6} {:^4}:  {:<6.2f} ({:>5.2f}), {:<6.2f} ({:>5.2f}), " \
        #           "{:<6.2f} ({:>5.2f}), {:<6.2f} ({:>5.2f}), " \
        #           "{:<6.2f} ({:>5.2f}), {:<6.2f} ({:>5.2f})   " \
        #           "{}".format('({})'.format(len(self.final_objects)), cons_pg,
        #                       np.median(uc_a), np.std(uc_a),
        #                       np.median(uc_b), np.std(uc_b),
        #                       np.median(uc_c), np.std(uc_c),
        #                       np.median(uc_alpha), np.std(uc_alpha),
        #                       np.median(uc_beta), np.std(uc_beta),
        #                       np.median(uc_gamma), np.std(uc_gamma), '')
        #
        # uc_table.append(uc_line)
        uc_table.append(uc_no_stdev)
        lattices = ', '.join(['{} ({})'.format(i[0], i[1]) for i in all_pgs])
        # uc_info = [len(self.final_objects), cons_pg, unit_cell, None,
        #            uc_line, lattices]
        uc_info = [len(self.info.cluster_iterable), best_sg, best_uc, None,
                   uc_no_stdev, lattices]
        uc_summary.append(uc_info)

      uc_table.append('\nMost common unit cell:\n')

      # select the most prevalent unit cell (most members in cluster)
      uc_freqs = [i[0] for i in uc_summary]
      uc_pick = uc_summary[np.argmax(uc_freqs)]
      uc_table.append(uc_pick[4])
      uc_table.append('\nBravais Lattices in Biggest Cluster: {}'
                      ''.format(uc_pick[5]))
      self.info.best_pg = str(uc_pick[1])
      self.info.best_uc = uc_pick[2]

      if uc_pick[3] is not None:
        self.prime_data_path = uc_pick[3]

      for item in uc_table:
        util.main_log(self.info.logfile, item, False)
      self.info.update(uc_table=uc_table)

      if self.gui_mode:
        return self.info.clusters

  def print_summary(self, write_files=True):
    """ Prints summary and appends to general log file. Also outputs some of it
        on stdout. Also writes out output list files.
    """

    assert self.info

    if not self.info.categories['integrated']:
      util.main_log(self.info.logfile, "NO IMAGES SUCCESSFULLY PROCESSSED!",
                    (not self.gui_mode))
      return

    summary = []
    summary.append("\n\n{:-^80}\n".format('SUMMARY'))
    categories = ['total', 'failed_triage', 'have_diffraction',
                  'failed_spotfinding', 'failed_indexing',
                  'failed_grid_search', 'failed_integration',
                  'failed_filter', 'integrated']
    for cat in categories:
      lst, fail, fn, _ = self.info.categories[cat]
      path = os.path.join(self.info.int_base, fn)
      if len(lst) > 0 or cat in ('integrated', 'diffraction'):
        summary.append('{: <20}: {}'.format('{} '.format(fail), len(lst)))
      with open(path, 'w') as cf:
        for item in lst:
          cf.write('{}\n'.format(item))
      if cat == 'integrated' and write_files:
        if not hasattr(self, 'prime_data_path'):
          self.prime_data_path = path

    summary.append('\n\nIOTA version {0}'.format(iota_version))
    summary.append("{}\n".format(now))

    for item in summary:
      util.main_log(self.info.logfile, "{}".format(item), False)
    self.info.update(summary=summary)

  def make_prime_input(self, filename='prime.phil', run_zero=False):
    """ Imports default PRIME input parameters, modifies correct entries and
        prints out a starting PHIL file to be used with PRIME
    """
    assert self.info

    pixel_size = self.info.pixel_size
    hres = self.info.stats['res']
    lres = self.info.stats['lres']

    # If symmetry / unit cell were not overridden from GUI, set from INFO
    if not self.best_pg:
      try:
        self.best_pg = self.info.best_pg.replace(" ", '')
      except AttributeError as e:
        print('PRIME INPUT ERROR, SPACE GROUP: ', e)
        self.best_pg = 'P1'

    if not self.best_uc:
      self.best_uc = self.info.best_uc

    # Determine crystal system from crystal symmetry
    sym = crystal.symmetry(space_group_symbol=self.best_pg)
    crystal_system = str(sym.space_group().crystal_system())

    # Determine number of images for indexing ambiguity resolution
    # My default: 1/2 of images or 300, whichever is smaller
    if len(self.info.categories['integrated']) >= 600:
      idx_ambiguity_sample = 300
      idx_ambiguity_selected = 100
    else:
      idx_ambiguity_sample = int(
        round(len(self.info.categories['integrated']) / 2))
      idx_ambiguity_selected = int(round(idx_ambiguity_sample / 3))

    # Set run number to 000 if running LivePRIME
    run_no = '000' if run_zero else '001'

    # Populate pertinent data parameters
    prime_params = mod_input.master_phil.extract()
    prime_params.run_no = os.path.join(os.path.dirname(self.prime_data_path),
                                       'prime/{}'.format(run_no))
    prime_params.data = [self.prime_data_path]
    prime_params.title = 'Auto-generated by IOTA v{} on {}' \
                         ''.format(iota_version, now)
    prime_params.scale.d_min = hres['mean']
    prime_params.scale.d_max = 8
    prime_params.postref.scale.d_min = hres['mean']
    prime_params.postref.scale.d_max = lres['max']
    prime_params.postref.crystal_orientation.d_min = hres['mean']
    prime_params.postref.crystal_orientation.d_max = lres['max']
    prime_params.postref.reflecting_range.d_min = hres['mean']
    prime_params.postref.reflecting_range.d_max = lres['max']
    prime_params.postref.unit_cell.d_min = hres['mean']
    prime_params.postref.unit_cell.d_max = lres['max']
    prime_params.postref.allparams.d_min = hres['mean']
    prime_params.postref.allparams.d_max = lres['max']
    prime_params.merge.d_min = hres['mean']
    prime_params.merge.d_max = lres['max']
    prime_params.target_unit_cell = uctbx.unit_cell(self.best_uc)
    prime_params.target_space_group = self.best_pg
    prime_params.target_crystal_system = crystal_system
    prime_params.pixel_size_mm = pixel_size
    prime_params.n_residues = 500
    prime_params.indexing_ambiguity.n_sample_frames = idx_ambiguity_sample
    prime_params.indexing_ambiguity.n_selected_frames = idx_ambiguity_selected

    # Determine which queue to run on (i.e. match IOTA queue)
    # Modify specific options based in IOTA settings
    # Queue options
    if (
            self.params.mp.method == 'lsf' and
            self.params.mp.queue is not None
    ):
      prime_params.queue.mode = 'bsub'
      prime_params.queue.qname = self.params.mp.queue

    # Number of processors (automatically, 1/2 of IOTA procs)
    prime_params.n_processors = int(self.params.mp.n_processors / 2)

    # Generate PRIME param PHIL
    prime_phil = mod_input.master_phil.format(python_object=prime_params)
    prime_file = os.path.join(self.info.int_base, filename)
    with open(prime_file, 'w') as pf:
      pf.write(prime_phil.as_str())

    return prime_phil

  def run_get_results(self, finished_objects=None):
    self.info.have_results = self.get_results(finished_objects=finished_objects)
    return self.info.have_results

  def run_all(self, get_results=True):
    if get_results:
      self.info.have_results = self.get_results()

    if self.info.have_results:
      self.print_results()
      self.unit_cell_analysis()
      self.print_summary()
      self.make_prime_input()

    return self.info
