# -*- coding: utf-8 -*-
from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 04/07/2015
Last Changed: 04/04/2018
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

import cPickle as pickle
from libtbx import easy_pickle as ep
from cctbx.uctbx import unit_cell

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D

# workaround to avoid unused import warning
assert Axes3D
assert cm
assert colors

import time
assert time

import iota.components.iota_misc as misc
from iota.components.iota_misc import Capturing
from prime.postrefine import mod_input

def isprop(v):
  ''' Test if attribute is a property '''
  return isinstance(v, property)

class AnalysisResult(object):
  pass

class Plotter(object):

  def __init__(self,
               params,
               final_objects,
               viz_dir=None):

    self.final_objects = final_objects
    self.params = params
    self.hm_file = os.path.join(viz_dir, 'heatmap.pdf')
    self.xy_file = os.path.join(viz_dir, 'beam_xy.pdf')
    self.hi_file = os.path.join(viz_dir, 'res_histogram.pdf')

    self.font = {'fontfamily':'sans-serif', 'fontsize':12}

  def plot_spotfinding_heatmap(self, write_files=False):

    hlist = [i.final['sph'] for i in self.final_objects]
    alist = [i.final['spa'] for i in self.final_objects]

    ch = max(hlist) - min(hlist) + 1
    ca = max(alist) - min(alist) + 1
    ints = [(i.final['sph'], i.final['spa']) for i in self.final_objects]
    ic = Counter(ints)

    hm_data = np.zeros((ch, ca))
    for i in ic.items():
      hm_data[i[0][0]-min(hlist), i[0][1]-min(alist)] = i[1]

    rows = range(min(hlist), max(hlist) + 1)
    cols = range(min(alist), max(alist) + 1)
    row_labels = [str(i) for i in rows]
    col_labels = [str(j) for j in cols]

    fig, ax = plt.subplots()
    fig.canvas.draw()
    heatmap = plt.pcolor(hm_data, cmap='Reds')
    ax.set_yticks(np.arange(len(rows))+.5, minor=False)
    ax.set_xticks(np.arange(len(cols))+.5, minor=False)
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
    ''' calculates beam xy and other parameters '''
    info = []

    if self.params.advanced.integrate_with == 'cctbx':
      img_pickle = self.final_objects[0].final['img']
      pixel_size = pickle.load(open(img_pickle, "rb"))['PIXEL_SIZE']
    elif self.params.advanced.integrate_with == 'dials':
      proc_pickle = self.final_objects[0].final['final']
      pixel_size = pickle.load(open(proc_pickle, 'rb'))['pixel_size']

    # Import relevant info
    for i in [j.final for j in self.final_objects]:
      try:
        info.append([i, i['beamX'], i['beamY'], i['wavelength'], i['distance'],
                    (i['a'], i['b'], i['c'], i['alpha'], i['beta'], i['gamma'])])
      except IOError, e:
        pass

    # Calculate beam center coordinates and distances
    beamX = [i[1] for i in info]
    beamY = [j[2] for j in info]
    beam_dist = [math.hypot(i[1] - np.median(beamX), i[2] - np.median(beamY))
                 for i in info]
    beam_dist_std = np.std(beam_dist)
    img_list = [[i[0], i[1], i[2], i[3], i[4], i[5], j] for i, j in zip(info,
                                                                beam_dist)]

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
    ''' Plot beam center coordinates and a histogram of distances from the median
        of beam center coordinates to each set of coordinates. Superpose a
        predicted mis-indexing shift by L +/- 1 (calculated for each axis).
    '''

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
      ax2_n, ax2_bins, ax2_patches = plt.hist(beam_dist, 20, normed=False,
                                              facecolor='b', alpha=0.75,
                                              histtype='stepfilled')
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

    # Get resolution values
    hres = [i.final['res'] for i in self.final_objects]
    lres = [i.final['lres'] for i in self.final_objects]

    # Plot figure
    fig = plt.figure(figsize=(9, 13))
    gsp = gridspec.GridSpec(2, 1)
    hr = fig.add_subplot(gsp[0, :])
    hr_n, hr_bins, hr_patches = plt.hist(hres, 20, normed=False, facecolor='b',
                                         alpha=0.75, histtype='stepfilled')
    hr_height = (np.max(hr_n) + 9) // 10 * 10
    hr.axis([np.min(hres), np.max(hres), 0, hr_height])
    reslim = 'High Resolution Limit ({})'.format(r'$\AA$')
    hr.set_xlabel(reslim, fontsize=15)
    hr.set_ylabel('No. of frames', fontsize=15)

    lr = fig.add_subplot(gsp[1, :])
    lr_n, lr_bins, lr_patches = plt.hist(lres, 20, normed=False, facecolor='b',
                                         alpha=0.75, histtype='stepfilled')
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

  def __init__(self,
               init,
               all_objects,
               gui_mode = False):

    self.ver = misc.iota_version
    self.now = init.now
    self.params = init.params
    self.output_dir = init.int_base
    self.viz_dir = init.viz_base
    self.gui_mode = gui_mode

    self.logfile = init.logfile
    self.prime_data_path = None

    self.analysis_result = AnalysisResult()

    self.cons_pg = None
    self.cons_uc = None
    self.clusters = []

    # Analyze image objects
    self.all_objects = all_objects
    self.final_objects = [i for i in all_objects if i.fail == None]

    if self.final_objects is not None:
      self.sorted_final_images = sorted(self.final_objects,
                                        key=lambda i: i.final['mos'])
      self.no_diff_objects = [i for i in self.all_objects if
                              i.fail == 'failed triage']
      self.diff_objects = [i for i in self.all_objects if
                           i.fail != 'failed triage']
      if self.params.advanced.integrate_with == 'cctbx':
        self.not_int_objects = [i for i in self.all_objects if
                                i.fail == 'failed grid search']
        self.filter_fail_objects = [i for i in self.all_objects if
                                    i.fail == 'failed prefilter']
      elif self.params.advanced.integrate_with == 'dials':
        self.not_spf_objects = [i for i in self.all_objects if
                                i.fail == 'failed spotfinding']
        self.not_idx_objects = [i for i in self.all_objects if
                                i.fail == 'failed indexing']
        self.filter_fail_objects = [i for i in self.all_objects if
                                    i.fail == 'failed filter']
        self.not_int_objects = [i for i in self.all_objects if
                                i.fail == 'failed integration']

      self.pickles = [i.final['final'] for i in self.final_objects]
      self.hres = [i.final['res'] for i in self.final_objects]
      self.lres = [i.final['lres'] for i in self.final_objects]
      self.spots = [i.final['strong'] for i in self.final_objects]
      self.mos = [i.final['mos'] for i in self.final_objects]
      if self.params.advanced.integrate_with == 'cctbx':
        self.h = [i.final['sph'] for i in self.final_objects]
        self.s = [i.final['sih'] for i in self.final_objects]
        self.a = [i.final['spa'] for i in self.final_objects]
      else:
        self.s = [0]
        self.h = [0]
        self.a = [0]

      self.analysis_result.__setattr__('all_objects', len(self.all_objects))
      self.analysis_result.__setattr__('diff_objects', len(self.diff_objects))
      self.analysis_result.__setattr__('filter_fail_objects', len(self.filter_fail_objects))
      self.analysis_result.__setattr__('final_objects', len(self.final_objects))
      self.analysis_result.__setattr__('no_diff_objects', len(self.no_diff_objects))
      self.analysis_result.__setattr__('not_int_objects', len(self.not_int_objects))
      if self.params.advanced.integrate_with == 'dials':
        self.analysis_result.__setattr__('not_spf_objects', len(self.not_spf_objects))
        self.analysis_result.__setattr__('not_idx_objects', len(self.not_idx_objects))
      self.analysis_result.__setattr__('lres', self.lres)
      self.analysis_result.__setattr__('hres', self.hres)
      self.analysis_result.__setattr__('mos', self.mos)
      self.analysis_result.__setattr__('h', self.h)
      self.analysis_result.__setattr__('s', self.s)
      self.analysis_result.__setattr__('a', self.a)

  def print_results(self, final_table=None):
    """ Prints diagnostics from the final integration run. """

    cons_s = Counter(self.s).most_common(1)[0][0]
    cons_h = Counter(self.h).most_common(1)[0][0]
    cons_a = Counter(self.a).most_common(1)[0][0]

    if final_table is None:
      final_table = []
      final_table.append("\n\n{:-^80}\n".format('ANALYSIS OF RESULTS'))

      # In case no images were integrated
      if self.final_objects is None:
        final_table.append('NO IMAGES INTEGRATED!')
      else:
        if self.params.advanced.integrate_with == 'cctbx':
          final_table.append("Avg. signal height:    {:<8.3f}  std. dev:   "
                             "{:<6.2f}  max: {:<3}  min: {:<3}  consensus: {:<3}"
                             "".format(np.mean(self.s), np.std(self.s),
                                       max(self.s), min(self.s), cons_s))
          final_table.append("Avg. spot height:      {:<8.3f}  std. dev:   "
                             "{:<6.2f}  max: {:<3}  min: {:<3}  consensus: {:<3}"
                             "".format(np.mean(self.h), np.std(self.h),
                                       max(self.h), min(self.h), cons_h))
          final_table.append("Avg. spot areas:       {:<8.3f}  std. dev:   "
                             "{:<6.2f}  max: {:<3}  min: {:<3}  consensus: {:<3}"
                            "".format(np.mean(self.a), np.std(self.a),
                                      max(self.a), min(self.a), cons_a))
        final_table.append("Avg. resolution:       {:<8.3f}  std. dev:   "
                           "{:<6.2f}  lowest: {:<6.3f}  highest: {:<6.3f}"
                          "".format(np.mean(self.hres), np.std(self.hres),
                                    max(self.hres), min(self.hres)))
        final_table.append("Avg. number of spots:  {:<8.3f}  std. dev:   {:<6.2f}"
                          "".format(np.mean(self.spots), np.std(self.spots)))
        final_table.append("Avg. mosaicity:        {:<8.3f}  std. dev:   {:<6.2f}"
                          "".format(np.mean(self.mos), np.std(self.mos)))

        # If more than one integrated image, plot various summary graphs
        if len(self.final_objects) > 1:
          plot = Plotter(self.params, self.final_objects, self.viz_dir)
          if self.params.analysis.summary_graphs:
            if ( self.params.advanced.integrate_with == 'cctbx' and
                   self.params.cctbx.grid_search.type != None
                ):
              plot.plot_spotfinding_heatmap(write_files=True)
            plot.plot_res_histogram(write_files=True)
            med_beamX, med_beamY, pixel_size = plot.plot_beam_xy(write_files=True,
                                                               return_values=True)
          else:
            beamXY_info = plot.calculate_beam_xy()
            beamX, beamY = beamXY_info[:2]
            med_beamX = np.median(beamX)
            med_beamY = np.median(beamY)
            pixel_size = beamXY_info[-1]

          final_table.append("Median Beam Center:    X = {:<4.2f}, Y = {:<4.2f}"
                             "".format(med_beamX, med_beamY))
          self.analysis_result.__setattr__('beamX_mm', med_beamX)
          self.analysis_result.__setattr__('beamY_mm', med_beamY)
          self.analysis_result.__setattr__('pixel_size', pixel_size)
      self.analysis_result.__setattr__('final_table', final_table)


    for item in final_table:
        misc.main_log(self.logfile, item, (not self.gui_mode))

  def unit_cell_analysis(self):
    """ Calls unit cell analysis module, which uses hierarchical clustering
        (Zeldin, et al, Acta D, 2015) to split integration results according to
        detected morphological groupings (if any). Most useful with preliminary
        integration without target unit cell specified. """

    # Will not run clustering if only one integration result found or if turned off
    if self.final_objects is None:
      self.cons_uc = None
      self.cons_pg = None
      misc.main_log(self.logfile,
                    "\n\n{:-^80}\n".format(' UNIT CELL ANALYSIS '), True)
      misc.main_log(self.logfile,
                    '\n UNIT CELL CANNOT BE DETERMINED!', True)

    elif len(self.final_objects) == 1:
      unit_cell = (self.final_objects[0].final['a'],
                   self.final_objects[0].final['b'],
                   self.final_objects[0].final['c'],
                   self.final_objects[0].final['alpha'],
                   self.final_objects[0].final['beta'],
                   self.final_objects[0].final['gamma'])
      point_group = self.final_objects[0].final['sg']
      misc.main_log(self.logfile,
                    "\n\n{:-^80}\n".format(' UNIT CELL ANALYSIS '), True)
      uc_line = "{:<6} {:^4}:  {:<6.2f}, {:<6.2f}, {:<6.2f}, {:<6.2f}, "\
                "{:<6.2f}, {:<6.2f}".format('(1)', point_group,
                      unit_cell[0], unit_cell[1], unit_cell[2],
                      unit_cell[3], unit_cell[4], unit_cell[5])
      misc.main_log(self.logfile, uc_line, True)

      self.cons_pg = point_group
      self.cons_uc = unit_cell

    else:
      uc_table = []
      uc_summary = []

      if self.params.analysis.run_clustering:
        # run hierarchical clustering analysis
        from xfel.clustering.cluster import Cluster
        counter = 0

        threshold = self.params.analysis.cluster_threshold
        cluster_limit = self.params.analysis.cluster_limit
        if self.params.analysis.cluster_n_images > 0:
          n_images = self.params.analysis.cluster_n_images
        else:
          n_images = len(self.final_objects)

        obj_list = []
        if n_images < len(self.final_objects):
          import random
          for i in range(n_images):
            random_number = random.randrange(0, len(self.final_objects))
            if self.final_objects[random_number] in obj_list:
              while self.final_objects[random_number] in obj_list:
                random_number = random.randrange(0, len(self.final_objects))
              obj_list.append(self.final_objects[random_number])
            else:
              obj_list.append(self.final_objects[random_number])
        if obj_list == []:
          obj_list = self.final_objects

        # Cluster from iterable (a teensy bit less precise, but somewhat
        # faster; will see if it's worth the loss of precision); don't know
        # why it's less precise!
        with Capturing() as suppressed_output:
          uc_iterable = []
          for obj in obj_list:
            unit_cell = (float(obj.final['a']),
                         float(obj.final['b']),
                         float(obj.final['c']),
                         float(obj.final['alpha']),
                         float(obj.final['beta']),
                         float(obj.final['gamma']),
                         obj.final['sg'])
            uc_iterable.append(unit_cell)
          ucs = Cluster.from_iterable(iterable=uc_iterable)

        # Cluster from files (slower, but a teensy bit more precise; will
        # keep here for now)
        # ucs = Cluster.from_files(pickle_list=self.pickles)

        # Do clustering
        clusters, _ = ucs.ab_cluster(threshold=threshold,
                                     log=False, write_file_lists=False,
                                     schnell=False, doplot=False)
        uc_table.append("\n\n{:-^80}\n"\
                        "".format(' UNIT CELL ANALYSIS '))

        # extract clustering info and add to summary output list
        if cluster_limit is None:
          if len(self.pickles) / 10 >= 10:
            cluster_limit = 10
          else:
            cluster_limit = len(self.pickles) / 10

        for cluster in clusters:
          sorted_pg_comp = sorted(cluster.pg_composition.items(),
                                    key=lambda x: -1 * x[1])
          pg_nums = [pg[1] for pg in sorted_pg_comp]
          cons_pg = sorted_pg_comp[np.argmax(pg_nums)]

          if len(cluster.members) > cluster_limit:
            counter += 1

            # Sort clustered images by mosaicity, lowest to highest
            cluster_filenames = [j.path for j in cluster.members]
            clustered_objects = [i for i in self.final_objects if \
                                 i.final['final'] in cluster_filenames]
            sorted_cluster = sorted(clustered_objects,
                                    key=lambda i: i.final['mos'])
            # Write to file
            if self.params.analysis.cluster_write_files:
              output_file = os.path.join(self.output_dir, "uc_cluster_{}.lst".format(counter))
              for obj in sorted_cluster:
                with open(output_file, 'a') as scf:
                  scf.write('{}\n'.format(obj.final['final']))

              mark_output = os.path.basename(output_file)
            else:
              mark_output = '*'
              output_file = None

            # Populate clustering info for GUI display
            uc_no_stdev = "{:<6.2f} {:<6.2f} {:<6.2f} " \
                          "{:<6.2f} {:<6.2f} {:<6.2f} " \
                          "".format(cluster.medians[0], cluster.medians[1],
                                    cluster.medians[2], cluster.medians[3],
                                    cluster.medians[4], cluster.medians[5])
            cluster_info = {'number':len(cluster.members),
                            'pg': cons_pg[0],
                            'uc':uc_no_stdev,
                            'filename':mark_output}
            self.clusters.append(cluster_info)

          else:
            mark_output = ''
            output_file = None

          # format and record output
          uc_line = "{:<6} {:^4}:  {:<6.2f} ({:>5.2f}), {:<6.2f} ({:>5.2f}), "\
                    "{:<6.2f} ({:>5.2f}), {:<6.2f} ({:>5.2f}), "\
                    "{:<6.2f} ({:>5.2f}), {:<6.2f} ({:>5.2f})   "\
                    "{}".format('({})'.format(len(cluster.members)), cons_pg[0],
                                          cluster.medians[0], cluster.stdevs[0],
                                          cluster.medians[1], cluster.stdevs[1],
                                          cluster.medians[2], cluster.stdevs[2],
                                          cluster.medians[3], cluster.stdevs[3],
                                          cluster.medians[4], cluster.stdevs[4],
                                          cluster.medians[5], cluster.stdevs[5],
                                          mark_output)
          uc_table.append(uc_line)
          lattices = ', '.join(['{} ({})'.format(i[0], i[1]) for
                                i in sorted_pg_comp])
          uc_info = [len(cluster.members), cons_pg[0], cluster.medians,
                     output_file, uc_line, lattices]
          uc_summary.append(uc_info)

      else:

        # generate average unit cell
        uc_table.append("\n\n{:-^80}\n" \
                        "".format(' UNIT CELL AVERAGING (no clustering) '))
        uc_a = [i.final['a'] for i in self.final_objects]
        uc_b = [i.final['b'] for i in self.final_objects]
        uc_c = [i.final['c'] for i in self.final_objects]
        uc_alpha = [i.final['alpha'] for i in self.final_objects]
        uc_beta = [i.final['beta'] for i in self.final_objects]
        uc_gamma = [i.final['gamma'] for i in self.final_objects]
        uc_sg = [i.final['sg'] for i in self.final_objects]
        cons_pg = Counter(uc_sg).most_common(1)[0][0]
        all_pgs = Counter(uc_sg).most_common()
        uc_line = "{:<6} {:^4}:  {:<6.2f} ({:>5.2f}), {:<6.2f} ({:>5.2f}), " \
                  "{:<6.2f} ({:>5.2f}), {:<6.2f} ({:>5.2f}), " \
                  "{:<6.2f} ({:>5.2f}), {:<6.2f} ({:>5.2f})   " \
                  "{}".format('({})'.format(len(self.final_objects)), cons_pg,
                              np.median(uc_a), np.std(uc_a),
                              np.median(uc_b), np.std(uc_b),
                              np.median(uc_c), np.std(uc_c),
                              np.median(uc_alpha), np.std(uc_alpha),
                              np.median(uc_beta), np.std(uc_beta),
                              np.median(uc_gamma), np.std(uc_gamma), '')
        unit_cell = (np.median(uc_a), np.median(uc_b), np.median(uc_c),
                     np.median(uc_alpha), np.median(uc_beta), np.median(uc_gamma))
        uc_table.append(uc_line)
        lattices = ', '.join(['{} ({})'.format(i[0], i[1]) for i in all_pgs])
        uc_info = [len(self.final_objects), cons_pg, unit_cell, None,
                   uc_line, lattices]
        uc_summary.append(uc_info)

      uc_table.append('\nMost common unit cell:\n')

      # select the most prevalent unit cell (most members in cluster)
      uc_freqs = [i[0] for i in uc_summary]
      uc_pick = uc_summary[np.argmax(uc_freqs)]
      uc_table.append(uc_pick[4])
      uc_table.append('\nBravais Lattices in Biggest Cluster: {}'
                      ''.format(uc_pick[5]))

      self.cons_pg = uc_pick[1]
      self.cons_uc = uc_pick[2]

      if uc_pick[3] != None:
        self.prime_data_path = uc_pick[3]

      for item in uc_table:
        misc.main_log(self.logfile, item, (not self.gui_mode))

      self.analysis_result.__setattr__('clusters', self.clusters)
      self.analysis_result.__setattr__('cons_pg', self.cons_pg)
      self.analysis_result.__setattr__('cons_uc', self.cons_uc)

      if self.gui_mode:
        return self.cons_pg, self.cons_uc, self.clusters


  def print_summary(self, write_files=True):
    """ Prints summary and appends to general log file. Also outputs some of it
        on stdout. Also writes out output list files.
    """

    summary = []

    misc.main_log(self.logfile,
                  "\n\n{:-^80}\n".format('SUMMARY'),
                  (not self.gui_mode))

    summary.append('raw images read in:                  {}'\
                   ''.format(len(self.all_objects)))
    summary.append('raw images with no diffraction:      {}'\
                   ''.format(len(self.no_diff_objects)))
    summary.append('raw images with diffraction:         {}'\
                   ''.format(len(self.diff_objects)))
    if self.params.advanced.integrate_with == 'cctbx':
      summary.append('failed indexing / integration:       {}'\
                     ''.format(len(self.not_int_objects)))
      summary.append('failed prefilter:                    {}'\
                     ''.format(len(self.filter_fail_objects)))
    elif self.params.advanced.integrate_with == 'dials':
      summary.append('failed spotfinding:                  {}'\
                     ''.format(len(self.not_spf_objects)))
      summary.append('failed indexing:                     {}'\
                     ''.format(len(self.not_idx_objects)))
      summary.append('failed filter:'
                     ''.format(len(self.filter_fail_objects)))
      summary.append('failed integration:                  {}'\
                     ''.format(len(self.not_int_objects)))
    summary.append('final integrated pickles:            {}'\
                   ''.format(len(self.sorted_final_images)))

    for item in summary:
      misc.main_log(self.logfile, "{}".format(item), (not self.gui_mode))

    misc.main_log(self.logfile, '\n\nIOTA version {0}'.format(self.ver))
    misc.main_log(self.logfile, "{}\n".format(self.now))


    # Write list files:
    if write_files:

      input_list_file = os.path.join(self.output_dir, 'input_images.lst')
      blank_images_file = os.path.join(self.output_dir, 'blank_images.lst')
      prefilter_fail_file = os.path.join(self.output_dir, 'failed_cctbx_prefilter.lst')
      spotfinding_fail_file = os.path.join(self.output_dir, 'failed_dials_spotfinding.lst')
      indexing_fail_file = os.path.join(self.output_dir, 'failed_dials_indexing.lst')
      not_integrated_file = os.path.join(self.output_dir, 'not_integrated.lst')
      integrated_file = os.path.join(self.output_dir, 'integrated.lst')
      int_images_file = os.path.join(self.output_dir, 'int_image_pickles.lst')
      analyzer_file = os.path.join(self.output_dir, 'analysis.pickle')

      if self.prime_data_path == None:
        self.prime_data_path = integrated_file

      if len(self.no_diff_objects) > 0:
        with open(blank_images_file, 'w') as bif:
          for obj in self.no_diff_objects:
            bif.write('{}\n'.format(obj.conv_img))

      if len(self.diff_objects) > 0:
        with open(input_list_file, 'w') as ilf:
          for obj in self.diff_objects:
            ilf.write('{}\n'.format(obj.conv_img))

      if len(self.not_int_objects) > 0:
        with open(not_integrated_file, 'w') as nif:
          for obj in self.not_int_objects:
            nif.write('{}\n'.format(obj.conv_img))

      if len(self.filter_fail_objects) > 0:
        with open(prefilter_fail_file, 'w') as pff:
          for obj in self.filter_fail_objects:
            pff.write('{}\n'.format(obj.conv_img))

      if self.params.advanced.integrate_with == 'dials':
        if len(self.not_spf_objects) > 0:
          with open(spotfinding_fail_file, 'w') as sff:
            for obj in self.not_spf_objects:
              sff.write('{}\n'.format(obj.conv_img))
        if len(self.not_idx_objects) > 0:
          with open(indexing_fail_file, 'w') as iff:
            for obj in self.not_idx_objects:
              iff.write('{}\n'.format(obj.conv_img))

      if len(self.final_objects) > 0:
        with open(integrated_file, 'w') as intf:
          for obj in self.sorted_final_images:
            intf.write('{}\n'.format(obj.final['final']))
        with open(int_images_file, 'w') as ipf:
          for obj in self.sorted_final_images:
            ipf.write('{}\n'.format(obj.final['img']))

      # Dump the Analyzer object into a pickle file for later fast recovery
      ep.dump(analyzer_file, self.analysis_result)




  def make_prime_input(self, filename='prime.phil'):
    """ Imports default PRIME input parameters, modifies correct entries and
        prints out a starting PHIL file to be used with PRIME
    """

    if self.params.advanced.integrate_with == 'cctbx':
      img_pickle = self.final_objects[0].final['img']
      pixel_size = pickle.load(open(img_pickle, "rb"))['PIXEL_SIZE']
    elif self.params.advanced.integrate_with == 'dials':
      proc_pickle = self.final_objects[0].final['final']
      pixel_size = pickle.load(open(proc_pickle, 'rb'))['pixel_size']

    triclinic = ['P1']
    monoclinic = ['C2', 'P2']
    orthorhombic = ['P222', 'C222', 'I222', 'F222']
    tetragonal = ['I4', 'I422', 'P4', 'P422']
    hexagonal = ['P3', 'P312', 'P321', 'P6', 'P622']
    rhombohedral = ['R3', 'R32']
    cubic = ['F23', 'F432', 'I23', 'I432', 'P23', 'P432']

    sg = self.cons_pg.replace(" ", "")
    uc = ['{:4.2f}'.format(i) for i in self.cons_uc]

    if sg in triclinic:
      crystal_system = 'Triclinic'
    elif sg in monoclinic:
      crystal_system = 'Monoclinic'
    elif sg in orthorhombic:
      crystal_system = 'Orthorhombic'
    elif sg in tetragonal:
      crystal_system = 'Tetragonal'
    elif sg in hexagonal:
      crystal_system = 'Hexagonal'
    elif sg in rhombohedral:
      crystal_system = 'Rhombohedral'
    elif sg in cubic:
      crystal_system = 'Cubic'
    else:
      crystal_system = 'None'

    # Determine number of images for indexing ambiguity resolution
    # My default: 1/2 of images or 300, whichever is smaller
    if len(self.final_objects) >= 600:
      idx_ambiguity_sample = 300
      idx_ambiguity_selected = 100
    else:
      idx_ambiguity_sample = int(len(self.final_objects) / 2)
      idx_ambiguity_selected = int(idx_ambiguity_sample / 3)

    prime_params = mod_input.master_phil.extract()

    prime_params.data = [self.prime_data_path]
    prime_params.run_no = os.path.join(os.path.dirname(self.prime_data_path),
                                       'prime/001')
    prime_params.title = 'Auto-generated by IOTA v{} on {}'.format(self.ver, self.now)
    prime_params.scale.d_min = np.mean(self.hres)
    prime_params.scale.d_max = 8
    prime_params.postref.scale.d_min = np.mean(self.hres)
    prime_params.postref.scale.d_max = np.max(self.lres)
    prime_params.postref.crystal_orientation.d_min = np.mean(self.hres)
    prime_params.postref.crystal_orientation.d_max = np.max(self.lres)
    prime_params.postref.reflecting_range.d_min = np.mean(self.hres)
    prime_params.postref.reflecting_range.d_max = np.max(self.lres)
    prime_params.postref.unit_cell.d_min = np.mean(self.hres)
    prime_params.postref.unit_cell.d_max = np.max(self.lres)
    prime_params.postref.allparams.d_min = np.mean(self.hres)
    prime_params.postref.allparams.d_max = np.max(self.lres)
    prime_params.merge.d_min = np.mean(self.hres)
    prime_params.merge.d_max = np.max(self.lres)
    prime_params.target_unit_cell = unit_cell(self.cons_uc)
    prime_params.target_space_group = sg
    prime_params.target_crystal_system = crystal_system
    prime_params.pixel_size_mm = pixel_size
    prime_params.n_residues = 500
    prime_params.indexing_ambiguity.n_sample_frames = idx_ambiguity_sample
    prime_params.indexing_ambiguity.n_selected_frames = idx_ambiguity_selected

    prime_phil = mod_input.master_phil.format(python_object=prime_params)

    with Capturing() as output:
      prime_phil.show()

    txt_out = ''
    for one_output in output:
      txt_out += one_output + '\n'

    prime_file = os.path.join(self.output_dir, filename)
    with open(prime_file, 'w') as pf:
      pf.write(txt_out)

    return prime_phil
