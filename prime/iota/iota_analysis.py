from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 04/07/2015
Last Changed: 02/24/2016
Description : Analyzes integration results and outputs them in an accessible
              format. Includes unit cell analysis by hierarchical clustering
              (Zeldin, et al., Acta Cryst D, 2013). In case of multiple clusters
              outputs a file with list of integrated pickles that comprise each
              cluster. Populates a PHIL file for PRIME with information from
              integration results (e.g. unit cell, resolution, data path, etc.)
'''

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import os
import numpy as np
from collections import Counter
import math

import cPickle as pickle
from cctbx.uctbx import unit_cell
from xfel.clustering.cluster import Cluster

import prime.iota.iota_misc as misc
from prime.iota.iota_misc import Capturing
from prime.postrefine import mod_input

class Analyzer(object):
  """ Class to analyze integration results """

  def __init__(self,
               init,
               all_objects,
               version):

    self.ver = version
    self.now = init.now
    self.params = init.params
    self.args = init.args
    self.output_dir = init.int_base

    self.prime_data_path = None
    self.all_objects = all_objects
    self.final_objects = [i for i in all_objects if i.fail == None]
    self.logfile = init.logfile

    self.cons_pg = None
    self.cons_uc = None

    self.pickles = [i.final['final'] for i in self.final_objects]
    self.res = [i.final['res'] for i in self.final_objects]
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

  def grid_search_heatmap(self):

    ch = max(self.h) - min(self.h) + 1
    ca = max(self.a) - min(self.a) + 1
    ints = [(i.final['sph'], i.final['spa']) for i in self.final_objects]
    ic = Counter(ints)

    hm_data = np.zeros((ch, ca))
    for i in ic.items():
      hm_data[i[0][0]-min(self.h), i[0][1]-min(self.a)] = i[1]

    rows = range(min(self.h), max(self.h) + 1)
    cols = range(min(self.a), max(self.a) + 1)
    row_labels = [str(i) for i in rows]
    col_labels = [str(j) for j in cols]

    fig, ax = plt.subplots()
    fig.canvas.draw()
    heatmap = plt.pcolor(hm_data, cmap='Reds')
    plt.colorbar()
    ax.set_yticks(np.arange(len(rows))+.5, minor=False)
    ax.set_xticks(np.arange(len(cols))+.5, minor=False)
    ax.set_yticklabels(row_labels, minor=False)
    ax.set_xticklabels(col_labels, minor=False)
    ax.set_xlabel('Spot area')
    ax.set_ylabel('Spot height')

    if self.args.analyze == None:
      fig.savefig(os.path.join(self.output_dir, 'visualization/heatmap.pdf'),
                               format='pdf', bbox_inches=0)

  def plot_beam_centers(self):
    ''' Plot beam center coordinates and a histogram of distances from the median
        of beam center coordinates to each set of coordinates. Superpose a
        predicted mis-indexing shift by L +/- 1 (calculated for each axis).
    '''
    from libtbx import easy_pickle as ep

    info = []
    wavelengths = []
    distances = []
    cells = []

    # Import relevant info
    for i in self.pickles:
      beam = ep.load(i)
      info.append([i, beam['xbeam'], beam['ybeam']])

    # Calculate beam center coordinates and distances
    beamX = [i[1] for i in info]
    beamY = [j[2] for j in info]
    beam_dist = [math.hypot(i[1] - np.median(beamX), i[2] - np.median(beamY)) for i in info]
    beam_dist_std = np.std(beam_dist)
    img_list = [[i[0], i[1], i[2], j] for i,j in zip(info, beam_dist)]

    # Separate out outliers
    outliers = [i for i in img_list if i[3] > 2 * beam_dist_std]
    clean = [i for i in img_list if i[3] <= 2 * beam_dist_std]
    cbeamX = [i[1] for i in clean]
    cbeamY = [j[2] for j in clean]
    obeamX = [i[1] for i in outliers]
    obeamY = [j[2] for j in outliers]

    # Calculate median wavelength, detector distance and unit cell params from
    # non-outliers only
    for i in clean:
      beam = ep.load(i[0])
      wavelengths.append(beam['wavelength'])
      distances.append(beam['distance'])
      cells.append(beam['observations'][0].unit_cell().parameters())

    wavelength = np.median(wavelengths)
    det_distance = np.median(distances)
    a = np.median([i[0] for i in cells])
    b = np.median([i[1] for i in cells])
    c = np.median([i[2] for i in cells])

    # Calculate predicted L +/- 1 misindexing distance for each cell edge
    aD = det_distance * math.tan(2 * math.asin(wavelength / (2 * a)))
    bD = det_distance * math.tan(2 * math.asin(wavelength / (2 * b)))
    cD = det_distance * math.tan(2 * math.asin(wavelength / (2 * c)))

    # Plot figure
    fig = plt.figure(figsize=(9,13))
    gsp = gridspec.GridSpec(2, 1, height_ratios=[3,1])

    # Calculate axis limits of beam center scatter plot
    ax1_delta = np.ceil(np.max(beam_dist))
    xmax = round(np.median(beamX) + ax1_delta)
    xmin = round(np.median(beamX) - ax1_delta)
    ymax = round(np.median(beamY) + ax1_delta)
    ymin = round(np.median(beamY) - ax1_delta)

    # Plot beam center scatter plot
    ax1 = fig.add_subplot(gsp[0,:], aspect='equal')
    ax1.axis([xmin, xmax, ymin, ymax])
    ax1.scatter(cbeamX, cbeamY, alpha=1, s=20, c='grey', lw=1)
    ax1.scatter(obeamX, obeamY, alpha=1, s=20, c='red', lw=1)
    ax1.plot(np.median(beamX), np.median(beamY), markersize=8, marker='o', c='yellow', lw=2)

    # Plot projected mis-indexing limits for all three axes
    circle_a = plt.Circle((np.median(beamX), np.median(beamY)), radius=aD, color='r', fill=False, clip_on=True)
    circle_b = plt.Circle((np.median(beamX), np.median(beamY)), radius=bD, color='g', fill=False, clip_on=True)
    circle_c = plt.Circle((np.median(beamX), np.median(beamY)), radius=cD, color='b', fill=False, clip_on=True)
    ax1.add_patch(circle_a)
    ax1.add_patch(circle_b)
    ax1.add_patch(circle_c)
    #ax1.annotate('a-shift ={:2.1f} mm'.format(aD), xy=(np.median(beamX)+aD, np.median(beamY)), ha='right', size=15)
    #ax1.annotate('b-shift ={:2.1f} mm'.format(bD), xy=(np.median(beamX)+bD, np.median(beamY)), ha='right', size=15)
    #ax1.annotate('c-shift ={:2.1f} mm'.format(cD), xy=(np.median(beamX)+cD, np.median(beamY)), ha='right', size=15)
    ax1.set_xlabel('BeamX (mm)', fontsize=15)
    ax1.set_ylabel('BeamY (mm)', fontsize=15)
    ax1.set_title('Beam Center Coordinates')

    # Plot histogram of distances to each beam center from median
    ax2_delta = np.max([abs(np.mean(beam_dist) - np.max(beam_dist)),
                       abs(np.mean(beam_dist) - np.min(beam_dist))])
    ax2 = fig.add_subplot(gsp[1, :])
    ax2_n, ax2_bins, ax2_patches = plt.hist(beam_dist, 20, normed=False, facecolor='b', alpha=0.75, histtype='stepfilled')
    ax2_height = (np.max(ax2_n) + 9) // 10 * 10
    ax2.axis([0, np.max(beam_dist), 0, ax2_height])
    ax2.set_xlabel('Distance from median (mm)', fontsize=15)
    ax2.set_ylabel('No. of images', fontsize=15)

 #   plt.tight_layout()

    if self.args.analyze == None:
      bc_file = os.path.join(self.output_dir, 'visualization/beam_center.pdf')
      fig.savefig(bc_file, format='pdf', bbox_inches=0)


    return np.median(beamX), np.median(beamY)

  def print_results(self):
    """ Prints diagnostics from the final integration run. """

    cons_s = Counter(self.s).most_common(1)[0][0]
    cons_h = Counter(self.h).most_common(1)[0][0]
    cons_a = Counter(self.a).most_common(1)[0][0]

    final_table = []
    final_table.append("\n\n{:-^80}\n".format('ANALYSIS OF RESULTS'))
    final_table.append("Total images:          {}".format(len(self.final_objects)))

    if self.params.advanced.integrate_with == 'cctbx':
      final_table.append("Avg. signal height:    {:<8.3f}  std. dev:    {:<6.2f}"\
                         "  max: {:<3}  min: {:<3}  consensus: {:<3}"\
                         "".format(np.mean(self.s), np.std(self.s),
                                   max(self.s), min(self.s), cons_s))
      final_table.append("Avg. spot height:      {:<8.3f}  std. dev:    {:<6.2f}"\
                         "  max: {:<3}  min: {:<3}  consensus: {:<3}"\
                         "".format(np.mean(self.h), np.std(self.h),
                                   max(self.h), min(self.h), cons_h))
      final_table.append("Avg. spot areas:       {:<8.3f}  std. dev:    {:<6.2f}"\
                        "  max: {:<3}  min: {:<3}  consensus: {:<3}"\
                        "".format(np.mean(self.a), np.std(self.a),
                                  max(self.a), min(self.a), cons_a))
    final_table.append("Avg. resolution:       {:<8.3f}  std. dev:    {:<6.2f}"\
                       "  lowest: {:<6.3f}  highest: {:<6.3f}"\
                      "".format(np.mean(self.res), np.std(self.res),
                                max(self.res), min(self.res)))
    final_table.append("Avg. number of spots:  {:<8.3f}  std. dev:    {:<6.2f}"\
                      "".format(np.mean(self.spots), np.std(self.spots)))
    final_table.append("Avg. mosaicity:        {:<8.3f}  std. dev:    {:<6.2f}"\
                      "".format(np.mean(self.mos), np.std(self.mos)))

    if len(self.final_objects) > 1:
      if ( self.params.advanced.integrate_with == 'cctbx' and
             self.params.cctbx.grid_search.type != None
          ):
        self.grid_search_heatmap()
      med_beamX, med_beamY = self.plot_beam_centers()
      final_table.append("Median Beam Center:    X = {:<4.2f}, Y = {:<4.2f}"\
                         "".format(med_beamX, med_beamY))

    for item in final_table:
        misc.main_log(self.logfile, item, True)

  def unit_cell_analysis(self,
                         write_files=True):
    """ Calls unit cell analysis module, which uses hierarchical clustering
        (Zeldin, et al, Acta D, 2015) to split integration results according to
        detected morphological groupings (if any). Most useful with preliminary
        integration without target unit cell specified. """

    # Will not run clustering if only one integration result found
    if len(self.final_objects) == 1:
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
      counter = 1

      # run hierarchical clustering analysis
      ucs = Cluster.from_files(self.pickles, use_b=True)
      clusters, _ = ucs.ab_cluster(self.params.analysis.cluster_threshold,
                                   log=False, write_file_lists=False,
                                   schnell=False, doplot=False)
      uc_table.append("\n\n{:-^80}\n"\
                      "".format(' UNIT CELL ANALYSIS '))

      # extract clustering info and add to summary output list
      for cluster in clusters:
        sorted_pg_comp = sorted(cluster.pg_composition.items(),
                                  key=lambda x: -1 * x[1])
        pg_nums = [pg[1] for pg in sorted_pg_comp]
        cons_pg = sorted_pg_comp[np.argmax(pg_nums)]

        # write out lists of output pickles that comprise clusters with > 1 members
        if len(cluster.members) > 1:
          counter += 1

          # Sort clustered images by mosaicity, lowest to highest
          cluster_filenames = [j.path for j in cluster.members]
          clustered_objects = [i for i in self.final_objects if \
                               i.final['final'] in cluster_filenames]
          sorted_cluster = sorted(clustered_objects,
                                  key=lambda i: i.final['mos'])
          # Write to file
          if write_files:
            output_file = os.path.join(self.output_dir, "uc_cluster_{}.lst".format(counter))
            for obj in sorted_cluster:
              with open(output_file, 'a') as scf:
                scf.write('{}\n'.format(obj.final['final']))

            mark_output = os.path.basename(output_file)
          else:
            mark_output = '*'
            output_file = None
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
        uc_info = [len(cluster.members), cons_pg[0], cluster.medians,
                   output_file, uc_line]
        uc_summary.append(uc_info)

      uc_table.append('\nMost common unit cell:\n')

      # select the most prevalent unit cell (most members in cluster)
      uc_freqs = [i[0] for i in uc_summary]
      uc_pick = uc_summary[np.argmax(uc_freqs)]
      uc_table.append(uc_pick[4])

      self.cons_pg = uc_pick[1]
      self.cons_uc = uc_pick[2]

      if uc_pick[3] != None:
        self.prime_data_path = uc_pick[3]

      for item in uc_table:
          misc.main_log(self.logfile, item, True)


  def print_summary(self, write_files=True):
    """ Prints summary and appends to general log file. Also outputs some of it
        on stdout. Also writes out output list files.
    """

    summary = []

    misc.main_log(self.logfile, "\n\n{:-^80}\n".format('SUMMARY'), True)

    summary.append('raw images read in:                  {}'\
                   ''.format(len(self.all_objects)))

    no_diff_objects = [i for i in self.all_objects if i.fail == 'failed triage']
    summary.append('raw images with no diffraction:      {}'\
                   ''.format(len(no_diff_objects)))

    diff_objects = [i for i in self.all_objects if i.fail != 'failed triage']
    summary.append('raw images with diffraction:         {}'\
                   ''.format(len(diff_objects)))

    not_int_objects = [i for i in self.all_objects if i.fail == 'failed grid search']
    summary.append('raw images not integrated:           {}'\
                   ''.format(len(not_int_objects)))

    if self.params.advanced.integrate_with == 'cctbx':
      prefilter_fail_objects = [i for i in self.all_objects if i.fail == 'failed prefilter']
      summary.append('images failed prefilter:             {}'\
                     ''.format(len(prefilter_fail_objects)))

    final_images = sorted(self.final_objects, key=lambda i: i.final['mos'])
    summary.append('final integrated pickles:            {}'\
                   ''.format(len(final_images)))

    for item in summary:
      misc.main_log(self.logfile, "{}".format(item), True)
    misc.main_log(self.logfile, '\n\nIOTA version {0}'.format(self.ver))
    misc.main_log(self.logfile, "{}\n".format(self.now))


    # Write list files:
    if write_files:

      input_list_file = os.path.join(self.output_dir, 'input_images.lst')
      blank_images_file = os.path.join(self.output_dir, 'blank_images.lst')
      prefilter_fail_file = os.path.join(self.output_dir, 'failed_prefilter.lst')
      not_integrated_file = os.path.join(self.output_dir, 'not_integrated.lst')
      integrated_file = os.path.join(self.output_dir, 'integrated.lst')
      int_images_file = os.path.join(self.output_dir, 'int_image_pickles.lst')

      if self.prime_data_path == None:
        self.prime_data_path = integrated_file

      if len(no_diff_objects) > 0:
        with open(blank_images_file, 'w') as bif:
          for obj in no_diff_objects:
            bif.write('{}\n'.format(obj.conv_img))

      if len(diff_objects) > 0:
        with open(input_list_file, 'w') as ilf:
          for obj in diff_objects:
            ilf.write('{}\n'.format(obj.conv_img))

      if len(not_int_objects) > 0:
        with open(not_integrated_file, 'w') as nif:
          for obj in not_int_objects:
            nif.write('{}\n'.format(obj.conv_img))

      if self.params.advanced.integrate_with == 'cctbx' and len(prefilter_fail_objects) > 0:
        with open(prefilter_fail_file, 'w') as pff:
          for obj in prefilter_fail_objects:
            pff.write('{}\n'.format(obj.conv_img))

      if len(self.final_objects) > 0:
        with open(integrated_file, 'w') as intf:
          for obj in final_images:
            intf.write('{}\n'.format(obj.final['final']))
        with open(int_images_file, 'w') as ipf:
          for obj in final_images:
            ipf.write('{}\n'.format(obj.final['img']))


  def make_prime_input(self):
    """ Imports default PRIME input parameters, modifies correct entries and
        prints out a starting PHIL file to be used with PRIME
    """

    img_pickle = self.final_objects[0].final['img']
    pixel_size = pickle.load(open(img_pickle, "rb"))['PIXEL_SIZE']

    triclinic = ['P1']
    monoclinic = ['C2', 'P2']
    orthorhombic = ['P222', 'C222', 'I222', 'F222']
    tetragonal = ['I4', 'I422', 'P4', 'P422']
    hexagonal = ['P3', 'P312', 'P321', 'P6', 'P622']
    rhombohedral = ['R3', 'R32']
    cubic = ['F23', 'F432', 'I23', 'I432', 'P23', 'P432']

    sg = self.cons_pg.replace(" ", "")

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

    prime_params = mod_input.master_phil.extract()

    prime_params.data = [self.prime_data_path]
    prime_params.run_no = '001'
    prime_params.title = 'Auto-generated by IOTA v{} on {}'.format(self.ver, self.now)
    prime_params.scale.d_min = np.mean(self.res)
    prime_params.postref.scale.d_min = np.mean(self.res)
    prime_params.postref.crystal_orientation.d_min = np.mean(self.res)
    prime_params.postref.reflecting_range.d_min = np.mean(self.res)
    prime_params.postref.unit_cell.d_min = np.mean(self.res)
    prime_params.postref.allparams.d_min = np.mean(self.res)
    prime_params.merge.d_min = np.mean(self.res)
    prime_params.target_unit_cell = unit_cell(self.cons_uc)
    prime_params.target_space_group = sg
    prime_params.target_crystal_system = crystal_system
    prime_params.pixel_size_mm = pixel_size

    prime_phil = mod_input.master_phil.format(python_object=prime_params)

    with Capturing() as output:
      prime_phil.show()

    txt_out = ''
    for one_output in output:
      txt_out += one_output + '\n'

    prime_file = os.path.join(self.output_dir, 'prime.phil')
    with open(prime_file, 'w') as pf:
      pf.write(txt_out)
