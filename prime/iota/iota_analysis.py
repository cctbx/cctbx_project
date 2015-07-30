from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 04/07/2015
Last Changed: 07/29/2015
Description : Analyzes integration results and outputs them in an accessible
              format. Includes unit cell analysis by hierarchical clustering
              (Zeldin, et al., Acta Cryst D, 2013). In case of multiple clusters
              outputs a file with list of integrated pickles that comprise each
              cluster. Populates a PHIL file for PRIME with information from
              integration results (e.g. unit cell, resolution, data path, etc.)
'''
import os
import numpy as np
from collections import Counter

import cPickle as pickle
from cctbx.uctbx import unit_cell
from xfel.clustering.cluster import Cluster

import prime.iota.iota_misc as misc
from prime.iota.iota_misc import Capturing
from prime.postrefine import mod_input




class Analyzer(object):
  """ Class to analyze integration results """

  def __init__(self,
               all_objects,
               logfile,
               version,
               now):

    self.ver = version
    self.now = now

    self.prime_data_path = None

    self.all_objects = all_objects
    self.final_objects = [i for i in all_objects if i.final['final'] != None]
    self.logfile = logfile

    self.cons_pg = None
    self.cons_uc = None

    self.pickles = [i.final['final'] for i in self.final_objects]
    self.h = [i.final['sph'] for i in self.final_objects]
    self.s = [i.final['sih'] for i in self.final_objects]
    self.a = [i.final['spa'] for i in self.final_objects]
    self.res = [i.final['res'] for i in self.final_objects]
    self.spots = [i.final['strong'] for i in self.final_objects]
    self.mos = [i.final['mos'] for i in self.final_objects]


  def print_results(self):
    """ Prints diagnostics from the final integration run. """

    cons_s = Counter(self.s).most_common(1)[0][0]
    cons_h = Counter(self.h).most_common(1)[0][0]
    cons_a = Counter(self.a).most_common(1)[0][0]

    final_table = []
    final_table.append("\n\n{:-^80}\n".format('ANALYSIS OF RESULTS'))
    final_table.append("Total images:          {}".format(len(self.final_objects)))
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

    for item in final_table:
        misc.main_log(self.logfile, item, True)


  def unit_cell_analysis(self, cluster_threshold, output_dir):
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
      clusters, _ = ucs.ab_cluster(cluster_threshold, log=False,
                                   write_file_lists=False, schnell=False,
                                   doplot=False)

      uc_table.append("\n\n{:-^80}\n"\
                      "".format(' UNIT CELL ANALYSIS '))

      # extract clustering info and add to summary output list
      for cluster in clusters:
        sorted_pg_comp = sorted(cluster.pg_composition.items(),
                                  key=lambda x: -1 * x[1])
        pg_nums = [pg[1] for pg in sorted_pg_comp]
        cons_pg = sorted_pg_comp[np.argmax(pg_nums)]

        output_file = os.path.join(output_dir, "uc_cluster_{}.lst".format(counter))

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
          for obj in sorted_cluster:
            with open(output_file, 'a') as scf:
              scf.write('{}\n'.format(obj.final['img']))

          mark_output = os.path.basename(output_file)
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


  def print_summary(self, int_base):
    """ Prints summary and appends to general log file. Also outputs some of it
        on stdout. Also writes out output list files.
    """

    summary = []
    input_list_file = os.path.join(int_base, 'input_images.lst')
    blank_images_file = os.path.join(int_base, 'blank_images.lst')
    prefilter_fail_file = os.path.join(int_base, 'failed_prefilter.lst')
    not_integrated_file = os.path.join(int_base, 'not_integrated.lst')
    integrated_file = os.path.join(int_base, 'integrated.lst')
    int_images_file = os.path.join(int_base, 'int_image_pickles.lst')

    if self.prime_data_path == None:
      self.prime_data_path = integrated_file


    misc.main_log(self.logfile, "\n\n{:-^80}\n".format('SUMMARY'), True)

    summary.append('raw images read in:                  {}'\
                   ''.format(len(self.all_objects)))

    no_diff_objects = [i for i in self.all_objects if i.triage == 'rejected']
    summary.append('raw images with no diffraction:      {}'\
                   ''.format(len(no_diff_objects)))
    if len(no_diff_objects) > 0:
      with open(blank_images_file, 'w') as bif:
        for obj in no_diff_objects:
          bif.write('{}\n'.format(obj.conv_img))

    diff_objects = [i for i in self.all_objects if i.triage == 'accepted']
    summary.append('raw images with diffraction:         {}'\
                   ''.format(len(diff_objects)))
    if len(diff_objects) > 0:
      with open(input_list_file, 'w') as ilf:
        for obj in diff_objects:
          ilf.write('{}\n'.format(obj.conv_img))

    not_int_objects = [i for i in diff_objects if len(i.grid) == 0]
    summary.append('raw images not integrated:           {}'\
                   ''.format(len(not_int_objects)))
    if len(not_int_objects) > 0:
      with open(not_integrated_file, 'w') as nif:
        for obj in not_int_objects:
          nif.write('{}\n'.format(obj.conv_img))

    int_objects = [i for i in diff_objects if len(i.grid) != 0]
    prefilter_fail_objects = [i for i in int_objects if not i.prefilter]
    summary.append('images failed prefilter:             {}'\
                   ''.format(len(prefilter_fail_objects)))
    if len(prefilter_fail_objects) > 0:
      with open(prefilter_fail_file, 'w') as pff:
        for obj in prefilter_fail_objects:
          pff.write('{}\n'.format(obj.conv_img))

    final_images = sorted(self.final_objects, key=lambda i: i.final['mos'])
    summary.append('final integrated pickles:            {}'\
                   ''.format(len(final_images)))
    if len(self.final_objects) > 0:
      with open(integrated_file, 'w') as intf:
        for obj in final_images:
          intf.write('{}\n'.format(obj.final['final']))
      with open(int_images_file, 'w') as ipf:
        for obj in final_images:
          ipf.write('{}\n'.format(obj.final['img']))

    for item in summary:
      misc.main_log(self.logfile, "{}".format(item), True)

    misc.main_log(self.logfile, '\n\nIOTA version {0}'.format(self.ver))
    misc.main_log(self.logfile, "{}\n".format(self.now))


  def make_prime_input(self, int_folder):
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

    prime_file = os.path.join(int_folder, 'prime.phil')
    with open(prime_file, 'w') as pf:
      pf.write(txt_out)
