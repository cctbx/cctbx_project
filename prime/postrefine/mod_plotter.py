from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 05/25/2016
Last Changed: 06/01/2016
Description : PRIME Result Plotter module
'''

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec

class Plotter(object):
  ''' Class with function to plot various PRIME charts (includes Table 1) '''

  def __init__(self, params, info, output_dir=None):
    self.params = params
    self.info = info
    self.output_dir = output_dir


  def table_one(self):
    ''' Constructs Table 1 for GUI or logging '''

    A = u'\N{ANGSTROM SIGN}'
    d = u'\N{DEGREE SIGN}'
    a = u'\N{GREEK SMALL LETTER ALPHA}'
    b = u'\N{GREEK SMALL LETTER BETA}'
    g = u'\N{GREEK SMALL LETTER GAMMA}'
    s = u'\N{GREEK SMALL LETTER SIGMA}'
    h = u'\u00BD'

    uc_edges = '{:4.2f}, {:4.2f}, {:4.2f}'.format(self.info['mean_a'][-1],
                                                  self.info['mean_b'][-1],
                                                  self.info['mean_c'][-1])
    uc_angles = '{:4.2f}, {:4.2f}, {:4.2f}'.format(self.info['mean_alpha'][-1],
                                                   self.info['mean_beta'][-1],
                                                   self.info['mean_gamma'][-1])
    res_total = '{:4.2f} - {:4.2f}'.format(self.info['total_res_max'][-1],
                                           self.info['total_res_min'][-1])
    res_last_shell = '{:4.2f} - {:4.2f}' \
                     ''.format(self.info['binned_resolution'][-1][-1],
                               self.info['binned_resolution'][-1][-2])
    t1_rlabels = [u'No. of images'.encode('utf-8'),
                  u'Space Group'.encode('utf-8'),
                  u'Cell dimensions'.encode('utf-8'),
                  u'  a, b, c ({})  '.format(A).encode('utf-8'),
                  u'  {}, {}, {} ({})    '.format(a, b, g, d).encode('utf-8'),
                  u'Resolution ({})  '.format(A).encode('utf-8'),
                  u'Completeness (%)'.encode('utf-8'),
                  u'Multiplicity'.encode('utf-8'),
                  u'I / {}(I) '.format(s).encode('utf-8'),
                  u'CC{} '.format(h).encode('utf-8'),
                  u'R_merge'.encode('utf-8')]

    t1_data = [['{}'.format(self.info['n_frames_good'][-1])],
               ['{}'.format(self.info['space_group_info'][
                              -1].symbol_and_number())],
               [''],
               ['{}'.format(uc_edges)],
               ['{}'.format(uc_angles)],
               ['{} ({})'.format(res_total, res_last_shell)],
               ['{:4.2f} ({:4.2f})'.format(self.info['total_completeness'][-1],
                                    self.info['binned_completeness'][-1][-1])],
               ['{:4.2f} ({:4.2f})'.format(self.info['total_n_obs'][-1],
                                          self.info['binned_n_obs'][-1][-1])],
               ['{:4.2f} ({:4.2f})'.format(self.info['total_i_o_sigi'][-1],
                                        self.info['binned_i_o_sigi'][-1][-1])],
               ['{:4.2f} ({:4.2f})'.format(self.info['total_cc12'][-1],
                                          self.info['binned_cc12'][-1][-1])],
               ['{:4.2f} ({:4.2f})'.format(self.info['total_rmerge'][-1],
                                          self.info['binned_rmerge'][-1][-1])]
               ]

    return t1_rlabels, t1_data


  def stat_charts(self):
    ''' Displays charts of CC1/2, Completeness, Multiplicity and I / sig(I)
        per resolution bin after the final cycle of post-refinement '''

    fig = plt.figure(figsize=(12, 12))
    gsp = gridspec.GridSpec(2, 2)

    fig.set_alpha(0)
    plt.rc('font', family='sans-serif', size=12)
    plt.rc('mathtext', default='regular')

    x = self.info['binned_resolution'][-1]
    bins = np.arange(len(x))
    xticks = bins[0::len(bins) // 6]
    xlabels = ["{:.2f}".format(i) for i in x]

    # Plot CC1/2 vs. resolution
    ax_cc12 = fig.add_subplot(gsp[0])
    reslabel = 'Resolution ({})'.format(r'$\AA$')
    ax_cc12.set_xlabel(reslabel, fontsize=15)
    ax_cc12.set_ylabel(r'$CC_{1/2}$ (%)', fontsize=15)
    ax_cc12.set_xticks(xticks)
    ax_cc12.set_xticklabels(xlabels)
    ax_cc12.grid(True)
    ax_cc12.plot(bins, self.info['binned_cc12'][0], c='#7fcdbb', lw=2)
    ax_cc12.plot(bins, self.info['binned_cc12'][-1], c='#2c7fb8', lw=3)

    # Plot Completeness vs. resolution
    ax_comp = fig.add_subplot(gsp[1])
    ax_comp.set_xlabel(reslabel, fontsize=15)
    ax_comp.set_ylabel('Completeness (%)', fontsize=15)
    ax_comp.set_xticks(xticks)
    ax_comp.set_xticklabels(xlabels)
    ax_comp.grid(True)
    ax_comp.plot(bins, self.info['binned_completeness'][0], c='#7fcdbb', lw=2)
    ax_comp.plot(bins, self.info['binned_completeness'][-1], c='#2c7fb8', lw=3)

    # Plot Multiplicity (no. of observations) vs. resolution
    ax_mult = fig.add_subplot(gsp[2])
    ax_mult.set_xlabel(reslabel, fontsize=15)
    ax_mult.set_ylabel('# of Observations', fontsize=15)
    ax_mult.set_xticks(xticks)
    ax_mult.set_xticklabels(xlabels)
    ax_mult.grid(True)
    ax_mult.plot(bins, self.info['binned_n_obs'][0], c='#7fcdbb', lw=2)
    ax_mult.plot(bins, self.info['binned_n_obs'][-1], c='#2c7fb8', lw=3)

    # Plot I / sig(I) vs. resolution
    ax_i_sigi = fig.add_subplot(gsp[3])
    ax_i_sigi.set_xlabel(reslabel, fontsize=15)
    ax_i_sigi.set_ylabel(r'I / $\sigma$(I)', fontsize=15)
    ax_i_sigi.set_xticks(xticks)
    ax_i_sigi.set_xticklabels(xlabels)
    ax_i_sigi.grid(True)
    ax_i_sigi.plot(bins, self.info['binned_i_o_sigi'][0], c='#7fcdbb', lw=2)
    ax_i_sigi.plot(bins, self.info['binned_i_o_sigi'][-1], c='#2c7fb8', lw=3)

    plt.tight_layout()
    plt.show()
