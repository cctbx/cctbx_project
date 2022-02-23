from __future__ import division, unicode_literals, print_function,\
  absolute_import

'''
Author      : Lyubimov, A.Y.
Created     : 05/25/2016
Last Changed: 07/17/2019
Description : PRIME Result Plotter module
'''

import numpy as np

from matplotlib import gridspec, rc

from iota.gui.plotter import Plotter as IOTAPlotter

class Plotter(IOTAPlotter):
  ''' Class with function to plot various PRIME charts (includes Table 1) '''

  def __init__(self, parent, info, output_dir=None, anomalous_flag=False,
               *args, **kwargs):
    super(Plotter, self).__init__(parent=parent, info=info, *args, **kwargs)
    self.target_anomalous_flag = anomalous_flag
    self.output_dir = output_dir

  def table_one(self, as_tex=False):
    ''' Constructs Table 1 for GUI or logging '''

    if as_tex:
      A = r'$\AA$'
      d = r'$\circ$'
      a = r'$\alpha$'
      b = r'$\beta$'
      g = r'$\gamma$'
      s = r'$\sigma$'
      h = r'$\frac{1}{2}$'
      rm = r'$\_merge$'
    else:
      A = '\N{ANGSTROM SIGN}'
      d = '\N{DEGREE SIGN}'
      a = '\N{GREEK SMALL LETTER ALPHA}'
      b = '\N{GREEK SMALL LETTER BETA}'
      g = '\N{GREEK SMALL LETTER GAMMA}'
      s = '\N{GREEK SMALL LETTER SIGMA}'
      h = '\u00BD'
      rm = '_merge'
    t1_rlabels = ['No. of accepted images',
                  'No. of rejected images',
                  'Space Group',
                  'Cell dimensions',
                  '  a, b, c ({})  '.format(A),
                  '  {}, {}, {} ({})    '.format(a, b, g, d),
                  'Resolution ({})  '.format(A),
                  'Completeness (%)',
                  'Multiplicity',
                  'I / {}(I) '.format(s),
                  'CC{} '.format(h),
                  'R{}'.format(rm)]

    uc_edges = '{:4.2f}, {:4.2f}, {:4.2f}'.format(self.info['mean_a'][-1],
                                                  self.info['mean_b'][-1],
                                                  self.info['mean_c'][-1])
    uc_angles = '{:4.2f}, {:4.2f}, {:4.2f}'.format(self.info['mean_alpha'][-1],
                                                   self.info['mean_beta'][-1],
                                                   self.info['mean_gamma'][-1])
    res_total = '{:4.2f} - {:4.2f}'.format(self.info['total_res_max'][-1],
                                           self.info['total_res_min'][-1])
    res_last_shell = '{:4.2f} - {:4.2f}' \
                     ''.format(self.info['binned_resolution'][-1][-2],
                               self.info['binned_resolution'][-1][-1])

    n_frames_bad = self.info['n_frames_bad_cc'][-1]      + \
                   self.info['n_frames_bad_G'][-1]       + \
                   self.info['n_frames_bad_uc'][-1]      + \
                   self.info['n_frames_bad_gamma_e'][-1] + \
                   self.info['n_frames_bad_SE'][-1]
    t1_data = [['{}'.format(self.info['n_frames_good'][-1])],
               ['{}'.format(n_frames_bad)],
               ['{}'.format(self.info['space_group_info'][-1])],
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
                                          self.info['binned_cc12'][-1][-1]*100)],
               ['{:4.2f} ({:4.2f})'.format(self.info['total_rmerge'][-1],
                                          self.info['binned_rmerge'][-1][-1])]
               ]

    return t1_rlabels, t1_data

  def flatten_table_one_data(self, as_tex=False):
    rlabels, tb1_data = self.table_one(as_tex=as_tex)

    # Flatten data list of lists (works this once, since each sub-list
    # contains a single item)
    tb1_data = list(zip(*tb1_data))[0]
    data = list(zip(rlabels, tb1_data))
    return data

  def table_one_text(self):
    data = self.flatten_table_one_data()
    return self.plot_table_text(data=data)

  def table_one_figure(self):
    data = self.flatten_table_one_data()
    self.plot_table(data=data)

  def stat_charts(self):
    ''' Displays charts of CC1/2, Completeness, Multiplicity and I / sig(I)
        per resolution bin after the final cycle of post-refinement '''

    gsp = gridspec.GridSpec(2, 2)
    gsp.update(wspace=0.1, hspace=0.1)

    self.figure.set_alpha(0)
    rc('font', family='sans-serif', size=12)
    rc('mathtext', default='regular')

    x = self.info['binned_resolution'][-1]
    bins = np.arange(len(x))
    xlabels = ["{:.2f}".format(i) for i in x]
    sel_bins = bins[0::len(bins) // 6]
    sel_xlabels = [xlabels[t] for t in sel_bins]

    # Plot CC1/2 vs. resolution
    ax_cc12 = self.figure.add_subplot(gsp[0])
    reslabel = 'Resolution ({})'.format(r'$\AA$')
    ax_cc12.set_xlabel(reslabel, fontsize=15)
    ax_cc12.ticklabel_format(axis='y', style='plain')
    ax_cc12.set_ylim(0, 100)

    if self.target_anomalous_flag:
      ax_cc12.set_ylabel(r'$CC_{1/2}$ anom (%)', fontsize=15)
    else:
      ax_cc12.set_ylabel(r'$CC_{1/2}$ (%)', fontsize=15)
    ax_cc12.set_xticks(sel_bins)
    ax_cc12.set_xticklabels(sel_xlabels)
    ax_cc12.grid(True)
    cc12_start_percent = [c * 100 for c in self.info['binned_cc12'][0]]
    cc12_end_percent = [c * 100 for c in self.info['binned_cc12'][-1]]
    start, = ax_cc12.plot(bins, cc12_start_percent, c='#7fcdbb', lw=2)
    end, = ax_cc12.plot(bins, cc12_end_percent, c='#2c7fb8', lw=3)
    labels = ['Initial', 'Final']
    ax_cc12.legend([start, end], labels, loc='lower left',
                          fontsize=9, fancybox=True)

    # Plot Completeness vs. resolution
    ax_comp = self.figure.add_subplot(gsp[1])
    ax_comp.set_xlabel(reslabel, fontsize=15)
    ax_comp.ticklabel_format(axis='y', style='plain')
    ax_comp.set_ylabel('Completeness (%)', fontsize=15)
    ax_comp.set_xticks(sel_bins)
    ax_comp.set_xticklabels(sel_xlabels)
    ax_comp.set_ylim(0, 100)
    ax_comp.grid(True)
    start, = ax_comp.plot(bins, self.info['binned_completeness'][0],
                         c='#7fcdbb', lw=2)
    end, = ax_comp.plot(bins, self.info['binned_completeness'][-1], c='#2c7fb8',
                   lw=3)
    labels = ['Initial', 'Final']
    ax_comp.legend([start, end], labels, loc='lower left',
                          fontsize=9, fancybox=True)

    # Plot Multiplicity (no. of observations) vs. resolution
    ax_mult = self.figure.add_subplot(gsp[2])
    ax_mult.set_xlabel(reslabel, fontsize=15)
    ax_mult.ticklabel_format(axis='y', style='plain')
    ax_mult.set_ylabel('# of Observations', fontsize=15)
    ax_mult.set_xticks(sel_bins)
    ax_mult.set_xticklabels(sel_xlabels)
    ax_mult.grid(True)
    start, = ax_mult.plot(bins, self.info['binned_n_obs'][0], c='#7fcdbb', lw=2)
    end, = ax_mult.plot(bins, self.info['binned_n_obs'][-1], c='#2c7fb8', lw=3)
    labels = ['Initial', 'Final']
    ax_mult.legend([start, end], labels, loc='lower left',
                          fontsize=9, fancybox=True)

    # Plot I / sig(I) vs. resolution
    ax_i_sigi = self.figure.add_subplot(gsp[3])
    ax_i_sigi.set_xlabel(reslabel, fontsize=15)
    ax_i_sigi.ticklabel_format(axis='y', style='plain')
    ax_i_sigi.set_ylabel(r'I / $\sigma$(I)', fontsize=15)
    ax_i_sigi.set_xticks(sel_bins)
    ax_i_sigi.set_xticklabels(sel_xlabels)
    ax_i_sigi.grid(True)
    start, = ax_i_sigi.plot(bins, self.info['binned_i_o_sigi'][0], c='#7fcdbb',
                        lw=2)
    end, = ax_i_sigi.plot(bins, self.info['binned_i_o_sigi'][-1], c='#2c7fb8',
                       lw=3)
    labels = ['Initial', 'Final']
    ax_i_sigi.legend([start, end], labels, loc='lower left',
                          fontsize=9, fancybox=True)

    # self.figure.set_tight_layout(True)
    self.draw(tight_layout=False)
