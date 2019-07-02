# LIBTBX_SET_DISPATCHER_NAME prime.plotstats
'''
Author      : Uervirojnangkoorn, M.
Created     : 7/11/2016
Description : Plot stats by post-refinement cycles.
'''
from __future__ import absolute_import, division, print_function
import matplotlib.pyplot as plt
import sys
from six.moves import cPickle as pickle
from six.moves import range

if len(sys.argv)==1:
  print('Use prime.plotstats to view different stats along post-refinement cycles.')
  print('Usage: prime.plotstats path/to/run/folder')
  exit()
run_no = sys.argv[1]
stat_pickle = pickle.load(open(run_no+"/pickle.stat","rb"))
total_i_o_sigi = stat_pickle['total_i_o_sigi']
total_completeness = stat_pickle['total_completeness']
total_rmerge = stat_pickle['total_rmerge']
total_n_obs = stat_pickle['total_n_obs']
total_cc12 = stat_pickle['total_cc12']
x = list(range(len(total_cc12)))
#plot
plt.subplot(1,2,1, title='CC1/2')
plt.plot(x, total_cc12, linewidth=2.0)
plt.grid(True)
plt.ylim([min(total_cc12)-5,100])
plt.xlabel('No. of cycles')
plt.subplot(1,2,2, title='Rmerge')
plt.plot(x, total_rmerge, linewidth=2.0)
plt.grid(True)
plt.xlabel('No. of cycles')
plt.show()
