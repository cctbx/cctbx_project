# LIBTBX_SET_DISPATCHER_NAME prime.comparestats
'''
Author      : Uervirojnangkoorn, M.
Created     : 6/25/2016
Description : Given list of log files, compare completeness, n_obs, rmerge, and cc12.
'''
from __future__ import absolute_import, division, print_function
import matplotlib.pyplot as plt
import sys

if len(sys.argv)==1:
  print('Use prime.comparestats to compare stats from different log files.')
  print('Usage: prime.comparestats log1.txt log2.txt log3.txt c=cycle_no n=n_bins')
  exit()
log_files = []
cycle_no = 0
for arg in sys.argv[1:]:
  param = arg.split('=')
  if len(param) == 1:
    log_files.append(param[0])
  elif len(param) == 2:
    if param[0]=='c':
      cycle_no = int(param[1])
    elif param[0] =='n':
      n_bins = int(param[1])

data_dict = {'onedsqr':[], 'complete':[], 'n_obs':[], 'r_merge':[], 'cc12':[]}
for log_file in log_files:
  pf = open(log_file,'r')
  data = pf.read().split('\n')
  if cycle_no == 0:
    search_word = '/mean_scaled_merge.mtz'
  else:
    search_word = '/postref_cycle_'+str(cycle_no)+'_merge.mtz'
  cn_row = 0
  for data_row in data:
    if data_row.endswith(search_word):
      found_row = cn_row + 3
      break
    cn_row += 1
  onedsqr = []
  complete = []
  n_obs = []
  r_merge = []
  cc12 = []
  for data_row in data[found_row:found_row+n_bins]:
    data_col = data_row.split()
    onedsqr.append(1/(float(data_col[3])**2))
    complete.append(float(data_col[4]))
    n_obs.append(float(data_col[8]))
    r_merge.append(float(data_col[9]))
    cc12.append(float(data_col[11]))
  data_dict['onedsqr'].append(onedsqr)
  data_dict['complete'].append(complete)
  data_dict['n_obs'].append(n_obs)
  data_dict['r_merge'].append(r_merge)
  data_dict['cc12'].append(cc12)

plt.subplot(221)
x = data_dict['onedsqr'][0]
for y in data_dict['complete']:
  plt.plot(x,y)
plt.legend(log_files, loc=[0,0])
plt.title('Completeness')
#plt.xlabel('1/d^2')
plt.grid(True)
#plt.ylim([95,100])
plt.subplot(222)
x = data_dict['onedsqr'][0]
for y in data_dict['n_obs']:
  plt.plot(x,y)
#plt.legend(log_files, loc=[0,0])
plt.title('No. observations')
#plt.xlabel('1/d^2')
plt.grid(True)
#plt.ylim([95,100])
plt.subplot(223)
x = data_dict['onedsqr'][0]
for y in data_dict['r_merge']:
  plt.plot(x,y)
#plt.legend(log_files, loc=[0,0])
plt.title('Rmerge')
plt.xlabel('1/d^2')
plt.grid(True)
#plt.ylim([95,100])
plt.subplot(224)
x = data_dict['onedsqr'][0]
for y in data_dict['cc12']:
  plt.plot(x,y)
#plt.legend(log_files, loc=[0,0])
plt.title('CC1/2')
plt.xlabel('1/d^2')
plt.grid(True)
#plt.ylim([95,100])
plt.show()
