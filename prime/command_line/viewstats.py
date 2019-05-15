# LIBTBX_SET_DISPATCHER_NAME prime.viewstats
'''
Author      : Uervirojnangkoorn, M.
Created     : 9/19/2014
Description : View convergence and other stats for post-refinement.
'''
from __future__ import absolute_import, division, print_function
import matplotlib.pyplot as plt
import sys
import numpy as np
from six.moves import range

if len(sys.argv)==1:
  print('Use prime.viewstats to view convergence of post-refined parameters.')
  print('Usage: prime.viewstats your_run_no')
  exit()
run_no = sys.argv[1]

#read .paramhist files and display refinement results
import os
cn_file = 0
for file_in in os.listdir(run_no):
  if file_in.endswith('.paramhist'):
    cn_file += 1

if cn_file == 0:
  print('Cannot find .paramhist file in your ', run_no, '.')
  print('To enable viewstats, rerun prime with flag_output_verbose=True in .phil file.')
  print('The .paramhist parameters will be recorded during the run.')
  exit()

param_file_list = []
for i in range(cn_file):
  param_file_list.append(run_no+'/'+str(i)+'.paramhist')

data_dict_list = []
for param_file in param_file_list:
  pf = open(param_file,'r')
  data = pf.read().split('\n')
  data_dict = {}
  n_data = 0
  for data_row in data:
    dc = data_row.split()
    #use row 1 to set n_col
    if n_data == 0:
      n_col = len(dc)
    if len(dc)==n_col:
      data_dict[dc[n_col-1]] = np.array([float(dc[i]) for i in range(n_col-1)])
      n_data += 1
  data_dict_list.append(data_dict)

#prepare test key
data_dict_0 = data_dict_list[0]
data_dict_1 = data_dict_list[1]
for i in range(n_data):
  test_key = list(data_dict_list[0].keys())[i]
  if (test_key in data_dict_0) and (test_key in data_dict_1):
    test_id = i
    break

#Fix Tpr and Txy for first data_dict
test_param_0_raw = np.array(data_dict_0[test_key])
test_param_1 = data_dict_1[test_key]
for key in data_dict_0.keys():
  if key in data_dict_1:
    data_param_0 = data_dict_0[key]
    data_param_1 = data_dict_1[key]
    data_param_0[1] = data_param_1[0]
    data_param_0[3] = data_param_1[2]
    data_dict_0[key] = data_param_0

test_param_0_update = data_dict_0[test_key]
test_delta_1_calc = np.absolute(test_param_1-test_param_0_update)
print('test id', test_id)
print('test key', test_key)
print('0th cycle (raw):', test_param_0_raw)
print('0th cycle (updated):', test_param_0_update)
print('1st cycle:', test_param_1)
print('delta (calc.):', test_delta_1_calc)

delta_dict_list = []
for i in range(len(data_dict_list)-1):
  data_dict = data_dict_list[i]
  data_dict_next = data_dict_list[i+1]
  delta_dict = {}
  for key in data_dict.keys():
    if (key in data_dict_next):
      delta_param = np.absolute(data_dict_next[key] - data_dict[key])
    else:
      delta_param = np.zeros(n_col-1)
    delta_dict[key] = delta_param
  delta_dict_list.append(delta_dict)

delta_dict_0 = delta_dict_list[0]
test_delta_1 = delta_dict_0[test_key]
print('delta (prog.):', test_delta_1)
print('delta diff.:', test_delta_1 - test_delta_1_calc)
print('sum of delta diff.', np.sum(np.absolute(test_delta_1 - test_delta_1_calc)))

x_range = list(range(1, len(delta_dict_list)+1))
x_label = []
for i in range(1, len(delta_dict_list)+1):
  x_label.append(str(i))
data_title = ['Tpr_i','Tpr','Txy_i','Txy','G','B','RotX','RotY','ry','rz','r0','re','voigt_nu','a','b','c','alpha','beta','gamma','CC1/2']
cn_plot = 1
for i in range(n_col-1):
  if i not in (0,2):
    data_series = []
    for delta_dict in delta_dict_list:
      narr = np.array([delta_dict[key][i] for key in delta_dict.keys()])
      data_series.append(narr)

      ax = plt.subplot(3, 6, cn_plot, title=data_title[i])
      plt.boxplot(data_series)
      plt.xticks(x_range, x_label)
      if data_title[i] in ('ry','rz','r0','re'):
        plt.ylim([0, 0.01])
      plt.grid(True)
    cn_plot += 1
plt.show()
