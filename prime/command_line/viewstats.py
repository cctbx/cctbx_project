from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME prime.viewstats
'''
Author      : Uervirojnangkoorn, M.
Created     : 9/19/2014
Description : View convergence and other stats for post-refinement.
'''

#read log file and summarize refinement results
import matplotlib.pyplot as plt
import sys
import numpy as np


file_a = open(sys.argv[1],'r')
data_a=file_a.read().split("\n")
search_word = 'refinement and merging'
line_no_now = 0
cc12_all=np.array([0.0]*10)
cciso_all=np.array([0.0]*10)
fpr_all=np.array([0.0]*10)
fxy_all=np.array([0.0]*10)

G_all=np.array([0.0]*10)
B_all=np.array([0.0]*10)
rotx_all=np.array([0.0]*10)
roty_all=np.array([0.0]*10)
ry_all=np.array([0.0]*10)
rx_all=np.array([0.0]*10)
re_all=np.array([0.0]*10)
a_all=np.array([0.0]*10)
b_all=np.array([0.0]*10)
c_all=np.array([0.0]*10)
alpha_all=np.array([0.0]*10)
beta_all=np.array([0.0]*10)
gamma_all=np.array([0.0]*10)

fpr_std=np.array([0.0]*10)
fxy_std=np.array([0.0]*10)
G_std=np.array([0.0]*10)
B_std=np.array([0.0]*10)
rotx_std=np.array([0.0]*10)
roty_std=np.array([0.0]*10)
ry_std=np.array([0.0]*10)
rx_std=np.array([0.0]*10)
re_std=np.array([0.0]*10)
a_std=np.array([0.0]*10)
b_std=np.array([0.0]*10)
c_std=np.array([0.0]*10)
alpha_std=np.array([0.0]*10)
beta_std=np.array([0.0]*10)
gamma_std=np.array([0.0]*10)

i_row = 0
fpr_init = 0
fpr_std_init = 0
fxy_init = 0
fxy_std_init = 0
for line in data_a:
  if line.find(search_word) > 0:
    data_col = data_a[line_no_now-2].split()
    cc12_all[i_row] = float(data_col[6])
    cciso_all[i_row] = float(data_col[8])

    if i_row == 1:
      #collect post-refinement target before
      data_col = data_a[line_no_now+9].split()
      fpr_init = float(data_col[2])
      fpr_std_init = float(data_col[len(data_col)-1].strip('(').strip(')'))

      data_col = data_a[line_no_now+10].split()
      fxy_init = float(data_col[3])
      fxy_std_init = float(data_col[len(data_col)-1].strip('(').strip(')'))

    data_col = data_a[line_no_now+12].split()
    fpr_all[i_row] = float(data_col[2])
    fpr_std[i_row] = float(data_col[len(data_col)-1].strip('(').strip(')'))

    data_col = data_a[line_no_now+13].split()
    fxy_all[i_row] = float(data_col[3])
    fxy_std[i_row] = float(data_col[len(data_col)-1].strip('(').strip(')'))

    data_col = data_a[line_no_now+14].split()
    G_all[i_row] = float(data_col[1])
    G_std[i_row] = float(data_col[len(data_col)-1].strip('(').strip(')'))

    data_col = data_a[line_no_now+15].split()
    B_all[i_row] = float(data_col[1])
    B_std[i_row] = float(data_col[len(data_col)-1].strip('(').strip(')'))

    data_col = data_a[line_no_now+16].split()
    rotx_all[i_row] = float(data_col[1])
    rotx_std[i_row] = float(data_col[len(data_col)-1].strip('(').strip(')'))

    data_col = data_a[line_no_now+17].split()
    roty_all[i_row] = float(data_col[1])
    roty_std[i_row] = float(data_col[len(data_col)-1].strip('(').strip(')'))

    data_col = data_a[line_no_now+18].split()
    ry_all[i_row] = float(data_col[1])
    ry_std[i_row] = float(data_col[len(data_col)-1].strip('(').strip(')'))

    data_col = data_a[line_no_now+19].split()
    rx_all[i_row] = float(data_col[1])
    rx_std[i_row] = float(data_col[len(data_col)-1].strip('(').strip(')'))

    data_col = data_a[line_no_now+19].split()
    re_all[i_row] = float(data_col[1])
    re_std[i_row] = float(data_col[len(data_col)-1].strip('(').strip(')'))

    data_col = data_a[line_no_now+22].split()
    a_all[i_row] = float(data_col[1])
    a_std[i_row] = float(data_col[len(data_col)-1].strip('(').strip(')'))

    data_col = data_a[line_no_now+22].split()
    b_all[i_row] = float(data_col[1])
    b_std[i_row] = float(data_col[len(data_col)-1].strip('(').strip(')'))

    data_col = data_a[line_no_now+23].split()
    c_all[i_row] = float(data_col[1])
    c_std[i_row] = float(data_col[len(data_col)-1].strip('(').strip(')'))

    data_col = data_a[line_no_now+25].split()
    alpha_all[i_row] = float(data_col[1])
    alpha_std[i_row] = float(data_col[len(data_col)-1].strip('(').strip(')'))

    data_col = data_a[line_no_now+26].split()
    beta_all[i_row] = float(data_col[1])
    beta_std[i_row] = float(data_col[len(data_col)-1].strip('(').strip(')'))

    data_col = data_a[line_no_now+27].split()
    gamma_all[i_row] = float(data_col[1])
    gamma_std[i_row] = float(data_col[len(data_col)-1].strip('(').strip(')'))

    i_row += 1


  line_no_now += 1

#reset fpr and fxy at 0th cycle
fpr_all[0] = fpr_init
fpr_std[0] = fpr_std_init
fxy_all[0] = fxy_init
fxy_std[0] = fxy_std_init

n_cycle = 0
for G, G_s, B, B_s, rotx, rotx_s, roty, roty_s, ry, ry_s, rx, rx_s, \
    re, re_s, a, a_s, b, b_s, c, c_s, alpha, alpha_s, beta, beta_s, gamma, gamma_s, \
    fpr, fpr_s, fxy, fxy_s, cc12, cciso in \
    zip(G_all, G_std, B_all, B_std, rotx_all, rotx_std, roty_all, roty_std, \
        ry_all, ry_std, rx_all, rx_std, re_all, re_std, \
        a_all, a_std, b_all, b_std, c_all, c_std, \
        alpha_all, alpha_std, beta_all, beta_std, gamma_all, gamma_std, \
        fpr_all, fpr_std, fxy_all, fxy_std, cc12_all, cciso_all):
  print '%2.0f %6.2f %6.2f %6.2f %6.2f %6.4f %6.4f %6.4f %6.4f %6.5f %6.5f %6.5f %6.5f %6.5f %6.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f'%(n_cycle, G, G_s, B, B_s, rotx, rotx_s, roty, roty_s, ry, ry_s, rx, rx_s, \
    re, re_s, a, a_s, b, b_s, c, c_s, alpha, alpha_s, beta, beta_s, gamma, gamma_s, \
    fpr, fpr_s, fxy, fxy_s, cc12, cciso)
  n_cycle += 1

  if n_cycle == i_row:
    break

#convert absolute parameters to delta
fpr_all_delta=np.array([0.0]*10)
fxy_all_delta=np.array([0.0]*10)
G_all_delta=np.array([0.0]*10)
B_all_delta=np.array([0.0]*10)
ry_all_delta=np.array([0.0]*10)
rx_all_delta=np.array([0.0]*10)
re_all_delta=np.array([0.0]*10)
a_all_delta=np.array([0.0]*10)
b_all_delta=np.array([0.0]*10)
c_all_delta=np.array([0.0]*10)
alpha_all_delta=np.array([0.0]*10)
beta_all_delta=np.array([0.0]*10)
gamma_all_delta=np.array([0.0]*10)

fpr_std_delta=np.array([0.0]*10)
fxy_std_delta=np.array([0.0]*10)
G_std_delta=np.array([0.0]*10)
B_std_delta=np.array([0.0]*10)
ry_std_delta=np.array([0.0]*10)
rx_std_delta=np.array([0.0]*10)
re_std_delta=np.array([0.0]*10)
a_std_delta=np.array([0.0]*10)
b_std_delta=np.array([0.0]*10)
c_std_delta=np.array([0.0]*10)
alpha_std_delta=np.array([0.0]*10)
beta_std_delta=np.array([0.0]*10)
gamma_std_delta=np.array([0.0]*10)
for i in range(i_row):
  if i == 0:
    fpr_all_delta[i] = fpr_all[i]
    fpr_std_delta[i] = fpr_std[i]
    fxy_all_delta[i] = fxy_all[i]
    fxy_std_delta[i] = fxy_std[i]
    G_all_delta[i] = G_all[i]
    G_std_delta[i] = G_std[i]
    B_all_delta[i] = B_all[i]
    B_std_delta[i] = B_std[i]
    ry_all_delta[i] = ry_all[i]
    ry_std_delta[i] = ry_std[i]
    rx_all_delta[i] = rx_all[i]
    rx_std_delta[i] = rx_std[i]
    re_all_delta[i] = re_all[i]
    re_std_delta[i] = re_std[i]
    a_all_delta[i] = a_all[i]
    a_std_delta[i] = a_std[i]
    b_all_delta[i] = b_all[i]
    b_std_delta[i] = b_std[i]
    c_all_delta[i] = c_all[i]
    c_std_delta[i] = c_std[i]
    alpha_all_delta[i] = alpha_all[i]
    alpha_std_delta[i] = alpha_std[i]
    beta_all_delta[i] = beta_all[i]
    beta_std_delta[i] = beta_std[i]
    gamma_all_delta[i] = gamma_all[i]
    gamma_std_delta[i] = gamma_std[i]
  else:
    fpr_all_delta[i] = np.absolute(fpr_all[i] - fpr_all[i-1])
    fpr_std_delta[i] = np.absolute(fpr_std[i] - fpr_std[i-1])
    fxy_all_delta[i] = np.absolute(fxy_all[i] - fxy_all[i-1])
    fxy_std_delta[i] = np.absolute(fxy_std[i] - fxy_std[i-1])
    G_all_delta[i] = np.absolute(G_all[i] - G_all[i-1])
    G_std_delta[i] = np.absolute(G_std[i] - G_std[i-1])
    B_all_delta[i] = np.absolute(B_all[i] - B_all[i-1])
    B_std_delta[i] = np.absolute(B_std[i] - B_std[i-1])
    ry_all_delta[i] = np.absolute(ry_all[i] - ry_all[i-1])
    ry_std_delta[i] = np.absolute(ry_std[i] - ry_std[i-1])
    rx_all_delta[i] = np.absolute(rx_all[i] - rx_all[i-1])
    rx_std_delta[i] = np.absolute(rx_std[i] - rx_std[i-1])
    re_all_delta[i] = np.absolute(re_all[i] - re_all[i-1])
    re_std_delta[i] = np.absolute(re_std[i] - re_std[i-1])
    a_all_delta[i] = np.absolute(a_all[i] - a_all[i-1])
    a_std_delta[i] = np.absolute(a_std[i] - a_std[i-1])
    b_all_delta[i] = np.absolute(b_all[i] - b_all[i-1])
    b_std_delta[i] = np.absolute(b_std[i] - b_std[i-1])
    c_all_delta[i] = np.absolute(c_all[i] - c_all[i-1])
    c_std_delta[i] = np.absolute(c_std[i] - c_std[i-1])
    alpha_all_delta[i] = np.absolute(alpha_all[i] - alpha_all[i-1])
    alpha_std_delta[i] = np.absolute(alpha_std[i] - alpha_std[i-1])
    beta_all_delta[i] = np.absolute(beta_all[i] - beta_all[i-1])
    beta_std_delta[i] = np.absolute(beta_std[i] - beta_std[i-1])
    gamma_all_delta[i] = np.absolute(gamma_all[i] - gamma_all[i-1])
    gamma_std_delta[i] = np.absolute(gamma_std[i] - gamma_std[i-1])


#prepare data for plotting
n_data = 200
G_series = []
B_series = []
rotx_series = []
roty_series = []
ry_series = []
rx_series = []
re_series = []
a_series = []
b_series = []
c_series = []
alpha_series = []
gamma_series = []
beta_series = []
fpr_series = []
fxy_series = []

for i in range(1,i_row):
  if G_std_delta[i] < 0.0000001:
    narr = np.zeros(n_data)
  else:
    narr = np.random.normal(G_all_delta[i], G_std_delta[i], n_data)

  G_series.append(narr)

  if B_std_delta[i] < 0.0000001:
    narr = np.zeros(n_data)
  else:
    narr = np.random.normal(B_all_delta[i], B_std_delta[i], n_data)
  B_series.append(narr)

  if rotx_std[i] < 0.0000001:
    narr = np.zeros(n_data)
  else:
    narr = np.random.normal(rotx_all[i], rotx_std[i], n_data)
  rotx_series.append(narr)

  if roty_std[i] < 0.0000001:
    narr = np.zeros(n_data)
  else:
    narr = np.random.normal(roty_all[i], roty_std[i], n_data)
  roty_series.append(narr)

  if ry_std_delta[i] < 0.0000001:
    narr = np.zeros(n_data)
  else:
    narr = np.random.normal(ry_all_delta[i], ry_std_delta[i], n_data)
  ry_series.append(narr)

  if rx_std_delta[i] < 0.0000001:
    narr = np.zeros(n_data)
  else:
    narr = np.random.normal(rx_all_delta[i], rx_std_delta[i], n_data)
  rx_series.append(narr)

  if re_std_delta[i] < 0.0000001:
    narr = np.random.normal(re_all_delta[i], 0.00001, n_data)
  else:
    narr = np.random.normal(re_all_delta[i], re_std_delta[i], n_data)
  re_series.append(narr)

  print a_all_delta[i], a_std_delta[i], c_all_delta[i], c_std_delta[i]
  if a_std_delta[i] < 0.0000001:
    narr = np.zeros(n_data)
  else:
    narr = np.random.normal(a_all_delta[i], a_std_delta[i], n_data)
  a_series.append(narr)

  if b_std_delta[i] < 0.0000001:
    narr = np.zeros(n_data)
  else:
    narr = np.random.normal(b_all_delta[i], b_std_delta[i], n_data)
  b_series.append(narr)

  if c_std_delta[i] < 0.0000001:
    narr = np.zeros(n_data)
  else:
    narr = np.random.normal(c_all_delta[i], c_std_delta[i], n_data)
  c_series.append(narr)

  if alpha_std_delta[i] < 0.0000001:
    narr = np.zeros(n_data)
  else:
    narr = np.random.normal(alpha_all_delta[i], alpha_std_delta[i], n_data)
  alpha_series.append(narr)

  if beta_std_delta[i] < 0.0000001:
    narr = np.zeros(n_data)
  else:
    narr = np.random.normal(beta_all_delta[i], beta_std_delta[i], n_data)
  beta_series.append(narr)

  if gamma_std_delta[i] < 0.0000001:
    narr = np.zeros(n_data)
  else:
    narr = np.random.normal(gamma_all_delta[i], gamma_std_delta[i], n_data)
  gamma_series.append(narr)

  if fpr_std_delta[i] < 0.0000001:
    narr = np.zeros(n_data)
  else:
    narr = np.random.normal(fpr_all_delta[i], fpr_std_delta[i], n_data)
  fpr_series.append(narr)

  if fxy_std_delta[i] < 0.0000001:
    narr = np.zeros(n_data)
  else:
    narr = np.random.normal(fxy_all_delta[i], fxy_std_delta[i], n_data)
  fxy_series.append(narr)


x_range = range(1,i_row)
x_label = []
for i in range(1,i_row):
  x_label.append(str(i))

font_plot = {'family':'Helvetica',
             'weight': 'normal',
             'size'  : 12,
            }
#begin plotting
ax = plt.subplot(331, title=r'$\Delta$G')
plt.boxplot(G_series)
#plt.ylabel(r'$\Delta$G', rotation=0)
plt.xticks(x_range, x_label)
ax.title.set_fontweight('bold')
ax.title.set_fontsize(12)

ax = plt.subplot(332, title=r'$\Delta$B-factor')
plt.boxplot(B_series)
plt.xticks(x_range, x_label)
ax.title.set_fontweight('bold')
ax.title.set_fontsize(12)

ax = plt.subplot(333, title=r'$\Delta$rotation on x-axis ($\degree$)')
plt.boxplot(rotx_series)
plt.xticks(x_range, x_label)
ax.title.set_fontweight('bold')
ax.title.set_fontsize(12)

ax = plt.subplot(334, title=r'$\Delta$rotation on y-axis ($\degree$)')
plt.boxplot(roty_series)
plt.xticks(x_range, x_label)
ax.title.set_fontweight('bold')
ax.title.set_fontsize(12)

ax = plt.subplot(335, title=r'$\Delta\gamma_y (\frac{1}{\AA})$')
plt.boxplot(ry_series)
plt.xticks(x_range, x_label)
ax.title.set_fontweight('bold')
ax.title.set_fontsize(14)

ax = plt.subplot(336, title=r'$\Delta\gamma_x (\frac{1}{\AA})$')
plt.boxplot(rx_series)
plt.xticks(x_range, x_label)
ax.title.set_fontweight('bold')
ax.title.set_fontsize(14)

ax = plt.subplot(337, title=r'$\Delta\gamma_e (\frac{1}{\AA})$')
plt.boxplot(re_series)
plt.xticks(x_range, x_label)
ax.title.set_fontweight('bold')
ax.title.set_fontsize(14)

ax = plt.subplot(338, title=r'$\Delta<F_{pr}>$')
plt.boxplot(fpr_series)
plt.xticks(x_range, x_label)
ax.title.set_fontweight('bold')
ax.title.set_fontsize(12)

ax = plt.subplot(339, title=r'$\Delta <F_{xy}> (mm^2)$')
plt.boxplot(fxy_series)
plt.xticks(x_range, x_label)
ax.title.set_fontweight('bold')
ax.title.set_fontsize(12)
plt.show()


ax = plt.subplot(231, title=r'$\Delta a (\AA)$')
plt.boxplot(a_series)
plt.xticks(x_range, x_label)
ax.title.set_fontweight('bold')
ax.title.set_fontsize(12)

ax = plt.subplot(232, title=r'$\Delta b (\AA)$')
plt.boxplot(b_series)
plt.xticks(x_range, x_label)
ax.title.set_fontweight('bold')
ax.title.set_fontsize(12)

ax = plt.subplot(233, title=r'$\Delta c (\AA)$')
plt.boxplot(c_series)
plt.xticks(x_range, x_label)
ax.title.set_fontweight('bold')
ax.title.set_fontsize(12)

ax = plt.subplot(234, title=r'$\Delta \alpha (\degree)$')
plt.boxplot(alpha_series)
plt.xticks(x_range, x_label)
ax.title.set_fontweight('bold')
ax.title.set_fontsize(12)

ax = plt.subplot(235, title=r'$\Delta \beta (\degree)$')
plt.boxplot(beta_series)
plt.xticks(x_range, x_label)
ax.title.set_fontweight('bold')
ax.title.set_fontsize(12)

ax = plt.subplot(236, title=r'$\Delta \gamma (\degree)$')
plt.boxplot(gamma_series)
plt.xticks(x_range, x_label)
ax.title.set_fontweight('bold')
ax.title.set_fontsize(12)

plt.show()
