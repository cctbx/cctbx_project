from __future__ import division
from io import StringIO

diamonds_short = '''"carat","cut","color","clarity","depth","table","price","x","y","z"
0.23,"Ideal","E","SI2",61.5,55,326,3.95,3.98,2.43
0.21,"Premium","E","SI1",59.8,61,326,3.89,3.84,2.31
0.23,"Good","E","VS1",56.9,65,327,4.05,4.07,2.31
0.29,"Premium","I","VS2",62.4,58,334,4.2,4.23,2.63
0.31,"Good","J","SI2",63.3,58,335,4.34,4.35,2.75
0.24,"Very Good","J","VVS2",62.8,57,336,3.94,3.96,2.48
0.24,"Very Good","I","VVS1",62.3,57,336,3.95,3.98,2.47
0.26,"Very Good","H","SI1",61.9,55,337,4.07,4.11,2.53
0.22,"Fair","E","VS2",65.1,61,337,3.87,3.78,2.49
0.23,"Very Good","H","VS1",59.4,61,338,4,4.05,2.39
0.3,"Good","J","SI1",64,55,339,4.25,4.28,2.73
0.23,"Ideal","J","VS1",62.8,56,340,3.93,3.9,2.46
0.22,"Premium","F","SI1",60.4,61,342,3.88,3.84,2.33
0.31,"Ideal","J","SI2",62.2,54,344,4.35,4.37,2.71
0.2,"Premium","E","SI2",60.2,62,345,3.79,3.75,2.27
0.32,"Premium","E","I1",60.9,58,345,4.38,4.42,2.68
0.3,"Ideal","I","SI2",62,54,348,4.31,4.34,2.68
0.3,"Good","J","SI1",63.4,54,351,4.23,4.29,2.7
0.3,"Good","J","SI1",63.8,56,351,4.23,4.26,2.71
'''

def mahalanobis_using_numpy_and_scipy(x=None, data=None, mu=None, cov=None):
  """
  Compute the Mahalanobis Distance between each row of x and the data
    x    : vector or matrix of data with, say, p columns.
    data : ndarray of the distribution from which Mahalanobis distance of each
           observation of x is to be computed.
    cov  : covariance matrix (p x p) of the distribution. If None, will be
           computed from data.
  """
  # error prone
  assert 0
  import numpy as np
  import scipy as sp
  if mu:
    x_minus_mu = x - mu
  else:
    print(x)
    print(np.mean(data))
    x_minus_mu = x - np.mean(data)
  print('x_minus_mu',x_minus_mu)
  if not cov: cov = np.cov(data)
  print(cov)
  inv_covmat = sp.linalg.inv(cov)
  print(inv_covmat)
  left_term = np.dot(x_minus_mu, inv_covmat)
  mahal = np.dot(left_term, x_minus_mu.T)
  return mahal.diagonal()

def main():
  from libtbx import math_utils
  import csv
  data=[]
  f=StringIO(diamonds_short)
  spamreader = csv.reader(f, delimiter=',', quotechar='|')
  for i, row in enumerate(spamreader):
    if not i: continue
    # print(', '.join(row))
    data.append([float(row[0]), float(row[4]), float(row[6])])

  zz=[(0.23, 61.5, 326),
      (0.23, 56.9, 329),
      (0.23, 56.9, 349),
     ]

  rc=math_utils.mahalanobis(zz, data)
  print(rc)
  rc=math_utils.mahalanobis_p_values(zz, data, verbose=True)
  print(rc)
  rc=math_utils.mahalanobis_p_values_outlier_indices(zz, data)
  assert rc==[2]
  for i, p in enumerate(math_utils.mahalanobis_p_values(zz, data)):
    outlier=''
    if i in rc:
      outlier='OUTLIER'
    print('  %s %0.4f %s' % (i,p,outlier))
  cov = math_utils.covariance_using_sklearn(data, verbose=True)
  rc=math_utils.mahalanobis(zz, cov=cov)
  print(rc)
  rc=math_utils.mahalanobis_p_values(zz, cov=cov, verbose=True)
  print(rc)
  rc=math_utils.mahalanobis_p_values_outlier_indices(zz, cov=cov)
  assert rc==[2]
  rc=math_utils.mahalanobis_p_values_outlier_indices(zz[:1], cov=cov)
  print(rc)

  # values = cov.covariance_.tolist()
  # cov = math_utils.covariance_using_sklearn(values=values, verbose=True)
  # assert 0
  from libtbx import easy_pickle
  easy_pickle.dump('tst_mahal.pickle', cov)
  cov=easy_pickle.load('tst_mahal.pickle')
  # import pickle
  # pf='tst_mahal.pickle'
  # f=open(pf, 'w')
  # pickle.dump(cov, f)
  # del f
  # f=open(pf, 'r')
  # cov=pickle.load(f)
  # del f
  rc=math_utils.mahalanobis_p_values_outlier_indices(zz, cov=cov)
  print(rc)
  assert rc==[2]

if __name__ == '__main__':
  main()
