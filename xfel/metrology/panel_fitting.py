# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

import math
from dials.array_family import flex
import numpy as np
from matplotlib import pyplot as plt
from sklearn.covariance import EmpiricalCovariance, MinCovDet

class Panel_MCD_Filter(object):
  """For analyzing CALC-OBS residuals on a detector panel, look at a population
of Bragg spots from a refl table.  Calculate a 3-feature fit using Minimum
Covariance Determinat, using Δx, Δy, ΔΨ(deg) as the 3 features.  Then filter
the refls using a Mahalanobis cutoff.
"""
  def __init__(self, lab_coords_x, lab_coords_y, data, i_panel, delta_scalar, params, verbose=True):
    training_data = []

    mean_x = flex.mean(lab_coords_x)
    mean_y = flex.mean(lab_coords_y)
    limit=delta_scalar * 10

    for ix in range(len(data)):
      if abs(lab_coords_x[ix] - mean_x) > limit: continue
      if abs(lab_coords_y[ix] - mean_y) > limit: continue
      if abs(data[ix])>1: continue
      training_data.append((lab_coords_x[ix],lab_coords_y[ix],data[ix]))
    if verbose: print("Training data is less",len(lab_coords_x) - len(training_data),end=" ")
    colorcode_set = []
    for ix in range(len(data)):
      colorcode_set.append((lab_coords_x[ix],lab_coords_y[ix],data[ix]))

    from sklearn.covariance import EmpiricalCovariance, MinCovDet
    # compare estimators learnt from the full data set with true parameters
    emp_cov = EmpiricalCovariance(assume_centered=False, store_precision=True).fit(X=training_data)
    # fit a Minimum Covariance Determinant (MCD) robust estimator to data
    robust_cov = MinCovDet(assume_centered=False, store_precision=True).fit(X=training_data)

    features = ["Δx","Δy","ΔΨ(deg)"]
    if verbose:
      print("%3d"%i_panel,end=" ")
      print("%4d items "%(len(training_data),),end=" ")
    for idx_report in range(len(features)):
      feature = features[idx_report]
      diag_elem = math.sqrt(emp_cov.covariance_[idx_report,idx_report])
      if verbose: print( "%s=%7.2f±%6.2f"%(feature, emp_cov.location_[idx_report], diag_elem),end=" ")

    if verbose: print("%4d items:"%(flex.bool(robust_cov.support_).count(True)),end=" ")
    for idx_report in range(len(features)):
      feature = features[idx_report]
      diag_elem = math.sqrt(robust_cov.covariance_[idx_report,idx_report])
      if verbose: print( "%s=%7.2f±%6.2f"%(feature, robust_cov.location_[idx_report], diag_elem),end=" ")

    disc = flex.double(robust_cov.mahalanobis(X=colorcode_set)) # this metric represents malahanobis ** 2
    disc_select = disc < (params.residuals.mcd_filter.mahalanobis_distance)**2
    if params.residuals.mcd_filter.keep == "outliers":
      disc_select = (disc_select==False)
    if verbose: print("OK %4.1f%%"%(100*(disc_select.count(True))/len(training_data)))
    self.lab_coords_x = lab_coords_x.select(disc_select)
    self.lab_coords_y = lab_coords_y.select(disc_select)
    self.data = data.select(disc_select)
    self.n_input = len(lab_coords_x)
    self.n_output = len(self.lab_coords_x)
    self.emp_cov = emp_cov
    self.rob_cov = robust_cov

  def scatter_coords(self): return self.lab_coords_x, self.lab_coords_y, self.data

  def plot_contours(self, ax, show=False):
    COV = self.emp_cov
    COV_slice = EmpiricalCovariance()
    COV_slice.location_ = np.array([ COV.location_[0], COV.location_[1] ])
    COV_slice.covariance_ = np.array([ COV.covariance_[ 0,0 ], COV.covariance_[ 0,1 ],
                                       COV.covariance_[ 1,0 ], COV.covariance_[ 1,1 ] ])
    COV_slice.covariance_ = COV_slice.covariance_.reshape((2,2))
    COV_slice.precision_ = np.array([ COV.precision_[ 0,0 ], COV.precision_[ 0,1 ],
                                      COV.precision_[ 1,0 ], COV.precision_[ 1,1 ] ])
    COV_slice.precision_ = COV_slice.precision_.reshape((2,2))

    # Show contours of the distance functions
    xx, yy = np.meshgrid(
          np.linspace(COV_slice.location_[0]-5*math.sqrt(COV_slice.covariance_[0,0]), COV_slice.location_[0]+5*math.sqrt(COV_slice.covariance_[0,0]), 100),
          np.linspace(COV_slice.location_[1]-5*math.sqrt(COV_slice.covariance_[1,1]), COV_slice.location_[1]+5*math.sqrt(COV_slice.covariance_[1,1]), 100),
    )
    zz = np.c_[xx.ravel(), yy.ravel()]

    # Empirical fit is not so good.  Don't plot this
    if False: # keep for debugging
      mahal_emp_cov = COV_slice.mahalanobis(zz)
      mahal_emp_cov = mahal_emp_cov.reshape(xx.shape)
      emp_cov_contour = ax.contour(xx, yy, np.sqrt(mahal_emp_cov),
                                  levels=[1.,2.,3.,4.,5.],
                                  #cmap=plt.cm.PuBu_r,
                                  cmap=plt.cm.cool_r,
                                  linestyles='dashed')

    COV = self.rob_cov
    COV_slice = EmpiricalCovariance()
    COV_slice.location_ = np.array([ COV.location_[0], COV.location_[1] ])
    COV_slice.covariance_ = np.array([ COV.covariance_[ 0,0 ], COV.covariance_[ 0,1 ],
                                       COV.covariance_[ 1,0 ], COV.covariance_[ 1,1 ] ])
    COV_slice.covariance_ = COV_slice.covariance_.reshape((2,2))
    COV_slice.precision_ = np.array([ COV.precision_[ 0,0 ], COV.precision_[ 0,1 ],
                                      COV.precision_[ 1,0 ], COV.precision_[ 1,1 ] ])
    COV_slice.precision_ = COV_slice.precision_.reshape((2,2))
    self.robust_model_XY = COV_slice

    # robust is better
    if show:
      mahal_robust_cov = COV_slice.mahalanobis(zz)
      mahal_robust_cov = mahal_robust_cov.reshape(xx.shape)
      robust_contour = ax.contour(xx, yy, np.sqrt(mahal_robust_cov),
                                 levels=[1.,2.,3.,4.,5.],
                                 #cmap=plt.cm.YlOrBr_r,
                                 cmap=plt.cm.spring_r,
                                 linestyles='dotted')

class three_feature_fit(Panel_MCD_Filter):
  """Analyze the relationship between radial, transverse, and deltapsi.
     Calculate a 3-feature fit using Minimum Covariance Determinant, using ΔR, ΔT, ΔΨ(deg)
     as the 3 features.
  """
  def __init__(self, delta_radial, delta_transverse, delta_psi, i_panel, params=None, verbose=True):
    training_data = []
    assert len(delta_radial) == len(delta_transverse) == len(delta_psi)
    for ix in range(len(delta_radial)):
      training_data.append((1000.*delta_radial[ix], 1000*delta_transverse[ix], delta_psi[ix]))

    from sklearn.covariance import EmpiricalCovariance, MinCovDet
    # compare estimators learnt from the full data set with true parameters
    self.emp_cov = EmpiricalCovariance(assume_centered=True, store_precision=True).fit(X=training_data)

    features = ["ΔR","ΔT","ΔΨ(deg)"]
    if verbose:
      print("%3d"%i_panel,end=" ")
      print("%4d items "%(len(training_data),),end=" ")
    diag_elem = flex.double(3)
    self.cross_correl = flex.double(3)
    for idx_report in range(len(features)):
      feature = features[idx_report]
      diag_elem[idx_report] = math.sqrt(self.emp_cov.covariance_[idx_report,idx_report])
      if verbose: print( "%s=%7.2f±%6.2f"%(feature, self.emp_cov.location_[idx_report], diag_elem[idx_report]),end=" ")
    for ic, cross_correlation in enumerate([(0,1),(1,2),(0,2)]):
      feature1 = features[cross_correlation[0]]
      feature2 = features[cross_correlation[1]]
      self.cross_correl[ic] = self.emp_cov.covariance_[cross_correlation[0],cross_correlation[1]] / (
                           diag_elem[cross_correlation[0]]*diag_elem[cross_correlation[1]]) # convert Cov to Corr
      if verbose: print( "%s,%s cc=%5.1f%%"%(feature1,feature2,100*self.cross_correl[ic]),end=" ")
    if verbose: print()
