from __future__ import absolute_import, division, print_function
from six.moves import range
'''
Author      : Uervirojnangkoorn, M.
Created     : 12/1/2014
Description : Utilitiy functions.
'''

import math
from cctbx.array_family import flex

class util_handler(object):
  '''
  classdocs
  '''
  def __init__(self):
    '''
    Constructor
    '''

  def calcmpe(self, phic_given, phiraw_given, is_deg):
    phi_error=[0]*len(phic_given);
    for i_phi in range(len(phic_given)):
      if is_deg:
        phi_error[i_phi]=math.acos(math.cos((phic_given[i_phi]*math.pi/180)-(phiraw_given[i_phi]*math.pi/180)));
      else:
        phi_error[i_phi]=math.acos(math.cos(phic_given[i_phi]-phiraw_given[i_phi]));

    mean_phi_error=sum(phi_error)/len(phi_error);

    return mean_phi_error;


  def calcphicc(self, magf, fomc, fomraw, phic_given, phiraw_given, is_deg):
    '''
    calculate mapcc as phase correlation
    sum((f*fom)^2*(cos(phic-phir)))/sum((f*fom)^2)
    '''
    phi_cc=[0]*len(phic_given)
    phi_error=[0]*len(phic_given)
    sum_magf_sq=0
    for i_phi in range(len(phic_given)):
      f_sq=math.pow(magf[i_phi], 2)*fomc[i_phi]*fomraw[i_phi]
      sum_magf_sq+=f_sq
      if is_deg:
        phi_cc[i_phi]=f_sq*math.cos((phic_given[i_phi]*math.pi/180)-(phiraw_given[i_phi]*math.pi/180))
        phi_error[i_phi]=math.acos(math.cos((phic_given[i_phi]*math.pi/180)-(phiraw_given[i_phi]*math.pi/180)))
      else:
        phi_cc[i_phi]=f_sq*math.cos(phic_given[i_phi]-phiraw_given[i_phi])
        phi_error[i_phi]=math.acos(math.cos(phic_given[i_phi]-phiraw_given[i_phi]))

    mean_phicc=sum(phi_cc)/sum_magf_sq
    mean_phierr=sum(phi_error)/len(phi_error)

    return mean_phicc, mean_phierr

  def calcphibar(self,phi_set):
    '''
    calculate centroid phases for given set
    input data structure
    [[phia_refl_1, phia_refl_2, phia_refl_3...],
    [phib_refl_1,phib_refl_2,phib_refl_3..],
    ...]
    average phia_refl_1, phib_refl_1,...
    '''
    import cmath

    n_pop_size=len(phi_set)
    n_refl=len(phi_set[0])

    flex_phi_bar=flex.double(n_refl)
    txt_phi_bar=""
    for i in range(n_refl):
      sum_phis=cmath.rect(0,0)

      for j in range(n_pop_size):
        sum_phis+=cmath.rect(1,phi_set[j][i]*math.pi/180)

      flex_phi_bar[i]=cmath.phase(sum_phis)*180/math.pi
      txt_phi_bar+=str(flex_phi_bar[i])+","

    txt_phi_bar+="\n"

    return flex_phi_bar,txt_phi_bar
