from __future__ import absolute_import, division, print_function
from scitbx.linalg import eigensystem
from scitbx.array_family import flex
from libtbx.utils import Sorry
from libtbx.utils import null_out
from scitbx import matrix
from math import acos,pi
from iotbx import pdb
import iotbx.pdb

class fab_elbow_angle(object):

  def __init__(self,
               pdb_hierarchy,
               chain_id_light='L',
               chain_id_heavy='H',
               limit_light=107,
               limit_heavy=113):
    '''
    Get elbow angle for Fragment antigen-binding (Fab)

    - Default heavy and light chains IDs are: H : heavy,  L : light
    - Default limit (cutoff) between variable and constant parts
      are residue number 107/113 for light/heavy chains
    - Variable domain is from residue 1 to limit.
      Constant domain form limit+1 to end.

    Reference:
    ----------
    Stanfield, et al., JMB 2006

    Usage example:
    --------------
    >>>fab = fab_elbow_angle(pdb_hierarchy=ph,limit_light=114,limit_heavy=118)
    >>>fab_angle = fab.fab_elbow_angle
    '''
    # create selection strings for the heavy/light var/const part of chains
    self.select_str(
      chain_ID_H=chain_id_heavy,
      limit_H=limit_heavy,
      chain_ID_L=chain_id_light,
      limit_L=limit_light)
    # get the hierarchy for and divide using selection strings
    self.pdb_hierarchy = pdb_hierarchy
    self.get_pdb_chains()
    # Get heavy to light reference vector before alignment !!!
    vh_end = self.pdb_var_H.atoms()[-1].xyz
    vl_end = self.pdb_var_L.atoms()[-1].xyz
    mid_H_to_L = self.norm_vec(start=vh_end,end=vl_end)
    # Get transformations objects
    tranformation_const= self.get_transformation(
      fixed_selection=self.pdb_const_H,
      moving_selection=self.pdb_const_L)
    tranformation_var = self.get_transformation(
      fixed_selection=self.pdb_var_H,
      moving_selection=self.pdb_var_L)
    # Get the angle and eigenvalues
    eigen_const = eigensystem.real_symmetric(tranformation_const.r.as_sym_mat3())
    eigen_var = eigensystem.real_symmetric(tranformation_var.r.as_sym_mat3())
    # c : consttant, v : variable
    eigenvectors_c = self.get_eigenvector(eigen_const)
    eigenvectors_v = self.get_eigenvector(eigen_var)
    # test eignevectors pointing in oposite directions
    if eigenvectors_c.dot(eigenvectors_v) > 0:
      eigenvectors_v = - eigenvectors_v
    # Calc Feb elbow angle
    angle = self.get_angle(vec1=eigenvectors_c, vec2=eigenvectors_v)
    # Test if elbow angle larger or smaller than 180
    zaxis = self.cross_product_as_unit_axis(eigenvectors_v, eigenvectors_c)
    xaxis = self.cross_product_as_unit_axis(eigenvectors_c, zaxis)
    if mid_H_to_L.dot(zaxis) <= 0:
        angle = 360 - angle
    self.fab_elbow_angle = angle
    # The cosine of the angles the vector from limit_heavy to limit_light
    # make with the axes x, y (eigenvectors_c), z
    # where the eigenvectors_v lies in the plane x - y
    self.cos_H_to_L_with_xaxis = mid_H_to_L.dot(xaxis)
    self.cos_H_to_L_with_yaxis = mid_H_to_L.dot(eigenvectors_c)
    self.cos_H_to_L_with_zaxis = mid_H_to_L.dot(zaxis)

  def norm_vec(self,start,end):
    '''retruns normalized vector that starts at "stat" and ends at "end"'''
    x = flex.double(end) - flex.double(start)
    l = x.norm()
    if l != 0:
      x = x/l
    return x

  def cross_product_as_unit_axis(self,a,b):
    a = tuple(a)
    b = tuple(b)
    x = flex.double(matrix.cross_product_matrix(a)*b)
    l = x.norm()
    if l != 0:
      x = x/l
    return x

  def get_angle(self,vec1,vec2,larger=True):
    '''retrun the larger angle between vec1 and vec2'''
    if vec1 and vec1:
      angle_cos = vec1.dot(vec2)
      acos_angle_cos = acos(angle_cos)
      assert acos_angle_cos != 0
      angle = 180/pi*acos_angle_cos
    else:
      angle = 0
    if (angle < 90) and larger: angle = 180 - angle
    if (angle > 90) and not larger: angle = 180 - angle
    return angle

  def get_eigenvector(self,eigen):
    '''
    Get the eigen vector for eigen value 1 and normalize it
    '''
    v = eigen.vectors()
    e = eigen.values()
    indx = None
    # select eigenvector that corespondes to a real egienvalue == 1
    for i,x in enumerate(e):
      if not isinstance(x,complex):
        if abs(1-x)<1e-6:
          indx = i
          break
    # make sure we have egienvalue == 1
    assert not indx
    eigenvector = v[indx:indx+3]
    # normalize
    eigenvector = eigenvector / eigenvector.dot(eigenvector)
    if e.all_eq(flex.double([1,1,1])):
      eigenvector = None
    return eigenvector

  def get_pdb_chains(self):
    '''Create seperate pdb hierarchy for each on the chains we want to align'''
    ph = self.pdb_hierarchy
    # test selection
    test = ph.atom_selection_cache().selection
    #
    self.pdb_var_H = ph.select(test(self.select_var_str_H))
    self.pdb_const_H = ph.select(test(self.select_const_str_H))
    self.pdb_var_L = ph.select(test(self.select_var_str_L))
    self.pdb_const_L = ph.select(test(self.select_const_str_L))
    if test(self.select_var_str_L).count(True) == 0:
      raise Sorry('No atoms in light chain. Check chain name')
    if test(self.select_var_str_H).count(True) == 0:
      raise Sorry('No atoms in heavy chain. Check chain name')

  def get_transformation(self,fixed_selection,moving_selection):
    from phenix.command_line import superpose_pdbs
    '''
    Align the moving pdb hierarchy on to the fixed one.
    Provides an object with rotation and translation info
    '''
    params = superpose_pdbs.master_params.extract()
    x = superpose_pdbs.manager(
      params,
      log=null_out(),
      write_output=False,
      save_lsq_fit_obj=True,
      pdb_hierarchy_fixed=fixed_selection,
      pdb_hierarchy_moving=moving_selection)
    return x.lsq_fit_obj

  def select_str(self,chain_ID_H,limit_H,chain_ID_L,limit_L):
    '''create selection strings for the heavy and light chains
    separating the variable and constant parts of the chains'''
    s1 = 'pepnames and (name ca or name n or name c) and altloc " "'
    s2 = 'chain {0} and resseq {1}:{2} and {3}'
    self.select_var_str_H = s2.format(chain_ID_H,1,limit_H,s1)
    self.select_const_str_H = s2.format(chain_ID_H,limit_H+1,'end',s1)
    self.select_var_str_L = s2.format(chain_ID_L,1,limit_L,s1)
    self.select_const_str_L = s2.format(chain_ID_L,limit_L+1,'end',s1)
