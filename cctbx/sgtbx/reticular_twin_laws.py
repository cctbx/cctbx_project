from __future__ import absolute_import, division, print_function
from cctbx import sgtbx
from cctbx import uctbx
from cctbx import crystal
from cctbx.array_family import flex
import cctbx.sgtbx.lattice_symmetry
from cctbx.sgtbx import reticular_pg_tools as rmpg
from cctbx.sgtbx import sub_lattice_tools as slt
from scitbx import matrix
import math, sys
import scitbx.math
from six.moves import zip

class symmetry_safe_sublattice_xs(object):
  def __init__(self,xsin,start_order=1,stop_order=5,max_delta=5.0):
    self.xs = xsin
    self.cb_op_to_niggli = xsin.change_of_basis_op_to_niggli_cell()
    self.xs_n = self.xs.change_basis( self.cb_op_to_niggli )
    self.basis = matrix.sqr( self.xs_n.unit_cell().orthogonalization_matrix() )

    self.max_delta = max_delta

    self.xs_expand = []
    self.matrices  = slt.generate_matrix_up_to_order( stop_order, start_order )
    self.add_to_niggli = []
    self.metric_r_values = []
    for mat in self.matrices:
      self.make_new_xs(mat)

  def metric_r(self,a,b):
    old = a.metrical_matrix()
    new = b.metrical_matrix()
    old = flex.double( old )
    new = flex.double( new )
    delta = flex.abs(old - new)
    result = 0
    top = 0
    bottom = 0
    for oi, ni in zip(old,new):
      tmp = math.fabs(oi-ni)
      top += tmp
      bottom += math.fabs(oi)
    assert (bottom>0)
    delta = 100.0*top/bottom
    return( delta )

  def make_new_xs(self, mat):
    new_basis = self.basis*mat.as_float()
    new_uc = uctbx.unit_cell( orthogonalization_matrix = new_basis )
    new_xs = crystal.symmetry( new_uc, "P1" )
    ecbn = new_xs.change_of_basis_op_to_niggli_cell()
    cb_new_xs = new_xs.change_basis( ecbn )

    new_xs = crystal.symmetry( unit_cell=cb_new_xs.unit_cell(),
                               space_group=sgtbx.lattice_symmetry.group(cb_new_xs.unit_cell(),self.max_delta),
                               assert_is_compatible_unit_cell=False,
                             )
    # gather some results please
    self.metric_r_values.append( self.metric_r(new_xs.unit_cell(),cb_new_xs.unit_cell() )   )
    self.xs_expand.append( new_xs )
    self.add_to_niggli.append( ecbn )

class ret_twin_law_info(object):
  def __init__(self, twin_laws, m, metric_r, base_xs, subl_xs ):
    self.m = m
    self.metric_r = metric_r
    self.twin_laws = twin_laws
    self.order = m.determinant()
    self.base_xs = base_xs
    self.subl_xs = subl_xs

  def show(self, out=None):
    if out is None:
      out = sys.stdout
    print(file=out)
    print("----- Reticular twin laws summary ----- ", file=out)
    print(file=out)
    print("Base crystal symmetry", file=out)
    self.base_xs.show_summary(f=out)
    print(file=out)
    print("Sublattice symmetry", file=out)
    self.subl_xs.show_summary(f=out)
    print(file=out)
    print("    Matrix M : ", file=out)
    print(self.m)
    print("    Order    : %i"%self.order, file=out)
    print("    Metric R : %3.2e"%self.metric_r, file=out)
    print("    Twin laws:", file=out)
    if self.twin_laws is not None:
      for law in self.twin_laws:
        print("           ", rmpg.as_hkl( law ), file=out)
    else:
      print("         No twin laws", file=out)
    print(file=out)

class reticular_twin_laws(object):
  def __init__(self, xs, max_delta=5.0, max_index=5):
    self.xs1 = xs
    cbop_prim = self.xs1.change_of_basis_op_to_niggli_cell()
    self.xs2 = self.xs1.change_basis(cbop_prim)
    self.cbop =  rmpg.cb_op_as_rational(cbop_prim)
    self.xs_sl_list = symmetry_safe_sublattice_xs(self.xs2,1,max_index,max_delta)
    self.base_to_niggli_inv = rmpg.cb_op_as_rational( self.xs_sl_list.cb_op_to_niggli).inverse()

    self.ori_sg = rmpg.construct_rational_point_group( self.xs2.space_group() )
    self.ori_sg.change_basis( self.cbop.inverse() )
    self.derived_laws = []

    for nxs, mat, add_cb_op, mr in zip(self.xs_sl_list.xs_expand,
                                       self.xs_sl_list.matrices,
                                       self.xs_sl_list.add_to_niggli,
                                       self.xs_sl_list.metric_r_values):
      these_laws=  self.construct_twin_laws( nxs, mat, add_cb_op, mr )
      if these_laws is not None:
        this_info =   self.construct_twin_laws( nxs, mat, add_cb_op, mr )
        self.derived_laws.append( this_info )


  def construct_twin_laws(self, xs, mat, additional_cb_op, r, show=False):
    # get the new symmetry object
    # A.cb(nig).cb(mat).cb(add_to_niggli)
    tmp_sg = xs.reflection_intensity_symmetry(anomalous_flag = self.xs2.space_group().is_chiral() ).space_group()
    new_sg = rmpg.construct_rational_point_group( tmp_sg )
    tmp_cb_op = rmpg.cb_op_as_rational( additional_cb_op ).inverse()

    new_sg.change_basis( tmp_cb_op )
    new_sg.change_basis( mat.inverse() )
    new_sg.change_basis( self.base_to_niggli_inv )
    new_sg.change_basis( self.cbop.inverse() )

    these_twin_laws = rmpg.build_reticular_twin_laws(self.ori_sg,new_sg )
    if these_twin_laws is not None:
      tlinfo = ret_twin_law_info(these_twin_laws,mat,r,self.xs1,xs )
      return tlinfo
    else:
      return None





  def show(self,out=None):
    if out is None:
      out=sys.stdout

    if len(self.derived_laws):
      for o in self.derived_laws:
        o.show(out=out)
    else:
      print("No reticular twin laws found.", file=out)



def exercise():
  """
  uc = uctbx.unit_cell( "10.079 10.079 48.409 90 90 120" )
  xs = crystal.symmetry( uc, "R32")
  reticular_twin_laws( xs , max_index=3, max_delta=1.5).show()

  uc = uctbx.unit_cell( "5.05 11.15 11.45 108.25 98.42 95.78" )
  xs = crystal.symmetry( uc, "P1" )
  reticular_twin_laws( xs, max_index=3 ).show()

  uc = uctbx.unit_cell( "127.6, 58.1, 51.2, 90, 97.2, 90" )
  xs = crystal.symmetry( uc, "C2" )
  reticular_twin_laws( xs, max_index=3, max_delta=1.5 ).show()

  uc = uctbx.unit_cell( "127.6, 152.1, 51.2, 90, 90.0, 90" )
  xs = crystal.symmetry( uc, "P222" )
  reticular_twin_laws( xs, max_index=3, max_delta=1.5 ).show()
  """

  uc  =  uctbx.unit_cell( "96 74 77 90 113 90")
  xs  = crystal.symmetry(uc, "C2")
  trlw = reticular_twin_laws( xs, max_index=2, max_delta=1.5 )
  #trlw.show()
  m = [2, 1, 0,0, 1, 0,0, 0, 1]
  for ii, jj in zip( trlw.derived_laws[0].m, m):
    assert(ii==jj)

  tl = [-1, 0, 0, 0, -1, 0,-1, 0, 1]
  dl = trlw.derived_laws[0].twin_laws[0]
  for ii, jj in zip(tl,dl):
    assert(ii==jj)

  uc = uctbx.unit_cell( "10.079 10.079 48.409 90 90 120" )
  xs = crystal.symmetry( uc, "R32")
  rtl = reticular_twin_laws( xs , max_index=8, max_delta=1.5)
  for tli in rtl.derived_laws:
    for ttl in tli.twin_laws:
      assert(ttl.determinant()==1) #check that these twin laws have det equal to 1



if (__name__ == "__main__"):
  exercise()
  print("OK")
