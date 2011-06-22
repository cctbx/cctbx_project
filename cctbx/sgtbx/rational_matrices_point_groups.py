from cctbx import sgtbx
from cctbx.sgtbx import sub_lattice_tools
from scitbx import matrix
from boost import rational
from libtbx.math_utils import ifloor

def as_hkl( op ):
  def row_as_hkl( row, txt=['h','k','l']):
      result = ""
      part = 0
      for n,j in zip(row,txt):
        nn=""
        if n != rational.int(0):
          part += 1
          if n >rational.int(0):
            if part==1:
              if n==rational.int(1):
                nn = j
              else:
                nn = str(n)+j
            if part > 1:
              if n==rational.int(1):
                nn = "+"+j
              else:
                nn = "+"+str(n)+j
          if n < rational.int(0) :
            if part==1:
              if n==rational.int(-1):
                nn = "-"+j
              else:
                nn = str(n)+j
            else:
              if n == rational.int(-1):
                nn+="-"+j
              else:
                nn = str(n)+j
        result += nn
      return result

  hklmat=op.as_mat3()
  hkl = row_as_hkl(hklmat[0:3]) + "," + row_as_hkl(hklmat[3:6]) +","+ row_as_hkl(hklmat[6:])
  return hkl

class rat_rot_group(object):
  def __init__(self):
    self.unit = matrix.sqr( [rational.int(1),
                               rational.int(0),
                               rational.int(0),
                               rational.int(0),
                               rational.int(1),
                               rational.int(0),
                               rational.int(0),
                               rational.int(0),
                               rational.int(1)
                               ]  )

    self.ops = [ self.unit ]
    self.max_ops = 100

  def expand(self, new_op=None,show=False):
    if new_op==None:
      new_op = self.ops[0]
    if not self.is_in_group( new_op ):
      self.ops.append( new_op )
    done=False
    while not done:
      new_ops = []
      for op in self.ops:
        pno = op*new_op
        if show:
          print pno, op, new_op
        if not self.is_in_group( pno ):
          new_ops.append( pno )
      if len(new_ops)==0:
        done=True
      else:
        self.ops = self.ops+new_ops
      assert len(self.ops)<= self.max_ops

  def is_in_group(self,this_op):
    inside=False
    for op in self.ops:
      if op==this_op:
        inside=True
    return inside

  def change_basis(self, cb_op):
    new_ops = []
    for op in self.ops:
      new_op = cb_op.inverse()*op*cb_op
      new_ops.append( new_op )
    self.ops = new_ops

  def show(self):
    for op in self.ops:
      print as_hkl( op.transpose() )


def rt_mx_as_rational(rot_mat):
  # make sure one provides an integer matrix!
  tmp_mat = rot_mat.num()
  rational_list = []
  for ij in tmp_mat:
    rational_list.append( rational.int( ifloor(ij) ) )
  return matrix.sqr( rational_list )

def cb_op_as_rational(cb_op):
   num = cb_op.c().r().num()
   den = cb_op.c().r().den()
   rational_list = []
   for ii in num:
     rational_list.append( rational.int(ii)/rational.int(den) )
   return matrix.sqr( rational_list ).inverse()

def construct_rational_point_group( space_group,  cb_op=None ):
  gr = rat_rot_group()
  if cb_op is None:
    cb_op = sgtbx.change_of_basis_op(  "a,b,c" )
    cb_op = cb_op_as_rational( cb_op )
  for s in space_group:
    tmp_r = rt_mx_as_rational( s.r() )
    gr.expand( tmp_r )
  gr.change_basis( cb_op )
  gr.expand(show=False)
  return gr


def compare_groups( sg1, sg2 ):
  if len(sg1.ops) == len(sg2.ops):
    equal=True
    for op1 in sg1.ops:
      this_one = False
      for op2 in sg2.ops:
        if op1 == op2:
          this_one = True
          break
      if not this_one:
        equal = False
    return equal
  else:
    return False

def tst_compare():
  sg1 = construct_rational_point_group( sgtbx.space_group_info( "P 2 2 2" ).group() )
  sg2 = construct_rational_point_group( sgtbx.space_group_info( "P 2 2 2" ).group() )
  sg3 = construct_rational_point_group( sgtbx.space_group_info( "P 4 2 2" ).group() )
  sg4 = construct_rational_point_group( sgtbx.space_group_info( "P 2 2 2 (a+b,a-b,2c)").group() )
  assert compare_groups( sg1, sg2 )
  assert not compare_groups( sg1, sg3 )
  assert not compare_groups( sg1, sg4 )

def tst_groups():
  cb_ops = sub_lattice_tools.generate_cb_op_up_to_order(7)
  mats = sub_lattice_tools.generate_matrix_up_to_order(7)
  base_group = sgtbx.space_group_info( "P 2 2 2" ).group()
  for cb_op, mat in zip(cb_ops, mats):
    rat_cb_op = mat
    extended_group=None
    try:
      extended_group = sgtbx.space_group_info( "P 2 2 2 (%s)"%cb_op ).group()
    except Exception: pass
    rbg = construct_rational_point_group( base_group, rat_cb_op )
    reg = None
    if extended_group is not None:
      reg = construct_rational_point_group( extended_group )
      assert compare_groups(reg, rbg)


def run():
  tst_groups()
  tst_compare()

if __name__ == "__main__":
  run()
  print "OK"
